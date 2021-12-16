import os
import sys
from copy import copy

from itertools import chain

import logging
logger = logging.getLogger(__name__) 

import bamboo
from bamboo.analysismodules import HistogramsModule, DataDrivenBackgroundHistogramsModule

from bamboo import treefunctions as op
from bamboo.plots import CutFlowReport, Plot, EquidistantBinning, SummedPlot
from bamboo.analysisutils import forceDefine

from bamboo.root import gbl 

sys.path.append(os.path.dirname(os.path.abspath(__file__))) # Add scripts in this directory
from BaseHHtobbWW import BaseNanoHHtobbWW
from plotDef import *
from selectionDef import *
from DDHelper import DataDrivenFake, DataDrivenDY
import mvaEvaluatorDL_nonres
import mvaEvaluatorDL_res

def switch_on_index(indexes, condition, contA, contB):
    if contA._base != contB._base:
        raise RuntimeError("The containers do not derive from the same base, this won't work")
    base = contA._base
    return [base[op.switch(condition, contA[index].idx, contB[index].idx)] for index in indexes]       
#===============================================================================================#
#                                       PlotterHHtobbWW                                         #
#===============================================================================================#
class PlotterNanoHHtobbWWDL(BaseNanoHHtobbWW,DataDrivenBackgroundHistogramsModule):
    """ Plotter module: HH->bbW(->e/µ nu)W(->e/µ nu) histograms from NanoAOD """
    def __init__(self, args):
        super(PlotterNanoHHtobbWWDL, self).__init__(args)

    def initialize(self):
        # Change the way the FakeExtrapolation is postProcessed (avoids overriding the `postProcess` method) 
        super(PlotterNanoHHtobbWWDL, self).initialize()
        if "FakeExtrapolation" in self.datadrivenContributions:
            contrib = self.datadrivenContributions["FakeExtrapolation"]
            self.datadrivenContributions["FakeExtrapolation"] = DataDrivenFake(contrib.name, contrib.config)
        if "DYEstimation" in self.datadrivenContributions: 
            contrib = self.datadrivenContributions["DYEstimation"]
            self.datadrivenContributions["DYEstimation"] = DataDrivenDY(contrib.name, contrib.config,"PseudoData" in self.datadrivenContributions)

    def definePlots(self, t, noSel, sample=None, sampleCfg=None): 
        noSel = super(PlotterNanoHHtobbWWDL,self).prepareObjects(t, noSel, sample, sampleCfg, 'DL')
        #---- Parameters -----#

        plots = []

        if hasattr(self,'base_plots'):
            plots.extend(self.base_plots)

        era = sampleCfg['era']

        #----- Machine Learning Model -----#                
        if self.args.analysis == 'nonres':
            model_num = "11"
            path_model = os.path.join(os.path.abspath(os.path.dirname(__file__)),'MachineLearning','ml-models','models','multi-classification','dnn',model_num,'model','model.pb')
            input_names = ["lep","jet","fat","met","hl","param","eventnr"]
            output_name = "Identity"

            print ("DNN model : %s"%path_model)
            if not os.path.exists(path_model):
                raise RuntimeError('Could not find model file %s'%path_model)
            try:
                DNN = op.mvaEvaluator(path_model,mvaType='Tensorflow',otherArgs=(input_names, output_name))
            except:
                raise RuntimeError('Could not load model %s'%path_model)

        if self.args.analysis == 'res':
            dnn_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)),'MachineLearning','ResonantModels')
            path_model_HighMass = os.path.join(dnn_dir,'Resonant_HighMass_Final_512x4_w0p1.pb')
            path_model_LowMass = os.path.join(dnn_dir,'Resonant_LowMass_Final_512x4_w1.pb')
            input_names_HighMass = []
            input_names_LowMass = []
            with open(os.path.join(dnn_dir,'Resonant_HighMass_Final_512x4_w0p1_inputs.txt'),'r') as handle:
                for line in handle:
                    input_names_HighMass.append(line.split()[0])
            with open(os.path.join(dnn_dir,'Resonant_LowMass_Final_512x4_w1_inputs.txt'),'r') as handle:
                for line in handle:
                    input_names_LowMass.append(line.split()[0])
            output_name = "Identity"

            print ("DNN model : %s"%path_model_HighMass)
            print ("DNN model : %s"%path_model_LowMass)
            if not os.path.exists(path_model_HighMass):
                raise RuntimeError('Could not find model file %s'%path_model_HighMass)
            if not os.path.exists(path_model_LowMass):
                raise RuntimeError('Could not find model file %s'%path_model_LowMass)
            try:
                DNN_HighMass = op.mvaEvaluator(path_model_HighMass,mvaType='Tensorflow',otherArgs=(input_names_HighMass, output_name))
            except:
                raise RuntimeError('Could not load model %s'%path_model_HighMass)
            try:
                DNN_LowMass = op.mvaEvaluator(path_model_LowMass,mvaType='Tensorflow',otherArgs=(input_names_LowMass, output_name))
            except:
                raise RuntimeError('Could not load model %s'%path_model_LowMass)

        #----- Dilepton selection -----#
        ElElSelObj,MuMuSelObj,ElMuSelObj = makeDoubleLeptonSelection(self,noSel)


        # Select the jets selections that will be done depending on user input #
        jet_level = ["Ak4","Ak8","Resolved0Btag","Resolved1Btag","Resolved2Btag","Boosted0Btag","Boosted1Btag"]
        jetplot_level = [arg for (arg,boolean) in self.args.__dict__.items() if arg in jet_level and boolean]
        if len(jetplot_level) == 0:  
            jetplot_level = jet_level # If nothing said, will do all
        jetsel_level = copy(jetplot_level)  # A plot level might need a previous selection that needs to be defined but not necessarily plotted
        if "Resolved0Btag" in jetsel_level or "Resolved1Btag" in jetsel_level or "Resolved2Btag" in jetsel_level:
            jetsel_level.append("Ak4") # Resolved needs the Ak4 selection
        if "Boosted0Btag" in jetsel_level or "Boosted1Btag" in jetsel_level:
            jetsel_level.append("Ak8") # Boosted needs the Ak8 selection

        #----- Select correct dilepton -----#
        OSElElDilepton = self.ElElFakeSel
        OSMuMuDilepton = self.MuMuFakeSel
        OSElMuDilepton = self.ElMuFakeSel

        self.beforeJetselection(ElElSelObj,'ElEl')
        self.beforeJetselection(MuMuSelObj,'MuMu')
        self.beforeJetselection(ElMuSelObj,'ElMu')

        #----- Boolean to know what datadriven is applied -----#
        DYCR   = 'DYEstimation' in self.datadrivenContributions.keys() and self.datadrivenContributions['DYEstimation'].usesSample(self.sample,self.sampleCfg)
        FakeCR = 'FakeExtrapolation'  in self.datadrivenContributions.keys() and self.datadrivenContributions['FakeExtrapolation'].usesSample(self.sample,self.sampleCfg) 

        #----- HME -----#
        # Must be done to branching of the RDF for jets (weight or cut) to avoid systematics being recomputed in the HME
        # But after the forceDefine on jet+met (in self.beforeJetselection)
        if self.args.analysis == 'res':
            #HME_resolved_per_channel = {'ElEl': self.computeResolvedHMEAfterLeptonSelections(sel     = ElElSelObj.sel,
            #                                                                                 l1      = OSElElDilepton[0][0],
            #                                                                                 l2      = OSElElDilepton[0][1],
            #                                                                                 bjets   = self.ak4JetsByBtagScore,
            #                                                                                 met     = self.corrMET)[0],
            #                            'MuMu': self.computeResolvedHMEAfterLeptonSelections(sel     = MuMuSelObj.sel,
            #                                                                                 l1      = OSMuMuDilepton[0][0],
            #                                                                                 l2      = OSMuMuDilepton[0][1],
            #                                                                                 bjets   = self.ak4JetsByBtagScore,
            #                                                                                 met     = self.corrMET)[0],
            #                            'ElMu': self.computeResolvedHMEAfterLeptonSelections(sel     = ElMuSelObj.sel,
            #                                                                                 l1      = OSElMuDilepton[0][0],
            #                                                                                 l2      = OSElMuDilepton[0][1],
            #                                                                                 bjets   = self.ak4JetsByBtagScore,
            #                                                                                 met     = self.corrMET)[0]}
            #HME_boosted_per_channel = {'ElEl': self.computeBoostedHMEAfterLeptonSelections(sel     = ElElSelObj.sel,
            #                                                                               l1      = OSElElDilepton[0][0],
            #                                                                               l2      = OSElElDilepton[0][1],
            #                                                                               fatjets = self.ak8BJets,
            #                                                                               met     = self.corrMET)[0],
            #                           'MuMu': self.computeBoostedHMEAfterLeptonSelections(sel     = MuMuSelObj.sel,
            #                                                                               l1      = OSMuMuDilepton[0][0],
            #                                                                               l2      = OSMuMuDilepton[0][1],
            #                                                                               fatjets = self.ak8BJets,
            #                                                                               met     = self.corrMET)[0],
            #                           'ElMu': self.computeBoostedHMEAfterLeptonSelections(sel     = ElMuSelObj.sel,
            #                                                                               l1      = OSElMuDilepton[0][0],
            #                                                                               l2      = OSElMuDilepton[0][1],
            #                                                                               fatjets = self.ak8BJets,
            #                                                                               met     = self.corrMET)[0]}
            treePath = '/home/ucl/cp3/fbury/scratch/SkimsHME/Skim{era}_{cat}_{channel}_{region}/results/{sample}.root'
            hmeReaders = {"SR":{},"FakeCR":{},"DYCR":{}}
            if self.args.Resolved1Btag:
                hmeReaders['SR']['Resolved1B']          = {'ElEl':None,'MuMu':None,'ElMu':None}
                if FakeCR:
                    hmeReaders['FakeCR']['Resolved1B']  = {'ElEl':None,'MuMu':None,'ElMu':None}
                if DYCR:
                    hmeReaders['DYCR']['Resolved1B']    = {'ElEl':None,'MuMu':None,'ElMu':None}
            if self.args.Resolved2Btag:
                hmeReaders['SR']['Resolved2B']          = {'ElEl':None,'MuMu':None,'ElMu':None}
                if FakeCR:
                    hmeReaders['FakeCR']['Resolved2B']  = {'ElEl':None,'MuMu':None,'ElMu':None}
                if DYCR:
                    hmeReaders['DYCR']['Resolved2B']    = {'ElEl':None,'MuMu':None,'ElMu':None}
            if self.args.Boosted1Btag:
                hmeReaders['SR']['Boosted1B']           = {'ElEl':None,'MuMu':None,'ElMu':None}
                if FakeCR:
                    hmeReaders['FakeCR']['Boosted1B']   = {'ElEl':None,'MuMu':None,'ElMu':None}
                if DYCR:
                    hmeReaders['DYCR']['Boosted1B']     = {'ElEl':None,'MuMu':None,'ElMu':None}

            for region in hmeReaders.keys():
                for cat in hmeReaders[region].keys():
                    for channel in hmeReaders[region][cat].keys():
                        if 'related-sample' in sampleCfg.keys():
                            pathToSkim = treePath.format(era     = self.era,
                                                         cat     = cat,
                                                         channel = channel, 
                                                         region  = region,
                                                         sample  = sampleCfg["related-sample"])
                        else:
                            pathToSkim = treePath.format(era     = self.era,
                                                         cat     = cat,
                                                         channel = channel, 
                                                         region  = region,
                                                         sample  = self.sample)
                        if not os.path.exists(pathToSkim):
                            raise RuntimeError(f'Could not find skim {pathToSkim}')
                        hmeReaders[region][cat][channel] = op.define("hme::HMEReader", 
                                                                     f'hme::HMEReader <<name>>{{"{pathToSkim}"}}; // for {self.sample.replace("-","")}',
                                                                     nameHint=f'bamboo_hmeReader{region}{cat}{channel}{self.sample.replace("-","")}')

        #----- Channel and trigger plots -----#

#        if not self.args.OnlyYield:
#            ChannelDictList = []
#            ChannelDictList.append({'channel':'ElEl','sel':ElElSelObj.sel,'suffix':ElElSelObj.selName})
#            ChannelDictList.append({'channel':'MuMu','sel':MuMuSelObj.sel,'suffix':MuMuSelObj.selName})
#            ChannelDictList.append({'channel':'ElMu','sel':ElMuSelObj.sel,'suffix':ElMuSelObj.selName})
#
#            for channelDict in ChannelDictList:
#                #----- Trigger plots -----#
#                plots.extend(doubleLeptonTriggerPlots(**channelDict,triggerDict=self.triggersPerPrimaryDataset))
#                #----- Dilepton plots -----#
#                #plots.extend(doubleLeptonChannelPlot(**channelDict,DilepElEl=OSElElDilepton,DilepMuMu=OSMuMuDilepton,DilepElMu=OSElMuDilepton))
#
#        LeptonKeys = ['channel','sel','dilepton','suffix','is_MC']
#        JetKeys    = ['channel','sel','leadjet','subleadjet','lead_is_b','sublead_is_b','suffix','is_MC']
#        commonItems = ['channel','sel','suffix']
            
        #----- Ak4 jets selection -----#
        if "Ak4" in jetsel_level:
            print ("...... Processing Ak4 jet selection")
            ElElSelObjAk4Jets = makeAtLeastTwoAk4JetSelection(self,ElElSelObj,copy_sel=True)
            MuMuSelObjAk4Jets = makeAtLeastTwoAk4JetSelection(self,MuMuSelObj,copy_sel=True)
            ElMuSelObjAk4Jets = makeAtLeastTwoAk4JetSelection(self,ElMuSelObj,copy_sel=True)

            if self.args.onlypost:
                ElElSelObjAk4Jets.record_yields = True
                MuMuSelObjAk4Jets.record_yields = True
                ElMuSelObjAk4Jets.record_yields = True
                ElElSelObjAk4Jets.yieldTitle = 'Channel $e^{+}e^{-}$'
                MuMuSelObjAk4Jets.yieldTitle = 'Channel $\mu^{+}\mu^{-}$'
                ElMuSelObjAk4Jets.yieldTitle = 'Channel $e^{\pm}\mu^{\mp}$'

#            # Jet and lepton plots #
#            ChannelDictList = []
#            if "Ak4" in jetplot_level:
#                # Cut flow report #
#                if not self.args.OnlyYield:
#                    ChannelDictList.append({'channel':'ElEl','sel':ElElSelObjAk4Jets.sel,'dilepton':OSElElDilepton[0],'leadjet':self.ak4Jets[0],'subleadjet':self.ak4Jets[1],'lead_is_b':False,'sublead_is_b':False,'suffix':ElElSelObjAk4Jets.selName,'is_MC':self.is_MC})
#                    ChannelDictList.append({'channel':'MuMu','sel':MuMuSelObjAk4Jets.sel,'dilepton':OSMuMuDilepton[0],'leadjet':self.ak4Jets[0],'subleadjet':self.ak4Jets[1],'lead_is_b':False,'sublead_is_b':False,'suffix':MuMuSelObjAk4Jets.selName,'is_MC':self.is_MC})
#                    ChannelDictList.append({'channel':'ElMu','sel':ElMuSelObjAk4Jets.sel,'dilepton':OSElMuDilepton[0],'leadjet':self.ak4Jets[0],'subleadjet':self.ak4Jets[1],'lead_is_b':False,'sublead_is_b':False,'suffix':ElMuSelObjAk4Jets.selName,'is_MC':self.is_MC})
#
#            JetsN = {'objName':'Ak4Jets','objCont':self.ak4Jets,'Nmax':15,'xTitle':'N(Ak4 jets)'}
#                                                        
#            for channelDict in ChannelDictList:
#                # Dilepton #
#                plots.extend(makeDileptonPlots(**{k:channelDict[k] for k in LeptonKeys}))
#                # Number of jets #
#                plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**JetsN))
#                # Ak4 Jets #
#                plots.extend(makeTwoAk4JetsPlots(**{k:channelDict[k] for k in JetKeys}))
#                # MET #
#                plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))

        ##### Ak8 jets selection #####
        if "Ak8" in jetsel_level:
            print ("...... Processing Ak8 jet selection")
            ElElSelObjAk8Jets = makeAtLeastOneAk8JetSelection(self,ElElSelObj,copy_sel=True)
            MuMuSelObjAk8Jets = makeAtLeastOneAk8JetSelection(self,MuMuSelObj,copy_sel=True)
            ElMuSelObjAk8Jets = makeAtLeastOneAk8JetSelection(self,ElMuSelObj,copy_sel=True)

            if self.args.onlypost:
                ElElSelObjAk8Jets.record_yields = True
                MuMuSelObjAk8Jets.record_yields = True
                ElMuSelObjAk8Jets.record_yields = True
                ElElSelObjAk8Jets.yieldTitle = 'Channel $e^{+}e^{-}$'
                MuMuSelObjAk8Jets.yieldTitle = 'Channel $\mu^{+}\mu^{-}$'
                ElMuSelObjAk8Jets.yieldTitle = 'Channel $e^{\pm}\mu^{\mp}$'
#            # Fatjets plots #
#            ChannelDictList = []
#            if "Ak8" in jetplot_level:
#                # Cut flow report #
#                if not self.args.OnlyYield:
#                    ChannelDictList.append({'channel':'ElEl','sel':ElElSelObjAk8Jets.sel,'dilepton':OSElElDilepton[0],'fatjet':self.ak8Jets[0],'suffix':ElElSelObjAk8Jets.selName,'is_MC':self.is_MC})
#                    ChannelDictList.append({'channel':'MuMu','sel':MuMuSelObjAk8Jets.sel,'dilepton':OSMuMuDilepton[0],'fatjet':self.ak8Jets[0],'suffix':MuMuSelObjAk8Jets.selName,'is_MC':self.is_MC})
#                    ChannelDictList.append({'channel':'ElMu','sel':ElMuSelObjAk8Jets.sel,'dilepton':OSElMuDilepton[0],'fatjet':self.ak8Jets[0],'suffix':ElMuSelObjAk8Jets.selName,'is_MC':self.is_MC})
#
#            FatJetKeys = ['channel','sel','fatjet','suffix']
#            FatJetsN = {'objName':'Ak8Jets','objCont':self.ak8Jets,'Nmax':5,'xTitle':'N(Ak8 jets)'}
#            
#            for channelDict in ChannelDictList:
#                # Dilepton #
#                plots.extend(makeDileptonPlots(**{k:channelDict[k] for k in LeptonKeys}))
#                # Number of jets #
#                plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**FatJetsN))
#                # Ak8 Jets #
#                plots.extend(makeDoubleLeptonAk8JetsPlots(**{k:channelDict[k] for k in FatJetKeys}))
#                # MET #
#                plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))
                
                     
        #----- Resolved selection -----#
        if any(item in ["Resolved0Btag","Resolved1Btag","Resolved2Btag"] for item in jetsel_level): # If any of resolved category is asked
            #----- select the jets -----#
            aka4JetsByBtagScore = op.sort(self.ak4Jets, lambda jet : -jet.btagDeepFlavB)

            container0b2j = [ t.Jet[aka4JetsByBtagScore[i].idx] for i in range(2) ]
            container1b1j = [ t.Jet[op.switch(op.rng_len(self.ak4BJets) == 1, self.ak4BJets[0].idx,  aka4JetsByBtagScore[0].idx)] ,
                              t.Jet[op.switch(op.rng_len(self.ak4BJets) == 1, self.ak4LightJetsByBtagScore[0].idx,  aka4JetsByBtagScore[1].idx)]]
            container2b0j = [ t.Jet[op.switch(op.rng_len(self.ak4BJets) >= 2, self.ak4BJets[i].idx, aka4JetsByBtagScore[i].idx)] for i in range(2) ]

#            ChannelDictList = []
#            #----- Resolved selection : 0 Btag -----#
            if "Resolved0Btag" in jetsel_level:
                print ("...... Processing Resolved jet (0 btag) selection")
                ElElSelObjAk4JetsExclusiveResolvedNoBtag = makeExclusiveResolvedNoBtagSelection(self,ElElSelObjAk4Jets,copy_sel=True)
                MuMuSelObjAk4JetsExclusiveResolvedNoBtag = makeExclusiveResolvedNoBtagSelection(self,MuMuSelObjAk4Jets,copy_sel=True)
                ElMuSelObjAk4JetsExclusiveResolvedNoBtag = makeExclusiveResolvedNoBtagSelection(self,ElMuSelObjAk4Jets,copy_sel=True)
#
#                if "Resolved0Btag" in jetplot_level:
#                    # Cut flow report #
#                    if not self.args.OnlyYield:
#                        ChannelDictList.append({'channel':'ElEl','sel':ElElSelObjAk4JetsExclusiveResolvedNoBtag.sel,'dilepton':OSElElDilepton[0],'leadjet':container0b2j[0],'subleadjet':container0b2j[1],'lead_is_b':False,'sublead_is_b':False,'suffix':ElElSelObjAk4JetsExclusiveResolvedNoBtag.selName,'is_MC':self.is_MC})
#                        ChannelDictList.append({'channel':'MuMu','sel':MuMuSelObjAk4JetsExclusiveResolvedNoBtag.sel,'dilepton':OSMuMuDilepton[0],'leadjet':container0b2j[0],'subleadjet':container0b2j[1],'lead_is_b':False,'sublead_is_b':False,'suffix':MuMuSelObjAk4JetsExclusiveResolvedNoBtag.selName,'is_MC':self.is_MC})
#                        ChannelDictList.append({'channel':'ElMu','sel':ElMuSelObjAk4JetsExclusiveResolvedNoBtag.sel,'dilepton':OSElMuDilepton[0],'leadjet':container0b2j[0],'subleadjet':container0b2j[1],'lead_is_b':False,'sublead_is_b':False,'suffix':ElMuSelObjAk4JetsExclusiveResolvedNoBtag.selName,'is_MC':self.is_MC})
#            #----  Resolved selection : 1 Btag  -----#
            if "Resolved1Btag" in jetsel_level:
                print ("...... Processing Resolved jet (1 btag) selection")
                ElElSelObjAk4JetsExclusiveResolvedOneBtag = makeExclusiveResolvedOneBtagSelection(self,ElElSelObjAk4Jets,copy_sel=True)
                MuMuSelObjAk4JetsExclusiveResolvedOneBtag = makeExclusiveResolvedOneBtagSelection(self,MuMuSelObjAk4Jets,copy_sel=True)
                ElMuSelObjAk4JetsExclusiveResolvedOneBtag = makeExclusiveResolvedOneBtagSelection(self,ElMuSelObjAk4Jets,copy_sel=True)
#
#                if "Resolved1Btag" in jetplot_level:
#                    # Cut flow report #
#                    if not self.args.OnlyYield:
#                        ChannelDictList.append({'channel':'ElEl','sel':ElElSelObjAk4JetsExclusiveResolvedOneBtag.sel,'dilepton':OSElElDilepton[0],'leadjet':container1b1j[0],'subleadjet':container1b1j[1],'lead_is_b':True,'sublead_is_b':False,'suffix':ElElSelObjAk4JetsExclusiveResolvedOneBtag.selName,'is_MC':self.is_MC})
#                        ChannelDictList.append({'channel':'MuMu','sel':MuMuSelObjAk4JetsExclusiveResolvedOneBtag.sel,'dilepton':OSMuMuDilepton[0],'leadjet':container1b1j[0],'subleadjet':container1b1j[1],'lead_is_b':True,'sublead_is_b':False,'suffix':MuMuSelObjAk4JetsExclusiveResolvedOneBtag.selName,'is_MC':self.is_MC})
#                        ChannelDictList.append({'channel':'ElMu','sel':ElMuSelObjAk4JetsExclusiveResolvedOneBtag.sel,'dilepton':OSElMuDilepton[0],'leadjet':container1b1j[0],'subleadjet':container1b1j[1],'lead_is_b':True,'sublead_is_b':False,'suffix':ElMuSelObjAk4JetsExclusiveResolvedOneBtag.selName,'is_MC':self.is_MC})
#            #----- Resolved selection : 2 Btags -----#
            if "Resolved2Btag" in jetsel_level:
                print ("...... Processing Resolved jet (2 btags) selection")
                ElElSelObjAk4JetsExclusiveResolvedTwoBtags = makeExclusiveResolvedTwoBtagsSelection(self,ElElSelObjAk4Jets,copy_sel=True)
                MuMuSelObjAk4JetsExclusiveResolvedTwoBtags = makeExclusiveResolvedTwoBtagsSelection(self,MuMuSelObjAk4Jets,copy_sel=True)
                ElMuSelObjAk4JetsExclusiveResolvedTwoBtags = makeExclusiveResolvedTwoBtagsSelection(self,ElMuSelObjAk4Jets,copy_sel=True)
#
#                if "Resolved2Btag" in jetplot_level:
#                    # Cut flow report #
#                    if not self.args.OnlyYield:
#                        ChannelDictList.append({'channel':'ElEl','sel':ElElSelObjAk4JetsExclusiveResolvedTwoBtags.sel,'dilepton':OSElElDilepton[0],'leadjet':container2b0j[0],'subleadjet':container2b0j[1],'lead_is_b':True,'sublead_is_b':True,'suffix':ElElSelObjAk4JetsExclusiveResolvedTwoBtags.selName,'is_MC':self.is_MC})
#                        ChannelDictList.append({'channel':'MuMu','sel':MuMuSelObjAk4JetsExclusiveResolvedTwoBtags.sel,'dilepton':OSMuMuDilepton[0],'leadjet':container2b0j[0],'subleadjet':container2b0j[1],'lead_is_b':True,'sublead_is_b':True,'suffix':MuMuSelObjAk4JetsExclusiveResolvedTwoBtags.selName,'is_MC':self.is_MC})
#                        ChannelDictList.append({'channel':'ElMu','sel':ElMuSelObjAk4JetsExclusiveResolvedTwoBtags.sel,'dilepton':OSElMuDilepton[0],'leadjet':container2b0j[0],'subleadjet':container2b0j[1],'lead_is_b':True,'sublead_is_b':True,'suffix':ElMuSelObjAk4JetsExclusiveResolvedTwoBtags.selName,'is_MC':self.is_MC})
#
#            # Lepton + jet Plots #
#
#            ResolvedJetsN = {'objName':'Ak4BJets','objCont':self.ak4BJets,'Nmax':5,'xTitle':'N(Ak4 Bjets)'}
#    
#            for channelDict in ChannelDictList:
#                # Dilepton #
#                plots.extend(makeDileptonPlots(**{k:channelDict[k] for k in LeptonKeys}))
#                # Number of jets #
#                plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**ResolvedJetsN))
#                plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**JetsN))
#                # Ak4 Jets #
#                plots.extend(makeTwoAk4JetsPlots(**{k:channelDict[k] for k in JetKeys}))
#                # MET #
#                plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))
        
        #----- Boosted selection -----#
        if any(item in ["Boosted0Btag","Boosted1Btag"] for item in jetsel_level): # If any of the boosted category is asked 
            container1fatb = [t.FatJet[op.switch(op.rng_len(self.ak8BJets) >= 1, self.ak8BJets[0].idx, self.ak8Jets[0].idx)]]

            ChannelDictList = []
            #----- Boosted selection : 0 Btag -----#
            if "Boosted0Btag" in jetsel_level:
                print ("...... Processing Boosted jet (0 btag) selection")
                ElElSelObjAk8JetsInclusiveBoostedNoBtag = makeInclusiveBoostedNoBtagSelection(self,ElElSelObjAk8Jets,copy_sel=True)
                MuMuSelObjAk8JetsInclusiveBoostedNoBtag = makeInclusiveBoostedNoBtagSelection(self,MuMuSelObjAk8Jets,copy_sel=True)
                ElMuSelObjAk8JetsInclusiveBoostedNoBtag = makeInclusiveBoostedNoBtagSelection(self,ElMuSelObjAk8Jets,copy_sel=True)
#                if "Boosted0Btag" in jetplot_level:
#                    # Cut flow report #
#                    if not self.args.OnlyYield:
#                        ChannelDictList.append({'channel':'ElEl','sel':ElElSelObjAk8JetsInclusiveBoostedNoBtag.sel,'dilepton':OSElElDilepton[0],'fatjet':container1fatb[0],'suffix':ElElSelObjAk8JetsInclusiveBoostedNoBtag.selName,'is_MC':self.is_MC})
#                        ChannelDictList.append({'channel':'MuMu','sel':MuMuSelObjAk8JetsInclusiveBoostedNoBtag.sel,'dilepton':OSMuMuDilepton[0],'fatjet':container1fatb[0],'suffix':MuMuSelObjAk8JetsInclusiveBoostedNoBtag.selName,'is_MC':self.is_MC})
#                        ChannelDictList.append({'channel':'ElMu','sel':ElMuSelObjAk8JetsInclusiveBoostedNoBtag.sel,'dilepton':OSElMuDilepton[0],'fatjet':container1fatb[0],'suffix':ElMuSelObjAk8JetsInclusiveBoostedNoBtag.selName,'is_MC':self.is_MC})
#            #----- Boosted selection : 1 Btag -----#
            if "Boosted1Btag" in jetsel_level:
                print ("...... Processing Boosted jet (1 btag) selection")
                ElElSelObjAk8JetsInclusiveBoostedOneBtag = makeInclusiveBoostedOneBtagSelection(self,ElElSelObjAk8Jets,copy_sel=True)
                MuMuSelObjAk8JetsInclusiveBoostedOneBtag = makeInclusiveBoostedOneBtagSelection(self,MuMuSelObjAk8Jets,copy_sel=True)
                ElMuSelObjAk8JetsInclusiveBoostedOneBtag = makeInclusiveBoostedOneBtagSelection(self,ElMuSelObjAk8Jets,copy_sel=True)
#                if "Boosted1Btag" in jetplot_level:
#                    # Cut flow report #
#                    if not self.args.OnlyYield:
#                        ChannelDictList.append({'channel':'ElEl','sel':ElElSelObjAk8JetsInclusiveBoostedOneBtag.sel,'dilepton':OSElElDilepton[0],'fatjet':container1fatb[0],'suffix':ElElSelObjAk8JetsInclusiveBoostedOneBtag.selName,'is_MC':self.is_MC})
#                        ChannelDictList.append({'channel':'MuMu','sel':MuMuSelObjAk8JetsInclusiveBoostedOneBtag.sel,'dilepton':OSMuMuDilepton[0],'fatjet':container1fatb[0],'suffix':MuMuSelObjAk8JetsInclusiveBoostedOneBtag.selName,'is_MC':self.is_MC})
#                        ChannelDictList.append({'channel':'ElMu','sel':ElMuSelObjAk8JetsInclusiveBoostedOneBtag.sel,'dilepton':OSElMuDilepton[0],'fatjet':container1fatb[0],'suffix':ElMuSelObjAk8JetsInclusiveBoostedOneBtag.selName,'is_MC':self.is_MC})
#
#            BoostedJetsN = {'objName':'Ak8BJets','objCont':self.ak8BJets,'Nmax':5,'xTitle':'N(Ak8 Bjets)'}
#            
#            for channelDict in ChannelDictList:
#                # Dilepton #
#                plots.extend(makeDileptonPlots(**{k:channelDict[k] for k in LeptonKeys}))
#                # Number of jets #
#                plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**BoostedJetsN))
#                plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**FatJetsN))
#                # Ak8 Jets #
#                plots.extend(makeDoubleLeptonAk8JetsPlots(**{k:channelDict[k] for k in FatJetKeys}))
#                # MET #
#                plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))
#        
#        #----- High-level combinations -----#
#        # NOTE : very time consuming 
#        ChannelDictList = []
#        if not self.args.OnlyYield:
#            # Resolved No Btag #
#            if "Resolved0Btag" in jetplot_level:
#                ChannelDictList.append({'channel': 'ElEl','met': self.corrMET,'l1':OSElElDilepton[0][0],'l2':OSElElDilepton[0][1],'j1':container0b2j[0],'j2':container0b2j[1],'sel':ElElSelObjAk4JetsExclusiveResolvedNoBtag.sel,'suffix':ElElSelObjAk4JetsExclusiveResolvedNoBtag.selName})
#                ChannelDictList.append({'channel': 'MuMu','met': self.corrMET,'l1':OSMuMuDilepton[0][0],'l2':OSMuMuDilepton[0][1],'j1':container0b2j[0],'j2':container0b2j[1],'sel':MuMuSelObjAk4JetsExclusiveResolvedNoBtag.sel,'suffix':MuMuSelObjAk4JetsExclusiveResolvedNoBtag.selName})
#                ChannelDictList.append({'channel': 'ElMu','met': self.corrMET,'l1':OSElMuDilepton[0][0],'l2':OSElMuDilepton[0][1],'j1':container0b2j[0],'j2':container0b2j[1],'sel':ElMuSelObjAk4JetsExclusiveResolvedNoBtag.sel,'suffix':ElMuSelObjAk4JetsExclusiveResolvedNoBtag.selName})
#            # Resolved One Btag #
#            if "Resolved1Btag" in jetplot_level:
#                ChannelDictList.append({'channel': 'ElEl','met': self.corrMET,'l1':OSElElDilepton[0][0],'l2':OSElElDilepton[0][1],'j1':container1b1j[0],'j2':container1b1j[1],'sel':ElElSelObjAk4JetsExclusiveResolvedOneBtag.sel,'suffix':ElElSelObjAk4JetsExclusiveResolvedOneBtag.selName})
#                ChannelDictList.append({'channel': 'MuMu','met': self.corrMET,'l1':OSMuMuDilepton[0][0],'l2':OSMuMuDilepton[0][1],'j1':container1b1j[0],'j2':container1b1j[1],'sel':MuMuSelObjAk4JetsExclusiveResolvedOneBtag.sel,'suffix':MuMuSelObjAk4JetsExclusiveResolvedOneBtag.selName})
#                ChannelDictList.append({'channel': 'ElMu','met': self.corrMET,'l1':OSElMuDilepton[0][0],'l2':OSElMuDilepton[0][1],'j1':container1b1j[0],'j2':container1b1j[1],'sel':ElMuSelObjAk4JetsExclusiveResolvedOneBtag.sel,'suffix':ElMuSelObjAk4JetsExclusiveResolvedOneBtag.selName})
#            # Resolved Two Btags  #
#            if "Resolved2Btag" in jetplot_level:
#                ChannelDictList.append({'channel': 'ElEl','met': self.corrMET,'l1':OSElElDilepton[0][0],'l2':OSElElDilepton[0][1],'j1':container2b0j[0],'j2':container2b0j[1],'sel':ElElSelObjAk4JetsExclusiveResolvedTwoBtags.sel,'suffix':ElElSelObjAk4JetsExclusiveResolvedTwoBtags.selName})
#                ChannelDictList.append({'channel': 'MuMu','met': self.corrMET,'l1':OSMuMuDilepton[0][0],'l2':OSMuMuDilepton[0][1],'j1':container2b0j[0],'j2':container2b0j[1],'sel':MuMuSelObjAk4JetsExclusiveResolvedTwoBtags.sel,'suffix':MuMuSelObjAk4JetsExclusiveResolvedTwoBtags.selName})
#                ChannelDictList.append({'channel': 'ElMu','met': self.corrMET,'l1':OSElMuDilepton[0][0],'l2':OSElMuDilepton[0][1],'j1':container2b0j[0],'j2':container2b0j[1],'sel':ElMuSelObjAk4JetsExclusiveResolvedTwoBtags.sel,'suffix':ElMuSelObjAk4JetsExclusiveResolvedTwoBtags.selName})
#            # Boosted No Btag #
#            if "Boosted0Btag" in jetplot_level:
#                ChannelDictList.append({'channel': 'ElEl','met': self.corrMET,'l1':OSElElDilepton[0][0],'l2':OSElElDilepton[0][1],'j1':container1fatb[0].subJet1,'j2':container1fatb[0].subJet2,'sel':ElElSelObjAk8JetsInclusiveBoostedNoBtag.sel,'suffix':ElElSelObjAk8JetsInclusiveBoostedNoBtag.selName})
#                ChannelDictList.append({'channel': 'MuMu','met': self.corrMET,'l1':OSMuMuDilepton[0][0],'l2':OSMuMuDilepton[0][1],'j1':container1fatb[0].subJet1,'j2':container1fatb[0].subJet2,'sel':MuMuSelObjAk8JetsInclusiveBoostedNoBtag.sel,'suffix':MuMuSelObjAk8JetsInclusiveBoostedNoBtag.selName})
#                ChannelDictList.append({'channel': 'ElMu','met': self.corrMET,'l1':OSElMuDilepton[0][0],'l2':OSElMuDilepton[0][1],'j1':container1fatb[0].subJet1,'j2':container1fatb[0].subJet2,'sel':ElMuSelObjAk8JetsInclusiveBoostedNoBtag.sel,'suffix':ElMuSelObjAk8JetsInclusiveBoostedNoBtag.selName})
#            # Boosted One Btag #
#            if "Boosted1Btag" in jetplot_level:
#                ChannelDictList.append({'channel': 'ElEl','met': self.corrMET,'l1':OSElElDilepton[0][0],'l2':OSElElDilepton[0][1],'j1':container1fatb[0].subJet1,'j2':container1fatb[0].subJet2,'sel':ElElSelObjAk8JetsInclusiveBoostedOneBtag.sel,'suffix':ElElSelObjAk8JetsInclusiveBoostedOneBtag.selName})
#                ChannelDictList.append({'channel': 'MuMu','met': self.corrMET,'l1':OSMuMuDilepton[0][0],'l2':OSMuMuDilepton[0][1],'j1':container1fatb[0].subJet1,'j2':container1fatb[0].subJet2,'sel':MuMuSelObjAk8JetsInclusiveBoostedOneBtag.sel,'suffix':MuMuSelObjAk8JetsInclusiveBoostedOneBtag.selName})
#                ChannelDictList.append({'channel': 'ElMu','met': self.corrMET,'l1':OSElMuDilepton[0][0],'l2':OSElMuDilepton[0][1],'j1':container1fatb[0].subJet1,'j2':container1fatb[0].subJet2,'sel':ElMuSelObjAk8JetsInclusiveBoostedOneBtag.sel,'suffix':ElMuSelObjAk8JetsInclusiveBoostedOneBtag.selName})
#
#        for channelDict in ChannelDictList:
#            plots.extend(makeDoubleLeptonHighLevelQuantities(**channelDict,HLL=self.HLL))

        #----- Selected variables : Resolved -----#
        ChannelDictList = []
        if not self.args.OnlyYield:
            # Resolved No Btag #
            if "Resolved0Btag" in jetplot_level:
                ChannelDictList.append({'channel': 'ElEl','jets':self.ak4Jets,'met': self.corrMET,'l1':OSElElDilepton[0][0],'l2':OSElElDilepton[0][1],'b1':container0b2j[0],'b2':container0b2j[1],'sel':ElElSelObjAk4JetsExclusiveResolvedNoBtag.sel,'suffix':ElElSelObjAk4JetsExclusiveResolvedNoBtag.selName})
                ChannelDictList.append({'channel': 'MuMu','jets':self.ak4Jets,'met': self.corrMET,'l1':OSMuMuDilepton[0][0],'l2':OSMuMuDilepton[0][1],'b1':container0b2j[0],'b2':container0b2j[1],'sel':MuMuSelObjAk4JetsExclusiveResolvedNoBtag.sel,'suffix':MuMuSelObjAk4JetsExclusiveResolvedNoBtag.selName})
                ChannelDictList.append({'channel': 'ElMu','jets':self.ak4Jets,'met': self.corrMET,'l1':OSElMuDilepton[0][0],'l2':OSElMuDilepton[0][1],'b1':container0b2j[0],'b2':container0b2j[1],'sel':ElMuSelObjAk4JetsExclusiveResolvedNoBtag.sel,'suffix':ElMuSelObjAk4JetsExclusiveResolvedNoBtag.selName})
            if "Resolved1Btag" in jetplot_level:
                ChannelDictList.append({'channel': 'ElEl','jets':self.ak4Jets,'met': self.corrMET,'l1':OSElElDilepton[0][0],'l2':OSElElDilepton[0][1],'b1':container1b1j[0],'b2':container1b1j[1],'sel':ElElSelObjAk4JetsExclusiveResolvedOneBtag.sel,'suffix':ElElSelObjAk4JetsExclusiveResolvedOneBtag.selName})
                ChannelDictList.append({'channel': 'MuMu','jets':self.ak4Jets,'met': self.corrMET,'l1':OSMuMuDilepton[0][0],'l2':OSMuMuDilepton[0][1],'b1':container1b1j[0],'b2':container1b1j[1],'sel':MuMuSelObjAk4JetsExclusiveResolvedOneBtag.sel,'suffix':MuMuSelObjAk4JetsExclusiveResolvedOneBtag.selName})
                ChannelDictList.append({'channel': 'ElMu','jets':self.ak4Jets,'met': self.corrMET,'l1':OSElMuDilepton[0][0],'l2':OSElMuDilepton[0][1],'b1':container1b1j[0],'b2':container1b1j[1],'sel':ElMuSelObjAk4JetsExclusiveResolvedOneBtag.sel,'suffix':ElMuSelObjAk4JetsExclusiveResolvedOneBtag.selName})
            # Resolved Two Btags  #
            if "Resolved2Btag" in jetplot_level:
                ChannelDictList.append({'channel': 'ElEl','jets':self.ak4Jets,'met': self.corrMET,'l1':OSElElDilepton[0][0],'l2':OSElElDilepton[0][1],'b1':container2b0j[0],'b2':container2b0j[1],'sel':ElElSelObjAk4JetsExclusiveResolvedTwoBtags.sel,'suffix':ElElSelObjAk4JetsExclusiveResolvedTwoBtags.selName})
                ChannelDictList.append({'channel': 'MuMu','jets':self.ak4Jets,'met': self.corrMET,'l1':OSMuMuDilepton[0][0],'l2':OSMuMuDilepton[0][1],'b1':container2b0j[0],'b2':container2b0j[1],'sel':MuMuSelObjAk4JetsExclusiveResolvedTwoBtags.sel,'suffix':MuMuSelObjAk4JetsExclusiveResolvedTwoBtags.selName})
                ChannelDictList.append({'channel': 'ElMu','jets':self.ak4Jets,'met': self.corrMET,'l1':OSElMuDilepton[0][0],'l2':OSElMuDilepton[0][1],'b1':container2b0j[0],'b2':container2b0j[1],'sel':ElMuSelObjAk4JetsExclusiveResolvedTwoBtags.sel,'suffix':ElMuSelObjAk4JetsExclusiveResolvedTwoBtags.selName})

        for channelDict in ChannelDictList:
            # Branch out the LO -> NLO reweighting #
            channelDict["sel"] = self.addSignalReweighting(channelDict["sel"])
            plots.extend(makeDoubleLeptonSelectedResolvedVariables(**channelDict,HLL=self.HLL))

        #----- Selected variables : Boosted -----#
        ChannelDictList = []
        if not self.args.OnlyYield:
            # Boosted No Btag #
            if "Boosted0Btag" in jetplot_level:
                ChannelDictList.append({'channel': 'ElEl','jets':self.ak4Jets,'met': self.corrMET,'l1':OSElElDilepton[0][0],'l2':OSElElDilepton[0][1],'B':container1fatb[0],'sel':ElElSelObjAk8JetsInclusiveBoostedNoBtag.sel,'suffix':ElElSelObjAk8JetsInclusiveBoostedNoBtag.selName})
                ChannelDictList.append({'channel': 'MuMu','jets':self.ak4Jets,'met': self.corrMET,'l1':OSMuMuDilepton[0][0],'l2':OSMuMuDilepton[0][1],'B':container1fatb[0],'sel':MuMuSelObjAk8JetsInclusiveBoostedNoBtag.sel,'suffix':MuMuSelObjAk8JetsInclusiveBoostedNoBtag.selName})
                ChannelDictList.append({'channel': 'ElMu','jets':self.ak4Jets,'met': self.corrMET,'l1':OSElMuDilepton[0][0],'l2':OSElMuDilepton[0][1],'B':container1fatb[0],'sel':ElMuSelObjAk8JetsInclusiveBoostedNoBtag.sel,'suffix':ElMuSelObjAk8JetsInclusiveBoostedNoBtag.selName})
            # Boosted One Btag #
            if "Boosted1Btag" in jetplot_level:
                ChannelDictList.append({'channel': 'ElEl','jets':self.ak4Jets,'met': self.corrMET,'l1':OSElElDilepton[0][0],'l2':OSElElDilepton[0][1],'B':container1fatb[0],'sel':ElElSelObjAk8JetsInclusiveBoostedOneBtag.sel,'suffix':ElElSelObjAk8JetsInclusiveBoostedOneBtag.selName})
                ChannelDictList.append({'channel': 'MuMu','jets':self.ak4Jets,'met': self.corrMET,'l1':OSMuMuDilepton[0][0],'l2':OSMuMuDilepton[0][1],'B':container1fatb[0],'sel':MuMuSelObjAk8JetsInclusiveBoostedOneBtag.sel,'suffix':MuMuSelObjAk8JetsInclusiveBoostedOneBtag.selName})
                ChannelDictList.append({'channel': 'ElMu','jets':self.ak4Jets,'met': self.corrMET,'l1':OSElMuDilepton[0][0],'l2':OSElMuDilepton[0][1],'B':container1fatb[0],'sel':ElMuSelObjAk8JetsInclusiveBoostedOneBtag.sel,'suffix':ElMuSelObjAk8JetsInclusiveBoostedOneBtag.selName})

        if not self.args.OnlyYield:
            for channelDict in ChannelDictList:
                channelDict["sel"] = self.addSignalReweighting(channelDict["sel"])
                plots.extend(makeDoubleLeptonSelectedBoostedVariables(**channelDict,HLL=self.HLL))

        #----- Machine Learning plots -----#
        selObjectDictList = []
        if not self.args.OnlyYield:
            if "Resolved0Btag" in jetplot_level:
                selObjectDictList.append({'channel':'ElEl','selObject':ElElSelObjAk4JetsExclusiveResolvedNoBtag})
                selObjectDictList.append({'channel':'MuMu','selObject':MuMuSelObjAk4JetsExclusiveResolvedNoBtag})
                selObjectDictList.append({'channel':'ElMu','selObject':ElMuSelObjAk4JetsExclusiveResolvedNoBtag})
            if "Resolved1Btag" in jetplot_level:
                selObjectDictList.append({'channel':'ElEl','selObject':ElElSelObjAk4JetsExclusiveResolvedOneBtag})
                selObjectDictList.append({'channel':'MuMu','selObject':MuMuSelObjAk4JetsExclusiveResolvedOneBtag})
                selObjectDictList.append({'channel':'ElMu','selObject':ElMuSelObjAk4JetsExclusiveResolvedOneBtag})
            if "Resolved2Btag" in jetplot_level:
                selObjectDictList.append({'channel':'ElEl','selObject':ElElSelObjAk4JetsExclusiveResolvedTwoBtags})
                selObjectDictList.append({'channel':'MuMu','selObject':MuMuSelObjAk4JetsExclusiveResolvedTwoBtags})
                selObjectDictList.append({'channel':'ElMu','selObject':ElMuSelObjAk4JetsExclusiveResolvedTwoBtags})
            if "Boosted0Btag" in jetplot_level:
                selObjectDictList.append({'channel':'ElEl','selObject':ElElSelObjAk8JetsInclusiveBoostedNoBtag})
                selObjectDictList.append({'channel':'MuMu','selObject':MuMuSelObjAk8JetsInclusiveBoostedNoBtag})
                selObjectDictList.append({'channel':'ElMu','selObject':ElMuSelObjAk8JetsInclusiveBoostedNoBtag})
            if "Boosted1Btag" in jetplot_level:
                selObjectDictList.append({'channel':'ElEl','selObject':ElElSelObjAk8JetsInclusiveBoostedOneBtag})
                selObjectDictList.append({'channel':'MuMu','selObject':MuMuSelObjAk8JetsInclusiveBoostedOneBtag})
                selObjectDictList.append({'channel':'ElMu','selObject':ElMuSelObjAk8JetsInclusiveBoostedOneBtag})

        dileptonsCont = {'ElEl':OSElElDilepton[0],'MuMu':OSMuMuDilepton[0],'ElMu':OSElMuDilepton[0]}
        

        for selObjectDict in selObjectDictList:
            channel = selObjectDict['channel']
            category = None
            if 'ResolvedOneBtag' in selObjectDict['selObject'].selName:
                category = "Resolved1B"
            if 'ResolvedTwoBtags' in selObjectDict['selObject'].selName:
                category = "Resolved2B"
            if 'BoostedOneBtag' in selObjectDict['selObject'].selName:
                category = "Boosted1B"
            assert category is not None
            dilepton = dileptonsCont[channel]

            if self.args.analysis == 'nonres':
                self.nodes = ['GGF','VBF','H', 'DY', 'ST', 'TT', 'TTVX', 'VVV', 'Rare']

                inputsLeps = mvaEvaluatorDL_nonres.returnLeptonsMVAInputs(
                                                self      = self,
                                                l1        = dilepton[0],
                                                l2        = dilepton[1],
                                                channel   = selObjectDict['channel'])
                inputsJets = mvaEvaluatorDL_nonres.returnJetsMVAInputs(
                                                self      = self,
                                                jets      = self.ak4Jets)
                inputsMET = mvaEvaluatorDL_nonres.returnMETMVAInputs(
                                                self      = self,
                                                met       = self.corrMET)     
                inputsFatjet = mvaEvaluatorDL_nonres.returnFatjetMVAInputs(
                                                self      = self,
                                                fatjets   = self.ak8Jets)
                inputsHL = mvaEvaluatorDL_nonres.returnHighLevelMVAInputs(
                                                self      = self,
                                                l1        = dilepton[0],
                                                l2        = dilepton[1],
                                                met       = self.corrMET,
                                                jets      = self.ak4Jets,
                                                bjets     = self.ak4JetsByBtagScore[:op.min(op.rng_len(self.ak4JetsByBtagScore),op.static_cast("std::size_t",op.c_int(2)))],
                                                electrons = self.electronsTightSel,
                                                muons     = self.muonsTightSel,
                                                channel   = selObjectDict['channel'])

                inputsParam = mvaEvaluatorDL_nonres.returnParamMVAInputs(self)
                inputsEventNr = mvaEvaluatorDL_nonres.returnEventNrMVAInputs(self,t)

               #print ("Lepton variables : %d"%len(inputsLeps))
               #print ("Jet variables    : %d"%len(inputsJets))
               #print ("Fatjet variables : %d"%len(inputsFatjet))
               #print ("MET variables    : %d"%len(inputsMET))
               #print ("HL variables     : %d"%len(inputsHL))
               #print ("Param variables  : %d"%len(inputsParam))
               #print ("Event variables  : %d"%len(inputsEventNr))

               #plots.extend(makeDoubleLeptonMachineLearningInputPlots(selObjectDict['selObject'].sel,selObjectDict['selObject'].selName,selObjectDict['channel'],inputsLeps))
               #plots.extend(makeDoubleLeptonMachineLearningInputPlots(selObjectDict['selObject'].sel,selObjectDict['selObject'].selName,selObjectDict['channel'],inputsJets))
               #plots.extend(makeDoubleLeptonMachineLearningInputPlots(selObjectDict['selObject'].sel,selObjectDict['selObject'].selName,selObjectDict['channel'],inputsFatjet))
               #plots.extend(makeDoubleLeptonMachineLearningInputPlots(selObjectDict['selObject'].sel,selObjectDict['selObject'].selName,selObjectDict['channel'],inputsMET))
               #plots.extend(makeDoubleLeptonMachineLearningInputPlots(selObjectDict['selObject'].sel,selObjectDict['selObject'].selName,selObjectDict['channel'],inputsHL))
               #plots.extend(makeDoubleLeptonMachineLearningInputPlots(selObjectDict['selObject'].sel,selObjectDict['selObject'].selName,selObjectDict['channel'],inputsParam))
               #plots.extend(makeDoubleLeptonMachineLearningInputPlots(selObjectDict['selObject'].sel,selObjectDict['selObject'].selName,selObjectDict['channel'],inputsEventNr))
                from  mvaEvaluatorDL_nonres import inputStaticCast
            
                inputs = [op.array("double",*inputStaticCast(inputsLeps,"float")),
                          op.array("double",*inputStaticCast(inputsJets,"float")),
                          op.array("double",*inputStaticCast(inputsFatjet,"float")),
                          op.array("double",*inputStaticCast(inputsMET,"float")),
                          op.array("double",*inputStaticCast(inputsHL,"float")),
                          op.array("double",*inputStaticCast(inputsParam,"float")),
                          op.array("long",*inputStaticCast(inputsEventNr,"long"))]
    
                output = DNN(*inputs)
                selObjNodesDict = makeDNNOutputNodesSelections(self,selObjectDict['selObject'],output,suffix=model_num)
    
                # Branch out the LO -> NLO reweighting #
                for node in selObjNodesDict.values():
                    node.sel = self.addSignalReweighting(node.sel)
    
                if not self.args.OnlyYield:
                    plots.extend(makeDoubleLeptonMachineLearningExclusiveOutputPlots(selObjNodesDict,output,self.nodes,channel=selObjectDict['channel']))
                    #plots.extend(makeDoubleLeptonMachineLearningInclusiveOutputPlots(selObjNodesDict,output,self.nodes,channel=selObjectDict['channel']))
                if self.args.PrintYield:
                    for selNode in selObjNodesDict.values():
                        self.yields.add(selNode.sel)

            if self.args.analysis == 'res':

                # HME computation #
#                if 'Resolved' in selObjectDict['selObject'].selName:
#                    HME = HME_resolved_per_channel[channel]
#                elif 'Boosted' in selObjectDict['selObject'].selName:
#                    HME = HME_boosted_per_channel[channel]

                # HME reader #
                event_info = (self.tree.run,self.tree.luminosityBlock,self.tree.event)

                # Determine what is the DY condition and associated HME Reader #
                if category == "Resolved1B":
                    hmeInputs = [dilepton[0].p4,                   # l1
                                 dilepton[1].p4,                   # l2
                                 self.ak4JetsByBtagScore[0].p4,    # b1
                                 self.ak4JetsByBtagScore[1].p4,    # b2
                                 self.corrMET.p4,                  # met
                                 op.c_bool(False)]                 # boosted_tag
                    # Only needed to recompute in case missing in TTree (eg when correction is applied and the event is selected while discarded in nominal)
                    if DYCR:
                        DYCond = (op.AND(op.rng_len(self.ak4BJets)==0,op.rng_len(self.ak8BJets)==0,op.rng_len(self.ak8Jets)==0,(self.tree.event//5)%2==0),
                                  hmeReaders['DYCR'][category][selObjectDict['channel']].getHME(*event_info,*hmeInputs))
                elif category == "Resolved2B":
                    hmeInputs = [dilepton[0].p4,                   # l1
                                 dilepton[1].p4,                   # l2
                                 self.ak4JetsByBtagScore[0].p4,    # b1
                                 self.ak4JetsByBtagScore[1].p4,    # b2
                                 self.corrMET.p4,                  # met
                                 op.c_bool(False)]                 # boosted_tag
                    if DYCR:
                        DYCond = (op.AND(op.rng_len(self.ak4BJets)==0,op.rng_len(self.ak8BJets)==0,op.rng_len(self.ak8Jets)==0,(self.tree.event//5)%2==1),
                                  hmeReaders['DYCR'][category][selObjectDict['channel']].getHME(*event_info,*hmeInputs))
                elif category == "Boosted1B":
                    hmeInputs = [dilepton[0].p4,                   # l1
                                 dilepton[1].p4,                   # l2
                                 op.switch(op.rng_len(self.ak8BJets)>0,self.ak8BJets[0].subJet1.p4,self.ak8Jets[0].subJet1.p4),      # b1
                                 op.switch(op.rng_len(self.ak8BJets)>0,self.ak8BJets[0].subJet2.p4,self.ak8Jets[0].subJet1.p4),      # b2
                                    # Because of DY we need the switch
                                 self.corrMET.p4,                  # met
                                 op.c_bool(True)]                  # boosted_tag
                    if DYCR:
                        DYCond = (op.AND(op.rng_len(self.ak8BJets) == 0,op.rng_len(self.ak4BJets) == 0),
                                  hmeReaders['DYCR'][category][selObjectDict['channel']].getHME(*event_info,*hmeInputs))
                else:
                    raise RuntimeError("Not understood category")

                # Determine what is the Fake condition and associated HME Reader #
                if channel == "ElEl" and FakeCR:
                    FakeCond = (op.AND(self.lambda_fakepair_ElEl(self.ElElFakeSel[0]),op.rng_len(self.electronsTightSel) + op.rng_len(self.muonsTightSel)<=2),
                                hmeReaders['FakeCR'][category][channel].getHME(*event_info,*hmeInputs))
                if channel == "MuMu" and FakeCR:
                    FakeCond = (op.AND(self.lambda_fakepair_MuMu(self.MuMuFakeSel[0]),op.rng_len(self.electronsTightSel) + op.rng_len(self.muonsTightSel)<=2),
                                hmeReaders['FakeCR'][category][channel].getHME(*event_info,*hmeInputs))
                if channel == "ElMu" and FakeCR:
                    FakeCond = (op.AND(self.lambda_fakepair_ElMu(self.ElMuFakeSel[0]),op.rng_len(self.electronsTightSel) + op.rng_len(self.muonsTightSel)<=2),
                                hmeReaders['FakeCR'][category][channel].getHME(*event_info,*hmeInputs))

                # Get pre-calculated value of HME #
                switchCond = [hmeReaders['SR'][category][selObjectDict['channel']].getHME(*event_info,*hmeInputs)] # SR condition (else)
                if DYCR:
                    switchCond.insert(0,DYCond)
                if FakeCR:
                    switchCond.insert(0,FakeCond)

                HME = op.defineOnFirstUse(op.multiSwitch(*switchCond))
                #HME = op.multiSwitch(*switchCond)
                #HME = op.forSystematicVariation(HME, "jet", "nominal")
                #forceDefine(HME,selObjectDict['selObject'].sel)

                # Add HME plots #
                if not self.args.OnlyYield:
                    plots.extend(makeDoubleLeptonHMEPlots(selObjectDict['selObject'].sel,selObjectDict['selObject'].selName,selObjectDict['channel'],HME))

                # Inputs of MVA #
                inputsAll = mvaEvaluatorDL_res.returnResonantMVAInputs(
                                                self      = self,
                                                l1        = dilepton[0],
                                                l2        = dilepton[1],
                                                channel   = selObjectDict['channel'],
                                                jets      = self.ak4Jets,
                                                fatjets   = self.ak8Jets,
                                                bjets     = self.ak4JetsByBtagScore[:op.min(op.rng_len(self.ak4JetsByBtagScore),op.static_cast("std::size_t",op.c_int(2)))],
                                                met       = self.corrMET,
                                                electrons = self.electronsTightSel,
                                                muons     = self.muonsTightSel)


                inputs = {inpName:val  for (inpName,_,_),val in inputsAll.items()}
                
                # Mass for resonance in parametric DNN #
                print ('Using parametric DNN with')
                if self.args.mass is not None:
                    masses = self.args.mass
                else:
                    print ('No mass requested, will run each sample with its mass')
                    masses = [sampleCfg['mass']]

                self.nodes = ['DY','GGF','H','Rare','ST','TT','TTVX','VVV']
                for mass in masses:
                    # Select correct DNN #
                    print ('... MH = {}'.format(mass))
                    if mass <= 500:
                        DNN = DNN_LowMass
                        input_names = input_names_LowMass
                    else:
                        DNN = DNN_HighMass
                        input_names = input_names_LowMass

                    # Define the inputs #
                    inputs['param'] = op.c_float(mass) 
                    inputsArr = []
                    for inpName in input_names:
                        if inpName not in inputs.keys():
                            for key in inputs.keys():
                                if inpName == key.replace('$','').replace(' ','').replace('_',''):
                                    inpName = key
                        
                        if inpName not in inputs.keys():
                            raise RuntimeError(f"Input node {inpName} not found in the inputs in bamboo")
                        inpVal = inputs[inpName]

                        if inpName == "eventnr":
                            inpType = "long"
                        else:
                            inpType = "float"
                        if isinstance(inpVal,list):
                            inpVal = [op.static_cast(inpType,inp) for inp in inpVal] 
                            inputsArr.append(op.array(inpType,*inpVal))
                        else:
                            inpVal = op.static_cast(inpType,inpVal)
                            inputsArr.append(op.array(inpType,inpVal))

                    output = DNN(*inputsArr)
                    selObjNodesDict = makeDNNOutputNodesSelections(self,selObjectDict['selObject'],output,suffix="M{}".format(int(mass)))
                    if not self.args.OnlyYield:
                        plots.extend(makeDoubleLeptonMachineLearningExclusiveOutputPlots(selObjNodesDict,output,self.nodes,channel=selObjectDict['channel']))
                        plots.extend(makeDoubleLeptonMachineLearningExclusiveOutputPlotsWithHME(selObjNodesDict,output,['GGF'],channel=selObjectDict['channel'],HME=HME))
                    if self.args.PrintYield:
                        for selNode in selObjNodesDict.values():
                            self.yields.add(selNode.sel)
                            
                            

                            
    
        #----- Add the Yield plots -----#
        if self.args.PrintYield or self.args.OnlyYield:
            plots.append(self.yields)

            
        #----- Return -----#
        return plots

    ### PostProcess ###
    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None):
        super(PlotterNanoHHtobbWWDL, self).postProcess(taskList, config, workdir, resultsdir, forSkimmer=False)
