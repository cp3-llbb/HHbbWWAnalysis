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

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)))) # Add scripts in this directory
from BaseHHtobbWW import BaseNanoHHtobbWW
from plotDef import *
from selectionDef import *
from DDHelper import DataDrivenFake, DataDrivenDY
from mvaEvaluatorDL import *

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

        era = sampleCfg['era']


        #----- Ratio reweighting variables (before lepton and jet selection) -----#
        if self.args.BtagReweightingOff or self.args.BtagReweightingOn:
            noSel = self.beforeJetselection(noSel)
            plots.append(objectsNumberPlot(channel="NoChannel",suffix='NoSelection',sel=noSel,objCont=self.ak4Jets,objName='Ak4Jets',Nmax=15,xTitle='N(Ak4 jets)'))
            plots.append(CutFlowReport("BtagReweightingCutFlowReport",noSel,printInLog=True,recursive=True))
            return plots

        #----- Stitching study -----#
        if self.args.DYStitchingPlots or self.args.WJetsStitchingPlots:
            if self.args.DYStitchingPlots and sampleCfg['group'] != 'DY':
                raise RuntimeError("Stitching is only done on DY MC samples")
            if self.args.WJetsStitchingPlots and sampleCfg['group'] != 'Wjets':
                raise RuntimeError("Stitching is only done on WJets MC samples")
            plots.extend(makeLHEPlots(noSel,t.LHE))
            plots.append(objectsNumberPlot(channel="NoChannel",suffix='NoSelection',sel=noSel,objCont=self.ak4Jets,objName='Ak4Jets',Nmax=15,xTitle='N(Ak4 jets)'))
            plots.append(CutFlowReport("DYStitchingCutFlowReport",noSel,printInLog=True,recursive=True))
            return plots

        #----- Machine Learning Model -----#                
        model_num = "09"
        if not self.args.OnlyYield:
            path_model = os.path.join(os.path.abspath(os.path.dirname(__file__)),'MachineLearning','ml-models','models','multi-classification','dnn',model_num,'model','model.pb')
            print ("DNN model : %s"%path_model)
            if not os.path.exists(path_model):
                raise RuntimeError('Could not find model file %s'%path_model)
            try:
                input_names = ["lep","jet","fat","met","hl","param","eventnr"]
                output_name = "Identity"
                DNN = op.mvaEvaluator(path_model,mvaType='Tensorflow',otherArgs=(input_names, output_name))
            except:
                raise RuntimeError('Could not load model %s'%path_model)

        #----- Dileptons -----#
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
        if self.args.POGID:
            OSElElDilepton = self.ElElTightSel
            OSMuMuDilepton = self.MuMuTightSel
            OSElMuDilepton = self.ElMuTightSel
        if self.args.TTHIDLoose or self.args.TTHIDTight:
            OSElElDilepton = self.ElElFakeSel
            OSMuMuDilepton = self.MuMuFakeSel
            OSElMuDilepton = self.ElMuFakeSel

        #----- DY reweighting -----#
        if 'DYEstimation' in self.datadrivenContributions:
            mode = 'mc' if "PseudoData" in self.datadrivenContributions else 'data'
            # Resolved : ElEl #
            self.ResolvedDYReweighting1bElEl = self.SF.get_scalefactor("lepton", ('DY_resolved_{}'.format(era),'ElEl_leadjetPt_{}_1b'.format(mode)), combine="weight", 
                                                                   systName="dy_ee_resolved_reweighting_1b", 
                                                                   additionalVariables={'Eta': lambda x : op.c_float(0.),'Pt': lambda x: self.HLL.getCorrBp4(x).Pt()})
            self.ResolvedDYReweighting2bElEl = self.SF.get_scalefactor("lepton", ('DY_resolved_{}'.format(era),'ElEl_leadjetPt_{}_2b'.format(mode)), combine="weight", 
                                                                   systName="dy_ee_resolved_reweighting_2b", 
                                                                   additionalVariables={'Eta': lambda x : op.c_float(0.),'Pt': lambda x: self.HLL.getCorrBp4(x).Pt()})
            # Resolved : MuMu #
            self.ResolvedDYReweighting1bMuMu = self.SF.get_scalefactor("lepton", ('DY_resolved_{}'.format(era),'MuMu_leadjetPt_{}_1b'.format(mode)), combine="weight", 
                                                                   systName="dy_mm_resolved_reweighting_1b", 
                                                                   additionalVariables={'Eta': lambda x : op.c_float(0.),'Pt': lambda x: self.HLL.getCorrBp4(x).Pt()})
            self.ResolvedDYReweighting2bMuMu = self.SF.get_scalefactor("lepton", ('DY_resolved_{}'.format(era),'MuMu_leadjetPt_{}_2b'.format(mode)), combine="weight", 
                                                                   systName="dy_mm_resolved_reweighting_1b", 
                                                                   additionalVariables={'Eta': lambda x : op.c_float(0.),'Pt': lambda x: self.HLL.getCorrBp4(x).Pt()})
            # Boosted : ElEl #
            self.BoostedDYReweighting1bElEl  =  self.SF.get_scalefactor("lepton", ('DY_boosted_{}'.format(era),'ElEl_fatjetsoftDropmass_{}_1b'.format(mode)), combine="weight", 
                                                                   systName="dy_ee_boosted_reweighting_1b", 
                                                                   additionalVariables={'Eta': lambda x : op.c_float(0.),'Pt': lambda x: x.msoftdrop})
            # Boosted : MuMu #
            self.BoostedDYReweighting1bMuMu  =  self.SF.get_scalefactor("lepton", ('DY_boosted_{}'.format(era),'MuMu_fatjetsoftDropmass_{}_1b'.format(mode)), combine="weight", 
                                                                   systName="dy_mm_boosted_reweighting_1b", 
                                                                   additionalVariables={'Eta': lambda x : op.c_float(0.),'Pt': lambda x: x.msoftdrop})
        else:
            self.ResolvedDYReweighting1bElEl = lambda dilep : None
            self.ResolvedDYReweighting2bElEl = lambda dilep : None
            self.ResolvedDYReweighting1bMuMu = lambda dilep : None
            self.ResolvedDYReweighting2bMuMu = lambda dilep : None
            self.BoostedDYReweighting1bElEl  = lambda dilep : None
            self.BoostedDYReweighting1bMuMu  = lambda dilep : None

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
                ChannelDictList.append({'channel': 'ElEl','met': self.corrMET,'l1':OSElElDilepton[0][0],'l2':OSElElDilepton[0][1],'b1':container0b2j[0],'b2':container0b2j[1],'sel':ElElSelObjAk4JetsExclusiveResolvedNoBtag.sel,'suffix':ElElSelObjAk4JetsExclusiveResolvedNoBtag.selName})
                ChannelDictList.append({'channel': 'MuMu','met': self.corrMET,'l1':OSMuMuDilepton[0][0],'l2':OSMuMuDilepton[0][1],'b1':container0b2j[0],'b2':container0b2j[1],'sel':MuMuSelObjAk4JetsExclusiveResolvedNoBtag.sel,'suffix':MuMuSelObjAk4JetsExclusiveResolvedNoBtag.selName})
                ChannelDictList.append({'channel': 'ElMu','met': self.corrMET,'l1':OSElMuDilepton[0][0],'l2':OSElMuDilepton[0][1],'b1':container0b2j[0],'b2':container0b2j[1],'sel':ElMuSelObjAk4JetsExclusiveResolvedNoBtag.sel,'suffix':ElMuSelObjAk4JetsExclusiveResolvedNoBtag.selName})
            if "Resolved1Btag" in jetplot_level:
                ChannelDictList.append({'channel': 'ElEl','met': self.corrMET,'l1':OSElElDilepton[0][0],'l2':OSElElDilepton[0][1],'b1':container1b1j[0],'b2':container1b1j[1],'sel':ElElSelObjAk4JetsExclusiveResolvedOneBtag.sel,'suffix':ElElSelObjAk4JetsExclusiveResolvedOneBtag.selName})
                ChannelDictList.append({'channel': 'MuMu','met': self.corrMET,'l1':OSMuMuDilepton[0][0],'l2':OSMuMuDilepton[0][1],'b1':container1b1j[0],'b2':container1b1j[1],'sel':MuMuSelObjAk4JetsExclusiveResolvedOneBtag.sel,'suffix':MuMuSelObjAk4JetsExclusiveResolvedOneBtag.selName})
                ChannelDictList.append({'channel': 'ElMu','met': self.corrMET,'l1':OSElMuDilepton[0][0],'l2':OSElMuDilepton[0][1],'b1':container1b1j[0],'b2':container1b1j[1],'sel':ElMuSelObjAk4JetsExclusiveResolvedOneBtag.sel,'suffix':ElMuSelObjAk4JetsExclusiveResolvedOneBtag.selName})
            # Resolved Two Btags  #
            if "Resolved2Btag" in jetplot_level:
                ChannelDictList.append({'channel': 'ElEl','met': self.corrMET,'l1':OSElElDilepton[0][0],'l2':OSElElDilepton[0][1],'b1':container2b0j[0],'b2':container2b0j[1],'sel':ElElSelObjAk4JetsExclusiveResolvedTwoBtags.sel,'suffix':ElElSelObjAk4JetsExclusiveResolvedTwoBtags.selName})
                ChannelDictList.append({'channel': 'MuMu','met': self.corrMET,'l1':OSMuMuDilepton[0][0],'l2':OSMuMuDilepton[0][1],'b1':container2b0j[0],'b2':container2b0j[1],'sel':MuMuSelObjAk4JetsExclusiveResolvedTwoBtags.sel,'suffix':MuMuSelObjAk4JetsExclusiveResolvedTwoBtags.selName})
                ChannelDictList.append({'channel': 'ElMu','met': self.corrMET,'l1':OSElMuDilepton[0][0],'l2':OSElMuDilepton[0][1],'b1':container2b0j[0],'b2':container2b0j[1],'sel':ElMuSelObjAk4JetsExclusiveResolvedTwoBtags.sel,'suffix':ElMuSelObjAk4JetsExclusiveResolvedTwoBtags.selName})

        for channelDict in ChannelDictList:
            plots.extend(makeDoubleLeptonSelectedResolvedVariables(**channelDict,HLL=self.HLL))

        #----- Selected variables : Boosted -----#
        ChannelDictList = []
        if not self.args.OnlyYield:
            # Boosted No Btag #
            if "Boosted0Btag" in jetplot_level:
                ChannelDictList.append({'channel': 'ElEl','met': self.corrMET,'l1':OSElElDilepton[0][0],'l2':OSElElDilepton[0][1],'B':container1fatb[0],'sel':ElElSelObjAk8JetsInclusiveBoostedNoBtag.sel,'suffix':ElElSelObjAk8JetsInclusiveBoostedNoBtag.selName})
                ChannelDictList.append({'channel': 'MuMu','met': self.corrMET,'l1':OSMuMuDilepton[0][0],'l2':OSMuMuDilepton[0][1],'B':container1fatb[0],'sel':MuMuSelObjAk8JetsInclusiveBoostedNoBtag.sel,'suffix':MuMuSelObjAk8JetsInclusiveBoostedNoBtag.selName})
                ChannelDictList.append({'channel': 'ElMu','met': self.corrMET,'l1':OSElMuDilepton[0][0],'l2':OSElMuDilepton[0][1],'B':container1fatb[0],'sel':ElMuSelObjAk8JetsInclusiveBoostedNoBtag.sel,'suffix':ElMuSelObjAk8JetsInclusiveBoostedNoBtag.selName})
            # Boosted One Btag #
            if "Boosted1Btag" in jetplot_level:
                ChannelDictList.append({'channel': 'ElEl','met': self.corrMET,'l1':OSElElDilepton[0][0],'l2':OSElElDilepton[0][1],'B':container1fatb[0],'sel':ElElSelObjAk8JetsInclusiveBoostedOneBtag.sel,'suffix':ElElSelObjAk8JetsInclusiveBoostedOneBtag.selName})
                ChannelDictList.append({'channel': 'MuMu','met': self.corrMET,'l1':OSMuMuDilepton[0][0],'l2':OSMuMuDilepton[0][1],'B':container1fatb[0],'sel':MuMuSelObjAk8JetsInclusiveBoostedOneBtag.sel,'suffix':MuMuSelObjAk8JetsInclusiveBoostedOneBtag.selName})
                ChannelDictList.append({'channel': 'ElMu','met': self.corrMET,'l1':OSElMuDilepton[0][0],'l2':OSElMuDilepton[0][1],'B':container1fatb[0],'sel':ElMuSelObjAk8JetsInclusiveBoostedOneBtag.sel,'suffix':ElMuSelObjAk8JetsInclusiveBoostedOneBtag.selName})

        for channelDict in ChannelDictList:
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
        self.nodes = ['GGF','VBF','H', 'DY', 'ST', 'TT', 'TTVX', 'VVV', 'Rare']
        
        for selObjectDict in selObjectDictList:
            dilepton = dileptonsCont[selObjectDict['channel']]

            #self.yields.add(selObjectDict['selObject'].sel,title=selObjectDict['selObject'].yieldTitle)

            inputsLeps = returnLeptonsMVAInputs(self      = self,
                                                l1        = dilepton[0],
                                                l2        = dilepton[1],
                                                channel   = selObjectDict['channel'])
            inputsJets =    returnJetsMVAInputs(self      = self,
                                                jets      = self.ak4Jets)
            inputsMET =      returnMETMVAInputs(self      = self,
                                                met       = self.corrMET)     
            inputsFatjet =  returnFatjetMVAInputs(self      = self,
                                                  fatjets   = self.ak8Jets)
            inputsHL = returnHighLevelMVAInputs(self      = self,
                                                l1        = dilepton[0],
                                                l2        = dilepton[1],
                                                met       = self.corrMET,
                                                jets      = self.ak4Jets,
                                                bjets     = self.ak4JetsByBtagScore[:op.min(op.rng_len(self.ak4JetsByBtagScore),op.static_cast("std::size_t",op.c_int(2)))],
                                                electrons = self.electronsTightSel,
                                                muons     = self.muonsTightSel,
                                                channel   = selObjectDict['channel'])

            inputsParam = returnParamMVAInputs(self)
            inputsEventNr = returnEventNrMVAInputs(self,t)

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
            
            inputs = [op.array("double",*inputStaticCast(inputsLeps,"float")),
                      op.array("double",*inputStaticCast(inputsJets,"float")),
                      op.array("double",*inputStaticCast(inputsFatjet,"float")),
                      op.array("double",*inputStaticCast(inputsMET,"float")),
                      op.array("double",*inputStaticCast(inputsHL,"float")),
                      op.array("double",*inputStaticCast(inputsParam,"float")),
                      op.array("long",*inputStaticCast(inputsEventNr,"long"))]

            output = DNN(*inputs)
            selObjNodesDict = makeDNNOutputNodesSelections(self,selObjectDict['selObject'],output,suffix=model_num)

            plots.extend(makeDoubleLeptonMachineLearningOutputPlots(selObjNodesDict,output,self.nodes,channel=selObjectDict['channel']))
            
    
        #----- Add the Yield plots -----#

        #plots.extend(cutFlowPlots)
        plots.append(self.yields)
            
        #----- Return -----#
        return plots
