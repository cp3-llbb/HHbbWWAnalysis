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
from mvaEvaluatorDL import returnLowLevelMVAInputs, returnHighLevelMVAInputs

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
        super(PlotterNanoHHtobbWWDL, self).initialize()
        # Change the way the FakeExtrapolation is postProcessed (avoids overriding the `postProcess` method) 
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

        cutFlowPlots = []

        era = sampleCfg['era']

        self.yieldPlots = makeYieldPlots(self.args.Synchronization)

        #----- Ratio reweighting variables (before lepton and jet selection) -----#
        if self.args.BtagReweightingOff or self.args.BtagReweightingOn:
            plots.append(objectsNumberPlot(channel="NoChannel",suffix='NoSelection',sel=noSel,objCont=self.ak4Jets,objName='Ak4Jets',Nmax=15,xTitle='N(Ak4 jets)'))
            plots.append(CutFlowReport("BtagReweightingCutFlowReport",noSel))
            return plots

        #----- Stitching study -----#
        if self.args.DYStitchingPlots or self.args.WJetsStitchingPlots:
            if self.args.DYStitchingPlots and sampleCfg['group'] != 'DY':
                raise RuntimeError("Stitching is only done on DY MC samples")
            if self.args.WJetsStitchingPlots and sampleCfg['group'] != 'Wjets':
                raise RuntimeError("Stitching is only done on WJets MC samples")
            plots.extend(makeLHEPlots(noSel,t.LHE))
            plots.append(objectsNumberPlot(channel="NoChannel",suffix='NoSelection',sel=noSel,objCont=self.ak4Jets,objName='Ak4Jets',Nmax=15,xTitle='N(Ak4 jets)'))
            plots.append(CutFlowReport("DYStitchingCutFlowReport",noSel))
            return plots

        #----- Machine Learning Model -----#                
        #path_model = os.path.join(os.path.abspath(os.path.dirname(__file__)),'MachineLearning','ml-models','models','multi-classification','dnn',self.args.Classifier,'model','model.pb')
        DNNs = {}
        for model_num in ["01","02","03","04"]:
        #for model_num in ["03"]:
            path_model = os.path.join(os.path.abspath(os.path.dirname(__file__)),'MachineLearning','ml-models','models','multi-classification','dnn',model_num,'model','model.pb')
            print ("DNN model : %s"%path_model)
            if not os.path.exists(path_model):
                raise RuntimeError('Could not find model file %s'%path_model)
            try:
                input_names = ["input_1", "input_2"]
                if model_num in ["03","04"]:
                    input_names.append("input_3")
                DNNs[model_num] = op.mvaEvaluator(path_model,mvaType='Tensorflow',otherArgs=(input_names, "model/output/Softmax"))
            except:
                raise RuntimeError('Could not load model %s'%path_model)


        #----- Dileptons -----#
        selObjectDict = makeDoubleLeptonSelection(self,noSel,plot_yield=True)
        # selObjectDict : keys -> level (str)
        #                 values -> [ElEl,MuMu,ElMu] x Selection object

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

        # Loop over lepton selection and start plotting #
        for selectionType, selectionList in selObjectDict.items():
            print ("... Processing %s lepton type"%selectionType)
            #----- Select correct dilepton -----#
            if selectionType == "Preselected":  
                OSElElDilepton = self.OSElElDileptonPreSel
                OSMuMuDilepton = self.OSMuMuDileptonPreSel
                OSElMuDilepton = self.OSElMuDileptonPreSel
            elif selectionType == "Fakeable":
                OSElElDilepton = self.ElElDileptonFakeSel
                OSMuMuDilepton = self.MuMuDileptonFakeSel
                OSElMuDilepton = self.ElMuDileptonFakeSel
            elif selectionType == "Tight":
                if self.args.POGID:
                    OSElElDilepton = self.ElElDileptonTightSel
                    OSMuMuDilepton = self.MuMuDileptonTightSel
                    OSElMuDilepton = self.ElMuDileptonTightSel
                if self.args.TTHIDLoose or self.args.TTHIDTight:
                    OSElElDilepton = switch_on_index([0],op.rng_len(self.ElElDileptonTightSel)>=1,self.ElElDileptonTightSel,self.ElElDileptonFakeExtrapolationSel)
                    OSMuMuDilepton = switch_on_index([0],op.rng_len(self.MuMuDileptonTightSel)>=1,self.MuMuDileptonTightSel,self.MuMuDileptonFakeExtrapolationSel)
                    OSElMuDilepton = switch_on_index([0],op.rng_len(self.ElMuDileptonTightSel)>=1,self.ElMuDileptonTightSel,self.ElMuDileptonFakeExtrapolationSel)
            elif selectionType == "FakeExtrapolation":
                OSElElDilepton = self.ElElDileptonFakeExtrapolationSel
                OSMuMuDilepton = self.MuMuDileptonFakeExtrapolationSel
                OSElMuDilepton = self.ElMuDileptonFakeExtrapolationSel


            #----- DY reweighting -----#
            if 'DYEstimation' in self.datadrivenContributions:
                mode = 'mc' if "PseudoData" in self.datadrivenContributions else 'data'
                # Resolved : ElEl #
                self.ResolvedDYReweighting1bElEl = self.SF.get_scalefactor("lepton", ('DY_resolved_{}'.format(era),'ElEl_leadjetPt_{}_1b'.format(mode)), combine="weight", 
                                                                       systName="dy_ee_resolved_reweighting_1b", 
                                                                       additionalVariables={'Eta': lambda x : op.c_float(0.),'Pt': lambda x: x.pt})
                self.ResolvedDYReweighting2bElEl = self.SF.get_scalefactor("lepton", ('DY_resolved_{}'.format(era),'ElEl_leadjetPt_{}_2b'.format(mode)), combine="weight", 
                                                                       systName="dy_ee_resolved_reweighting_2b", 
                                                                       additionalVariables={'Eta': lambda x : op.c_float(0.),'Pt': lambda x: x.pt})
                # Resolved : MuMu #
                self.ResolvedDYReweighting1bMuMu = self.SF.get_scalefactor("lepton", ('DY_resolved_{}'.format(era),'MuMu_leadjetPt_{}_1b'.format(mode)), combine="weight", 
                                                                       systName="dy_mm_resolved_reweighting_1b", 
                                                                       additionalVariables={'Eta': lambda x : op.c_float(0.),'Pt': lambda x: x.pt})
                self.ResolvedDYReweighting2bMuMu = self.SF.get_scalefactor("lepton", ('DY_resolved_{}'.format(era),'MuMu_leadjetPt_{}_2b'.format(mode)), combine="weight", 
                                                                       systName="dy_mm_resolved_reweighting_1b", 
                                                                       additionalVariables={'Eta': lambda x : op.c_float(0.),'Pt': lambda x: x.pt})
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


            #----- Separate selections ------#
            ElElSelObjectDilepton = selectionList[0]
            MuMuSelObjectDilepton = selectionList[1]
            ElMuSelObjectDilepton = selectionList[2]

            #----- Channel and trigger plots -----#

            if not self.args.OnlyYield:
                ChannelDictList = []
                ChannelDictList.append({'channel':'ElEl','sel':ElElSelObjectDilepton.sel,'suffix':ElElSelObjectDilepton.selName})
                ChannelDictList.append({'channel':'MuMu','sel':MuMuSelObjectDilepton.sel,'suffix':MuMuSelObjectDilepton.selName})
                ChannelDictList.append({'channel':'ElMu','sel':ElMuSelObjectDilepton.sel,'suffix':ElMuSelObjectDilepton.selName})

                for channelDict in ChannelDictList:
                    #----- Trigger plots -----#
                    plots.extend(doubleLeptonTriggerPlots(**channelDict,triggerDict=self.triggersPerPrimaryDataset))
                    #----- Dilepton plots -----#
                    #plots.extend(doubleLeptonChannelPlot(**channelDict,DilepElEl=OSElElDilepton,DilepMuMu=OSMuMuDilepton,DilepElMu=OSElMuDilepton))

            LeptonKeys = ['channel','sel','dilepton','suffix','is_MC']
            JetKeys    = ['channel','sel','leadjet','subleadjet','lead_is_b','sublead_is_b','suffix','is_MC']
            commonItems = ['channel','sel','suffix']
                
            #----- Ak4 jets selection -----#
            if "Ak4" in jetsel_level:
                print ("...... Processing Ak4 jet selection")
                ElElSelObjectDileptonAk4Jets = makeAtLeastTwoAk4JetSelection(self,ElElSelObjectDilepton,copy_sel=True,plot_yield=True)
                MuMuSelObjectDileptonAk4Jets = makeAtLeastTwoAk4JetSelection(self,MuMuSelObjectDilepton,copy_sel=True,plot_yield=True)
                ElMuSelObjectDileptonAk4Jets = makeAtLeastTwoAk4JetSelection(self,ElMuSelObjectDilepton,copy_sel=True,plot_yield=True)

                # Jet and lepton plots #
                ChannelDictList = []
                if "Ak4" in jetplot_level:
                    # Cut flow report #
                    cutFlowPlots.append(CutFlowReport(ElElSelObjectDileptonAk4Jets.selName,ElElSelObjectDileptonAk4Jets.sel))
                    cutFlowPlots.append(CutFlowReport(MuMuSelObjectDileptonAk4Jets.selName,MuMuSelObjectDileptonAk4Jets.sel))
                    cutFlowPlots.append(CutFlowReport(ElMuSelObjectDileptonAk4Jets.selName,ElMuSelObjectDileptonAk4Jets.sel))
                    if not self.args.OnlyYield:
                        ChannelDictList.append({'channel':'ElEl','sel':ElElSelObjectDileptonAk4Jets.sel,'dilepton':OSElElDilepton[0],'leadjet':self.ak4Jets[0],'subleadjet':self.ak4Jets[1],'lead_is_b':False,'sublead_is_b':False,'suffix':ElElSelObjectDileptonAk4Jets.selName,'is_MC':self.is_MC})
                        ChannelDictList.append({'channel':'MuMu','sel':MuMuSelObjectDileptonAk4Jets.sel,'dilepton':OSMuMuDilepton[0],'leadjet':self.ak4Jets[0],'subleadjet':self.ak4Jets[1],'lead_is_b':False,'sublead_is_b':False,'suffix':MuMuSelObjectDileptonAk4Jets.selName,'is_MC':self.is_MC})
                        ChannelDictList.append({'channel':'ElMu','sel':ElMuSelObjectDileptonAk4Jets.sel,'dilepton':OSElMuDilepton[0],'leadjet':self.ak4Jets[0],'subleadjet':self.ak4Jets[1],'lead_is_b':False,'sublead_is_b':False,'suffix':ElMuSelObjectDileptonAk4Jets.selName,'is_MC':self.is_MC})

                JetsN = {'objName':'Ak4Jets','objCont':self.ak4Jets,'Nmax':15,'xTitle':'N(Ak4 jets)'}
                                                            
                for channelDict in ChannelDictList:
                    # Dilepton #
                    plots.extend(makeDileptonPlots(**{k:channelDict[k] for k in LeptonKeys}))
                    # Number of jets #
                    plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**JetsN))
                    # Ak4 Jets #
                    plots.extend(makeTwoAk4JetsPlots(**{k:channelDict[k] for k in JetKeys}))
                    # MET #
                    plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))

            ##### Ak8 jets selection #####
            if "Ak8" in jetsel_level:
                print ("...... Processing Ak8 jet selection")
                ElElSelObjectDileptonAk8Jets = makeAtLeastOneAk8JetSelection(self,ElElSelObjectDilepton,copy_sel=True,plot_yield=True)
                MuMuSelObjectDileptonAk8Jets = makeAtLeastOneAk8JetSelection(self,MuMuSelObjectDilepton,copy_sel=True,plot_yield=True)
                ElMuSelObjectDileptonAk8Jets = makeAtLeastOneAk8JetSelection(self,ElMuSelObjectDilepton,copy_sel=True,plot_yield=True)

                # Fatjets plots #
                ChannelDictList = []
                if "Ak8" in jetplot_level:
                    # Cut flow report #
                    cutFlowPlots.append(CutFlowReport(ElElSelObjectDileptonAk8Jets.selName,ElElSelObjectDileptonAk8Jets.sel))
                    cutFlowPlots.append(CutFlowReport(MuMuSelObjectDileptonAk8Jets.selName,MuMuSelObjectDileptonAk8Jets.sel))
                    cutFlowPlots.append(CutFlowReport(ElMuSelObjectDileptonAk8Jets.selName,ElMuSelObjectDileptonAk8Jets.sel))
                    if not self.args.OnlyYield:
                        ChannelDictList.append({'channel':'ElEl','sel':ElElSelObjectDileptonAk8Jets.sel,'dilepton':OSElElDilepton[0],'fatjet':self.ak8Jets[0],'suffix':ElElSelObjectDileptonAk8Jets.selName,'is_MC':self.is_MC})
                        ChannelDictList.append({'channel':'MuMu','sel':MuMuSelObjectDileptonAk8Jets.sel,'dilepton':OSMuMuDilepton[0],'fatjet':self.ak8Jets[0],'suffix':MuMuSelObjectDileptonAk8Jets.selName,'is_MC':self.is_MC})
                        ChannelDictList.append({'channel':'ElMu','sel':ElMuSelObjectDileptonAk8Jets.sel,'dilepton':OSElMuDilepton[0],'fatjet':self.ak8Jets[0],'suffix':ElMuSelObjectDileptonAk8Jets.selName,'is_MC':self.is_MC})

                FatJetKeys = ['channel','sel','fatjet','suffix']
                FatJetsN = {'objName':'Ak8Jets','objCont':self.ak8Jets,'Nmax':5,'xTitle':'N(Ak8 jets)'}
                
                for channelDict in ChannelDictList:
                    # Dilepton #
                    plots.extend(makeDileptonPlots(**{k:channelDict[k] for k in LeptonKeys}))
                    # Number of jets #
                    plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**FatJetsN))
                    # Ak8 Jets #
                    plots.extend(makeDoubleLeptonAk8JetsPlots(**{k:channelDict[k] for k in FatJetKeys}))
                    # MET #
                    plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))

                         
            #----- Resolved selection -----#
            if any(item in ["Resolved0Btag","Resolved1Btag","Resolved2Btag"] for item in jetsel_level): # If any of resolved category is asked
                #----- select the jets -----#
                aka4JetsByBtagScore = op.sort(self.ak4Jets, lambda jet : -jet.btagDeepFlavB)

                container0b2j = [ t.Jet[self.ak4LightJetsByBtagScore[i].idx] for i in range(2) ]
                container1b1j = [ t.Jet[op.switch(op.rng_len(self.ak4BJets) == 1, self.ak4BJets[0].idx,  aka4JetsByBtagScore[0].idx)] ,
                                  t.Jet[op.switch(op.rng_len(self.ak4BJets) == 1, self.ak4LightJetsByBtagScore[0].idx,  aka4JetsByBtagScore[1].idx)]]
                container2b0j = [ t.Jet[op.switch(op.rng_len(self.ak4BJets) >= 2, self.ak4BJets[i].idx, aka4JetsByBtagScore[i].idx)] for i in range(2) ]

                ChannelDictList = []
                #----- Resolved selection : 0 Btag -----#
                if "Resolved0Btag" in jetsel_level:
                    print ("...... Processing Resolved jet (0 btag) selection")
                    ElElSelObjectDileptonAk4JetsExclusiveResolvedNoBtag = makeExclusiveResolvedNoBtagSelection(self,ElElSelObjectDileptonAk4Jets,copy_sel=True,plot_yield=True)
                    MuMuSelObjectDileptonAk4JetsExclusiveResolvedNoBtag = makeExclusiveResolvedNoBtagSelection(self,MuMuSelObjectDileptonAk4Jets,copy_sel=True,plot_yield=True)
                    ElMuSelObjectDileptonAk4JetsExclusiveResolvedNoBtag = makeExclusiveResolvedNoBtagSelection(self,ElMuSelObjectDileptonAk4Jets,copy_sel=True,plot_yield=True)

                    if "Resolved0Btag" in jetplot_level:
                        # Cut flow report #
                        cutFlowPlots.append(CutFlowReport(ElElSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.selName,ElElSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.sel))
                        cutFlowPlots.append(CutFlowReport(MuMuSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.selName,MuMuSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.sel))
                        cutFlowPlots.append(CutFlowReport(ElMuSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.selName,ElMuSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.sel))
                        if not self.args.OnlyYield:
                            ChannelDictList.append({'channel':'ElEl','sel':ElElSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.sel,'dilepton':OSElElDilepton[0],'leadjet':container0b2j[0],'subleadjet':container0b2j[1],'lead_is_b':False,'sublead_is_b':False,'suffix':ElElSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.selName,'is_MC':self.is_MC})
                            ChannelDictList.append({'channel':'MuMu','sel':MuMuSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.sel,'dilepton':OSMuMuDilepton[0],'leadjet':container0b2j[0],'subleadjet':container0b2j[1],'lead_is_b':False,'sublead_is_b':False,'suffix':MuMuSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.selName,'is_MC':self.is_MC})
                            ChannelDictList.append({'channel':'ElMu','sel':ElMuSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.sel,'dilepton':OSElMuDilepton[0],'leadjet':container0b2j[0],'subleadjet':container0b2j[1],'lead_is_b':False,'sublead_is_b':False,'suffix':ElMuSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.selName,'is_MC':self.is_MC})
                #----  Resolved selection : 1 Btag  -----#
                if "Resolved1Btag" in jetsel_level:
                    print ("...... Processing Resolved jet (1 btag) selection")
                    ElElSelObjectDileptonAk4JetsExclusiveResolvedOneBtag = makeExclusiveResolvedOneBtagSelection(self,ElElSelObjectDileptonAk4Jets,copy_sel=True,plot_yield=True)
                    MuMuSelObjectDileptonAk4JetsExclusiveResolvedOneBtag = makeExclusiveResolvedOneBtagSelection(self,MuMuSelObjectDileptonAk4Jets,copy_sel=True,plot_yield=True)
                    ElMuSelObjectDileptonAk4JetsExclusiveResolvedOneBtag = makeExclusiveResolvedOneBtagSelection(self,ElMuSelObjectDileptonAk4Jets,copy_sel=True,plot_yield=True)

                    if "Resolved1Btag" in jetplot_level:
                        # Cut flow report #
                        cutFlowPlots.append(CutFlowReport(ElElSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.selName,ElElSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.sel))
                        cutFlowPlots.append(CutFlowReport(MuMuSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.selName,MuMuSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.sel))
                        cutFlowPlots.append(CutFlowReport(ElMuSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.selName,ElMuSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.sel))
                        if not self.args.OnlyYield:
                            ChannelDictList.append({'channel':'ElEl','sel':ElElSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.sel,'dilepton':OSElElDilepton[0],'leadjet':container1b1j[0],'subleadjet':container1b1j[1],'lead_is_b':True,'sublead_is_b':False,'suffix':ElElSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.selName,'is_MC':self.is_MC})
                            ChannelDictList.append({'channel':'MuMu','sel':MuMuSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.sel,'dilepton':OSMuMuDilepton[0],'leadjet':container1b1j[0],'subleadjet':container1b1j[1],'lead_is_b':True,'sublead_is_b':False,'suffix':MuMuSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.selName,'is_MC':self.is_MC})
                            ChannelDictList.append({'channel':'ElMu','sel':ElMuSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.sel,'dilepton':OSElMuDilepton[0],'leadjet':container1b1j[0],'subleadjet':container1b1j[1],'lead_is_b':True,'sublead_is_b':False,'suffix':ElMuSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.selName,'is_MC':self.is_MC})
                #----- Resolved selection : 2 Btags -----#
                if "Resolved2Btag" in jetsel_level:
                    print ("...... Processing Resolved jet (2 btags) selection")
                    ElElSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags = makeExclusiveResolvedTwoBtagsSelection(self,ElElSelObjectDileptonAk4Jets,copy_sel=True,plot_yield=True)
                    MuMuSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags = makeExclusiveResolvedTwoBtagsSelection(self,MuMuSelObjectDileptonAk4Jets,copy_sel=True,plot_yield=True)
                    ElMuSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags = makeExclusiveResolvedTwoBtagsSelection(self,ElMuSelObjectDileptonAk4Jets,copy_sel=True,plot_yield=True)

                    if "Resolved2Btag" in jetplot_level:
                        # Cut flow report #
                        cutFlowPlots.append(CutFlowReport(ElElSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.selName,ElElSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.sel))
                        cutFlowPlots.append(CutFlowReport(MuMuSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.selName,MuMuSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.sel))
                        cutFlowPlots.append(CutFlowReport(ElMuSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.selName,ElMuSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.sel))
                        if not self.args.OnlyYield:
                            ChannelDictList.append({'channel':'ElEl','sel':ElElSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.sel,'dilepton':OSElElDilepton[0],'leadjet':container2b0j[0],'subleadjet':container2b0j[1],'lead_is_b':True,'sublead_is_b':True,'suffix':ElElSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.selName,'is_MC':self.is_MC})
                            ChannelDictList.append({'channel':'MuMu','sel':MuMuSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.sel,'dilepton':OSMuMuDilepton[0],'leadjet':container2b0j[0],'subleadjet':container2b0j[1],'lead_is_b':True,'sublead_is_b':True,'suffix':MuMuSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.selName,'is_MC':self.is_MC})
                            ChannelDictList.append({'channel':'ElMu','sel':ElMuSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.sel,'dilepton':OSElMuDilepton[0],'leadjet':container2b0j[0],'subleadjet':container2b0j[1],'lead_is_b':True,'sublead_is_b':True,'suffix':ElMuSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.selName,'is_MC':self.is_MC})

                # Lepton + jet Plots #

                ResolvedJetsN = {'objName':'Ak4BJets','objCont':self.ak4BJets,'Nmax':5,'xTitle':'N(Ak4 Bjets)'}
        
                for channelDict in ChannelDictList:
                    # Dilepton #
                    plots.extend(makeDileptonPlots(**{k:channelDict[k] for k in LeptonKeys}))
                    # Number of jets #
                    plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**ResolvedJetsN))
                    plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**JetsN))
                    # Ak4 Jets #
                    plots.extend(makeTwoAk4JetsPlots(**{k:channelDict[k] for k in JetKeys}))
                    # MET #
                    plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))
            
            #----- Boosted selection -----#
            if any(item in ["Boosted0Btag","Boosted1Btag"] for item in jetsel_level): # If any of the boosted category is asked 
                container1fatb = [t.FatJet[op.switch(op.rng_len(self.ak8BJets) >= 1, self.ak8BJets[0].idx, self.ak8Jets[0].idx)]]
                ChannelDictList = []
                #----- Boosted selection : 0 Btag -----#
                if "Boosted0Btag" in jetsel_level:
                    print ("...... Processing Boosted jet (0 btag) selection")
                    ElElSelObjectDileptonAk8JetsInclusiveBoostedNoBtag = makeInclusiveBoostedNoBtagSelection(self,ElElSelObjectDileptonAk8Jets,copy_sel=True,plot_yield=True)
                    MuMuSelObjectDileptonAk8JetsInclusiveBoostedNoBtag = makeInclusiveBoostedNoBtagSelection(self,MuMuSelObjectDileptonAk8Jets,copy_sel=True,plot_yield=True)
                    ElMuSelObjectDileptonAk8JetsInclusiveBoostedNoBtag = makeInclusiveBoostedNoBtagSelection(self,ElMuSelObjectDileptonAk8Jets,copy_sel=True,plot_yield=True)
                    if "Boosted0Btag" in jetplot_level:
                        # Cut flow report #
                        cutFlowPlots.append(CutFlowReport(ElElSelObjectDileptonAk8JetsInclusiveBoostedNoBtag.selName,ElElSelObjectDileptonAk8JetsInclusiveBoostedNoBtag.sel))
                        cutFlowPlots.append(CutFlowReport(MuMuSelObjectDileptonAk8JetsInclusiveBoostedNoBtag.selName,MuMuSelObjectDileptonAk8JetsInclusiveBoostedNoBtag.sel))
                        cutFlowPlots.append(CutFlowReport(ElMuSelObjectDileptonAk8JetsInclusiveBoostedNoBtag.selName,ElMuSelObjectDileptonAk8JetsInclusiveBoostedNoBtag.sel))
                        if not self.args.OnlyYield:
                            ChannelDictList.append({'channel':'ElEl','sel':ElElSelObjectDileptonAk8JetsInclusiveBoostedNoBtag.sel,'dilepton':OSElElDilepton[0],'fatjet':container1fatb[0],'suffix':ElElSelObjectDileptonAk8JetsInclusiveBoostedNoBtag.selName,'is_MC':self.is_MC})
                            ChannelDictList.append({'channel':'MuMu','sel':MuMuSelObjectDileptonAk8JetsInclusiveBoostedNoBtag.sel,'dilepton':OSMuMuDilepton[0],'fatjet':container1fatb[0],'suffix':MuMuSelObjectDileptonAk8JetsInclusiveBoostedNoBtag.selName,'is_MC':self.is_MC})
                            ChannelDictList.append({'channel':'ElMu','sel':ElMuSelObjectDileptonAk8JetsInclusiveBoostedNoBtag.sel,'dilepton':OSElMuDilepton[0],'fatjet':container1fatb[0],'suffix':ElMuSelObjectDileptonAk8JetsInclusiveBoostedNoBtag.selName,'is_MC':self.is_MC})
                #----- Boosted selection : 1 Btag -----#
                if "Boosted1Btag" in jetsel_level:
                    print ("...... Processing Boosted jet (1 btag) selection")
                    ElElSelObjectDileptonAk8JetsInclusiveBoostedOneBtag = makeInclusiveBoostedOneBtagSelection(self,ElElSelObjectDileptonAk8Jets,copy_sel=True,plot_yield=True)
                    MuMuSelObjectDileptonAk8JetsInclusiveBoostedOneBtag = makeInclusiveBoostedOneBtagSelection(self,MuMuSelObjectDileptonAk8Jets,copy_sel=True,plot_yield=True)
                    ElMuSelObjectDileptonAk8JetsInclusiveBoostedOneBtag = makeInclusiveBoostedOneBtagSelection(self,ElMuSelObjectDileptonAk8Jets,copy_sel=True,plot_yield=True)
                    if "Boosted1Btag" in jetplot_level:
                        # Cut flow report #
                        cutFlowPlots.append(CutFlowReport(ElElSelObjectDileptonAk8JetsInclusiveBoostedOneBtag.selName,ElElSelObjectDileptonAk8JetsInclusiveBoostedOneBtag.sel))
                        cutFlowPlots.append(CutFlowReport(MuMuSelObjectDileptonAk8JetsInclusiveBoostedOneBtag.selName,MuMuSelObjectDileptonAk8JetsInclusiveBoostedOneBtag.sel))
                        cutFlowPlots.append(CutFlowReport(ElMuSelObjectDileptonAk8JetsInclusiveBoostedOneBtag.selName,ElMuSelObjectDileptonAk8JetsInclusiveBoostedOneBtag.sel))
                        if not self.args.OnlyYield:
                            ChannelDictList.append({'channel':'ElEl','sel':ElElSelObjectDileptonAk8JetsInclusiveBoostedOneBtag.sel,'dilepton':OSElElDilepton[0],'fatjet':container1fatb[0],'suffix':ElElSelObjectDileptonAk8JetsInclusiveBoostedOneBtag.selName,'is_MC':self.is_MC})
                            ChannelDictList.append({'channel':'MuMu','sel':MuMuSelObjectDileptonAk8JetsInclusiveBoostedOneBtag.sel,'dilepton':OSMuMuDilepton[0],'fatjet':container1fatb[0],'suffix':MuMuSelObjectDileptonAk8JetsInclusiveBoostedOneBtag.selName,'is_MC':self.is_MC})
                            ChannelDictList.append({'channel':'ElMu','sel':ElMuSelObjectDileptonAk8JetsInclusiveBoostedOneBtag.sel,'dilepton':OSElMuDilepton[0],'fatjet':container1fatb[0],'suffix':ElMuSelObjectDileptonAk8JetsInclusiveBoostedOneBtag.selName,'is_MC':self.is_MC})

                BoostedJetsN = {'objName':'Ak8BJets','objCont':self.ak8BJets,'Nmax':5,'xTitle':'N(Ak8 Bjets)'}
                
                for channelDict in ChannelDictList:
                    # Dilepton #
                    plots.extend(makeDileptonPlots(**{k:channelDict[k] for k in LeptonKeys}))
                    # Number of jets #
                    plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**BoostedJetsN))
                    plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**FatJetsN))
                    # Ak8 Jets #
                    plots.extend(makeDoubleLeptonAk8JetsPlots(**{k:channelDict[k] for k in FatJetKeys}))
                    # MET #
                    plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))
            
            #----- High-level combinations -----#
            # NOTE : very time consuming 
            ChannelDictList = []
            if not self.args.OnlyYield:
                # Resolved No Btag #
                if "Resolved0Btag" in jetplot_level:
                    ChannelDictList.append({'channel': 'ElEl','met': self.corrMET,'l1':OSElElDilepton[0][0],'l2':OSElElDilepton[0][1],'j1':container0b2j[0],'j2':container0b2j[1],'sel':ElElSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.sel,'suffix':ElElSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.selName})
                    ChannelDictList.append({'channel': 'MuMu','met': self.corrMET,'l1':OSMuMuDilepton[0][0],'l2':OSMuMuDilepton[0][1],'j1':container0b2j[0],'j2':container0b2j[1],'sel':MuMuSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.sel,'suffix':MuMuSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.selName})
                    ChannelDictList.append({'channel': 'ElMu','met': self.corrMET,'l1':OSElMuDilepton[0][0],'l2':OSElMuDilepton[0][1],'j1':container0b2j[0],'j2':container0b2j[1],'sel':ElMuSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.sel,'suffix':ElMuSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.selName})
                # Resolved One Btag #
                if "Resolved1Btag" in jetplot_level:
                    ChannelDictList.append({'channel': 'ElEl','met': self.corrMET,'l1':OSElElDilepton[0][0],'l2':OSElElDilepton[0][1],'j1':container1b1j[0],'j2':container1b1j[1],'sel':ElElSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.sel,'suffix':ElElSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.selName})
                    ChannelDictList.append({'channel': 'MuMu','met': self.corrMET,'l1':OSMuMuDilepton[0][0],'l2':OSMuMuDilepton[0][1],'j1':container1b1j[0],'j2':container1b1j[1],'sel':MuMuSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.sel,'suffix':MuMuSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.selName})
                    ChannelDictList.append({'channel': 'ElMu','met': self.corrMET,'l1':OSElMuDilepton[0][0],'l2':OSElMuDilepton[0][1],'j1':container1b1j[0],'j2':container1b1j[1],'sel':ElMuSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.sel,'suffix':ElMuSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.selName})
                # Resolved Two Btags  #
                if "Resolved2Btag" in jetplot_level:
                    ChannelDictList.append({'channel': 'ElEl','met': self.corrMET,'l1':OSElElDilepton[0][0],'l2':OSElElDilepton[0][1],'j1':container2b0j[0],'j2':container2b0j[1],'sel':ElElSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.sel,'suffix':ElElSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.selName})
                    ChannelDictList.append({'channel': 'MuMu','met': self.corrMET,'l1':OSMuMuDilepton[0][0],'l2':OSMuMuDilepton[0][1],'j1':container2b0j[0],'j2':container2b0j[1],'sel':MuMuSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.sel,'suffix':MuMuSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.selName})
                    ChannelDictList.append({'channel': 'ElMu','met': self.corrMET,'l1':OSElMuDilepton[0][0],'l2':OSElMuDilepton[0][1],'j1':container2b0j[0],'j2':container2b0j[1],'sel':ElMuSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.sel,'suffix':ElMuSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.selName})
                # Boosted No Btag #
                if "Boosted0Btag" in jetplot_level:
                    ChannelDictList.append({'channel': 'ElEl','met': self.corrMET,'l1':OSElElDilepton[0][0],'l2':OSElElDilepton[0][1],'j1':container1fatb[0].subJet1,'j2':container1fatb[0].subJet2,'sel':ElElSelObjectDileptonAk8JetsInclusiveBoostedNoBtag.sel,'suffix':ElElSelObjectDileptonAk8JetsInclusiveBoostedNoBtag.selName})
                    ChannelDictList.append({'channel': 'MuMu','met': self.corrMET,'l1':OSMuMuDilepton[0][0],'l2':OSMuMuDilepton[0][1],'j1':container1fatb[0].subJet1,'j2':container1fatb[0].subJet2,'sel':MuMuSelObjectDileptonAk8JetsInclusiveBoostedNoBtag.sel,'suffix':MuMuSelObjectDileptonAk8JetsInclusiveBoostedNoBtag.selName})
                    ChannelDictList.append({'channel': 'ElMu','met': self.corrMET,'l1':OSElMuDilepton[0][0],'l2':OSElMuDilepton[0][1],'j1':container1fatb[0].subJet1,'j2':container1fatb[0].subJet2,'sel':ElMuSelObjectDileptonAk8JetsInclusiveBoostedNoBtag.sel,'suffix':ElMuSelObjectDileptonAk8JetsInclusiveBoostedNoBtag.selName})
                # Boosted One Btag #
                if "Boosted1Btag" in jetplot_level:
                    ChannelDictList.append({'channel': 'ElEl','met': self.corrMET,'l1':OSElElDilepton[0][0],'l2':OSElElDilepton[0][1],'j1':container1fatb[0].subJet1,'j2':container1fatb[0].subJet2,'sel':ElElSelObjectDileptonAk8JetsInclusiveBoostedOneBtag.sel,'suffix':ElElSelObjectDileptonAk8JetsInclusiveBoostedOneBtag.selName})
                    ChannelDictList.append({'channel': 'MuMu','met': self.corrMET,'l1':OSMuMuDilepton[0][0],'l2':OSMuMuDilepton[0][1],'j1':container1fatb[0].subJet1,'j2':container1fatb[0].subJet2,'sel':MuMuSelObjectDileptonAk8JetsInclusiveBoostedOneBtag.sel,'suffix':MuMuSelObjectDileptonAk8JetsInclusiveBoostedOneBtag.selName})
                    ChannelDictList.append({'channel': 'ElMu','met': self.corrMET,'l1':OSElMuDilepton[0][0],'l2':OSElMuDilepton[0][1],'j1':container1fatb[0].subJet1,'j2':container1fatb[0].subJet2,'sel':ElMuSelObjectDileptonAk8JetsInclusiveBoostedOneBtag.sel,'suffix':ElMuSelObjectDileptonAk8JetsInclusiveBoostedOneBtag.selName})

            for channelDict in ChannelDictList:
                plots.extend(makeDoubleLeptonHighLevelQuantities(**channelDict,HLL=self.HLL))

            #----- Machine Learning plots -----#
            selObjectDictList = []
            bjetsCont = {}
            nbtagconv = {}
            if not self.args.OnlyYield:
                if "Ak4" in jetplot_level:
                    selObjectDictList.append({'channel':'ElEl','selObject':ElElSelObjectDileptonAk4Jets})
                    selObjectDictList.append({'channel':'MuMu','selObject':MuMuSelObjectDileptonAk4Jets})
                    selObjectDictList.append({'channel':'ElMu','selObject':ElMuSelObjectDileptonAk4Jets})
                if "Resolved0Btag" in jetplot_level:
                    selObjectDictList.append({'channel':'ElEl','selObject':ElElSelObjectDileptonAk4JetsExclusiveResolvedNoBtag})
                    selObjectDictList.append({'channel':'MuMu','selObject':MuMuSelObjectDileptonAk4JetsExclusiveResolvedNoBtag})
                    selObjectDictList.append({'channel':'ElMu','selObject':ElMuSelObjectDileptonAk4JetsExclusiveResolvedNoBtag})
                    bjetsCont['NoBtag'] = container0b2j
                    nbtagconv['NoBtag'] = 0
                if "Resolved1Btag" in jetplot_level:
                    selObjectDictList.append({'channel':'ElEl','selObject':ElElSelObjectDileptonAk4JetsExclusiveResolvedOneBtag})
                    selObjectDictList.append({'channel':'MuMu','selObject':MuMuSelObjectDileptonAk4JetsExclusiveResolvedOneBtag})
                    selObjectDictList.append({'channel':'ElMu','selObject':ElMuSelObjectDileptonAk4JetsExclusiveResolvedOneBtag})
                    bjetsCont['OneBtag'] = container1b1j
                    nbtagconv['OneBtag'] = 1
                if "Resolved2Btag" in jetplot_level:
                    selObjectDictList.append({'channel':'ElEl','selObject':ElElSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags})
                    selObjectDictList.append({'channel':'MuMu','selObject':MuMuSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags})
                    selObjectDictList.append({'channel':'ElMu','selObject':ElMuSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags})
                    bjetsCont['TwoBtags'] = container2b0j
                    nbtagconv['TwoBtags'] = 2

            dileptonsCont = {'ElEl':OSElElDilepton[0],'MuMu':OSMuMuDilepton[0],'ElMu':OSElMuDilepton[0]}
            self.nodes = ['HH', 'H', 'DY', 'ST', 'TT', 'TTVX', 'VVV', 'Rare']
            
            for selObjectDict in selObjectDictList:
                dilepton = dileptonsCont[selObjectDict['channel']]
                resolved_str = [res_str for res_str in bjetsCont.keys() if res_str in selObjectDict["selObject"].selName][0]
                dibjet = bjetsCont[resolved_str]
                inputsLL = returnLowLevelMVAInputs(self     = self,
                                                   dilepton = dilepton,
                                                   jets     = self.ak4Jets,
                                                   channel  = selObjectDict['channel'])
                inputsHL = returnHighLevelMVAInputs(self      = self,
                                                    dilepton  = dilepton,
                                                    b1        = dibjet[0],
                                                    b2        = dibjet[1],
                                                    met       = self.corrMET,
                                                    jets      = self.ak4Jets,
                                                    nbtags    = nbtagconv[resolved_str],
                                                    channel  = selObjectDict['channel'])
                plots.extend(makeDoubleLeptonMachineLearningInputPlots(selObjectDict['selObject'].sel,selObjectDict['selObject'].selName,selObjectDict['channel'],inputsLL))
                plots.extend(makeDoubleLeptonMachineLearningInputPlots(selObjectDict['selObject'].sel,selObjectDict['selObject'].selName,selObjectDict['channel'],inputsHL))
                
                for model_num,DNN in DNNs.items():
                    if model_num in ["01","02"]:
                        inputs = inputsLL
                    elif model_num in ["03","04"]:
                        inputs = {**inputsLL,**inputsHL}
                    else:
                        raise RuntimeError("Failed to understand model number")

                    output = DNN(*inputs.values())
                    selObjNodesDict = makeDNNOutputNodesSelections(self,selObjectDict['selObject'],output,plot_yield=True,suffix=model_num)
                    for selObjNode in selObjNodesDict.values():
                        cutFlowPlots.append(CutFlowReport(selObjNode.selName,selObjNode.sel))
                    if not self.args.OnlyYield:
                        plots.extend(makeDoubleLeptonMachineLearningOutputPlots(selObjNodesDict,output,self.nodes,channel=selObjectDict['channel']))
            #        plots.append(plots.append(Plot.make1D("DNN_argmax_"+selObjectDict['selObject'].selName,
            #                                              op.rng_max_element_index(output),
            #                                              selObjectDict['selObject'].sel,
            #                                              EquidistantBinning(10,0.,10.),
            #                                              xTitle = 'DNN output')))
                
    
        #----- Add the Yield plots -----#
        plots.extend(self.yieldPlots.returnPlots())

        #plots.extend(cutFlowPlots)

        #----- Return -----#
        return plots
