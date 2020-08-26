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

def switch_on_index(indexes, condition, contA, contB):
    if contA._base != contB._base:
        raise RuntimeError("The containers do not derive from the same base, this won't work")
    base = contA._base
    return [base[op.switch(condition, contA[index].idx, contB[index].idx)] for index in indexes]       


#===============================================================================================#
#                                       PlotterHHtobbWW                                         #
#===============================================================================================#
class PlotterNanoHHtobbWWSL(BaseNanoHHtobbWW,DataDrivenBackgroundHistogramsModule):
    """ Plotter module: HH->bbW(->e/µ nu)W(->e/µ nu) histograms from NanoAOD """
    def __init__(self, args):
        super(PlotterNanoHHtobbWWSL, self).__init__(args)

    def initialize(self):
        super(PlotterNanoHHtobbWWSL, self).initialize()
        # Change the way the FakeExtrapolation is postProcesses (avoids overriding the `postProcess` method) 
        if "FakeExtrapolation" in self.datadrivenContributions:
            contrib = self.datadrivenContributions["FakeExtrapolation"]
            self.datadrivenContributions["FakeExtrapolation"] = DataDrivenFake(contrib.name, contrib.config)
        if "DYEstimation" in self.datadrivenContributions: 
            contrib = self.datadrivenContributions["DYEstimation"]
            self.datadrivenContributions["DYEstimation"] = DataDrivenDY(contrib.name, contrib.config)


    def definePlots(self, t, noSel, sample=None, sampleCfg=None): 
        noSel = super(PlotterNanoHHtobbWWSL,self).prepareObjects(t, noSel, sample, sampleCfg, 'SL')
        
        #----- Machine Learning Model -----#
        # Whad Tagger #
        
        # Event Level #
        # SM
        path_fullRecoSM_even_simple_model   = os.path.join(os.path.abspath(os.path.dirname(__file__)),'MachineLearning','hh_bb1l_SM_Wjj_simple_full_reco_only_noIndPt_even.xml')
        path_missRecoSM_even_simple_model   = os.path.join(os.path.abspath(os.path.dirname(__file__)),'MachineLearning','hh_bb1l_SM_Wj1_even.xml')
        path_fullRecoSM_odd_simple_model    = os.path.join(os.path.abspath(os.path.dirname(__file__)),'MachineLearning','hh_bb1l_SM_Wjj_simple_full_reco_only_noIndPt_odd.xml')
        path_missRecoSM_odd_simple_model    = os.path.join(os.path.abspath(os.path.dirname(__file__)),'MachineLearning','hh_bb1l_SM_Wj1_odd.xml')
        # 900_Radion
        path_fullReco900R_even_simple_model = os.path.join(os.path.abspath(os.path.dirname(__file__)),'MachineLearning','hh_bb1l_X900GeV_Wjj_simple_full_reco_even.xml')
        path_missReco900R_even_simple_model = os.path.join(os.path.abspath(os.path.dirname(__file__)),'MachineLearning','hh_bb1l_X900GeV_Wj1_even.xml')
        path_fullReco900R_odd_simple_model  = os.path.join(os.path.abspath(os.path.dirname(__file__)),'MachineLearning','hh_bb1l_X900GeV_Wjj_simple_full_reco_odd.xml')
        path_missReco900R_odd_simple_model  = os.path.join(os.path.abspath(os.path.dirname(__file__)),'MachineLearning','hh_bb1l_X900GeV_Wj1_odd.xml')
        #if not os.path.exists(path_model):
            #raise RuntimeError('Could not find model file %s'%path_model)
        try:
            BDT_fullRecoSM_simple_even   = op.mvaEvaluator(path_fullRecoSM_even_simple_model, mvaType='TMVA')
            BDT_missRecoSM_simple_even   = op.mvaEvaluator(path_missRecoSM_even_simple_model, mvaType='TMVA')
            BDT_fullRecoSM_simple_odd    = op.mvaEvaluator(path_fullRecoSM_odd_simple_model, mvaType='TMVA')
            BDT_missRecoSM_simple_odd    = op.mvaEvaluator(path_missRecoSM_odd_simple_model, mvaType='TMVA')
            BDT_fullReco900R_simple_even = op.mvaEvaluator(path_fullReco900R_even_simple_model, mvaType='TMVA')
            BDT_missReco900R_simple_even = op.mvaEvaluator(path_missReco900R_even_simple_model, mvaType='TMVA')
            BDT_fullReco900R_simple_odd  = op.mvaEvaluator(path_fullReco900R_odd_simple_model, mvaType='TMVA')
            BDT_missReco900R_simple_odd  = op.mvaEvaluator(path_missReco900R_odd_simple_model, mvaType='TMVA')
        except:
            raise RuntimeError('Could not load model %s'%path_model)
        
        plots = []

        cutFlowPlots = []

        era = sampleCfg['era']

        self.sample = sample
        self.sampleCfg = sampleCfg
        self.era = era

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


        #----- Singleleptons -----#
        selObjectDict = makeSingleLeptonSelection(self,noSel,plot_yield=True)
        # selObjectDict : keys -> level (str)
        #                 values -> [ElEl,MuMu,ElMu] x Selection object
        # Select the jets selections that will be done depending on user input #
        jet_level = ["Ak4","Ak8",
                     "LooseResolved0b3j","LooseResolved1b2j","LooseResolved2b1j",
                     "TightResolved0b4j","TightResolved1b3j","TightResolved2b2j",
                     "SemiBoostedHbbWtoJ","SemiBoostedHbbWtoJJ","SemiBoostedWjj",
                     "Boosted"]
        jetplot_level = [arg for (arg,boolean) in self.args.__dict__.items() if arg in jet_level and boolean]
        if len(jetplot_level) == 0:  
            jetplot_level = jet_level # If nothing said, will do all
        jetsel_level = copy(jetplot_level)  # A plot level might need a previous selection that needs to be defined but not necessarily plotted
        if any("Resolved" in item for item in jetsel_level):     
            jetsel_level.append("Ak4") # Resolved needs the Ak4 selection
        if any("Boosted" in item for item in jetsel_level):     
            jetsel_level.append("Ak8") # SemiBoosted & Boosted needs the Ak8 selection

        # Selections:    
        # Loop over lepton selection and start plotting #
        for selectionType, selectionList in selObjectDict.items():
            print ("... Processing %s lepton type"%selectionType)
            #----- Select correct dilepton -----#
            if selectionType == "Preselected":  
                ElColl = self.electronsPreSel
                MuColl = self.muonsPreSel
            elif selectionType == "Fakeable":
                ElColl = self.leadElectronFakeSel
                MuColl = self.leadMuonFakeSel
            elif selectionType == "Tight":
                ElColl = [t.Electron[op.switch(op.rng_len(self.leadElectronTightSel) == 1, self.leadElectronTightSel[0].idx, 
                                                         self.leadElectronFakeExtrapolationSel[0].idx)]]
                MuColl = [t.Muon[op.switch(op.rng_len(self.leadMuonTightSel) == 1, self.leadMuonTightSel[0].idx, 
                                                     self.leadMuonFakeExtrapolationSel[0].idx)]]
            elif selectionType == "FakeExtrapolation":
                ElColl = self.leadElectronFakeExtrapolationSel
                MuColl = self.leadMuonFakeExtrapolationSel

            #----- Separate selections ------#
            ElSelObj = selectionList[0]
            MuSelObj = selectionList[1]

            if not self.args.OnlyYield:
                ChannelDictList = []
                ChannelDictList.append({'channel':'El','sel':ElSelObj.sel,'suffix':ElSelObj.selName})
                ChannelDictList.append({'channel':'Mu','sel':MuSelObj.sel,'suffix':MuSelObj.selName})
                
                for channelDict in ChannelDictList:
                    #----- Trigger plots -----#
                    plots.extend(singleLeptonTriggerPlots(**channelDict, triggerDict=self.triggersPerPrimaryDataset))
                    #----- Lepton plots -----#
                    # Singlelepton channel plots #
                    #plots.extend(singleLeptonChannelPlot(**channelDict, SinlepEl=ElColl, SinlepMu=MuColl, suffix=ElSelObj.selName))
        

            # ------------- test ------------- #
            testDictList=[]
            testKeys = ['channel','sel','jet','suffix']
            _ElSelObj = copy(ElSelObj)
            _MuSelObj = copy(MuSelObj)
            _ElSelObj.selName    += 'hasMinOneAk8Jet'
            _ElSelObj.yieldTitle += " + at least one Ak8"
            _ElSelObj.refine(cut = [op.rng_len(self.ak8Jets) >= 1], weight = None)
            _MuSelObj.selName    += 'hasMinOneAk8Jet'
            _MuSelObj.yieldTitle += " + at least one Ak8"
            _MuSelObj.refine(cut = [op.rng_len(self.ak8Jets) >= 1], weight = None)
            testDictList.append({'channel':'El','sel':_ElSelObj.sel,'jet':self.ak8Jets[0],'suffix':_ElSelObj.selName})
            testDictList.append({'channel':'Mu','sel':_MuSelObj.sel,'jet':self.ak8Jets[0],'suffix':_MuSelObj.selName})
            for testDict in testDictList:
                plots.extend(makeBtagDDplots(**{k:testDict[k] for k in testKeys}))
            #--------------------------------- #

            #----- Ak4 jets selection -----#
            LeptonKeys  = ['channel','sel','lep','suffix','is_MC']
            JetKeys     = ['channel','sel','j1','j2','j3','j4','suffix','nJet','nbJet','is_MC']
            commonItems = ['channel','sel','suffix']
            
            if "Ak4" in jetsel_level:
                print("... Processing Ak4Jets Selection for Resolved category")
                ChannelDictList = []
                JetsN    = {'objName':'Ak4Jets','objCont':self.ak4Jets,'Nmax':10,'xTitle':'N(Ak4 jets)'}
                FatJetsN = {'objName':'Ak8Jets','objCont':self.ak8Jets,'Nmax':5,'xTitle':'N(Ak8 jets)'}

                if any("LooseResolved" in key for key in jetsel_level):
                    print ("...... Processing Ak4 jet selection Loose (nAk4Jets = 3)")
                    ElSelObjAk4JetsLoose = makeCoarseResolvedSelection(self,ElSelObj,nJet=3,copy_sel=True,plot_yield=True)
                    MuSelObjAk4JetsLoose = makeCoarseResolvedSelection(self,MuSelObj,nJet=3,copy_sel=True,plot_yield=True)
                    if not self.args.OnlyYield:
                        # cutFlow Report #
                        cutFlowPlots.append(CutFlowReport(ElSelObjAk4JetsLoose.selName, ElSelObjAk4JetsLoose.sel))
                        cutFlowPlots.append(CutFlowReport(MuSelObjAk4JetsLoose.selName, MuSelObjAk4JetsLoose.sel))

                        if "Ak4" in jetplot_level and any("LooseResolved" in key for key in jetplot_level):
                            ChannelDictList.append({'channel':'El','sel':ElSelObjAk4JetsLoose.sel,'sinlepton':ElColl[0],
                                                    'j1':self.ak4Jets[0],'j2':self.ak4Jets[1],'j3':self.ak4Jets[2],'j4':None,
                                                    'nJet':3,'nbJet':0,
                                                    'suffix':ElSelObjAk4JetsLoose.selName,
                                                    'is_MC':self.is_MC})
                            ChannelDictList.append({'channel':'Mu','sel':MuSelObjAk4JetsLoose.sel,'sinlepton':MuColl[0],
                                                    'j1':self.ak4Jets[0],'j2':self.ak4Jets[1],'j3':self.ak4Jets[2],'j4':None,
                                                    'nJet':3,'nbJet':0,
                                                    'suffix':MuSelObjAk4JetsLoose.selName,
                                                    'is_MC':self.is_MC})

                if any("TightResolved" in key for key in jetsel_level):        
                    print ("...... Processing Ak4 jet selection Tight (nAk4Jets >= 4)")
                    ElSelObjAk4JetsTight = makeCoarseResolvedSelection(self,ElSelObj,nJet=4,copy_sel=True,plot_yield=True)
                    MuSelObjAk4JetsTight = makeCoarseResolvedSelection(self,MuSelObj,nJet=4,copy_sel=True,plot_yield=True)

                    # Jet and lepton plots #
                    if not self.args.OnlyYield:
                        # cutFlow Report #
                        cutFlowPlots.append(CutFlowReport(ElSelObjAk4JetsTight.selName, ElSelObjAk4JetsTight.sel))
                        cutFlowPlots.append(CutFlowReport(MuSelObjAk4JetsTight.selName, MuSelObjAk4JetsTight.sel))
                        if "Ak4" in jetplot_level and any("TightResolved" in key for key in jetplot_level):
                            ChannelDictList.append({'channel':'El','sel':ElSelObjAk4JetsTight.sel,'sinlepton':ElColl[0],
                                                    'j1':self.ak4Jets[0],'j2':self.ak4Jets[1],'j3':self.ak4Jets[2],'j4':self.ak4Jets[3],
                                                    'nJet':4,'nbJet':0,
                                                    'suffix':ElSelObjAk4JetsTight.selName,
                                                    'is_MC':self.is_MC})
                            ChannelDictList.append({'channel':'Mu','sel':MuSelObjAk4JetsTight.sel,'sinlepton':MuColl[0],
                                                    'j1':self.ak4Jets[0],'j2':self.ak4Jets[1],'j3':self.ak4Jets[2],'j4':self.ak4Jets[3],
                                                    'nJet':4,'nbJet':0,
                                                    'suffix':MuSelObjAk4JetsTight.selName,
                                                    'is_MC':self.is_MC})

                for channelDict in ChannelDictList:
                    # Dilepton #
                    plots.extend(makeSinleptonPlots(**{k:channelDict[k] for k in LeptonKeys}))
                    # Number of jets #
                    plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**JetsN))
                    plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**FatJetsN))
                    # Ak4 Jets #
                    plots.extend(makeAk4JetsPlots(**{k:channelDict[k] for k in JetKeys},HLL=self.HLL))
                    # MET #
                    plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))

            ##### Ak8-b jets selection #####
            if "Ak8" in jetsel_level:
                print ("...... Processing Ak8b jet selection for SemiBoosted & Boosted Category")
                ElSelObjAk8bJets = makeCoarseBoostedSelection(self,ElSelObj,copy_sel=True,plot_yield=True)
                MuSelObjAk8bJets = makeCoarseBoostedSelection(self,MuSelObj,copy_sel=True,plot_yield=True)

                FatJetKeys = ['channel','sel','j1','j2','j3','has1fat','suffix']
                FatJetsN   = {'objName':'Ak8Jets','objCont':self.ak8Jets,'Nmax':5,'xTitle':'N(Ak8 jets)'}
                SlimJetsN  = {'objName':'Ak4Jets','objCont':self.ak4Jets,'Nmax':10,'xTitle':'N(Ak4 jets)'}

                # Fatjets plots #
                ChannelDictList = []
                if not self.args.OnlyYield:
                    # cutFlow Report #
                    cutFlowPlots.append(CutFlowReport(ElSelObjAk8bJets.selName, ElSelObjAk8bJets.sel))
                    cutFlowPlots.append(CutFlowReport(MuSelObjAk8bJets.selName, MuSelObjAk8bJets.sel))
                    if "Ak8" in jetplot_level:
                        ChannelDictList.append({'channel':'El','sel':ElSelObjAk8bJets.sel,'sinlepton':ElColl[0],
                                                'j1':self.ak8BJets[0],'j2':None,'j3':None,'has1fat':True,
                                                'suffix':ElSelObjAk8bJets.selName,'is_MC':self.is_MC})
                        ChannelDictList.append({'channel':'Mu','sel':MuSelObjAk8bJets.sel,'sinlepton':MuColl[0],
                                                'j1':self.ak8BJets[0],'j2':None,'j3':None,'has1fat':True,
                                                'suffix':MuSelObjAk8bJets.selName,'is_MC':self.is_MC})

                for channelDict in ChannelDictList:
                    # Dilepton #
                    plots.extend(makeSinleptonPlots(**{k:channelDict[k] for k in LeptonKeys}))
                    # Number of jets #
                    plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**FatJetsN))
                    plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**SlimJetsN))
                    # Ak8 Jets #
                    plots.extend(makeSingleLeptonAk8JetsPlots(**{k:channelDict[k] for k in FatJetKeys}))
                    # MET #
                    plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))

                         
            #-----------------------------|||||||||||||||| Resolved selection ||||||||||||||-----------------------------------#
            if any("Resolved" in item for item in jetsel_level):
                ResolvedKeys = ['channel','sel','met','lep','j1','j2','j3','j4','suffix','nJet','nbJet']
                ChannelDictList = []
                ChannelDictListML = []

                # Resolved Selection (Loose) #
                #----- Resolved selection : 0 Btag -----#
                if "LooseResolved0b3j" in jetsel_level:
                    print ("......... Processing Loose Resolved jet (0 btag i.e. bTaggedJets = 0 & nLightJets = 3) selection")
                    ElSelObjAk4JetsLooseExclusiveResolved0b3j = makeExclusiveLooseResolvedJetComboSelection(self,ElSelObjAk4JetsLoose,nbJet=0,copy_sel=True,plot_yield=True)
                    MuSelObjAk4JetsLooseExclusiveResolved0b3j = makeExclusiveLooseResolvedJetComboSelection(self,MuSelObjAk4JetsLoose,nbJet=0,copy_sel=True,plot_yield=True)

                    if not self.args.OnlyYield and "LooseResolved0b3j" in jetplot_level:
                        # Cut flow report #
                        cutFlowPlots.append(CutFlowReport(ElSelObjAk4JetsLooseExclusiveResolved0b3j.selName,ElSelObjAk4JetsLooseExclusiveResolved0b3j.sel))
                        cutFlowPlots.append(CutFlowReport(MuSelObjAk4JetsLooseExclusiveResolved0b3j.selName,MuSelObjAk4JetsLooseExclusiveResolved0b3j.sel))

                        ChannelDictList.append({'channel':'El','sel':ElSelObjAk4JetsLooseExclusiveResolved0b3j.sel,'lep':ElColl[0],'met':self.corrMET,
                                                'j1':self.ak4LightJetsByPt[0],'j2':self.ak4LightJetsByPt[1],'j3':self.ak4LightJetsByPt[2],'j4':None,
                                                'nJet':3,'nbJet':0,
                                                'suffix':ElSelObjAk4JetsLooseExclusiveResolved0b3j.selName,
                                                'is_MC':self.is_MC})
                        ChannelDictList.append({'channel':'Mu','sel':MuSelObjAk4JetsLooseExclusiveResolved0b3j.sel,'lep':MuColl[0],'met':self.corrMET,
                                                'j1':self.ak4LightJetsByPt[0],'j2':self.ak4LightJetsByPt[1],'j3':self.ak4LightJetsByPt[2],'j4':None,
                                                'nJet':3,'nbJet':0,
                                                'suffix':MuSelObjAk4JetsLooseExclusiveResolved0b3j.selName,
                                                'is_MC':self.is_MC})


                #----  Resolved selection : 1 Btag  -----#
                if "LooseResolved1b2j" in jetsel_level:
                    print ("......... Processing Resolved jet (1 btag i.e. bTaggedJets = 1 & nLightJets = 2) selection")
                    ElSelObjAk4JetsLooseExclusiveResolved1b2j = makeExclusiveLooseResolvedJetComboSelection(self,ElSelObjAk4JetsLoose,nbJet=1,copy_sel=True,plot_yield=True)
                    MuSelObjAk4JetsLooseExclusiveResolved1b2j = makeExclusiveLooseResolvedJetComboSelection(self,MuSelObjAk4JetsLoose,nbJet=1,copy_sel=True,plot_yield=True)

                    if not self.args.OnlyYield and "LooseResolved1b2j" in jetplot_level:
                        # Cut flow report #
                        cutFlowPlots.append(CutFlowReport(ElSelObjAk4JetsLooseExclusiveResolved1b2j.selName,ElSelObjAk4JetsLooseExclusiveResolved1b2j.sel))
                        cutFlowPlots.append(CutFlowReport(MuSelObjAk4JetsLooseExclusiveResolved1b2j.selName,MuSelObjAk4JetsLooseExclusiveResolved1b2j.sel))

                        ChannelDictList.append({'channel':'El','sel':ElSelObjAk4JetsLooseExclusiveResolved1b2j.sel,'lep':ElColl[0],'met':self.corrMET,
                                                'j1':self.ak4BJets[0],'j2':self.ak4LightJetsByBtagScore[0],'j3':self.remainingJets[0],'j4':None,
                                                'nJet':3,'nbJet':1,
                                                'suffix':ElSelObjAk4JetsLooseExclusiveResolved1b2j.selName,
                                                'is_MC':self.is_MC})
                        ChannelDictList.append({'channel':'Mu','sel':MuSelObjAk4JetsLooseExclusiveResolved1b2j.sel,'lep':MuColl[0],'met':self.corrMET,
                                                'j1':self.ak4BJets[0],'j2':self.ak4LightJetsByBtagScore[0],'j3':self.remainingJets[0],'j4':None,
                                                'nJet':3,'nbJet':1,
                                                'suffix':MuSelObjAk4JetsLooseExclusiveResolved1b2j.selName,
                                                'is_MC':self.is_MC})

                #----- Resolved selection : 2 Btags -----#
                if "LooseResolved2b1j" in jetsel_level:
                    print ("......... Processing Resolved jet (2 btag i.e. bTaggedJets = 2 & nLightJets = 1) selection")
                    ElSelObjAk4JetsLooseExclusiveResolved2b1j = makeExclusiveLooseResolvedJetComboSelection(self,ElSelObjAk4JetsLoose,nbJet=2,copy_sel=True,plot_yield=True)
                    MuSelObjAk4JetsLooseExclusiveResolved2b1j = makeExclusiveLooseResolvedJetComboSelection(self,MuSelObjAk4JetsLoose,nbJet=2,copy_sel=True,plot_yield=True)

                    if not self.args.OnlyYield and "LooseResolved2b1j" in jetplot_level:
                        # Cut flow report #
                        cutFlowPlots.append(CutFlowReport(ElSelObjAk4JetsLooseExclusiveResolved2b1j.selName,ElSelObjAk4JetsLooseExclusiveResolved2b1j.sel))
                        cutFlowPlots.append(CutFlowReport(MuSelObjAk4JetsLooseExclusiveResolved2b1j.selName,MuSelObjAk4JetsLooseExclusiveResolved2b1j.sel))

                        ChannelDictList.append({'channel':'El',
                                                'sel':ElSelObjAk4JetsLooseExclusiveResolved2b1j.sel,'lep':ElColl[0],'met':self.corrMET,
                                                'j1':self.ak4BJets[0],'j2':self.ak4BJets[1],'j3':self.ak4LightJetsByPt[0],'j4':None,
                                                'nJet':3,'nbJet':2,
                                                'suffix':ElSelObjAk4JetsLooseExclusiveResolved2b1j.selName,
                                                'is_MC':self.is_MC})
                        ChannelDictList.append({'channel':'Mu','sel':MuSelObjAk4JetsLooseExclusiveResolved2b1j.sel,'lep':MuColl[0],'met':self.corrMET,
                                                'j1':self.ak4BJets[0],'j2':self.ak4BJets[1],'j3':self.ak4LightJetsByPt[0],'j4':None,
                                                'nJet':3,'nbJet':2,
                                                'suffix':MuSelObjAk4JetsLooseExclusiveResolved2b1j.selName,
                                                'is_MC':self.is_MC})

                #------------------ Resolved selection (Tight) ---------------------#
                #----- Resolved selection : 0 Btag -----#
                if "TightResolved0b4j" in jetsel_level:
                    print ("......... Processing Tight Resolved jet (0 btag i.e. bTaggedJets = 0 & nLightJets >= 4) selection")
                    ElSelObjAk4JetsTightExclusiveResolved0b4j = makeExclusiveTightResolvedJetComboSelection(self,ElSelObjAk4JetsTight,nbJet=0,copy_sel=True,plot_yield=True)
                    MuSelObjAk4JetsTightExclusiveResolved0b4j = makeExclusiveTightResolvedJetComboSelection(self,MuSelObjAk4JetsTight,nbJet=0,copy_sel=True,plot_yield=True)

                    if not self.args.OnlyYield and "TightResolved0b4j" in jetplot_level:
                        ChannelDictList.append({'channel':'El','sel':ElSelObjAk4JetsTightExclusiveResolved0b4j.sel,'lep':ElColl[0],'met':self.corrMET,
                                                'j1':self.ak4LightJetsByPt[0],'j2':self.ak4LightJetsByPt[1],'j3':self.ak4LightJetsByPt[2],'j4':self.ak4LightJetsByPt[3],
                                                'nJet':4,'nbJet':0,
                                                'suffix':ElSelObjAk4JetsTightExclusiveResolved0b4j.selName,
                                                'is_MC':self.is_MC})
                        ChannelDictList.append({'channel':'Mu','sel':MuSelObjAk4JetsTightExclusiveResolved0b4j.sel,'lep':MuColl[0],'met':self.corrMET,
                                                'j1':self.ak4LightJetsByPt[0],'j2':self.ak4LightJetsByPt[1],'j3':self.ak4LightJetsByPt[2],'j4':self.ak4LightJetsByPt[3],
                                                'nJet':4,'nbJet':0,
                                                'suffix':MuSelObjAk4JetsTightExclusiveResolved0b4j.selName,
                                                'is_MC':self.is_MC})


                #----  Resolved selection : 1 Btag  -----#
                if "TightResolved1b3j" in jetsel_level:
                    print ("......... Processing Resolved jet (1 btag i.e. bTaggedJets = 1 & nLightJets >= 3) selection")
                    ElSelObjAk4JetsTightExclusiveResolved1b3j = makeExclusiveTightResolvedJetComboSelection(self,ElSelObjAk4JetsTight,nbJet=1,copy_sel=True,plot_yield=True)
                    MuSelObjAk4JetsTightExclusiveResolved1b3j = makeExclusiveTightResolvedJetComboSelection(self,MuSelObjAk4JetsTight,nbJet=1,copy_sel=True,plot_yield=True)

                    if not self.args.OnlyYield and "TightResolved1b3j" in jetplot_level:

                        cutFlowPlots.append(CutFlowReport(ElSelObjAk4JetsTightExclusiveResolved1b3j.selName,ElSelObjAk4JetsTightExclusiveResolved1b3j.sel))
                        cutFlowPlots.append(CutFlowReport(MuSelObjAk4JetsTightExclusiveResolved1b3j.selName,MuSelObjAk4JetsTightExclusiveResolved1b3j.sel))

                        # Simple reco method to find the jet paits coming from W
                        # https://gitlab.cern.ch/cms-hh-bbww/cms-hh-to-bbww/-/blob/master/Legacy/signal_extraction.md
                        self.remainingJetPairs_mu1b3j = op.combine(self.remainingJets, N=2)
                        self.lambda_chooseWjj_mu1b3j  = lambda dijet: op.abs((dijet[0].p4+dijet[1].p4+MuColl[0].p4+self.corrMET.p4).M() - 
                                                                             (bJetCorrP4(self.ak4BJets[0]) + bJetCorrP4(self.ak4LightJetsByBtagScore[0])).M()) 
                        self.WjjPairs_mu1b3j = op.sort(self.remainingJetPairs_mu1b3j, self.lambda_chooseWjj_mu1b3j)
                        
                        self.remainingJetPairs_el1b3j = op.combine(self.remainingJets, N=2)
                        self.lambda_chooseWjj_el1b3j  = lambda dijet: op.abs((dijet[0].p4+dijet[1].p4+ElColl[0].p4+self.corrMET.p4).M() - 
                                                                             (bJetCorrP4(self.ak4BJets[0]) + bJetCorrP4(self.ak4LightJetsByBtagScore[0])).M()) 
                        self.WjjPairs_el1b3j = op.sort(self.remainingJetPairs_el1b3j, self.lambda_chooseWjj_el1b3j)

                        
                        ChannelDictList.append({'channel':'El','sel':ElSelObjAk4JetsTightExclusiveResolved1b3j.sel,'lep':ElColl[0],'met':self.corrMET,
                                                'j1':self.ak4BJets[0],'j2':self.ak4LightJetsByBtagScore[0],'j3':self.WjjPairs_el1b3j[0][0],'j4':self.WjjPairs_el1b3j[0][1],
                                                'nJet':4,'nbJet':1,
                                                'suffix':ElSelObjAk4JetsTightExclusiveResolved1b3j.selName,
                                                'is_MC':self.is_MC})
                        ChannelDictList.append({'channel':'Mu','sel':MuSelObjAk4JetsTightExclusiveResolved1b3j.sel,'lep':MuColl[0],'met':self.corrMET,
                                                'j1':self.ak4BJets[0],'j2':self.ak4LightJetsByBtagScore[0],'j3':self.WjjPairs_mu1b3j[0][0],'j4':self.WjjPairs_mu1b3j[0][1],
                                                'nJet':4,'nbJet':1,
                                                'suffix':MuSelObjAk4JetsTightExclusiveResolved1b3j.selName,
                                                'is_MC':self.is_MC})

                #----- Resolved selection : 2 Btags -----#
                if "TightResolved2b2j" in jetsel_level:
                    print ("......... Processing Resolved jet (2 btag i.e. bTaggedJets >= 2 & nLightJets >= 2) selection")    
                    ElSelObjAk4JetsTightExclusiveResolved2b2j = makeExclusiveTightResolvedJetComboSelection(self,ElSelObjAk4JetsTight,nbJet=2,copy_sel=True,plot_yield=True)
                    MuSelObjAk4JetsTightExclusiveResolved2b2j = makeExclusiveTightResolvedJetComboSelection(self,MuSelObjAk4JetsTight,nbJet=2,copy_sel=True,plot_yield=True)

                    if not self.args.OnlyYield and "TightResolved2b2j" in jetplot_level:
                        # Cut flow report #
                        cutFlowPlots.append(CutFlowReport(ElSelObjAk4JetsTightExclusiveResolved2b2j.selName,ElSelObjAk4JetsTightExclusiveResolved2b2j.sel))
                        cutFlowPlots.append(CutFlowReport(MuSelObjAk4JetsTightExclusiveResolved2b2j.selName,MuSelObjAk4JetsTightExclusiveResolved2b2j.sel))

                        # Simple reco method to find the jet paits coming from W
                        # https://gitlab.cern.ch/cms-hh-bbww/cms-hh-to-bbww/-/blob/master/Legacy/signal_extraction.md
                        self.lightJetPairs_mu2b2j     = op.combine(self.ak4LightJetsByPt, N=2)
                        self.lambda_chooseWjj_mu2b2j  = lambda dijet: op.abs((dijet[0].p4+dijet[1].p4+MuColl[0].p4+self.corrMET.p4).M() - 
                                                                             (self.HLL.bJetCorrP4(self.ak4BJets[0]) + self.HLL.bJetCorrP4(self.ak4LightJetsByBtagScore[0])).M()) 
                        self.WjjPairs_mu2b2j          = op.sort(self.lightJetPairs_mu2b2j, self.lambda_chooseWjj_mu2b2j)
                        
                        self.lightJetPairs_el2b2j     = op.combine(self.ak4LightJetsByPt, N=2)
                        self.lambda_chooseWjj_el2b2j  = lambda dijet: op.abs((dijet[0].p4+dijet[1].p4+ElColl[0].p4+self.corrMET.p4).M() - 
                                                                             (self.HLL.bJetCorrP4(self.ak4BJets[0]) + self.HLL.bJetCorrP4(self.ak4LightJetsByBtagScore[0])).M())
                        self.WjjPairs_el2b2j          = op.sort(self.lightJetPairs_el2b2j, self.lambda_chooseWjj_el2b2j)
                        

                        ChannelDictList.append({'channel':'El','sel':ElSelObjAk4JetsTightExclusiveResolved2b2j.sel,'lep':ElColl[0],'met':self.corrMET,
                                                'j1':self.ak4BJets[0],'j2':self.ak4BJets[1],'j3':self.WjjPairs_el2b2j[0][0],'j4':self.WjjPairs_el2b2j[0][1],
                                                'nJet':4,'nbJet':2,
                                                'suffix':ElSelObjAk4JetsTightExclusiveResolved2b2j.selName,
                                                'is_MC':self.is_MC})
                        ChannelDictList.append({'channel':'Mu','sel':MuSelObjAk4JetsTightExclusiveResolved2b2j.sel,'lep':MuColl[0],'met':self.corrMET,
                                                'j1':self.ak4BJets[0],'j2':self.ak4BJets[1],'j3':self.WjjPairs_mu2b2j[0][0],'j4':self.WjjPairs_mu2b2j[0][1],
                                                'nJet':4,'nbJet':2,
                                                'suffix':MuSelObjAk4JetsTightExclusiveResolved2b2j.selName,
                                                'is_MC':self.is_MC})

                # Lepton + jet Plots #
                ResolvedBTaggedJetsN = {'objName':'Ak4BJets','objCont':self.ak4BJets,'Nmax':5,'xTitle':'N(Ak4 Bjets)'}
                ResolvedLightJetsN   = {'objName':'Ak4LightJets','objCont':self.ak4LightJetsByBtagScore,'Nmax':10,'xTitle':'N(Ak4 Lightjets)'}
        
                for channelDict in ChannelDictList:
                    # Dilepton #
                    plots.extend(makeSinleptonPlots(**{k:channelDict[k] for k in LeptonKeys}))
                    # Number of jets #
                    plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**ResolvedBTaggedJetsN))
                    plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**ResolvedLightJetsN))
                    # Ak4 Jets #
                    plots.extend(makeAk4JetsPlots(**{k:channelDict[k] for k in JetKeys},HLL=self.HLL))
                    # MET #
                    plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))
                    # High level #
                    plots.extend(makeHighLevelPlotsResolved(**{k:channelDict[k] for k in ResolvedKeys},HLL=self.HLL))

            if any("SemiBoosted" in key for key in jetsel_level):
                print ("......... processing Semi-Boosted Category")
                ChannelDictList   = []
                FatJetKeys     = ['channel','sel','j1','j2','j3','has1fat1slim','has1fat2slim','suffix']
                FatJetsN       = {'objName':'Ak8Jets','objCont':self.ak8Jets,'Nmax':5,'xTitle':'N(Ak8 jets)'}
                FatBJetsN      = {'objName':'Ak8BJets','objCont':self.ak8BJets,'Nmax':5,'xTitle':'N(Ak8 b-jets)'}
                BoostedKeys    = ['channel','sel','met','lep','j1','j2','suffix','bothAreFat']
                if "SemiBoostedHbbWtoJ" in jetsel_level:
                    print ("............ Processing Semi-Boosted category (Htobb:Ak8 + WtoJ)")
                    ElSelObjAk8bJetsHbbBoostedWtoJ = makeSemiBoostedHbbSelection(self,ElSelObjAk8bJets,nNonb=1,copy_sel=True,plot_yield=True)
                    MuSelObjAk8bJetsHbbBoostedWtoJ = makeSemiBoostedHbbSelection(self,MuSelObjAk8bJets,nNonb=1,copy_sel=True,plot_yield=True)
                    if not self.args.OnlyYield and "SemiBoostedHbbWtoJ" in jetplot_level:
                        # Cut flow report #
                        cutFlowPlots.append(CutFlowReport(ElSelObjAk8bJetsHbbBoostedWtoJ.selName,ElSelObjAk8bJetsHbbBoostedWtoJ.sel))
                        cutFlowPlots.append(CutFlowReport(MuSelObjAk8bJetsHbbBoostedWtoJ.selName,MuSelObjAk8bJetsHbbBoostedWtoJ.sel))
 
                        ChannelDictList.append({'channel':'El','sel':ElSelObjAk8bJetsHbbBoostedWtoJ.sel,'lep':ElColl[0],'met':self.corrMET,
                                                'j1':self.ak8BJets[0],'j2':self.ak4JetsCleanedFromAk8b[0],'j3':None,
                                                'has1fat1slim':True,'has1fat2slim':False,'bothAreFat':False,
                                                'suffix':ElSelObjAk8bJetsHbbBoostedWtoJ.selName,
                                                'is_MC':self.is_MC})
                        ChannelDictList.append({'channel':'Mu','sel':MuSelObjAk8bJetsHbbBoostedWtoJ.sel,'lep':MuColl[0],'met':self.corrMET,
                                                'j1':self.ak8BJets[0],'j2':self.ak4JetsCleanedFromAk8b[0],'j3':None,
                                                'has1fat1slim':True,'has1fat2slim':False,'bothAreFat':False,
                                                'suffix':MuSelObjAk8bJetsHbbBoostedWtoJ.selName,
                                                'is_MC':self.is_MC})

                if "SemiBoostedHbbWtoJJ" in jetsel_level:
                    print ("............ Processing Semi-Boosted category (Htobb:Ak8 + WtoJJ)")
                    ElSelObjAk8bJetsHbbBoostedWtoJJ = makeSemiBoostedHbbSelection(self,ElSelObjAk8bJets,nNonb=2,copy_sel=True,plot_yield=True)
                    MuSelObjAk8bJetsHbbBoostedWtoJJ = makeSemiBoostedHbbSelection(self,MuSelObjAk8bJets,nNonb=2,copy_sel=True,plot_yield=True)
                    if not self.args.OnlyYield and "SemiBoostedHbbWtoJJ" in jetplot_level:

                        self.nonbJetPairs_muSB        = op.combine(self.ak4JetsCleanedFromAk8b, N=2)
                        self.lambda_chooseWjj_muSB    = lambda dijet: op.abs((dijet[0].p4+dijet[1].p4+MuColl[0].p4+self.corrMET.p4).M() - self.ak8BJets[0].p4.M()) 
                        self.WjjPairs_muSB            = op.sort(self.nonbJetPairs_muSB, self.lambda_chooseWjj_muSB)
                        self.nonbJetPairs_elSB        = op.combine(self.ak4JetsCleanedFromAk8b, N=2)
                        self.lambda_chooseWjj_elSB    = lambda dijet: op.abs((dijet[0].p4+dijet[1].p4+ElColl[0].p4+self.corrMET.p4).M() - self.ak8BJets[0].p4.M()) 
                        self.WjjPairs_elSB            = op.sort(self.nonbJetPairs_elSB, self.lambda_chooseWjj_elSB)
                        

                        # cut flow report
                        cutFlowPlots.append(CutFlowReport(ElSelObjAk8bJetsHbbBoostedWtoJJ.selName,ElSelObjAk8bJetsHbbBoostedWtoJJ.sel))
                        cutFlowPlots.append(CutFlowReport(MuSelObjAk8bJetsHbbBoostedWtoJJ.selName,MuSelObjAk8bJetsHbbBoostedWtoJJ.sel))


                        ChannelDictList.append({'channel':'El','sel':ElSelObjAk8bJetsHbbBoostedWtoJJ.sel,'lep':ElColl[0],'met':self.corrMET,
                                                'j1':self.ak8BJets[0],'j2':self.WjjPairs_elSB[0][0],'j3':self.WjjPairs_elSB[0][1],
                                                'has1fat1slim':False,'has1fat2slim':True,'bothAreFat':False,
                                                'suffix':ElSelObjAk8bJetsHbbBoostedWtoJJ.selName,
                                                'is_MC':self.is_MC})
                        ChannelDictList.append({'channel':'Mu','sel':MuSelObjAk8bJetsHbbBoostedWtoJJ.sel,'lep':MuColl[0],'met':self.corrMET,
                                                'j1':self.ak8BJets[0],'j2':self.WjjPairs_muSB[0][0],'j3':self.WjjPairs_muSB[0][1],
                                                'has1fat1slim':False,'has1fat2slim':True,'bothAreFat':False,
                                                'suffix':MuSelObjAk8bJetsHbbBoostedWtoJJ.selName,
                                                'is_MC':self.is_MC})
                for channelDict in ChannelDictList:
                    # Dilepton #
                    plots.extend(makeSinleptonPlots(**{k:channelDict[k] for k in LeptonKeys}))
                    # Number of jets #
                    plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**FatJetsN))
                    plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**FatBJetsN))       
                    # Ak8 Jets #
                    plots.extend(makeSingleLeptonAk8JetsPlots(**{k:channelDict[k] for k in FatJetKeys}))
                    # MET #
                    plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))
                    # HighLevel #
                    plots.extend(makeHighLevelPlotsBoosted(**{k:channelDict[k] for k in BoostedKeys},HLL=self.HLL))
            
            if "Boosted" in jetsel_level:
                print ("......... Processing Boosted category selection")        
                ChannelDictList= []
                FatJetKeys     = ['channel','sel','j1','j2','has2fat','suffix']
                FatJetsN       = {'objName':'Ak8Jets','objCont':self.ak8Jets,'Nmax':5,'xTitle':'N(Ak8 jets)'}
                FatBJetsN      = {'objName':'Ak8BJets','objCont':self.ak8BJets,'Nmax':5,'xTitle':'N(Ak8 b-jets)'}
                ElSelObjAk8JetsInclusiveBoosted = makeInclusiveBoostedSelection(self,ElSelObjAk8Jets,copy_sel=True,plot_yield=True)
                MuSelObjAk8JetsInclusiveBoosted = makeInclusiveBoostedSelection(self,MuSelObjAk8Jets,copy_sel=True,plot_yield=True)
                if not self.args.OnlyYield and "Boosted" in jetplot_level:
                    # Cut flow report #
                    cutFlowlots.append(CutFlowReport(ElSelObjAk8JetsInclusiveBoosted.selName,ElSelObjAk8JetsInclusiveBoosted.sel))
                    cutFlowlots.append(CutFlowReport(MuSelObjAk8JetsInclusiveBoosted.selName,MuSelObjAk8JetsInclusiveBoosted.sel))
                    
                    ChannelDictList.append({'channel':'El','sel':ElSelObjAk8JetsInclusiveBoosted.sel,'sinlepton':ElColl[0],
                                            'j1':self.ak8BJets[0],'j2':self.ak8nonBJets[0],'has2fat':True,
                                            'suffix':ElSelObjAk8JetsInclusiveBoosted.selName,
                                            'is_MC':self.is_MC})
                    ChannelDictList.append({'channel':'Mu','sel':MuSelObjAk8JetsInclusiveBoosted.sel,'sinlepton':MuColl[0],
                                            'j1':self.ak8Jets[0],'j2':self.ak8nonBJets[0],'has2fat':True,
                                            'suffix':MuSelObjAk8JetsInclusiveBoosted.selName,
                                            'is_MC':self.is_MC})
                    
                for channelDict in ChannelDictList:
                    # Dilepton #
                    plots.extend(makeSinleptonPlots(**{k:channelDict[k] for k in LeptonKeys}))
                    # Number of jets #
                    plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**FatJetsN))
                    plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**FatBJetsN))       
                    # Ak8 Jets #
                    plots.extend(makeSingleLeptonAk8JetsPlots(**{k:channelDict[k] for k in FatJetKeys}))
                    # MET #
                    plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))

            
            #----- Machine Learning plots -----#
            channelDictList = []
            if not self.args.OnlyYield:
                '''
                if "LooseResolved2b1j" in jetplot_level:
                    if self.args.Classifier == 'BDT-SM':
                        channelDictList.append({'channel':'El','lep':ElColl[0],
                                                'jets':self.ak4Jets,'bJets':self.ak4BJets,'lJets':self.ak4LightJetsByPt,
                                                'j1':self.ak4BJets[0],'j2':self.ak4BJets[1],'j3':self.ak4LightJetsByPt[0],'j4':None,
                                                'sel':ElSelObjAk4JetsLooseExclusiveResolved2b1j.sel,
                                                'suffix':ElSelObjAk4JetsLooseExclusiveResolved2b1j.selName+'_simple_SMBDT',
                                                'model': BDT_missRecoSM_simple})
                    if self.args.Classifier == 'BDT-Rad900':
                        channelDictList.append({'channel':'El','lep':ElColl[0],
                                                'jets':self.ak4Jets,'bJets':self.ak4BJets,'lJets':self.ak4LightJetsByPt,
                                                'j1':self.ak4BJets[0],'j2':self.ak4BJets[1],'j3':self.ak4LightJetsByPt[0],'j4':None,
                                                'sel':ElSelObjAk4JetsLooseExclusiveResolved2b1j.sel,
                                                'suffix':ElSelObjAk4JetsLooseExclusiveResolved2b1j.selName+'_simple_900RBDT',
                                                'model': BDT_missReco900R_simple})
                '''
                if "TightResolved2b2j" in jetplot_level:
                    if self.args.Classifier == 'BDT-SM':
                        channelDictList.append({'channel':'El','fakeLepColl':self.electronsFakeSel,'lep':ElColl[0],'met':self.corrMET,
                                                'jets':self.ak4Jets,'bJets':self.ak4BJets,'lJets':self.ak4LightJetsByPt,
                                                'j1':self.ak4BJets[0],'j2':self.ak4BJets[1],'j3':self.WjjPairs_el2b2j[0][0],'j4':self.WjjPairs_el2b2j[0][1],
                                                'sel':ElSelObjAk4JetsTightExclusiveResolved2b2j.sel,
                                                'suffix':ElSelObjAk4JetsTightExclusiveResolved2b2j.selName+'_simple_SMBDT',
                                                'model_even':BDT_fullRecoSM_simple_even, 'model_odd':BDT_fullRecoSM_simple_odd})
                        channelDictList.append({'channel':'Mu','fakeLepColl':self.muonsFakeSel,'lep':MuColl[0],'met':self.corrMET,
                                                'jets':self.ak4Jets,'bJets':self.ak4BJets,'lJets':self.ak4LightJetsByPt,
                                                'j1':self.ak4BJets[0],'j2':self.ak4BJets[1],'j3':self.WjjPairs_el2b2j[0][0],'j4':self.WjjPairs_el2b2j[0][1],
                                                'sel':MuSelObjAk4JetsTightExclusiveResolved2b2j.sel,
                                                'suffix':MuSelObjAk4JetsTightExclusiveResolved2b2j.selName+'_simple_SMBDT',
                                                'model_even':BDT_fullRecoSM_simple_even, 'model_odd':BDT_fullRecoSM_simple_odd})

                    elif self.args.Classifier == 'BDT-Rad900':
                        channelDictList.append({'channel':'El','fakeLepColl':self.electronsFakeSel,'lep':ElColl[0],'met':self.corrMET,
                                                'jets':self.ak4Jets,'bJets':self.ak4BJets,'lJets':self.ak4LightJetsByPt,
                                                'j1':self.ak4BJets[0],'j2':self.ak4BJets[1],'j3':self.WjjPairs_mu2b2j[0][0],'j4':self.WjjPairs_mu2b2j[0][1],
                                                'sel':ElSelObjAk4JetsTightExclusiveResolved2b2j.sel,
                                                'suffix':ElSelObjAk4JetsTightExclusiveResolved2b2j.selName+'_simple_900RBDT',
                                                'model_even':BDT_fullReco900R_simple_even, 'model_odd':BDT_fullReco900R_simple_odd})
                        channelDictList.append({'channel':'Mu','fakeLepColl':self.muonsFakeSel,'lep':MuColl[0],'met':self.corrMET,
                                                'jets':self.ak4Jets,'bJets':self.ak4BJets,'lJets':self.ak4LightJetsByPt,
                                                'j1':self.ak4BJets[0],'j2':self.ak4BJets[1],'j3':self.WjjPairs_mu2b2j[0][0],'j4':self.WjjPairs_mu2b2j[0][1],
                                                'sel':MuSelObjAk4JetsTightExclusiveResolved2b2j.sel,
                                                'suffix':MuSelObjAk4JetsTightExclusiveResolved2b2j.selName+'_simple_900RBDT',
                                                'model_even':BDT_fullReco900R_simple_even, 'model_odd':BDT_fullReco900R_simple_odd})

            for channelDict in channelDictList:
                plots.extend(makeSingleLeptonMachineLearningPlotsBDT(**channelDict,event=t.event,HLL=self.HLL))
            
        #----- Add the Yield plots -----#
        plots.extend(self.yieldPlots.returnPlots())
        #plots.extend(cutFlowPlots)
        return plots

