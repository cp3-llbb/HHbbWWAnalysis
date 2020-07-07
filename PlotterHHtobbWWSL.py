import os
import sys
from copy import copy

from bamboo.analysismodules import HistogramsModule
from bamboo import treefunctions as op
from bamboo.plots import Plot, EquidistantBinning, SummedPlot

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)))) # Add scripts in this directory
from BaseHHtobbWW import BaseNanoHHtobbWW
from plotDefSL import *
from selectionDefSL import *


#===============================================================================================#
#                                       PlotterHHtobbWW                                         #
#===============================================================================================#
class PlotterNanoHHtobbWWSL(BaseNanoHHtobbWW,HistogramsModule):
    """ Plotter module: HH->bbW(->e/µ nu)W(->e/µ nu) histograms from NanoAOD """
    def __init__(self, args):
        super(PlotterNanoHHtobbWWSL, self).__init__(args)

    def definePlots(self, t, noSel, sample=None, sampleCfg=None): 
        noSel = super(PlotterNanoHHtobbWWSL,self).prepareObjects(t, noSel, sample, sampleCfg, 'SL')
        plots = []
        era = sampleCfg['era']
        self.yieldPlots = makeYieldPlots()

        #----- Ratio reweighting variables (before lepton and jet selection) -----#
        if self.args.BtagReweightingOff or self.args.BtagReweightingOn:
            plots.append(objectsNumberPlot(channel="NoChannel",suffix='NoSelection',sel=noSel,objCont=self.ak4Jets,objName='Ak4Jets',Nmax=15,xTitle='N(Ak4 jets)'))
            return plots

        #----- Singleleptons -----#
        selObjectDict = makeSingleLeptonSelection(self,noSel,plot_yield=True)
        # selObjectDict : keys -> level (str)
        #                 values -> [ElEl,MuMu,ElMu] x Selection object
        # Select the jets selections that will be done depending on user input #
        jet_level = ["Ak4","Ak8",
                     "LooseResolved0b3j","LooseResolved1b2j","LooseResolved2b1j",
                     "TightResolved0b4j","TightResolved1b3j","TightResolved2b2j",
                     "SemiBoostedHbb","SemiBoostedWjj",
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
                ElColl = self.leadElectronTightSel
                MuColl = self.leadMuonTightSel
            elif selectionType == "FakeExtrapolation":
                ElColl = self.leadElectronFakeExtrapolationSel
                MuColl = self.leadMuonFakeExtrapolationSel

            #----- Separate selections ------#
            ElSelObj = selectionList[0]
            MuSelObj = selectionList[1]

            if not self.args.OnlyYield:
                #----- Trigger plots -----#
                plots.extend(triggerPlots(sel=ElSelObj.sel, triggerDict=self.triggersPerPrimaryDataset, suffix=ElSelObj.selName, channel="El"))
                plots.extend(triggerPlots(sel=MuSelObj.sel, triggerDict=self.triggersPerPrimaryDataset, suffix=MuSelObj.selName, channel="Mu"))

                #----- Lepton plots -----#
                # Singlelepton channel plots #
                plots.extend(channelPlot(sel=ElSelObj.sel, SinlepEl=ElColl, SinlepMu=MuColl, suffix=ElSelObj.selName, channel = "El"))
                plots.extend(channelPlot(sel=MuSelObj.sel, SinlepEl=ElColl, SinlepMu=MuColl, suffix=MuSelObj.selName, channel = "Mu"))
        

            #----- Ak4 jets selection -----#
            LeptonKeys  = ['channel','sel','sinlepton','suffix','is_MC']
            JetKeys     = ['channel','sel','j1','j2','j3','j4','suffix','nJet','nbJet','is_MC']
            commonItems = ['channel','sel','suffix']

            if "Ak4" in jetsel_level:
                print("... Processing Ak4Jets Selection for Resolved category")
                ChannelDictList = []
                JetsN    = {'objName':'Ak4Jets','objCont':self.ak4Jets,'Nmax':10,'xTitle':'N(Ak4 jets)'}
                FatJetsN = {'objName':'Ak8Jets','objCont':self.ak8Jets,'Nmax':5,'xTitle':'N(Ak8 jets)'}

                if any("LooseResolved" in key for key in jetsel_level):
                    print ("...... Processing Ak4 jet selection Loose (nAk4Jets >= 3)")
                    ElSelObjAk4JetsLoose = makeAk4JetSelection(self,ElSelObj,nJet=3,copy_sel=True,plot_yield=True)
                    MuSelObjAk4JetsLoose = makeAk4JetSelection(self,MuSelObj,nJet=3,copy_sel=True,plot_yield=True)
                    if not self.args.OnlyYield:
                        # cutFlow Report #
                        plots.append(CutFlowReport(ElSelObjAk4JetsLoose.selName, ElSelObjAk4JetsLoose.sel))
                        plots.append(CutFlowReport(MuSelObjAk4JetsLoose.selName, MuSelObjAk4JetsLoose.sel))

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
                    ElSelObjAk4JetsTight = makeAk4JetSelection(self,ElSelObj,nJet=4,copy_sel=True,plot_yield=True)
                    MuSelObjAk4JetsTight = makeAk4JetSelection(self,MuSelObj,nJet=4,copy_sel=True,plot_yield=True)

                    # Jet and lepton plots #
                    if not self.args.OnlyYield:
                        # cutFlow Report #
                        plots.append(CutFlowReport(ElSelObjAk4JetsTight.selName, ElSelObjAk4JetsTight.sel))
                        plots.append(CutFlowReport(MuSelObjAk4JetsTight.selName, MuSelObjAk4JetsTight.sel))
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
                    plots.extend(makeAk4JetsPlots(**{k:channelDict[k] for k in JetKeys}))
                    # MET #
                    plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))

            ##### Ak8 jets selection #####
            if "Ak8" in jetsel_level:
                print ("...... Processing Ak8 jet selection for SemiBoosted & Boosted Category")
                ElSelObjAk8Jets = makeAk8JetSelection(self,ElSelObj,copy_sel=True,plot_yield=True)
                MuSelObjAk8Jets = makeAk8JetSelection(self,MuSelObj,copy_sel=True,plot_yield=True)

                FatJetKeys = ['channel','sel','j1','j2','has1fat','suffix']
                FatJetsN   = {'objName':'Ak8Jets','objCont':self.ak8Jets,'Nmax':5,'xTitle':'N(Ak8 jets)'}
                SlimJetsN  = {'objName':'Ak4Jets','objCont':self.ak4Jets,'Nmax':10,'xTitle':'N(Ak4 jets)'}

                # Fatjets plots #
                ChannelDictList = []
                if not self.args.OnlyYield:
                    # cutFlow Report #
                    plots.append(CutFlowReport(ElSelObjAk8Jets.selName, ElSelObjAk8Jets.sel))
                    plots.append(CutFlowReport(MuSelObjAk8Jets.selName, MuSelObjAk8Jets.sel))
                    if "Ak8" in jetplot_level:
                        ChannelDictList.append({'channel':'El','sel':ElSelObjAk8Jets.sel,'sinlepton':ElColl[0],
                                                'j1':self.ak8Jets[0],'j2':self.ak8Jets[0],'has1fat':True,
                                                'suffix':ElSelObjAk8Jets.selName,'is_MC':self.is_MC})
                        ChannelDictList.append({'channel':'Mu','sel':MuSelObjAk8Jets.sel,'sinlepton':MuColl[0],
                                                'fatjet':self.ak8Jets[0],'j2':self.ak8Jets[0],'has1fat':True,
                                                'suffix':MuSelObjAk8Jets.selName,'is_MC':self.is_MC})

                for channelDict in ChannelDictList:
                    # Dilepton #
                    plots.extend(makeSinleptonPlots(**{k:channelDict[k] for k in LeptonKeys}))
                    # Number of jets #
                    plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**FatJetsN))
                    plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**SlimJetsN))
                    # Ak8 Jets #
                    plots.extend(makeAk8JetsPlots(**{k:channelDict[k] for k in FatJetKeys}))
                    # MET #
                    plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))

                         
            #-----------------------------|||||||||||||||| Resolved selection ||||||||||||||-----------------------------------#
            if any("Resolved" in item for item in jetsel_level):
                ChannelDictList = []
                # Resolved Selection (Loose) #
                #----- Resolved selection : 0 Btag -----#
                if "LooseResolved0b3j" in jetsel_level:
                    print ("......... Processing Loose Resolved jet (0 btag i.e. bTaggedJets = 0 & nLightJets >= 3) selection")
                    ElSelObjAk4JetsLooseExclusiveResolved0b3j = makeExclusiveLooseResolvedJetComboSelection(self,ElSelObjAk4JetsLoose,nbJet=0,copy_sel=True,plot_yield=True)
                    MuSelObjAk4JetsLooseExclusiveResolved0b3j = makeExclusiveLooseResolvedJetComboSelection(self,MuSelObjAk4JetsLoose,nbJet=0,copy_sel=True,plot_yield=True)

                    if not self.args.OnlyYield:
                        # Cut flow report #
                        plots.append(CutFlowReport(ElSelObjAk4JetsLooseExclusiveResolved0b3j.selName,ElSelObjAk4JetsLooseExclusiveResolved0b3j.sel))
                        plots.append(CutFlowReport(MuSelObjAk4JetsLooseExclusiveResolved0b3j.selName,MuSelObjAk4JetsLooseExclusiveResolved0b3j.sel))
                        if "LooseResolved0b3j" in jetplot_level:
                            ChannelDictList.append({'channel':'El','sel':ElSelObjAk4JetsLooseExclusiveResolved0b3j.sel,'sinlepton':ElColl[0],
                                                    'j1':self.ak4LightJetsByPt[0],'j2':self.ak4LightJetsByPt[1],'j3':self.ak4LightJetsByPt[2],'j4':None,
                                                    'nJet':3,'nbJet':0,
                                                    'suffix':ElSelObjAk4JetsLooseExclusiveResolved0b3j.selName,
                                                    'is_MC':self.is_MC})
                            ChannelDictList.append({'channel':'Mu','sel':MuSelObjAk4JetsLooseExclusiveResolved0b3j.sel,'sinlepton':MuColl[0],
                                                    'j1':self.ak4LightJetsByPt[0],'j2':self.ak4LightJetsByPt[1],'j3':self.ak4LightJetsByPt[2],'j4':None,
                                                    'nJet':3,'nbJet':0,
                                                    'suffix':MuSelObjAk4JetsLooseExclusiveResolved0b3j.selName,
                                                    'is_MC':self.is_MC})


                #----  Resolved selection : 1 Btag  -----#
                if "LooseResolved1b2j" in jetsel_level:
                    print ("......... Processing Resolved jet (1 btag i.e. bTaggedJets = 1 & nLightJets >= 2) selection")
                    ElSelObjAk4JetsLooseExclusiveResolved1b2j = makeExclusiveLooseResolvedJetComboSelection(self,ElSelObjAk4JetsLoose,nbJet=1,copy_sel=True,plot_yield=True)
                    MuSelObjAk4JetsLooseExclusiveResolved1b2j = makeExclusiveLooseResolvedJetComboSelection(self,MuSelObjAk4JetsLoose,nbJet=1,copy_sel=True,plot_yield=True)

                    if not self.args.OnlyYield:
                        # Cut flow report #
                        plots.append(CutFlowReport(ElSelObjAk4JetsLooseExclusiveResolved1b2j.selName,ElSelObjAk4JetsLooseExclusiveResolved1b2j.sel))
                        plots.append(CutFlowReport(MuSelObjAk4JetsLooseExclusiveResolved1b2j.selName,MuSelObjAk4JetsLooseExclusiveResolved1b2j.sel))
                        if "LooseResolved1b2j" in jetplot_level:
                            ChannelDictList.append({'channel':'El','sel':ElSelObjAk4JetsLooseExclusiveResolved1b2j.sel,'sinlepton':ElColl[0],
                                                    'j1':self.ak4BJets[0],'j2':self.ak4LightJets[0],'j3':self.ak4LightJets[1],'j4':None,
                                                    'nJet':3,'nbJet':1,
                                                    'suffix':ElSelObjAk4JetsLooseExclusiveResolved1b2j.selName,
                                                    'is_MC':self.is_MC})
                            ChannelDictList.append({'channel':'Mu','sel':MuSelObjAk4JetsLooseExclusiveResolved1b2j.sel,'sinlepton':MuColl[0],
                                                    'j1':self.ak4BJets[0],'j2':self.ak4LightJets[0],'j3':self.ak4LightJets[1],'j4':None,
                                                    'nJet':3,'nbJet':1,
                                                    'suffix':MuSelObjAk4JetsLooseExclusiveResolved1b2j.selName,
                                                    'is_MC':self.is_MC})


                #----- Resolved selection : 2 Btags -----#
                if "LooseResolved2b1j" in jetsel_level:
                    print ("......... Processing Resolved jet (2 btag i.e. bTaggedJets = 2 & nLightJets >= 1) selection")
                    ElSelObjAk4JetsLooseExclusiveResolved2b1j = makeExclusiveLooseResolvedJetComboSelection(self,ElSelObjAk4JetsLoose,nbJet=2,copy_sel=True,plot_yield=True)
                    MuSelObjAk4JetsLooseExclusiveResolved2b1j = makeExclusiveLooseResolvedJetComboSelection(self,MuSelObjAk4JetsLoose,nbJet=2,copy_sel=True,plot_yield=True)

                    if not self.args.OnlyYield:
                        # Cut flow report #
                        plots.append(CutFlowReport(ElSelObjAk4JetsLooseExclusiveResolved2b1j.selName,ElSelObjAk4JetsLooseExclusiveResolved2b1j.sel))
                        plots.append(CutFlowReport(MuSelObjAk4JetsLooseExclusiveResolved2b1j.selName,MuSelObjAk4JetsLooseExclusiveResolved2b1j.sel))
                        if "LooseResolved2b1j" in jetplot_level:
                            ChannelDictList.append({'channel':'El','sel':ElSelObjAk4JetsLooseExclusiveResolved2b1j.sel,'sinlepton':ElColl[0],
                                                    'j1':self.ak4BJets[0],'j2':self.ak4BJets[1],'j3':self.ak4LightJetsByPt[0],'j4':None,
                                                    'nJet':3,'nbJet':2,
                                                    'suffix':ElSelObjAk4JetsLooseExclusiveResolved2b1j.selName,
                                                    'is_MC':self.is_MC})
                            ChannelDictList.append({'channel':'Mu','sel':MuSelObjAk4JetsLooseExclusiveResolved2b1j.sel,'sinlepton':MuColl[0],
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
                        ChannelDictList.append({'channel':'El','sel':ElSelObjAk4JetsTightExclusiveResolved0b4j.sel,'sinlepton':ElColl[0],
                                                'j1':self.ak4LightJetsByPt[0],'j2':self.ak4LightJetsByPt[1],'j3':self.ak4LightJetsByPt[2],'j4':self.ak4LightJetsByPt[3],
                                                'nJet':4,'nbJet':0,
                                                'suffix':ElSelObjAk4JetsTightExclusiveResolved0b4j.selName,
                                                'is_MC':self.is_MC})
                        ChannelDictList.append({'channel':'Mu','sel':MuSelObjAk4JetsTightExclusiveResolved0b4j.sel,'sinlepton':MuColl[0],
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
                        ChannelDictList.append({'channel':'El','sel':ElSelObjAk4JetsTightExclusiveResolved1b3j.sel,'sinlepton':ElColl[0],
                                                'j1':self.ak4BJets[0],'j2':self.ak4LightJets[0],'j3':self.ak4LightJets[1],'j4':self.ak4LightJets[2],
                                                'nJet':4,'nbJet':1,
                                                'suffix':ElSelObjAk4JetsTightExclusiveResolved1b3j.selName,
                                                'is_MC':self.is_MC})
                        ChannelDictList.append({'channel':'Mu','sel':MuSelObjAk4JetsTightExclusiveResolved1b3j.sel,'sinlepton':MuColl[0],
                                                'j1':self.ak4BJets[0],'j2':self.ak4LightJets[0],'j3':self.ak4LightJets[1],'j4':self.ak4LightJets[2],
                                                'nJet':4,'nbJet':1,
                                                'suffix':MuSelObjAk4JetsTightExclusiveResolved1b3j.selName,
                                                'is_MC':self.is_MC})


                #----- Resolved selection : 2 Btags -----#
                if "TightResolved2b2j" in jetsel_level:
                    print ("......... Processing Resolved jet (2 btag i.e. bTaggedJets = 2 & nLightJets >= 2) selection")    
                    ElSelObjAk4JetsTightExclusiveResolved2b2j = makeExclusiveTightResolvedJetComboSelection(self,ElSelObjAk4JetsTight,nbJet=2,copy_sel=True,plot_yield=True)
                    MuSelObjAk4JetsTightExclusiveResolved2b2j = makeExclusiveTightResolvedJetComboSelection(self,MuSelObjAk4JetsTight,nbJet=2,copy_sel=True,plot_yield=True)

                    if not self.args.OnlyYield:
                        # Cut flow report #
                        plots.append(CutFlowReport(ElSelObjAk4JetsTightExclusiveResolved2b2j.selName,ElSelObjAk4JetsTightExclusiveResolved2b2j.sel))
                        plots.append(CutFlowReport(MuSelObjAk4JetsTightExclusiveResolved2b2j.selName,MuSelObjAk4JetsTightExclusiveResolved2b2j.sel))
                        if "TightResolved2b2j" in jetplot_level:
                            ChannelDictList.append({'channel':'El','sel':ElSelObjAk4JetsTightExclusiveResolved2b2j.sel,'sinlepton':ElColl[0],
                                                    'j1':self.ak4BJets[0],'j2':self.ak4BJets[1],'j3':self.ak4LightJetsByPt[0],'j4':self.ak4LightJetsByPt[1],
                                                    'nJet':4,'nbJet':2,
                                                    'suffix':ElSelObjAk4JetsTightExclusiveResolved2b2j.selName,
                                                    'is_MC':self.is_MC})
                            ChannelDictList.append({'channel':'Mu','sel':MuSelObjAk4JetsTightExclusiveResolved2b2j.sel,'sinlepton':MuColl[0],
                                                    'j1':self.ak4BJets[0],'j2':self.ak4BJets[1],'j3':self.ak4LightJetsByPt[0],'j4':self.ak4LightJetsByPt[1],
                                                    'nJet':4,'nbJet':2,
                                                    'suffix':MuSelObjAk4JetsTightExclusiveResolved2b2j.selName,
                                                    'is_MC':self.is_MC})

                # Lepton + jet Plots #
                ResolvedBTaggedJetsN = {'objName':'Ak4BJets','objCont':self.ak4BJets,'Nmax':5,'xTitle':'N(Ak4 Bjets)'}
                ResolvedLightJetsN   = {'objName':'Ak4LightJets','objCont':self.ak4LightJets,'Nmax':10,'xTitle':'N(Ak4 Lightjets)'}
        
                for channelDict in ChannelDictList:
                    # Dilepton #
                    plots.extend(makeSinleptonPlots(**{k:channelDict[k] for k in LeptonKeys}))
                    # Number of jets #
                    plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**ResolvedBTaggedJetsN))
                    plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**ResolvedLightJetsN))
                    # Ak4 Jets #
                    plots.extend(makeAk4JetsPlots(**{k:channelDict[k] for k in JetKeys}))
                    # MET #
                    plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))


            if any("SemiBoosted" in key for key in jetsel_level):
                print ("......... processing Semi-Boosted Category")
                ChannelDictList= []
                FatJetKeys     = ['channel','sel','j1','j2','has1fat1slim','suffix']
                FatJetsN       = {'objName':'Ak8Jets','objCont':self.ak8Jets,'Nmax':5,'xTitle':'N(Ak8 jets)'}
                FatBJetsN      = {'objName':'Ak8BJets','objCont':self.ak8BJets,'Nmax':5,'xTitle':'N(Ak8 b-jets)'}
                if "SemiBoostedHbb" in jetsel_level:
                    print ("............ Processing Semi-Boosted category (Htobb:Ak8 + Wtojj:Ak4 >= 1)")
                    ElSelObjAk8JetsHbbBoosted = makeSemiBoostedHbbSelection(self,ElSelObjAk8Jets,copy_sel=True,plot_yield=True)
                    MuSelObjAk8JetsHbbBoosted = makeSemiBoostedHbbSelection(self,MuSelObjAk8Jets,copy_sel=True,plot_yield=True)
                    if not self.args.OnlyYield:
                        # Cut flow report #
                        plots.append(CutFlowReport(ElSelObjAk8JetsHbbBoosted.selName,ElSelObjAk8JetsHbbBoosted.sel))
                        plots.append(CutFlowReport(MuSelObjAk8JetsHbbBoosted.selName,MuSelObjAk8JetsHbbBoosted.sel))
                        if"SemiBoostedHbb" in jetplot_level:
                            ChannelDictList.append({'channel':'El','sel':ElSelObjAk8JetsHbbBoosted.sel,'sinlepton':ElColl[0],
                                                    'j1':self.ak8BJets[0],'j2':self.ak4JetsCleanedFromAk8b[0],'has1fat1slim':True,
                                                    'suffix':ElSelObjAk8JetsHbbBoosted.selName,
                                                    'is_MC':self.is_MC})
                            ChannelDictList.append({'channel':'Mu','sel':MuSelObjAk8JetsHbbBoosted.sel,'sinlepton':MuColl[0],
                                                    'j1':self.ak8Jets[0],'j2':self.ak4JetsCleanedFromAk8b[0],'has1fat1slim':True,
                                                    'suffix':MuSelObjAk8JetsHbbBoosted.selName,
                                                    'is_MC':self.is_MC})

                for channelDict in ChannelDictList:
                    # Dilepton #
                    plots.extend(makeSinleptonPlots(**{k:channelDict[k] for k in LeptonKeys}))
                    # Number of jets #
                    plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**FatJetsN))
                    plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**FatBJetsN))       
                    # Ak8 Jets #
                    plots.extend(makeAk8JetsPlots(**{k:channelDict[k] for k in FatJetKeys}))
                    # MET #
                    plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))


            if "Boosted" in jetsel_level:
                print ("......... Processing Boosted category selection")        
                ChannelDictList= []
                FatJetKeys     = ['channel','sel','j1','j2','has2fat','suffix']
                FatJetsN       = {'objName':'Ak8Jets','objCont':self.ak8Jets,'Nmax':5,'xTitle':'N(Ak8 jets)'}
                FatBJetsN      = {'objName':'Ak8BJets','objCont':self.ak8BJets,'Nmax':5,'xTitle':'N(Ak8 b-jets)'}
                ElSelObjAk8JetsInclusiveBoosted = makeInclusiveBoostedSelection(self,ElSelObjAk8Jets,copy_sel=True,plot_yield=True)
                MuSelObjAk8JetsInclusiveBoosted = makeInclusiveBoostedSelection(self,MuSelObjAk8Jets,copy_sel=True,plot_yield=True)
                if not self.args.OnlyYield:
                    # Cut flow report #
                    plots.append(CutFlowReport(ElSelObjAk8JetsInclusiveBoosted.selName,ElSelObjAk8JetsInclusiveBoosted.sel))
                    plots.append(CutFlowReport(MuSelObjAk8JetsInclusiveBoosted.selName,MuSelObjAk8JetsInclusiveBoosted.sel))
                    if "Boosted" in jetplot_level:
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
                    plots.extend(makeAk8JetsPlots(**{k:channelDict[k] for k in FatJetKeys}))
                    # MET #
                    plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))


            #----- High-level combinations -----#
            if not self.args.OnlyYield:
                ResolvedKeys = ['channel','sel','met','lep','j1','j2','j3','j4','suffix','nJet','nbJet']
                BoostedKeys  = ['channel','sel','met','lep','j1','j2','suffix','bothAreFat']
                if any("Resolved" in item for item in jetplot_level): 
                    ChannelDictList = []

                    # LooseResolved
                    # Resolved No Btag #
                    if "LooseResolved0b3j" in jetplot_level:
                        ChannelDictList.append({'channel': 'El','met': self.corrMET,'lep':ElColl[0],
                                                'j1':self.ak4LightJets[0],'j2':self.ak4LightJets[1],'j3':self.ak4LightJets[2],'j4':None,
                                                'nJet':3,'nbJet':0,
                                                'sel':ElSelObjAk4JetsLooseExclusiveResolved0b3j.sel,
                                                'suffix':ElSelObjAk4JetsLooseExclusiveResolved0b3j.selName})
                        ChannelDictList.append({'channel': 'Mu','met': self.corrMET,'lep':MuColl[0],
                                                'j1':self.ak4LightJets[0],'j2':self.ak4LightJets[1],'j3':self.ak4LightJets[2],'j4':None,
                                                'nJet':3,'nbJet':0,
                                                'sel':MuSelObjAk4JetsLooseExclusiveResolved0b3j.sel,
                                                'suffix':MuSelObjAk4JetsLooseExclusiveResolved0b3j.selName})
                    # Resolved One Btag Two jet #
                    if "LooseResolved1b2j" in jetplot_level:
                        ChannelDictList.append({'channel': 'El','met': self.corrMET,'lep':ElColl[0],
                                                'j1':self.ak4BJets[0],'j2':self.ak4LightJets[0],'j3':self.ak4LightJets[1],'j4':None,
                                                'nJet':3,'nbJet':1,
                                                'sel':ElSelObjAk4JetsLooseExclusiveResolved1b2j.sel,
                                                'suffix':ElSelObjAk4JetsLooseExclusiveResolved1b2j.selName})
                        ChannelDictList.append({'channel': 'Mu','met': self.corrMET,'lep':MuColl[0],
                                                'j1':self.ak4BJets[0],'j2':self.ak4LightJets[0],'j3':self.ak4LightJets[1],'j4':None,
                                                'nJet':3,'nbJet':1,
                                                'sel':MuSelObjAk4JetsLooseExclusiveResolved1b2j.sel,
                                                'suffix':MuSelObjAk4JetsLooseExclusiveResolved1b2j.selName})
                    # Resolved Two Btag One jet #      
                    if "LooseResolved2b1j" in jetplot_level:
                        ChannelDictList.append({'channel': 'El','met': self.corrMET,'lep':ElColl[0],
                                                'j1':self.ak4BJets[0],'j2':self.ak4BJets[1],'j3':self.ak4LightJets[0],'j4':None,
                                                'nJet':3,'nbJet':2,
                                                'sel':ElSelObjAk4JetsLooseExclusiveResolved2b1j.sel,
                                                'suffix':ElSelObjAk4JetsLooseExclusiveResolved2b1j.selName})
                        ChannelDictList.append({'channel': 'Mu','met': self.corrMET,'lep':MuColl[0],
                                                'j1':self.ak4BJets[0],'j2':self.ak4BJets[1],'j3':self.ak4LightJets[0],'j4':None,
                                                'nJet':3,'nbJet':2,
                                                'sel':MuSelObjAk4JetsLooseExclusiveResolved2b1j.sel,
                                                'suffix':MuSelObjAk4JetsLooseExclusiveResolved2b1j.selName})
                     # Resolved Two Btag Two jet #      
                    if "TightResolved2b2j" in jetplot_level:
                        ChannelDictList.append({'channel': 'El','met': self.corrMET,'lep':ElColl[0],
                                                'j1':self.ak4BJets[0],'j2':self.ak4BJets[1],'j3':self.ak4LightJets[0],'j4':self.ak4LightJets[1],
                                                'nJet':4,'nbJet':2,
                                                'sel':ElSelObjAk4JetsTightExclusiveResolved2b2j.sel,
                                                'suffix':ElSelObjAk4JetsTightExclusiveResolved2b2j.selName})
                        ChannelDictList.append({'channel': 'Mu','met': self.corrMET,'lep':MuColl[0],
                                                'j1':self.ak4BJets[0],'j2':self.ak4BJets[1],'j3':self.ak4LightJets[0],'j4':self.ak4LightJets[1],
                                                'nJet':4,'nbJet':2,
                                                'sel':MuSelObjAk4JetsTightExclusiveResolved2b2j.sel,
                                                'suffix':MuSelObjAk4JetsTightExclusiveResolved2b2j.selName})                 

                    for channelDict in ChannelDictList:
                        plots.extend(makeHighLevelPlotsResolved(**{k:channelDict[k] for k in ResolvedKeys}))

                # To plot the highlevel distributions for both of the semiBoosted & boosted scenario #    
                if any("Boosted" in item for item in jetplot_level):
                    hasEle  = True if op.rng_len(ElColl) > 0 else False 
                    hasMuon = True if op.rng_len(MuColl) > 0 else False
                    ChannelDictList=[]
                    # SemiBoosted #
                    if "SemiBoostedHbb" in jetplot_level:
                        ChannelDictList.append({'channel': 'El','met': self.corrMET,'lep':ElColl[0],
                                                'j1':self.ak8BJets[0],'j2':self.ak4JetsCleanedFromAk8b[0],
                                                'bothAreFat':False,'sel':ElSelObjAk8JetsHbbBoosted.sel,
                                                'suffix':ElSelObjAk8JetsHbbBoosted.selName})
                        ChannelDictList.append({'channel': 'Mu','met': self.corrMET,'lep':MuColl[0],
                                                'j1':self.ak8BJets[0],'j2':self.ak4JetsCleanedFromAk8b[0],
                                                'bothAreFat':False,'sel':MuSelObjAk8JetsHbbBoosted.sel,
                                                'suffix':MuSelObjAk8JetsHbbBoosted.selName})                 
                    # Boosted #
                    if "Boosted" in jetplot_level:
                        ChannelDictList.append({'channel': 'El','met': self.corrMET,'lep':ElColl[0],
                                                'j1':self.ak8BJets[0],'j2':self.ak8nonBJets[0],
                                                'bothAreFat':True,'sel':ElSelObjAk8JetsInclusiveBoosted.sel,
                                                'suffix':ElSelObjAk8JetsInclusiveBoosted.selName})
                        ChannelDictList.append({'channel': 'Mu','met': self.corrMET,'lep':MuColl[0],
                                                'j1':self.ak8BJets[0],'j2':self.ak8nonBJets[0],
                                                'bothAreFat':True,'sel':MuSelObjAk8JetsInclusiveBoosted.sel,
                                                'suffix':MuSelObjAk8JetsInclusiveBoosted.selName})                 

                    for channelDict in ChannelDictList:
                        plots.extend(makeHighLevelPlotsBoosted(**{k:channelDict[k] for k in BoostedKeys}))

        #----- Add the Yield plots -----#
        plots.extend(self.yieldPlots.returnPlots())

        return plots

