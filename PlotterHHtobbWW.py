import os
import sys
from copy import copy

from bamboo.analysismodules import HistogramsModule
from bamboo import treefunctions as op
from bamboo.plots import Plot, EquidistantBinning, SummedPlot

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)))) # Add scripts in this directory
from BaseHHtobbWW import BaseNanoHHtobbWW
from plotDef import *
from selectionDef import *


#===============================================================================================#
#                                       PlotterHHtobbWW                                         #
#===============================================================================================#
class PlotterNanoHHtobbWW(BaseNanoHHtobbWW,HistogramsModule):
    """ Plotter module: HH->bbW(->e/µ nu)W(->e/µ nu) histograms from NanoAOD """
    def __init__(self, args):
        super(PlotterNanoHHtobbWW, self).__init__(args)

    def definePlots(self, t, noSel, sample=None, sampleCfg=None): 
        noSel = super(PlotterNanoHHtobbWW,self).prepareObjects(t, noSel, sample, sampleCfg)

        plots = []

        era = sampleCfg['era']

        self.yieldPlots = makeYieldPlots()

        #----- Dileptons -----#
        selObjectDict = makeLeptonSelection(self,noSel,plot_yield=True)
        # selObjectDict : keys -> level (str)
        #                 values -> [ElEl,MuMu,ElMu] x Selection object

        # Select the jets selections that will be done depending on user input #
        jet_level = ["Ak4","Ak8","Resolved0Btag","Resolved1Btag","Resolved2Btag","Boosted"]
        jetplot_level = [arg for (arg,boolean) in self.args.__dict__.items() if arg in jet_level and boolean]
        if len(jetplot_level) == 0:  
            jetplot_level = jet_level # If nothing said, will do all
        jetsel_level = copy(jetplot_level)  # A plot level might need a previous selection that needs to be defined but not necessarily plotted
        if "Resolved0Btag" in jetsel_level or "Resolved1Btag" in jetsel_level or "Resolved2Btag" in jetsel_level:
            jetsel_level.append("Ak4") # Resolved needs the Ak4 selection
        if "Boosted" in jetsel_level:
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
                OSElElDilepton = self.OSElElDileptonFakeSel
                OSMuMuDilepton = self.OSMuMuDileptonFakeSel
                OSElMuDilepton = self.OSElMuDileptonFakeSel
            elif selectionType == "Tight":
                OSElElDilepton = self.OSElElDileptonTightSel
                OSMuMuDilepton = self.OSMuMuDileptonTightSel
                OSElMuDilepton = self.OSElMuDileptonTightSel
            elif selectionType == "FakeExtrapolation":
                OSElElDilepton = self.OSElElDileptonFakeExtrapolationSel
                OSMuMuDilepton = self.OSMuMuDileptonFakeExtrapolationSel
                OSElMuDilepton = self.OSElMuDileptonFakeExtrapolationSel

            #----- Separate selections ------#
            ElElSelObjectDilepton = selectionList[0]
            MuMuSelObjectDilepton = selectionList[1]
            ElMuSelObjectDilepton = selectionList[2]

            if not self.args.OnlyYield:
                #----- Trigger plots -----#
                plots.extend(triggerPlots(sel         = ElElSelObjectDilepton.sel,
                                          triggerDict = self.triggersPerPrimaryDataset,
                                          suffix      = ElElSelObjectDilepton.selName,
                                          channel     = "ElEl"))
                plots.extend(triggerPlots(sel         = MuMuSelObjectDilepton.sel,
                                          triggerDict = self.triggersPerPrimaryDataset,
                                          suffix      = MuMuSelObjectDilepton.selName,
                                          channel     = "MuMu"))
                plots.extend(triggerPlots(sel         = ElMuSelObjectDilepton.sel,
                                          triggerDict = self.triggersPerPrimaryDataset,
                                          suffix      = ElMuSelObjectDilepton.selName,
                                          channel     = "ElMu"))

                #----- Lepton plots -----#
                # Dilepton channel plots #
                plots.extend(channelPlot(sel        = ElElSelObjectDilepton.sel,
                                         DilepElEl  = OSElElDilepton,
                                         DilepMuMu  = OSMuMuDilepton,
                                         DilepElMu  = OSElMuDilepton,
                                         suffix     = ElElSelObjectDilepton.selName,
                                         channel    = "ElEl"))
                plots.extend(channelPlot(sel        = MuMuSelObjectDilepton.sel,
                                         DilepElEl  = OSElElDilepton,
                                         DilepMuMu  = OSMuMuDilepton,
                                         DilepElMu  = OSElMuDilepton,
                                         suffix     = MuMuSelObjectDilepton.selName,
                                         channel    = "MuMu"))
                plots.extend(channelPlot(sel        = ElMuSelObjectDilepton.sel,
                                         DilepElEl  = OSElElDilepton,
                                         DilepMuMu  = OSMuMuDilepton,
                                         DilepElMu  = OSElMuDilepton,
                                         suffix     = ElMuSelObjectDilepton.selName,
                                         channel    = "ElMu"))
        
            #----- Ak4 jets selection -----#
            LeptonKeys = ['channel','sel','dilepton','suffix','is_MC']
            JetKeys    = ['channel','sel','leadjet','subleadjet','lead_is_b','sublead_is_b','suffix','is_MC']
            commonItems = ['channel','sel','suffix']
            if "Ak4" in jetsel_level:
                print ("...... Processing Ak4 jet selection")
                ElElSelObjectDileptonAk4Jets = makeAk4JetSelection(self,ElElSelObjectDilepton,copy_sel=True,plot_yield=True)
                MuMuSelObjectDileptonAk4Jets = makeAk4JetSelection(self,MuMuSelObjectDilepton,copy_sel=True,plot_yield=True)
                ElMuSelObjectDileptonAk4Jets = makeAk4JetSelection(self,ElMuSelObjectDilepton,copy_sel=True,plot_yield=True)

                # Jet and lepton plots #
                ChannelDictList = []
                if not self.args.OnlyYield and "Ak4" in jetplot_level:
                    ChannelDictList.append({'channel':'ElEl','sel':ElElSelObjectDileptonAk4Jets.sel,'dilepton':OSElElDilepton[0],'leadjet':self.ak4Jets[0],'subleadjet':self.ak4Jets[1],'lead_is_b':False,'sublead_is_b':False,'suffix':ElElSelObjectDileptonAk4Jets.selName,'is_MC':self.is_MC})
                    ChannelDictList.append({'channel':'MuMu','sel':MuMuSelObjectDileptonAk4Jets.sel,'dilepton':OSMuMuDilepton[0],'leadjet':self.ak4Jets[0],'subleadjet':self.ak4Jets[1],'lead_is_b':False,'sublead_is_b':False,'suffix':MuMuSelObjectDileptonAk4Jets.selName,'is_MC':self.is_MC})
                    ChannelDictList.append({'channel':'ElMu','sel':ElMuSelObjectDileptonAk4Jets.sel,'dilepton':OSElMuDilepton[0],'leadjet':self.ak4Jets[0],'subleadjet':self.ak4Jets[1],'lead_is_b':False,'sublead_is_b':False,'suffix':ElMuSelObjectDileptonAk4Jets.selName,'is_MC':self.is_MC})

                JetsN = {'objName':'Ak4Jets','objCont':self.ak4Jets,'Nmax':10,'xTitle':'N(Ak4 jets)'}
                                                            
                for channelDict in ChannelDictList:
                    # Dilepton #
                    plots.extend(makeDileptonPlots(**{k:channelDict[k] for k in LeptonKeys}))
                    # Number of jets #
                    plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**JetsN))
                    # Ak4 Jets #
                    plots.extend(makeAk4JetsPlots(**{k:channelDict[k] for k in JetKeys}))
                    # MET #
                    plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))

            ##### Ak8 jets selection #####
            if "Ak8" in jetsel_level:
                print ("...... Processing Ak8 jet selection")
                ElElSelObjectDileptonAk8Jets = makeAk8JetSelection(self,ElElSelObjectDilepton,copy_sel=True,plot_yield=True)
                MuMuSelObjectDileptonAk8Jets = makeAk8JetSelection(self,MuMuSelObjectDilepton,copy_sel=True,plot_yield=True)
                ElMuSelObjectDileptonAk8Jets = makeAk8JetSelection(self,ElMuSelObjectDilepton,copy_sel=True,plot_yield=True)

                # Fatjets plots #
                ChannelDictList = []
                if not self.args.OnlyYield and "Ak8" in jetplot_level:
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
                    plots.extend(makeAk8JetsPlots(**{k:channelDict[k] for k in FatJetKeys}))
                    # MET #
                    plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))

                         
            #----- Resolved selection -----#
            if any(item in ["Resolved0Btag","Resolved1Btag","Resolved2Btag"] for item in jetsel_level): # If any of resolved category is asked
                ChannelDictList = []
                #----- Resolved selection : 0 Btag -----#
                if "Resolved0Btag" in jetsel_level:
                    print ("...... Processing Resolved jet (0 btag) selection")
                    ElElSelObjectDileptonAk4JetsExclusiveResolvedNoBtag = makeExclusiveResolvedNoBtagSelection(self,ElElSelObjectDileptonAk4Jets,copy_sel=True,plot_yield=True)
                    MuMuSelObjectDileptonAk4JetsExclusiveResolvedNoBtag = makeExclusiveResolvedNoBtagSelection(self,MuMuSelObjectDileptonAk4Jets,copy_sel=True,plot_yield=True)
                    ElMuSelObjectDileptonAk4JetsExclusiveResolvedNoBtag = makeExclusiveResolvedNoBtagSelection(self,ElMuSelObjectDileptonAk4Jets,copy_sel=True,plot_yield=True)

                    if not self.args.OnlyYield and "Resolved0Btag" in jetplot_level:
                        ChannelDictList.append({'channel':'ElEl','sel':ElElSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.sel,'dilepton':OSElElDilepton[0],'leadjet':self.ak4LightJets[0],'subleadjet':self.ak4LightJets[1],'lead_is_b':False,'sublead_is_b':False,'suffix':ElElSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.selName,'is_MC':self.is_MC})
                        ChannelDictList.append({'channel':'MuMu','sel':MuMuSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.sel,'dilepton':OSMuMuDilepton[0],'leadjet':self.ak4LightJets[0],'subleadjet':self.ak4LightJets[1],'lead_is_b':False,'sublead_is_b':False,'suffix':MuMuSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.selName,'is_MC':self.is_MC})
                        ChannelDictList.append({'channel':'ElMu','sel':ElMuSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.sel,'dilepton':OSElMuDilepton[0],'leadjet':self.ak4LightJets[0],'subleadjet':self.ak4LightJets[1],'lead_is_b':False,'sublead_is_b':False,'suffix':ElMuSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.selName,'is_MC':self.is_MC})
                #----  Resolved selection : 1 Btag  -----#
                if "Resolved1Btag" in jetsel_level:
                    print ("...... Processing Resolved jet (1 btag) selection")
                    ElElSelObjectDileptonAk4JetsExclusiveResolvedOneBtag = makeExclusiveResolvedOneBtagSelection(self,ElElSelObjectDileptonAk4Jets,copy_sel=True,plot_yield=True)
                    MuMuSelObjectDileptonAk4JetsExclusiveResolvedOneBtag = makeExclusiveResolvedOneBtagSelection(self,MuMuSelObjectDileptonAk4Jets,copy_sel=True,plot_yield=True)
                    ElMuSelObjectDileptonAk4JetsExclusiveResolvedOneBtag = makeExclusiveResolvedOneBtagSelection(self,ElMuSelObjectDileptonAk4Jets,copy_sel=True,plot_yield=True)

                    if not self.args.OnlyYield and "Resolved1Btag" in jetplot_level:
                        ChannelDictList.append({'channel':'ElEl','sel':ElElSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.sel,'dilepton':OSElElDilepton[0],'leadjet':self.ak4BJets[0],'subleadjet':self.ak4LightJets[0],'lead_is_b':True,'sublead_is_b':False,'suffix':ElElSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.selName,'is_MC':self.is_MC})
                        ChannelDictList.append({'channel':'MuMu','sel':MuMuSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.sel,'dilepton':OSMuMuDilepton[0],'leadjet':self.ak4BJets[0],'subleadjet':self.ak4LightJets[0],'lead_is_b':True,'sublead_is_b':False,'suffix':MuMuSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.selName,'is_MC':self.is_MC})
                        ChannelDictList.append({'channel':'ElMu','sel':ElMuSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.sel,'dilepton':OSElMuDilepton[0],'leadjet':self.ak4BJets[0],'subleadjet':self.ak4LightJets[0],'lead_is_b':True,'sublead_is_b':False,'suffix':ElMuSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.selName,'is_MC':self.is_MC})
                #----- Resolved selection : 2 Btags -----#
                if "Resolved2Btag" in jetsel_level:
                    print ("...... Processing Resolved jet (2 btags) selection")
                    ElElSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags = makeExclusiveResolvedTwoBtagsSelection(self,ElElSelObjectDileptonAk4Jets,copy_sel=True,plot_yield=True)
                    MuMuSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags = makeExclusiveResolvedTwoBtagsSelection(self,MuMuSelObjectDileptonAk4Jets,copy_sel=True,plot_yield=True)
                    ElMuSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags = makeExclusiveResolvedTwoBtagsSelection(self,ElMuSelObjectDileptonAk4Jets,copy_sel=True,plot_yield=True)

                    if not self.args.OnlyYield and "Resolved2Btag" in jetplot_level:
                        ChannelDictList.append({'channel':'ElEl','sel':ElElSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.sel,'dilepton':OSElElDilepton[0],'leadjet':self.ak4BJets[0],'subleadjet':self.ak4BJets[1],'lead_is_b':True,'sublead_is_b':True,'suffix':ElElSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.selName,'is_MC':self.is_MC})
                        ChannelDictList.append({'channel':'MuMu','sel':MuMuSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.sel,'dilepton':OSMuMuDilepton[0],'leadjet':self.ak4BJets[0],'subleadjet':self.ak4BJets[1],'lead_is_b':True,'sublead_is_b':True,'suffix':MuMuSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.selName,'is_MC':self.is_MC})
                        ChannelDictList.append({'channel':'ElMu','sel':ElMuSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.sel,'dilepton':OSElMuDilepton[0],'leadjet':self.ak4BJets[0],'subleadjet':self.ak4BJets[1],'lead_is_b':True,'sublead_is_b':True,'suffix':ElMuSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.selName,'is_MC':self.is_MC})

                # Lepton + jet Plots #

                ResolvedJetsN = {'objName':'Ak4BJets','objCont':self.ak4BJets,'Nmax':5,'xTitle':'N(Ak4 Bjets)'}
        
                for channelDict in ChannelDictList:
                    # Dilepton #
                    plots.extend(makeDileptonPlots(**{k:channelDict[k] for k in LeptonKeys}))
                    # Number of jets #
                    plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**ResolvedJetsN))
                    # Ak4 Jets #
                    plots.extend(makeAk4JetsPlots(**{k:channelDict[k] for k in JetKeys}))
                    # MET #
                    plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))
            
            #----- Boosted selection -----#
            if "Boosted" in jetsel_level:
                print ("...... Processing Boosted jet selection")
                ElElSelObjectDileptonAk8JetsInclusiveBoosted = makeInclusiveBoostedSelection(self,ElElSelObjectDileptonAk8Jets,copy_sel=True,plot_yield=True)
                MuMuSelObjectDileptonAk8JetsInclusiveBoosted = makeInclusiveBoostedSelection(self,MuMuSelObjectDileptonAk8Jets,copy_sel=True,plot_yield=True)
                ElMuSelObjectDileptonAk8JetsInclusiveBoosted = makeInclusiveBoostedSelection(self,ElMuSelObjectDileptonAk8Jets,copy_sel=True,plot_yield=True)

                # Lepton + jet Plots #
                ChannelDictList = []
                if not self.args.OnlyYield and "Boosted" in jetplot_level:
                    ChannelDictList.append({'channel':'ElEl','sel':ElElSelObjectDileptonAk8JetsInclusiveBoosted.sel,'dilepton':OSElElDilepton[0],'fatjet':self.ak8BJets[0],'suffix':ElElSelObjectDileptonAk8JetsInclusiveBoosted.selName,'is_MC':self.is_MC})
                    ChannelDictList.append({'channel':'MuMu','sel':MuMuSelObjectDileptonAk8JetsInclusiveBoosted.sel,'dilepton':OSMuMuDilepton[0],'fatjet':self.ak8BJets[0],'suffix':MuMuSelObjectDileptonAk8JetsInclusiveBoosted.selName,'is_MC':self.is_MC})
                    ChannelDictList.append({'channel':'ElMu','sel':ElMuSelObjectDileptonAk8JetsInclusiveBoosted.sel,'dilepton':OSElMuDilepton[0],'fatjet':self.ak8BJets[0],'suffix':ElMuSelObjectDileptonAk8JetsInclusiveBoosted.selName,'is_MC':self.is_MC})

                BoostedJetsN = {'objName':'Ak8BJets','objCont':self.ak8BJets,'Nmax':5,'xTitle':'N(Ak8 Bjets)'}
                
                for channelDict in ChannelDictList:
                    # Dilepton #
                    plots.extend(makeDileptonPlots(**{k:channelDict[k] for k in LeptonKeys}))
                    # Number of jets #
                    plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**BoostedJetsN))
                    # Ak8 Jets #
                    plots.extend(makeAk8JetsPlots(**{k:channelDict[k] for k in FatJetKeys}))
                    # MET #
                    plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))
            
            #----- High-level combinations -----#
            ChannelDictList = []
            if not self.args.OnlyYield:
                # Resolved No Btag #
                if "Resolved0Btag" in jetplot_level:
                    ChannelDictList.append({'channel': 'ElEl','met': self.corrMET,'l1':OSElElDilepton[0][0],'l2':OSElElDilepton[0][1],'j1':self.ak4LightJets[0],'j2':self.ak4LightJets[1],'sel':ElElSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.sel,'suffix':ElElSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.selName})
                    ChannelDictList.append({'channel': 'MuMu','met': self.corrMET,'l1':OSMuMuDilepton[0][0],'l2':OSMuMuDilepton[0][1],'j1':self.ak4LightJets[0],'j2':self.ak4LightJets[1],'sel':MuMuSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.sel,'suffix':MuMuSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.selName})
                    ChannelDictList.append({'channel': 'ElMu','met': self.corrMET,'l1':OSElMuDilepton[0][0],'l2':OSElMuDilepton[0][1],'j1':self.ak4LightJets[0],'j2':self.ak4LightJets[1],'sel':ElMuSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.sel,'suffix':ElMuSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.selName})
                # Resolved One Btag #
                if "Resolved1Btag" in jetplot_level:
                    ChannelDictList.append({'channel': 'ElEl','met': self.corrMET,'l1':OSElElDilepton[0][0],'l2':OSElElDilepton[0][1],'j1':self.ak4BJets[0],'j2':self.ak4LightJets[0],'sel':ElElSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.sel,'suffix':ElElSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.selName})
                    ChannelDictList.append({'channel': 'MuMu','met': self.corrMET,'l1':OSMuMuDilepton[0][0],'l2':OSMuMuDilepton[0][1],'j1':self.ak4BJets[0],'j2':self.ak4LightJets[0],'sel':MuMuSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.sel,'suffix':MuMuSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.selName})
                    ChannelDictList.append({'channel': 'ElMu','met': self.corrMET,'l1':OSElMuDilepton[0][0],'l2':OSElMuDilepton[0][1],'j1':self.ak4BJets[0],'j2':self.ak4LightJets[0],'sel':ElMuSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.sel,'suffix':ElMuSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.selName})
                # Resolved Two Btags  #
                if "Resolved2Btag" in jetplot_level:
                    ChannelDictList.append({'channel': 'ElEl','met': self.corrMET,'l1':OSElElDilepton[0][0],'l2':OSElElDilepton[0][1],'j1':self.ak4BJets[0],'j2':self.ak4BJets[1],'sel':ElElSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.sel,'suffix':ElElSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.selName})
                    ChannelDictList.append({'channel': 'MuMu','met': self.corrMET,'l1':OSMuMuDilepton[0][0],'l2':OSMuMuDilepton[0][1],'j1':self.ak4BJets[0],'j2':self.ak4BJets[1],'sel':MuMuSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.sel,'suffix':MuMuSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.selName})
                    ChannelDictList.append({'channel': 'ElMu','met': self.corrMET,'l1':OSElMuDilepton[0][0],'l2':OSElMuDilepton[0][1],'j1':self.ak4BJets[0],'j2':self.ak4BJets[1],'sel':ElMuSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.sel,'suffix':ElMuSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.selName})
                # Boosted #
                if "Boosted" in jetplot_level:
                    ChannelDictList.append({'channel': 'ElEl','met': self.corrMET,'l1':OSElElDilepton[0][0],'l2':OSElElDilepton[0][1],'j1':self.ak8BJets[0].subJet1,'j2':self.ak8BJets[0].subJet2,'sel':ElElSelObjectDileptonAk8JetsInclusiveBoosted.sel,'suffix':ElElSelObjectDileptonAk8JetsInclusiveBoosted.selName})
                    ChannelDictList.append({'channel': 'MuMu','met': self.corrMET,'l1':OSMuMuDilepton[0][0],'l2':OSMuMuDilepton[0][1],'j1':self.ak8BJets[0].subJet1,'j2':self.ak8BJets[0].subJet2,'sel':MuMuSelObjectDileptonAk8JetsInclusiveBoosted.sel,'suffix':MuMuSelObjectDileptonAk8JetsInclusiveBoosted.selName})
                    ChannelDictList.append({'channel': 'ElMu','met': self.corrMET,'l1':OSElMuDilepton[0][0],'l2':OSElMuDilepton[0][1],'j1':self.ak8BJets[0].subJet1,'j2':self.ak8BJets[0].subJet2,'sel':ElMuSelObjectDileptonAk8JetsInclusiveBoosted.sel,'suffix':ElMuSelObjectDileptonAk8JetsInclusiveBoosted.selName})

            for channelDict in ChannelDictList:
                plots.extend(makeHighLevelQuantities(**channelDict))

        #----- Add the Yield plots -----#
        plots.extend(self.yieldPlots.returnPlots())

        return plots

