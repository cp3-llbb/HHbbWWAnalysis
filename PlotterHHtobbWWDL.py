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
from fakeHelper import DataDrivenFake

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
        # Change the way the FakeExtrapolation is postProcesses (avoids overriding the `postProcess` method) 
        if "FakeExtrapolation" in self.datadrivenContributions:
            contrib = self.datadrivenContributions["FakeExtrapolation"]
            self.datadrivenContributions["FakeExtrapolation"] = DataDrivenFake(contrib.name, contrib.config)


    def definePlots(self, t, noSel, sample=None, sampleCfg=None): 
        noSel = super(PlotterNanoHHtobbWWDL,self).prepareObjects(t, noSel, sample, sampleCfg, 'DL')

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
            #plots.append(CutFlowReport("BtagReweightingCutFlowReport",noSel))
            return plots

        #----- DY stitching study -----#
        if self.args.DYStitching:
            if sampleCfg['group'] != 'DY':
                raise RuntimeError("Stitching is only done on DY MC samples")
            plots.append(Plot.make1D("LHE_Njets",
                                     t.LHE.Njets,
                                     noSel,
                                     EquidistantBinning(5,0.,5),
                                     xTitle='LHE Njets'))
            #plots.append(CutFlowReport("DYStitchingCutFlowReport",noSel))
            return plots


        #----- Dileptons -----#
        selObjectDict = makeDoubleLeptonSelection(self,noSel,plot_yield=True)
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
                OSElElDilepton = self.ElElDileptonFakeSel
                OSMuMuDilepton = self.MuMuDileptonFakeSel
                OSElMuDilepton = self.ElMuDileptonFakeSel
            elif selectionType == "Tight":
                OSElElDilepton = switch_on_index([0],op.rng_len(self.ElElDileptonTightSel)>=1,self.ElElDileptonTightSel,self.ElElDileptonFakeExtrapolationSel)
                OSMuMuDilepton = switch_on_index([0],op.rng_len(self.MuMuDileptonTightSel)>=1,self.MuMuDileptonTightSel,self.MuMuDileptonFakeExtrapolationSel)
                OSElMuDilepton = switch_on_index([0],op.rng_len(self.ElMuDileptonTightSel)>=1,self.ElMuDileptonTightSel,self.ElMuDileptonFakeExtrapolationSel)
#                OSElElDilepton = [self.ElElDileptonFakeSel[op.switch(op.rng_len(self.ElElDileptonTightSel)>=1,self.ElElDileptonTightSel[0].index,self.ElElDileptonFakeExtrapolationSel[0].index)]]
#                OSMuMuDilepton = [self.MuMuDileptonFakeSel[op.switch(op.rng_len(self.MuMuDileptonTightSel)>=1,self.MuMuDileptonTightSel[0].index,self.MuMuDileptonFakeExtrapolationSel[0].index)]]
#                OSElMuDilepton = [self.ElMuDileptonFakeSel[op.switch(op.rng_len(self.ElMuDileptonTightSel)>=1,self.ElMuDileptonTightSel[0].index,self.ElMuDileptonFakeExtrapolationSel[0].index)]]
            elif selectionType == "FakeExtrapolation":
                OSElElDilepton = self.ElElDileptonFakeExtrapolationSel
                OSMuMuDilepton = self.MuMuDileptonFakeExtrapolationSel
                OSElMuDilepton = self.ElMuDileptonFakeExtrapolationSel

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
                    plots.extend(triggerPlots(**channelDict,triggerDict=self.triggersPerPrimaryDataset))
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
                    plots.extend(makeAk8JetsPlots(**{k:channelDict[k] for k in FatJetKeys}))
                    # MET #
                    plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))

                         
            #----- Resolved selection -----#
            if any(item in ["Resolved0Btag","Resolved1Btag","Resolved2Btag"] for item in jetsel_level): # If any of resolved category is asked
                #----- select the jets -----#
                aka4JetsByBtagScore = op.sort(self.ak4Jets, lambda jet : -jet.btagDeepFlavB)

                container0b2j = [ t.Jet[self.ak4LightJets[i].idx] for i in range(2) ]
                container1b1j = [ t.Jet[op.switch(op.rng_len(self.ak4BJets) == 1, self.ak4BJets[0].idx,  aka4JetsByBtagScore[0].idx)] ,
                                  t.Jet[op.switch(op.rng_len(self.ak4BJets) == 1, self.ak4LightJets[0].idx,  aka4JetsByBtagScore[1].idx)]]
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
            if "Boosted" in jetsel_level:
                print ("...... Processing Boosted jet selection")
                ElElSelObjectDileptonAk8JetsInclusiveBoosted = makeInclusiveBoostedSelection(self,ElElSelObjectDileptonAk8Jets,copy_sel=True,plot_yield=True)
                MuMuSelObjectDileptonAk8JetsInclusiveBoosted = makeInclusiveBoostedSelection(self,MuMuSelObjectDileptonAk8Jets,copy_sel=True,plot_yield=True)
                ElMuSelObjectDileptonAk8JetsInclusiveBoosted = makeInclusiveBoostedSelection(self,ElMuSelObjectDileptonAk8Jets,copy_sel=True,plot_yield=True)

                # Lepton + jet Plots #
                ChannelDictList = []
                if "Boosted" in jetplot_level:
                    # Cut flow report #
                    cutFlowPlots.append(CutFlowReport(ElElSelObjectDileptonAk8JetsInclusiveBoosted.selName,ElElSelObjectDileptonAk8JetsInclusiveBoosted.sel))
                    cutFlowPlots.append(CutFlowReport(MuMuSelObjectDileptonAk8JetsInclusiveBoosted.selName,MuMuSelObjectDileptonAk8JetsInclusiveBoosted.sel))
                    cutFlowPlots.append(CutFlowReport(ElMuSelObjectDileptonAk8JetsInclusiveBoosted.selName,ElMuSelObjectDileptonAk8JetsInclusiveBoosted.sel))
                    if not self.args.OnlyYield:
                        ChannelDictList.append({'channel':'ElEl','sel':ElElSelObjectDileptonAk8JetsInclusiveBoosted.sel,'dilepton':OSElElDilepton[0],'fatjet':self.ak8BJets[0],'suffix':ElElSelObjectDileptonAk8JetsInclusiveBoosted.selName,'is_MC':self.is_MC})
                        ChannelDictList.append({'channel':'MuMu','sel':MuMuSelObjectDileptonAk8JetsInclusiveBoosted.sel,'dilepton':OSMuMuDilepton[0],'fatjet':self.ak8BJets[0],'suffix':MuMuSelObjectDileptonAk8JetsInclusiveBoosted.selName,'is_MC':self.is_MC})
                        ChannelDictList.append({'channel':'ElMu','sel':ElMuSelObjectDileptonAk8JetsInclusiveBoosted.sel,'dilepton':OSElMuDilepton[0],'fatjet':self.ak8BJets[0],'suffix':ElMuSelObjectDileptonAk8JetsInclusiveBoosted.selName,'is_MC':self.is_MC})

                BoostedJetsN = {'objName':'Ak8BJets','objCont':self.ak8BJets,'Nmax':5,'xTitle':'N(Ak8 Bjets)'}
                
                for channelDict in ChannelDictList:
                    # Dilepton #
                    plots.extend(makeDileptonPlots(**{k:channelDict[k] for k in LeptonKeys}))
                    # Number of jets #
                    plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**BoostedJetsN))
                    plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**FatJetsN))
                    # Ak8 Jets #
                    plots.extend(makeAk8JetsPlots(**{k:channelDict[k] for k in FatJetKeys}))
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
                # Boosted #
                if "Boosted" in jetplot_level:
                    ChannelDictList.append({'channel': 'ElEl','met': self.corrMET,'l1':OSElElDilepton[0][0],'l2':OSElElDilepton[0][1],'j1':self.ak8BJets[0].subJet1,'j2':self.ak8BJets[0].subJet2,'sel':ElElSelObjectDileptonAk8JetsInclusiveBoosted.sel,'suffix':ElElSelObjectDileptonAk8JetsInclusiveBoosted.selName})
                    ChannelDictList.append({'channel': 'MuMu','met': self.corrMET,'l1':OSMuMuDilepton[0][0],'l2':OSMuMuDilepton[0][1],'j1':self.ak8BJets[0].subJet1,'j2':self.ak8BJets[0].subJet2,'sel':MuMuSelObjectDileptonAk8JetsInclusiveBoosted.sel,'suffix':MuMuSelObjectDileptonAk8JetsInclusiveBoosted.selName})
                    ChannelDictList.append({'channel': 'ElMu','met': self.corrMET,'l1':OSElMuDilepton[0][0],'l2':OSElMuDilepton[0][1],'j1':self.ak8BJets[0].subJet1,'j2':self.ak8BJets[0].subJet2,'sel':ElMuSelObjectDileptonAk8JetsInclusiveBoosted.sel,'suffix':ElMuSelObjectDileptonAk8JetsInclusiveBoosted.selName})

            for channelDict in ChannelDictList:
                plots.extend(makeHighLevelQuantities(**channelDict))

        #----- Add the Yield plots -----#
        plots.extend(self.yieldPlots.returnPlots())

        #plots.extend(cutFlowPlots)

        #----- Return -----#
        return plots
