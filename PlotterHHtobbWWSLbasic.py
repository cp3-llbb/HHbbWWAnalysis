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
from bamboo.analysisutils import printCutFlowReports, addPrintout

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)))) # Add scripts in this directory
from BaseHHtobbWW import BaseNanoHHtobbWW
from plotDef import *
from selectionDef import *
from JPA import *
import mvaEvaluatorSL_nonres
from DDHelper import DataDrivenFake, DataDrivenDY
from bamboo.root import gbl
import ROOT
from functools import partial

if not hasattr(gbl, "bamboo_printEntry"):
    gbl.gInterpreter.Declare("""
    bool bamboo_printEntry(long entry, float val) {
      std::cout << "Processing entry #" << entry << ": val " << val << std::endl;
      return true;
    }""")
    gbl.gInterpreter.Declare("""
    bool bamboo_printComb3_(long entry, const rdfhelpers::Combination<3>& comb, std::size_t idx, double value, UInt_t nJets, UInt_t nComb, const ROOT::VecOps::RVec<rdfhelpers::Combination<3>>& combs, const ROOT::VecOps::RVec<std::size_t>& rng) {
      std::cout << "Combination #" << idx << " (value " << value << ") for " << entry << ": " << comb.get(0) << ", " << comb.get(1) << ", " << comb.get(2) << ", nJets=" << nJets << ", nSelJets=" << rng.size() << ", nCombinations=" << nComb << std::endl;
      if ( ( comb.get(0) >= nJets ) || ( comb.get(1) >= nJets ) || ( comb.get(2) >= nJets ) ) {
        std::cout << "All combinations: " << std::endl;
        for ( const auto& icomb : combs ) {
          std::cout << " - " << icomb.get(0) << ", " << icomb.get(1) << ", " << icomb.get(2) << std::endl;
        }
        std::cout << "From range with size " << rng.size() << ": ";
        for ( auto idx : rng) {
          std::cout << " " << idx;
        }
        std::cout << std::endl;
      }
      return true;
    }""")
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
            self.datadrivenContributions["DYEstimation"] = DataDrivenDY(contrib.name, contrib.config,"PseudoData" in self.datadrivenContributions)

    def definePlots(self, t, noSel, sample=None, sampleCfg=None): 
        noSel = super(PlotterNanoHHtobbWWSL,self).prepareObjects(t, noSel, sample, sampleCfg, 'SL')
        era = sampleCfg['era']
        plots = []

        if hasattr(self,'base_plots'):
            plots.extend(self.base_plots)

        cutFlowPlots = []
        
        #----- Machine Learning Model -----#                
        model_num = "01"
        path_model = os.path.join(os.path.abspath(os.path.dirname(__file__)),'MachineLearning','ml-models','models','multi-classification','dnn','SL',model_num,'model','model.pb')
        input_names = ["lep","jet","fat","met","hl","param","eventnr"]
        output_name = "Identity"
        
        if not self.args.OnlyYield:
            print ("DNN model : %s"%path_model)
            if not os.path.exists(path_model):
                raise RuntimeError('Could not find model file %s'%path_model)
            try:
                DNN = op.mvaEvaluator(path_model,mvaType='Tensorflow',otherArgs=(input_names, output_name))
            except:
                raise RuntimeError('Could not load model %s'%path_model)

        self.sample = sample
        self.sampleCfg = sampleCfg
        self.era = era
        
        self.yieldPlots = makeYieldPlots(self.args.Synchronization)
        
        #yields = CutFlowReport("yields")
        #plots.append(yields)

        
        #----- Singleleptons -----#
        #ElSelObj,MuSelObj = makeSingleLeptonSelection(self,noSel,plot_yield=True)
        ElSelObj,MuSelObj = makeSingleLeptonSelection(self,noSel,plot_yield=True,use_dd=True,fake_selection=self.args.FakeCR)

        #----- Apply jet corrections -----#
        self.beforeJetselection(ElSelObj.sel,'El')
        self.beforeJetselection(MuSelObj.sel,'Mu')

        # selObjectDict : keys -> level (str)
        #                 values -> [El,Mu] x Selection object
        # Select the jets selections that will be done depending on user input #
        resolved_args = ["Resolved2Btag","Resolved1Btag"]
        boosted_args  = ["Boosted"]
        jet_level     = resolved_args + boosted_args
        jet_level.append("Ak4")  # to call all resolved categories
        jet_level.append("Ak8")  # to call all boosted categories
        jetplot_level = [arg for (arg,boolean) in self.args.__dict__.items() if arg in jet_level and boolean]
        if len(jetplot_level) == 0:  
            jetplot_level = jet_level # If nothing said, will do all
        jetsel_level = copy(jetplot_level)  # A plot level might need a previous selection that needs to be defined but not necessarily plotted

        if any(item in boosted_args for item in jetsel_level):
            jetsel_level.append("Ak8") # SemiBoosted & Boosted needs the Ak8 selection
        if any(item in resolved_args for item in jetsel_level):
            jetsel_level.append("Ak4") # Resolved needs the Ak4 selection 

        logger.info ('jetSel_Level: {}'.format(jetsel_level))
            
        # Selections:    
        #---- Lepton selection ----#
        ElColl = [t.Electron[op.switch(op.rng_len(self.electronsTightSel) == 1, self.electronsTightSel[0].idx, self.electronsFakeSel[0].idx)]]
        MuColl = [t.Muon[op.switch(op.rng_len(self.muonsTightSel) == 1, self.muonsTightSel[0].idx, self.muonsFakeSel[0].idx)]]

        '''
        if not self.args.OnlyYield:
            ChannelDictList = []
            ChannelDictList.append({'channel':'El','sel':ElSelObj.sel,'suffix':ElSelObj.selName})
            ChannelDictList.append({'channel':'Mu','sel':MuSelObj.sel,'suffix':MuSelObj.selName})
            for channelDict in ChannelDictList:
                #----- Trigger plots -----#
                plots.extend(singleLeptonTriggerPlots(**channelDict, triggerDict=self.triggersPerPrimaryDataset))
        '''

        LeptonKeys   = ['channel','sel','lepton','suffix','is_MC']
        JetKeys      = ['channel','sel','jet1','jet2','jet3','jet4','suffix','nJet','nbJet','is_MC']
        commonItems  = ['channel','sel','suffix']
            
        selObjectDictList = []
        # ========================== Resolved Categories ========================= #
        if any(item in resolved_args for item in jetsel_level):
            ChannelDictListR = []
            logger.info ("... Processing Ak4Jets Selection for Resolved category : nAk4Jets >= 3 + nAk4BJets >= 1 + nAk8BJets == 0")
            ElSelObjResolved = makeResolvedSelection(self,ElSelObj,copy_sel=True)
            MuSelObjResolved = makeResolvedSelection(self,MuSelObj,copy_sel=True)

            if "Resolved2Btag" in jetplot_level:
                logger.info ('Resolved2Btag Node Selection')
                ElSelObjResolved2b        = makeExclusiveResolvedSelection(self, ElSelObjResolved, nbJet=2, copy_sel=True)
                MuSelObjResolved2b        = makeExclusiveResolvedSelection(self, MuSelObjResolved, nbJet=2, copy_sel=True)

                selObjectDictList.append({'channel':'El','selObject':ElSelObjResolved2b,'category':'Resolved','VBFJets':self.VBFJetPairsResolved})
                selObjectDictList.append({'channel':'Mu','selObject':MuSelObjResolved2b,'category':'Resolved','VBFJets':self.VBFJetPairsResolved})

                #if self.args.onlypost:
                ElSelObjResolved2b.record_yields = True
                MuSelObjResolved2b.record_yields = True
                ElSelObjResolved2b.yieldTitle = 'Resolved2b Channel $e^{\pm}$'
                MuSelObjResolved2b.yieldTitle = 'Resolved2b Channel $\mu^{\pm}$'
                if self.args.PrintYield:
                    self.yields.add(ElSelObjResolved2b.sel)
                    self.yields.add(MuSelObjResolved2b.sel)

                #addPrintout(ElSelObjResolved2b.sel, "bamboo_printEntry", op.extVar("ULong_t", "rdfentry_"), t.event)
                #addPrintout(MuSelObjResolved2b.sel, "bamboo_printEntry", op.extVar("ULong_t", "rdfentry_"), t.event)
                
                if not self.args.OnlyYield:
                    ChannelDictListR.append({'channel':'El','selObj':ElSelObjResolved2b,'sel':ElSelObjResolved2b.sel,
                                             'lepton':ElColl[0],'jets':self.ak4Jets,
                                             'bjets':self.bJetsByScore,'wjets':self.wJetsByPt,
                                             'nJet':op.static_cast("UInt_t",op.rng_len(self.ak4Jets)),'nbJet':2,
                                             'suffix':ElSelObjResolved2b.selName,
                                             'is_MC':self.is_MC})
                    ChannelDictListR.append({'channel':'Mu','selObj':MuSelObjResolved2b,'sel':MuSelObjResolved2b.sel,
                                             'lepton':MuColl[0],'jets':self.ak4Jets,
                                             'bjets':self.bJetsByScore,'wjets':self.wJetsByPt,
                                             'nJet':op.static_cast("UInt_t",op.rng_len(self.ak4Jets)),'nbJet':2,
                                             'suffix':MuSelObjResolved2b.selName,
                                             'is_MC':self.is_MC})

            if "Resolved1Btag" in jetplot_level:
                logger.info ('Resolved1Btag Node Selection')
                ElSelObjResolved1b        = makeExclusiveResolvedSelection(self, ElSelObjResolved, nbJet=1, copy_sel=True)
                MuSelObjResolved1b        = makeExclusiveResolvedSelection(self, MuSelObjResolved, nbJet=1, copy_sel=True)

                plots.append(CutFlowReport("ElSelObjResolved1b", ElSelObjResolved1b.sel))
                plots.append(CutFlowReport("MuSelObjResolved1b", MuSelObjResolved1b.sel))

                selObjectDictList.append({'channel':'El','selObject':ElSelObjResolved1b,'category':'Resolved','VBFJets':self.VBFJetPairsResolved})
                selObjectDictList.append({'channel':'Mu','selObject':MuSelObjResolved1b,'category':'Resolved','VBFJets':self.VBFJetPairsResolved})
                self.yields.add(ElSelObjResolved1b.sel, "ElSelObjResolved1b")
                self.yields.add(MuSelObjResolved1b.sel, "MuSelObjResolved1b")

                #if self.args.onlypost:
                ElSelObjResolved1b.record_yields = True
                MuSelObjResolved1b.record_yields = True
                ElSelObjResolved1b.yieldTitle = 'Resolved1b Channel $e^{\pm}$'
                MuSelObjResolved1b.yieldTitle = 'Resolved1b Channel $\mu^{\pm}$'
                if self.args.PrintYield:
                    self.yields.add(ElSelObjResolved1b.sel)
                    self.yields.add(MuSelObjResolved1b.sel)

                #addPrintout(ElSelObjResolved1b.sel, "bamboo_printEntry", op.extVar("ULong_t", "rdfentry_"), t.event)
                #addPrintout(MuSelObjResolved1b.sel, "bamboo_printEntry", op.extVar("ULong_t", "rdfentry_"), t.event)
                
                if not self.args.OnlyYield:
                    ChannelDictListR.append({'channel':'El','selObj':ElSelObjResolved1b,'sel':ElSelObjResolved1b.sel,
                                             'lepton':ElColl[0],'jets':self.ak4Jets,
                                             'bjets':self.bJetsByScore,'wjets':self.wJetsByPt,
                                             'nJet':op.static_cast("UInt_t",op.rng_len(self.ak4Jets)),'nbJet':1,
                                             'suffix':ElSelObjResolved1b.selName,
                                             'is_MC':self.is_MC})
                    ChannelDictListR.append({'channel':'Mu','selObj':MuSelObjResolved1b,'sel':MuSelObjResolved1b.sel,
                                             'lepton':MuColl[0],'jets':self.ak4Jets,
                                             'bjets':self.bJetsByScore,'wjets':self.wJetsByPt,
                                             'nJet':op.static_cast("UInt_t",op.rng_len(self.ak4Jets)),'nbJet':1,
                                             'suffix':MuSelObjResolved1b.selName,
                                             'is_MC':self.is_MC})
            '''
            for channelDict in ChannelDictListR:
                # Singlelepton #
                #plots.extend(makeSinleptonPlots(**{k:channelDict[k] for k in LeptonKeys}))
                # Number of jets #
                #plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**JetsN))
                #plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**FatJetsN))
                # Ak4 Jets #
                #plots.extend(makeAk4JetsPlots(**{k:channelDict[k] for k in JetKeys},HLL=self.HLL))
                # MET #
                #plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))
                # High level #
                ##plots.extend(makeHighLevelPlotsResolved(**{k:channelDict[k] for k in ResolvedKeys},HLL=self.HLL))
            '''

        # ========================== Boosted Categories ========================= #
        if any(item in boosted_args for item in jetsel_level):
            ChannelDictListB = []
            FatJetKeys  = ['channel','sel','jet1','jet2','jet3','jet4','has1fat1slim','has1fat2slim','suffix']

            logger.info ("...... Processing Boosted Category : nAk8BJets >= 1, nAk4JetsCleanedFromAk8b >= 1")
            ElSelObjBoosted = makeBoostedSelection(self,ElSelObj,copy_sel=True)
            MuSelObjBoosted = makeBoostedSelection(self,MuSelObj,copy_sel=True)

            selObjectDictList.append({'channel':'El','selObject':ElSelObjBoosted,'category':'Boosted','VBFJets':self.VBFJetPairsBoosted})
            selObjectDictList.append({'channel':'Mu','selObject':MuSelObjBoosted,'category':'Boosted','VBFJets':self.VBFJetPairsBoosted})

            ElSelObjBoosted.record_yields = True
            MuSelObjBoosted.record_yields = True
            ElSelObjBoosted.yieldTitle = 'Boosted Channel $e^{\pm}$' 
            MuSelObjBoosted.yieldTitle = 'Boosted Channel $\mu^{\pm}$'
            if self.args.PrintYield:
                self.yields.add(ElSelObjBoosted.sel)
                self.yields.add(MuSelObjBoosted.sel)

            #addPrintout(ElSelObjBoosted.sel, "bamboo_printEntry", op.extVar("ULong_t", "rdfentry_"), t.event)
            #addPrintout(MuSelObjBoosted.sel, "bamboo_printEntry", op.extVar("ULong_t", "rdfentry_"), t.event)

            if not self.args.OnlyYield:
                ChannelDictListB.append({'channel':'El','selObj':ElSelObjBoosted, 'sel':ElSelObjBoosted.sel,
                                         'lep':ElColl[0],'met':self.corrMET,'jets':self.ak4JetsCleanedFromAk8b,
                                         'jet1':self.ak8BJets[0],'jet2':None,'jet3':self.ak4JetsCleanedFromAk8b[0],'jet4':self.ak4JetsCleanedFromAk8b[1],
                                         'has1fat1slim':False,'has1fat2slim':True,'bothAreFat':False,
                                         'suffix':ElSelObjBoosted.selName,
                                         'is_MC':self.is_MC, 'jpaArg':'Hbb2Wj'})
                ChannelDictListB.append({'channel':'Mu','selObj':MuSelObjBoosted,'sel':MuSelObjBoosted.sel,
                                         'lep':MuColl[0],'met':self.corrMET,'jets':self.ak4JetsCleanedFromAk8b,
                                         'jet1':self.ak8BJets[0],'jet2':None,'jet3':self.ak4JetsCleanedFromAk8b[0],'jet4':self.ak4JetsCleanedFromAk8b[1],
                                         'has1fat1slim':False,'has1fat2slim':True,'bothAreFat':False,
                                         'suffix':MuSelObjBoosted.selName,
                                         'is_MC':self.is_MC, 'jpaArg':'Hbb2Wj'})
                
            '''
            for channelDict in ChannelDictListB:
                # Dilepton #
                #plots.extend(makeSinleptonPlots(**{k:channelDict[k] for k in LeptonKeys}))
                # Number of jets #
                ##plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**FatJetsN))
                ##plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**SlimJetsN))
                # Ak8 Jets #
                #plots.extend(makeSingleLeptonAk8JetsPlots(**{k:channelDict[k] for k in FatJetKeys},nMedBJets=self.nMediumBTaggedSubJets, HLL=self.HLL))
                # MET #
                #plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))
                # HighLevel #
                ##plots.extend(makeHighLevelPlotsBoosted(**{k:channelDict[k] for k in BoostedKeys}, HLL=self.HLL))
            '''

        # ML
        self.nodes = ['GGF','VBF','TT','ST','WJets','H','Other']

        leptonCont = {'El':ElColl[0],'Mu':MuColl[0]}
        for selObjectDict in selObjectDictList:
            channel      = selObjectDict['channel']
            lepton       = leptonCont[channel]
            vbf          = selObjectDict['VBFJets']
            inputsLeps   = mvaEvaluatorSL_nonres.returnLeptonsMVAInputs   (self  = self, lep  = lepton)
            inputsJets   = mvaEvaluatorSL_nonres.returnJetsMVAInputs      (self  = self, bjets = self.bJetsByScore, jets = self.probableWJets)
            inputsMET    = mvaEvaluatorSL_nonres.returnMETMVAInputs       (self  = self, met  = self.corrMET)
            inputsFatjet = mvaEvaluatorSL_nonres.returnFatjetMVAInputs    (self  = self, fatjets = self.ak8Jets)
            inputsHL     = mvaEvaluatorSL_nonres.returnHighLevelMVAInputs (self  = self,
                                                                           lep   = lepton,
                                                                           bjets = self.bJetsByScore,
                                                                           wjets = self.wJetsByPt,
                                                                           VBFJetPairs = vbf,
                                                                           channel   = selObjectDict['channel'])
            inputsParam   = mvaEvaluatorSL_nonres.returnParamMVAInputs    (self)
            inputsEventNr = mvaEvaluatorSL_nonres.returnEventNrMVAInputs  (self,t)

            print ("Lepton variables : %d"%len(inputsLeps))                                 
            print ("Jet variables    : %d"%len(inputsJets))
            print ("Fatjet variables : %d"%len(inputsFatjet))                       
            print ("MET variables    : %d"%len(inputsMET)) 
            print ("HL variables     : %d"%len(inputsHL))
            print ("Param variables  : %d"%len(inputsParam)) 
            print ("Event variables  : %d"%len(inputsEventNr))
            
            from  mvaEvaluatorSL_nonres import inputStaticCast

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
                
            plots.extend(makeDoubleLeptonMachineLearningExclusiveOutputPlots(selObjNodesDict,output,self.nodes,channel=selObjectDict['channel']))
            if self.args.PrintYield:
                for selNode in selObjNodesDict.values():
                    self.yields.add(selNode.sel)
            
        #----- Add the Yield plots -----#
        if self.args.PrintYield or self.args.OnlyYield:
            plots.append(self.yields)
            #plots.append(yields)

        return plots

    ### PostProcess ###
    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None):
        super(PlotterNanoHHtobbWWSL, self).postProcess(taskList, config, workdir, resultsdir, forSkimmer=False)
