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
from JPA import *
import mvaEvaluatorSL_nonres
from DDHelper import DataDrivenFake, DataDrivenDY
from bamboo.root import gbl
import ROOT
from functools import partial

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
        cutFlowPlots = []
        
        #----- Machine Learning Model -----#                
        model_num = "12"
        path_model = os.path.join(os.path.abspath(os.path.dirname(__file__)),'MachineLearning','ml-models','multi-classification','dnn',model_num,'model','model.pb')
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
        
        #----- Ratio reweighting variables (before lepton and jet selection) -----#
        #if self.args.BtagReweightingOff or self.args.BtagReweightingOn:
        #    #plots.append(objectsNumberPlot(channel="NoChannel",suffix='NoSelection',sel=noSel,objCont=self.ak4Jets,objName='Ak4Jets',Nmax=15,xTitle='N(Ak4 jets)'))
        #    #plots.append(CutFlowReport("BtagReweightingCutFlowReport",noSel))
        #    return plots
            
        #----- Stitching study -----#
        if self.args.DYStitchingPlots or self.args.WJetsStitchingPlots:
            if self.args.DYStitchingPlots and sampleCfg['group'] != 'DY':
                raise RuntimeError("Stitching is only done on DY MC samples")
            if self.args.WJetsStitchingPlots and sampleCfg['group'] != 'Wjets':
                raise RuntimeError("Stitching is only done on WJets MC samples")
            #plots.extend(makeLHEPlots(noSel,t.LHE))
            #plots.append(objectsNumberPlot(channel="NoChannel",suffix='NoSelection',sel=noSel,objCont=self.ak4Jets,objName='Ak4Jets',Nmax=15,xTitle='N(Ak4 jets)'))
            #plots.append(CutFlowReport("DYStitchingCutFlowReport",noSel))
            return plots
            
            
        #----- Singleleptons -----#
        ElSelObj,MuSelObj = makeSingleLeptonSelection(self,noSel,plot_yield=True)

        #----- Apply jet corrections -----#
        ElSelObj.sel = self.beforeJetselection(ElSelObj.sel,'El')
        MuSelObj.sel = self.beforeJetselection(MuSelObj.sel,'Mu')

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
            
        #inputsEventNr    = returnEventNr(self, t)

        selObjectDictList = []
        if any(item in resolved_args for item in jetsel_level):
            ChannelDictListR = []
            logger.info ("... Processing Ak4Jets Selection for Resolved category : nAk4Jets >= 3 + nAk4BJets >= 1 + nAk8BJets == 0")
            ElSelObjResolved = makeResolvedSelection(self,ElSelObj,copy_sel=True,plot_yield=True)
            MuSelObjResolved = makeResolvedSelection(self,MuSelObj,copy_sel=True,plot_yield=True)

            VBFJetPairsResolved      = mvaEvaluatorSL_nonres.VBFJetPairs_Resolved(self)

            if "Resolved2Btag" in jetplot_level:
                logger.info ('Resolved2Btag Node Selection')
                ElSelObjResolved2b        = makeExclusiveResolvedSelection(self, ElSelObjResolved, nbJet=2, copy_sel=True)
                MuSelObjResolved2b        = makeExclusiveResolvedSelection(self, MuSelObjResolved, nbJet=2, copy_sel=True)

                selObjectDictList.append({'channel':'El','selObject':ElSelObjResolved2b,'category':'Resolved','VBFJets':VBFJetPairsResolved})
                selObjectDictList.append({'channel':'Mu','selObject':MuSelObjResolved2b,'category':'Resolved','VBFJets':VBFJetPairsResolved})

                if self.args.onlypost:
                    ElSelObjResolved2b.record_yields = True
                    MuSelObjResolved2b.record_yields = True
                    ElSelObjResolved2b.yieldTitle = 'Resolved2b Channel $e^{\pm}$'
                    MuSelObjResolved2b.yieldTitle = 'Resolved2b Channel $\mu^{\pm}$'
                
                if not self.args.OnlyYield:
                    ChannelDictListR.append({'channel':'El','selObj':ElSelObjResolved2b,'sel':ElSelObjResolved2b.sel,
                                             'lepton':ElColl[0],'jets':self.ak4Jets,
                                             'bjets':self.bJetsByScore,'wjets':self.wJetsByPt,
                                             'nbJet':2,'suffix':ElSelObjResolved2b.selName,
                                             'is_MC':self.is_MC})
                    ChannelDictListR.append({'channel':'Mu','selObj':MuSelObjResolved2b,'sel':MuSelObjResolved2b.sel,
                                             'lepton':MuColl[0],'jets':self.ak4Jets,
                                             'bjets':self.bJetsByScore,'wjets':self.wJetsByPt,
                                             'nbJet':2,'suffix':MuSelObjResolved2b.selName,
                                             'is_MC':self.is_MC})

            if "Resolved1Btag" in jetplot_level:
                logger.info ('Resolved1Btag Node Selection')
                ElSelObjResolved1b        = makeExclusiveResolvedSelection(self, ElSelObjResolved, nbJet=1, copy_sel=True)
                MuSelObjResolved1b        = makeExclusiveResolvedSelection(self, MuSelObjResolved, nbJet=1, copy_sel=True)

                selObjectDictList.append({'channel':'El','selObject':ElSelObjResolved1b,'category':'Resolved','VBFJets':VBFJetPairsResolved})
                selObjectDictList.append({'channel':'Mu','selObject':MuSelObjResolved1b,'category':'Resolved','VBFJets':VBFJetPairsResolved})

                if self.args.onlypost:
                    ElSelObjResolved1b.record_yields = True
                    MuSelObjResolved1b.record_yields = True
                    ElSelObjResolved1b.yieldTitle = 'Resolved1b Channel $e^{\pm}$'
                    MuSelObjResolved1b.yieldTitle = 'Resolved1b Channel $\mu^{\pm}$'
                
                if not self.args.OnlyYield:
                    ChannelDictListR.append({'channel':'El','selObj':ElSelObjResolved1b,'sel':ElSelObjResolved1b.sel,
                                             'lepton':ElColl[0],'jets':self.ak4Jets,
                                             'bjets':self.bJetsByScore,'wjets':self.wJetsByPt,
                                             'nbJet':1,'suffix':ElSelObjResolved1b.selName,
                                             'is_MC':self.is_MC})
                    ChannelDictListR.append({'channel':'Mu','selObj':MuSelObjResolved1b,'sel':MuSelObjResolved1b.sel,
                                             'lepton':MuColl[0],'jets':self.ak4Jets,
                                             'bjets':self.bJetsByScore,'wjets':self.wJetsByPt,
                                             'nbJet':1,'suffix':MuSelObjResolved1b.selName,
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

        # ========================== JPA Boosted Categories ========================= #
        if any(item in boosted_args for item in jetsel_level):
            ChannelDictListB = []
            FatJetKeys  = ['channel','sel','jet1','jet2','jet3','jet4','has1fat1slim','has1fat2slim','suffix']

            logger.info ("...... Processing Boosted Category : nAk8BJets >= 1, nAk4JetsCleanedFromAk8b >= 1")
            ElSelObjBoosted = makeBoostedSelection(self,ElSelObj,copy_sel=True,plot_yield=True)
            MuSelObjBoosted = makeBoostedSelection(self,MuSelObj,copy_sel=True,plot_yield=True)

            VBFJetPairsBoosted      = mvaEvaluatorSL_nonres.VBFJetPairs_Boosted(self)

            selObjectDictList.append({'channel':'El','selObject':ElSelObjBoosted,'category':'Boosted','VBFJets':VBFJetPairsBoosted})
            selObjectDictList.append({'channel':'Mu','selObject':MuSelObjBoosted,'category':'Boosted','VBFJets':VBFJetPairsBoosted})

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
                

                inputsCommon  = returnCommonInputs_Boosted(self)
                inputsClassic = returnClassicInputs_Boosted(self, **{k:channelDict[k] for k in ClassicInputKeys})
                inputsLBN     = returnLBNInputs_Boosted(self, **{k:channelDict[k] for k in LBNInputKeys})
                input_names   = [key[0] for key in inputsClassic.keys()] + ['LBN_inputs','eventnr']
                inputs_array  = [op.array("double",val) for val in inputStaticCast(inputsCommon,"float")]
                inputs_array.append(op.array("double",*inputStaticCast(inputsClassic,"float")))
                inputs_array.append(op.array("double",*inputStaticCast(inputsLBN,"float")))
                inputs_array.append(op.array("long",*inputStaticCast(inputsEventNr,"long")))
                
                DNN_SM  = op.mvaEvaluator (path_model_boosted_SM, mvaType='Tensorflow',otherArgs=(input_names, output_name))
                DNN_BSM = op.mvaEvaluator(path_model_boosted_BSM, mvaType='Tensorflow',otherArgs=(input_names, output_name))
                DNN_Score_SM  = DNN_SM (*inputs_array)
                DNN_Score_BSM = DNN_BSM (*inputs_array)
                
                selObjNodesDict_SM     = makeDNNOutputNodesSelections(self,channelDict['selObj'],DNN_Score_SM,suffix='_SM_')
                selObjNodesDict_BSM    = makeDNNOutputNodesSelections(self,channelDict['selObj'],DNN_Score_BSM,suffix='_BSM_')
                logger.info('Filling DNN responses Boosted')
                plots.extend(makeDoubleLeptonMachineLearningOutputPlots(selObjNodesDict_SM,DNN_Score_SM,self.nodes,channel=selObjectDNNDict['channel']))
                plots.extend(makeDoubleLeptonMachineLearningOutputPlots(selObjNodesDict_BSM,DNN_Score_BSM,self.nodes,channel=selObjectDNNDict['channel']))
            '''

        # ML
        self.nodes = ['GGF','VBF','TT','ST','WJets','H','Other']

        leptonCont = {'El':ElColl[0],'Mu':MuColl[0]}
        for selObjectDict in selObjectDictList:
            channel = selObjectDict['channel']
            lepton  = leptonCont[channel]
            vbf     = selObjectDict['VBFJets']
            inputsLeps = mvaEvaluatorSL_nonres.returnLeptonsMVAInputs (self  = self, lep  = lepton)
            inputsJets = mvaEvaluatorSL_nonres.returnJetsMVAInputs    (self  = self, jets = self.ak4Jets)
            inputsMET = mvaEvaluatorSL_nonres.returnMETMVAInputs      (self  = self, met  = self.corrMET)
            inputsFatjet = mvaEvaluatorSL_nonres.returnFatjetMVAInputs(self  = self, fatjets = self.ak8Jets)
            inputsHL = mvaEvaluatorSL_nonres.returnHighLevelMVAInputs (self  = self,
                                                                       lep   = lepton,
                                                                       bjets = self.bJetsByScore,
                                                                       wjets = self.wJetsByPt,
                                                                       VBFJetPairs = vbf,
                                                                       channel   = selObjectDict['channel'])
            inputsParam   = mvaEvaluatorSL_nonres.returnParamMVAInputs(self)
            inputsEventNr = mvaEvaluatorSL_nonres.returnEventNrMVAInputs(self,t)
            
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
            

        #----- Add the Yield plots -----#
        plots.append(self.yields)
        #plots.extend(cutFlowPlots)
        return plots

