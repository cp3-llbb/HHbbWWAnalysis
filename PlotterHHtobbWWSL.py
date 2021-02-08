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
from mvaEvaluatorSL import *
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
        # --------------------------- Machine Learning Model --------------------------- #
        # -------- for JPA --------- #
        era = sampleCfg['era']
        basepath = os.path.join(os.path.abspath(os.path.dirname(__file__)),'MachineLearning','ml-models','JPA_Loose_ttH')
        resolvedModelDict = dict()
        resolvedModelDict['2b2Wj'] = [os.path.join(basepath, 'bb1l_jpa_4jet_resolved_even_all_mass_all_node_'+str(era)+'.xml'), 
                                      os.path.join(basepath, 'bb1l_jpa_4jet_resolved_odd_all_mass_all_node_'+str(era)+'.xml')]
        resolvedModelDict['2b1Wj'] = [os.path.join(basepath, 'bb1l_jpa_missingWJet_resolved_even_all_mass_all_node_'+str(era)+'.xml'), 
                                      os.path.join(basepath, 'bb1l_jpa_missingWJet_resolved_odd_all_mass_all_node_'+str(era)+'.xml')]
        resolvedModelDict['2b0Wj'] = [os.path.join(basepath, 'bb1l_jpa_missingAllWJet_even_all_mass_all_node_'+str(era)+'.xml'), 
                                      os.path.join(basepath, 'bb1l_jpa_missingAllWJet_odd_all_mass_all_node_'+str(era)+'.xml')]
        resolvedModelDict['1b2Wj'] = [os.path.join(basepath, 'bb1l_jpa_missingBJet_even_all_mass_all_node_'+str(era)+'.xml'), 
                                      os.path.join(basepath, 'bb1l_jpa_missingBJet_odd_all_mass_all_node_'+str(era)+'.xml')]
        resolvedModelDict['1b1Wj'] = [os.path.join(basepath, 'bb1l_jpa_missingBJet_missingWJet_even_all_mass_all_node_'+str(era)+'.xml'), 
                                      os.path.join(basepath, 'bb1l_jpa_missingBJet_missingWJet_odd_all_mass_all_node_'+str(era)+'.xml')]
        resolvedModelDict['1b0Wj'] = [os.path.join(basepath, 'bb1l_jpa_missingBJet_missingAllWJet_even_all_mass_all_node_'+str(era)+'.xml'), 
                                      os.path.join(basepath, 'bb1l_jpa_missingBJet_missingAllWJet_odd_all_mass_all_node_'+str(era)+'.xml')]
        resolvedModelDict['evCat'] = [os.path.join(basepath, 'bb1l_jpa_evtCat_resolved_even_all_mass_all_node_'+str(era)+'.xml'), 
                                      os.path.join(basepath, 'bb1l_jpa_evtCat_resolved_odd_all_mass_all_node_'+str(era)+'.xml')]
        
        boostedModelDict  = dict()
        boostedModelDict['Hbb2Wj'] = [os.path.join(basepath, 'bb1l_jpa_4jet_boosted_even_all_mass_all_node_'+str(era)+'.xml'), 
                                      os.path.join(basepath, 'bb1l_jpa_4jet_boosted_odd_all_mass_all_node_'+str(era)+'.xml')]
        boostedModelDict['Hbb1Wj'] = [os.path.join(basepath, 'bb1l_jpa_missingWJet_boosted_even_all_mass_all_node_'+str(era)+'.xml'), 
                                      os.path.join(basepath, 'bb1l_jpa_missingWJet_boosted_odd_all_mass_all_node_'+str(era)+'.xml')]
        boostedModelDict['evCat']  = [os.path.join(basepath, 'bb1l_jpa_evtCat_boosted_even_all_mass_all_node_'+str(era)+'.xml'), 
                                      os.path.join(basepath, 'bb1l_jpa_evtCat_boosted_odd_all_mass_all_node_'+str(era)+'.xml')]
        
        # ----------------------------------- NodeList ------------------------------------- #
        # keep the exact same order of nodes as mentioned in respective xml files
        ResolvedJPANodeList = ['2b2Wj','2b1Wj','1b2Wj','2b0Wj','1b1Wj','1b0Wj','0b']
        BoostedJPANodeList  = ['Hbb2Wj','Hbb1Wj','Hbb0Wj']

        
        # ---------- LBN+DNN models ----------- #
        #path_model_resolved = '/home/ucl/cp3/gsaha/bamboodev/HHbbWWAnalysis/MachineLearning/ml-models/DNN/Resolved/ModelOneEpoch.pb'
        #path_model_resolved = os.path.join('/home/ucl/cp3/gsaha/bamboodev/HHbbWWAnalysis/MachineLearning/ml-models/DNN/Resolved','Model'+str(era)+'_5x256_crossval.pb')
        path_model_resolved = os.path.join('/home/ucl/cp3/gsaha/bamboodev/HHbbWWAnalysis/MachineLearning/ml-models/DNN/Resolved','Resolved'+str(era)+'.pb')
        path_model_boosted  = os.path.join('/home/ucl/cp3/gsaha/bamboodev/HHbbWWAnalysis/MachineLearning/ml-models/DNN/Resolved','Boosted_allEras_v2.pb')
        logger.info('DNN_Model_Resolved {}'.format(path_model_resolved))
        logger.info('DNN_Model_Boosted {}'.format(path_model_boosted))
        plots = []
        cutFlowPlots = []
        
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

        # selObjectDict : keys -> level (str)
        #                 values -> [El,Mu] x Selection object
        # Select the jets selections that will be done depending on user input #
        resolved_args = ["Res2b2Wj","Res2b1Wj","Res2b0Wj","Res1b2Wj","Res1b1Wj","Res1b0Wj","Res0b","Resolved"]
        boosted_args  = ["Hbb2Wj","Hbb1Wj","Hbb0Wj","Boosted"]
        jet_level = resolved_args + boosted_args
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
            
        def cleanVBFwithJPA_Resolved(jpaJets, nJpaJets):
            #return lambda j : op.OR(*(op.deltaR(jpaJets[i].p4, j.p4) > 0.8 for i in range(nJpaJets)))
            return lambda j : op.AND(*(op.deltaR(jpaJets[i].p4, j.p4) > 0.8 for i in range(nJpaJets)))
        def cleanVBFwithJPA_Boosted(jpaJets, nJpaJets):
            #return lambda j : op.AND(op.rng_len(self.ak8BJets) >= 1, op.AND(op.OR(*(op.deltaR(jpaJets[i].p4, j.p4) > 0.8 for i in range(nJpaJets))), op.deltaR(self.ak8Jets[0].p4, j.p4) > 1.2))
            return lambda j : op.AND(op.rng_len(self.ak8BJets) >= 1, op.AND(op.AND(*(op.deltaR(jpaJets[i].p4, j.p4) > 0.8 for i in range(nJpaJets))), op.deltaR(self.ak8Jets[0].p4, j.p4) > 1.2))
        def cleanVBFwithJPA_Boosted_rest():
            return lambda j : op.AND(op.rng_len(self.ak8BJets) >= 1, op.deltaR(self.ak8Jets[0].p4, j.p4) > 1.2)

        # Selections:    
        #---- Lepton selection ----#
        ElColl = [t.Electron[op.switch(op.rng_len(self.electronsTightSel) == 1, 
                                       self.electronsTightSel[0].idx, 
                                       self.electronsFakeSel[0].idx)]]
        MuColl = [t.Muon[op.switch(op.rng_len(self.muonsTightSel) == 1, 
                                   self.muonsTightSel[0].idx, 
                                   self.muonsFakeSel[0].idx)]]

        if not self.args.OnlyYield:
            ChannelDictList = []
            ChannelDictList.append({'channel':'El','sel':ElSelObj.sel,'suffix':ElSelObj.selName})
            ChannelDictList.append({'channel':'Mu','sel':MuSelObj.sel,'suffix':MuSelObj.selName})
                
            #for channelDict in ChannelDictList:
            #    #----- Trigger plots -----#
            #    plots.extend(singleLeptonTriggerPlots(**channelDict, triggerDict=self.triggersPerPrimaryDataset))
            
        LeptonKeys   = ['channel','sel','lep','suffix','is_MC']
        JetKeys      = ['channel','sel','j1','j2','j3','j4','suffix','nJet','nbJet','is_MC']
        commonItems  = ['channel','sel','suffix']
            
        #----- Ak4 jets selection -----#
        if "Ak4" in jetsel_level:
            logger.info ("... Processing Ak4Jets Selection for Resolved category : nAk4Jets >= 3")
            ResolvedKeys = ['channel','sel','met','lep','j1','j2','j3','j4','suffix','nJet','nbJet']
            JetsN        = {'objName':'Ak4Jets','objCont':self.ak4Jets,'Nmax':10,'xTitle':'nAk4Jets'}
            FatJetsN     = {'objName':'Ak8Jets','objCont':self.ak8Jets,'Nmax':5,'xTitle':'nAk8Jets'}
            
            ElSelObjResolved = makeResolvedSelection(self,ElSelObj,copy_sel=True,plot_yield=True)
            MuSelObjResolved = makeResolvedSelection(self,MuSelObj,copy_sel=True,plot_yield=True)
        
            if self.args.onlypost:
                ElSelObjResolved.record_yields = True
                MuSelObjResolved.record_yields = True
                ElSelObjResolved.yieldTitle = 'Resolved Channel $e^{\pm}$'
                MuSelObjResolved.yieldTitle = 'Resolved Channel $\mu^{\pm}$'
                
            ChannelDictListAk4 = []
            if "Ak4" in jetplot_level and not self.args.OnlyYield:
                logger.info ('...... Base Resolved Selection : nAk4Jets >= 3')
                ChannelDictListAk4.append({'channel':'El','sel':ElSelObjResolved.sel,'lep':ElColl[0],'met':self.corrMET,
                                           'j1':self.ak4Jets[0],'j2':self.ak4Jets[1],'j3':self.ak4Jets[2],'j4':None,
                                           'nJet':3,'nbJet':0,'suffix':ElSelObjResolved.selName,'is_MC':self.is_MC})
                ChannelDictListAk4.append({'channel':'Mu','sel':MuSelObjResolved.sel,'lep':MuColl[0],'met':self.corrMET,
                                           'j1':self.ak4Jets[0],'j2':self.ak4Jets[1],'j3':self.ak4Jets[2],'j4':None,
                                           'nJet':3,'nbJet':0,'suffix':MuSelObjResolved.selName,'is_MC':self.is_MC})

                for channelDict in ChannelDictListAk4:
                    # Singlelepton #
                    plots.extend(makeSinleptonPlots(**{k:channelDict[k] for k in LeptonKeys}))
                    # Number of jets #
                    plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**JetsN))
                    plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**FatJetsN))
                    # Ak4 Jets #
                    plots.extend(makeAk4JetsPlots(**{k:channelDict[k] for k in JetKeys},HLL=self.HLL))
                    # MET #
                    plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))
                    # High level #
                    plots.extend(makeHighLevelPlotsResolved(**{k:channelDict[k] for k in ResolvedKeys},HLL=self.HLL))

                
        ##### Ak8-b jets selection #####
        if "Ak8" in jetsel_level:
            logger.info ("...... Processing Ak8b jet selection for SemiBoosted & Boosted Category")
            FatJetKeys  = ['channel','sel','j1','j2','j3','has1fat1slim','has1fat2slim','suffix']
            FatJetsN    = {'objName':'Ak8Jets','objCont':self.ak8Jets,'Nmax':5,'xTitle':'N(Ak8 jets)'}
            SlimJetsN   = {'objName':'Ak4Jets','objCont':self.ak4JetsCleanedFromAk8b,'Nmax':10,'xTitle':'N(Ak4 jets)'}
            BoostedKeys = ['channel','sel','met','lep','j1','j2','j3','suffix','bothAreFat','has1fat2slim']
            
            ElSelObjBoosted = makeBoostedSelection(self,ElSelObj,copy_sel=True,plot_yield=True)
            MuSelObjBoosted = makeBoostedSelection(self,MuSelObj,copy_sel=True,plot_yield=True)
            
            if self.args.onlypost:
                ElSelObjBoosted.record_yields = True
                MuSelObjBoosted.record_yields = True
                ElSelObjBoosted.yieldTitle = 'Boosted Channel $e^{\pm}$'
                MuSelObjBoosted.yieldTitle = 'Boosted Channel $\mu^{\pm}$'
                
            # Fatjets plots #
            ChannelDictListAk8 = []
            if "Ak8" in jetplot_level and not self.args.OnlyYield:
                logger.info ('boosted selection...')
                ChannelDictListAk8.append({'channel':'El','sel':ElSelObjBoosted.sel,'lep':ElColl[0],'met':self.corrMET,
                                           'j1':self.ak8BJets[0],'j2':self.ak4JetsCleanedFromAk8b[0],'j3':None,'has1fat1slim':True,'has1fat2slim':False,'bothAreFat':False,
                                           'suffix':ElSelObjBoosted.selName,'is_MC':self.is_MC})
                ChannelDictListAk8.append({'channel':'Mu','sel':MuSelObjBoosted.sel,'lep':MuColl[0],'met':self.corrMET,
                                           'j1':self.ak8BJets[0],'j2':self.ak4JetsCleanedFromAk8b[0],'j3':None,'has1fat1slim':True,'has1fat2slim':False,'bothAreFat':False,
                                           'suffix':MuSelObjBoosted.selName,'is_MC':self.is_MC})

                for channelDict in ChannelDictListAk8:
                    # Dilepton #
                    plots.extend(makeSinleptonPlots(**{k:channelDict[k] for k in LeptonKeys}))
                    # Number of jets #
                    plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**FatJetsN))
                    plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**SlimJetsN))
                    # Ak8 Jets #
                    plots.extend(makeSingleLeptonAk8JetsPlots(**{k:channelDict[k] for k in FatJetKeys},nMedBJets=self.nMediumBTaggedSubJets, HLL=self.HLL))
                    # MET #
                    plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))
                    # HighLevel #
                    plots.extend(makeHighLevelPlotsBoosted(**{k:channelDict[k] for k in BoostedKeys}, HLL=self.HLL))

        
        leptonCont = {'El':ElColl[0],'Mu':MuColl[0]}
        #self.nodes = ['DY','GGF','H','Rare','ST','TT','VBF','WJets']
        #self.nodes = ['GGF','H','Rare','ST','TT','VBF','WJets']
        self.nodes = ['Ewk','GGF','H','Top','VBF','WJets']
        inputsEventNr    = returnEventNr(self, t)
        output_name = "Identity"
        # ========================== JPA Resolved Categories ========================= #
        if any(item in resolved_args for item in jetsel_level):
            ChannelDictList = []
            selObjectDNNDictListResolved = []
            # dict = {'key':'Node', 'value' : [refined selObj, [JPAjetIndices]]}
            elL1OutList, elL2OutList, ElResolvedSelObjJetsIdxPerJpaNodeDict = findJPACategoryResolved (self, ElSelObjResolved, ElColl[0],self.muonsPreSel, self.electronsPreSel, 
                                                                                                       self.ak4Jets, self.ak4BJetsLoose,self.ak4BJets, self.corrMET, 
                                                                                                       resolvedModelDict, t.event,self.HLL, ResolvedJPANodeList, 
                                                                                                       plot_yield=True)
            muL1OutList, muL2OutList, MuResolvedSelObjJetsIdxPerJpaNodeDict = findJPACategoryResolved (self, MuSelObjResolved, MuColl[0],self.muonsPreSel, self.electronsPreSel, 
                                                                                                       self.ak4Jets, self.ak4BJetsLoose,self.ak4BJets, self.corrMET, 
                                                                                                       resolvedModelDict, t.event,self.HLL, ResolvedJPANodeList,
                                                                                                       plot_yield=True)
            if "Res2b2Wj" in jetplot_level or "Resolved" in jetplot_level:
                logger.info ('...... JPA : 2b2Wj Node Selection')
                ElSelObjResolved2b2Wj        = ElResolvedSelObjJetsIdxPerJpaNodeDict.get('2b2Wj')[0]
                ElSelObjResolved2b2WjJets    = ElResolvedSelObjJetsIdxPerJpaNodeDict.get('2b2Wj')[1]
                MuSelObjResolved2b2Wj        = MuResolvedSelObjJetsIdxPerJpaNodeDict.get('2b2Wj')[0]
                MuSelObjResolved2b2WjJets    = MuResolvedSelObjJetsIdxPerJpaNodeDict.get('2b2Wj')[1]

                ElSelObjResolved2b2WjVBFJets = op.sort(op.combine(op.select(self.VBFJets, cleanVBFwithJPA_Resolved(ElSelObjResolved2b2WjJets, 4)), N=2, pred=self.lambda_VBFPair),
                                                       lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))
                MuSelObjResolved2b2WjVBFJets = op.sort(op.combine(op.select(self.VBFJets, cleanVBFwithJPA_Resolved(MuSelObjResolved2b2WjJets, 4)), N=2, pred=self.lambda_VBFPair),
                                                       lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))

                if self.args.onlypost:
                    ElSelObjResolved2b2Wj.record_yields = True
                    MuSelObjResolved2b2Wj.record_yields = True
                    ElSelObjResolved2b2Wj.yieldTitle = 'Resolved2b2Wj Channel $e^{\pm}$'
                    MuSelObjResolved2b2Wj.yieldTitle = 'Resolved2b2Wj Channel $\mu^{\pm}$'
                

                inputsClassic_2b2Wj_El = returnClassicInputs_2b2Wj(self, ElColl[0], ElSelObjResolved2b2WjJets[0], ElSelObjResolved2b2WjJets[1], ElSelObjResolved2b2WjJets[2], 
                                                                   ElSelObjResolved2b2WjJets[3], ElSelObjResolved2b2WjVBFJets, elL2OutList)
                inputsLBN_2b2Wj_El     = returnLBNInputs_2b2Wj(self, ElColl[0], ElSelObjResolved2b2WjJets[0], ElSelObjResolved2b2WjJets[1], ElSelObjResolved2b2WjJets[2], 
                                                               ElSelObjResolved2b2WjJets[3])
                inputsClassic_2b2Wj_Mu = returnClassicInputs_2b2Wj(self, MuColl[0], MuSelObjResolved2b2WjJets[0], MuSelObjResolved2b2WjJets[1], MuSelObjResolved2b2WjJets[2], 
                                                                   MuSelObjResolved2b2WjJets[3], MuSelObjResolved2b2WjVBFJets, muL2OutList)
                inputsLBN_2b2Wj_Mu     = returnLBNInputs_2b2Wj(self, MuColl[0], MuSelObjResolved2b2WjJets[0], MuSelObjResolved2b2WjJets[1], MuSelObjResolved2b2WjJets[2], 
                                                               MuSelObjResolved2b2WjJets[3])
                input_names_2b2Wj      = [key[0] for key in inputsClassic_2b2Wj_El.keys()] + ['LBN_inputs','eventnr']
                DNN_2b2Wj   = op.mvaEvaluator(path_model_resolved,mvaType='Tensorflow',otherArgs=(input_names_2b2Wj, output_name))
                inputs_array_2b2Wj_El = [op.array("double",val) for val in inputStaticCast(inputsClassic_2b2Wj_El,"float")]
                inputs_array_2b2Wj_El.append(op.array("double",*inputStaticCast(inputsLBN_2b2Wj_El,"float")))
                inputs_array_2b2Wj_El.append(op.array("long",*inputStaticCast(inputsEventNr,"long")))
                inputs_array_2b2Wj_Mu = [op.array("double",val) for val in inputStaticCast(inputsClassic_2b2Wj_Mu,"float")]
                inputs_array_2b2Wj_Mu.append(op.array("double",*inputStaticCast(inputsLBN_2b2Wj_Mu,"float")))
                inputs_array_2b2Wj_Mu.append(op.array("long",*inputStaticCast(inputsEventNr,"long")))
                DNNScore_2b2Wj_El = DNN_2b2Wj(*inputs_array_2b2Wj_El)
                DNNScore_2b2Wj_Mu = DNN_2b2Wj(*inputs_array_2b2Wj_Mu)

                selObjectDNNDictListResolved.append({'channel':'El','selObject':ElSelObjResolved2b2Wj,'DNN_Inputs':inputs_array_2b2Wj_El,'DNN_Score':DNNScore_2b2Wj_El})
                selObjectDNNDictListResolved.append({'channel':'Mu','selObject':MuSelObjResolved2b2Wj,'DNN_Inputs':inputs_array_2b2Wj_Mu,'DNN_Score':DNNScore_2b2Wj_Mu})

                #selObjectDNNDictListResolved.append({'channel':'El','selObject':ElSelObjResolved2b2Wj,'DNN_Inputs':inputs_array_2b2Wj_El,'DNN_Score':DNNScore_2b2Wj_El})
                #selObjectDNNDictListResolved.append({'channel':'Mu','selObject':MuSelObjResolved2b2Wj,'DNN_Inputs':inputs_array_2b2Wj_Mu,'DNN_Score':DNNScore_2b2Wj_Mu})

                plots.append(Plot.make1D("El_2b2WjDNNMaxScore_%s"%ElSelObjResolved2b2Wj.selName, op.rng_max_element_by(DNNScore_2b2Wj_El), ElSelObjResolved2b2Wj.sel, EquidistantBinning(400, 0.0, 1.0)))
                plots.append(Plot.make1D("Mu_2b2WjDNNMaxScore_%s"%MuSelObjResolved2b2Wj.selName, op.rng_max_element_by(DNNScore_2b2Wj_Mu), MuSelObjResolved2b2Wj.sel, EquidistantBinning(400, 0.0, 1.0)))
                plots.append(Plot.make1D("El_2b2WjnDNNMaxScores_%s"%ElSelObjResolved2b2Wj.selName, op.rng_len(DNNScore_2b2Wj_El), ElSelObjResolved2b2Wj.sel, EquidistantBinning(10, 0.0, 10.)))
                plots.append(Plot.make1D("Mu_2b2WjnDNNMaxScores_%s"%MuSelObjResolved2b2Wj.selName, op.rng_len(DNNScore_2b2Wj_Mu), MuSelObjResolved2b2Wj.sel, EquidistantBinning(10, 0.0, 10.)))
                
                if not self.args.OnlyYield:
                    ChannelDictList.append({'channel':'El','sel':ElSelObjResolved2b2Wj.sel,'lep':ElColl[0],'met':self.corrMET,
                                            'j1':ElSelObjResolved2b2WjJets[0],'j2':ElSelObjResolved2b2WjJets[1],
                                            'j3':ElSelObjResolved2b2WjJets[2],'j4':ElSelObjResolved2b2WjJets[3],
                                            'nJet':4,'nbJet':2,'suffix':ElSelObjResolved2b2Wj.selName,
                                            'is_MC':self.is_MC})
                    ChannelDictList.append({'channel':'Mu','sel':MuSelObjResolved2b2Wj.sel,'lep':MuColl[0],'met':self.corrMET,
                                            'j1':MuSelObjResolved2b2WjJets[0],'j2':MuSelObjResolved2b2WjJets[1],
                                            'j3':MuSelObjResolved2b2WjJets[2],'j4':MuSelObjResolved2b2WjJets[3],
                                            'nJet':4,'nbJet':2,'suffix':MuSelObjResolved2b2Wj.selName,
                                            'is_MC':self.is_MC})

            if "Res2b1Wj" in jetplot_level or "Resolved" in jetplot_level:
                logger.info ('...... JPA : 2b1Wj Node Selection')
                ElSelObjResolved2b1Wj        = ElResolvedSelObjJetsIdxPerJpaNodeDict.get('2b1Wj')[0]
                ElSelObjResolved2b1WjJets    = ElResolvedSelObjJetsIdxPerJpaNodeDict.get('2b1Wj')[1]
                MuSelObjResolved2b1Wj        = MuResolvedSelObjJetsIdxPerJpaNodeDict.get('2b1Wj')[0]
                MuSelObjResolved2b1WjJets    = MuResolvedSelObjJetsIdxPerJpaNodeDict.get('2b1Wj')[1]

                if self.args.onlypost:
                    ElSelObjResolved2b1Wj.record_yields = True
                    MuSelObjResolved2b1Wj.record_yields = True
                    ElSelObjResolved2b1Wj.yieldTitle = 'Resolved2b1Wj Channel $e^{\pm}$'
                    MuSelObjResolved2b1Wj.yieldTitle = 'Resolved2b1Wj Channel $\mu^{\pm}$'
                
                ElSelObjResolved2b1WjVBFJets = op.sort(op.combine(op.select(self.VBFJets, cleanVBFwithJPA_Resolved(ElSelObjResolved2b1WjJets, 3)), N=2, pred=self.lambda_VBFPair),
                                                       lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))
                MuSelObjResolved2b1WjVBFJets = op.sort(op.combine(op.select(self.VBFJets, cleanVBFwithJPA_Resolved(MuSelObjResolved2b1WjJets, 3)), N=2, pred=self.lambda_VBFPair),
                                                       lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))


                inputsClassic_2b1Wj_El = returnClassicInputs_2b1Wj(self, ElColl[0], ElSelObjResolved2b1WjJets[0], ElSelObjResolved2b1WjJets[1], ElSelObjResolved2b1WjJets[2], None, ElSelObjResolved2b1WjVBFJets, elL2OutList)
                inputsLBN_2b1Wj_El     = returnLBNInputs_2b1Wj(self, ElColl[0], ElSelObjResolved2b1WjJets[0], ElSelObjResolved2b1WjJets[1], ElSelObjResolved2b1WjJets[2], None)
                inputsClassic_2b1Wj_Mu = returnClassicInputs_2b1Wj(self, MuColl[0], MuSelObjResolved2b1WjJets[0], MuSelObjResolved2b1WjJets[1], MuSelObjResolved2b1WjJets[2], None, MuSelObjResolved2b1WjVBFJets, muL2OutList)
                inputsLBN_2b1Wj_Mu     = returnLBNInputs_2b1Wj(self, MuColl[0], MuSelObjResolved2b1WjJets[0], MuSelObjResolved2b1WjJets[1], MuSelObjResolved2b1WjJets[2], None)
                input_names_2b1Wj = [key[0] for key in inputsClassic_2b1Wj_El.keys()] + ['LBN_inputs','eventnr']
                DNN_2b1Wj = op.mvaEvaluator(path_model_resolved,mvaType='Tensorflow',otherArgs=(input_names_2b1Wj, output_name))
                inputs_array_2b1Wj_El = [op.array("double",val) for val in inputStaticCast(inputsClassic_2b1Wj_El,"float")]
                inputs_array_2b1Wj_El.append(op.array("double",*inputStaticCast(inputsLBN_2b1Wj_El,"float")))
                inputs_array_2b1Wj_El.append(op.array("long",*inputStaticCast(inputsEventNr,"long")))
                inputs_array_2b1Wj_Mu = [op.array("double",val) for val in inputStaticCast(inputsClassic_2b1Wj_Mu,"float")]
                inputs_array_2b1Wj_Mu.append(op.array("double",*inputStaticCast(inputsLBN_2b1Wj_Mu,"float")))
                inputs_array_2b1Wj_Mu.append(op.array("long",*inputStaticCast(inputsEventNr,"long")))
                DNNScore_2b1Wj_El = DNN_2b1Wj(*inputs_array_2b1Wj_El)
                DNNScore_2b1Wj_Mu = DNN_2b1Wj(*inputs_array_2b1Wj_Mu)

                selObjectDNNDictListResolved.append({'channel':'El','selObject':ElSelObjResolved2b1Wj,'DNN_Inputs':inputs_array_2b1Wj_El,'DNN_Score':DNNScore_2b1Wj_El})
                selObjectDNNDictListResolved.append({'channel':'Mu','selObject':MuSelObjResolved2b1Wj,'DNN_Inputs':inputs_array_2b1Wj_Mu,'DNN_Score':DNNScore_2b1Wj_Mu})


                if not self.args.OnlyYield:
                    ChannelDictList.append({'channel':'El','sel':ElSelObjResolved2b1Wj.sel,'lep':ElColl[0],'met':self.corrMET,
                                            'j1':ElSelObjResolved2b1WjJets[0],'j2':ElSelObjResolved2b1WjJets[1],
                                            'j3':ElSelObjResolved2b1WjJets[2],'j4':None,
                                            'nJet':3,'nbJet':2,'suffix':ElSelObjResolved2b1Wj.selName,
                                            'is_MC':self.is_MC})
                    ChannelDictList.append({'channel':'Mu','sel':MuSelObjResolved2b1Wj.sel,'lep':MuColl[0],'met':self.corrMET,
                                            'j1':MuSelObjResolved2b1WjJets[0],'j2':MuSelObjResolved2b1WjJets[1],
                                            'j3':MuSelObjResolved2b1WjJets[2],'j4':None,
                                            'nJet':3,'nbJet':2,'suffix':MuSelObjResolved2b1Wj.selName,
                                            'is_MC':self.is_MC})
                    
            if "Res2b0Wj" in jetplot_level or "Resolved" in jetplot_level:
                logger.info ('...... JPA : 2b0Wj Node Selection')
                ElSelObjResolved2b0Wj          = ElResolvedSelObjJetsIdxPerJpaNodeDict.get('2b0Wj')[0]
                ElSelObjResolved2b0WjJets      = ElResolvedSelObjJetsIdxPerJpaNodeDict.get('2b0Wj')[1]
                MuSelObjResolved2b0Wj          = MuResolvedSelObjJetsIdxPerJpaNodeDict.get('2b0Wj')[0]
                MuSelObjResolved2b0WjJets      = MuResolvedSelObjJetsIdxPerJpaNodeDict.get('2b0Wj')[1]

                if self.args.onlypost:
                    ElSelObjResolved2b0Wj.record_yields = True
                    MuSelObjResolved2b0Wj.record_yields = True
                    ElSelObjResolved2b0Wj.yieldTitle = 'Resolved2b0Wj Channel $e^{\pm}$'
                    MuSelObjResolved2b0Wj.yieldTitle = 'Resolved2b0Wj Channel $\mu^{\pm}$'                

                ElSelObjResolved2b0WjVBFJets = op.sort(op.combine(op.select(self.VBFJets, cleanVBFwithJPA_Resolved(ElSelObjResolved2b0WjJets, 2)), N=2, pred=self.lambda_VBFPair),
                                                       lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))
                MuSelObjResolved2b0WjVBFJets = op.sort(op.combine(op.select(self.VBFJets, cleanVBFwithJPA_Resolved(MuSelObjResolved2b0WjJets, 2)), N=2, pred=self.lambda_VBFPair),
                                                       lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))


                inputsClassic_2b0Wj_El = returnClassicInputs_2b0Wj(self, ElColl[0], ElSelObjResolved2b0WjJets[0], ElSelObjResolved2b0WjJets[1], None, None, ElSelObjResolved2b0WjVBFJets, elL2OutList)
                inputsLBN_2b0Wj_El     = returnLBNInputs_2b0Wj(self, ElColl[0], ElSelObjResolved2b0WjJets[0], ElSelObjResolved2b0WjJets[1], None, None)
                inputsClassic_2b0Wj_Mu = returnClassicInputs_2b0Wj(self, MuColl[0], MuSelObjResolved2b0WjJets[0], MuSelObjResolved2b0WjJets[1], None, None, MuSelObjResolved2b0WjVBFJets, muL2OutList)
                inputsLBN_2b0Wj_Mu     = returnLBNInputs_2b0Wj(self, MuColl[0], MuSelObjResolved2b0WjJets[0], MuSelObjResolved2b0WjJets[1], None, None)
                input_names_2b0Wj = [key[0] for key in inputsClassic_2b0Wj_El.keys()] + ['LBN_inputs','eventnr']
                DNN_2b0Wj = op.mvaEvaluator(path_model_resolved,mvaType='Tensorflow',otherArgs=(input_names_2b0Wj, output_name))
                inputs_array_2b0Wj_El = [op.array("double",val) for val in inputStaticCast(inputsClassic_2b0Wj_El,"float")]
                inputs_array_2b0Wj_El.append(op.array("double",*inputStaticCast(inputsLBN_2b0Wj_El,"float")))
                inputs_array_2b0Wj_El.append(op.array("long",*inputStaticCast(inputsEventNr,"long")))
                inputs_array_2b0Wj_Mu = [op.array("double",val) for val in inputStaticCast(inputsClassic_2b0Wj_Mu,"float")]
                inputs_array_2b0Wj_Mu.append(op.array("double",*inputStaticCast(inputsLBN_2b0Wj_Mu,"float")))
                inputs_array_2b0Wj_Mu.append(op.array("long",*inputStaticCast(inputsEventNr,"long")))
                DNNScore_2b0Wj_El = DNN_2b0Wj(*inputs_array_2b0Wj_El)
                DNNScore_2b0Wj_Mu = DNN_2b0Wj(*inputs_array_2b0Wj_Mu)

                selObjectDNNDictListResolved.append({'channel':'El','selObject':ElSelObjResolved2b0Wj,'DNN_Inputs':inputs_array_2b0Wj_El,'DNN_Score':DNNScore_2b0Wj_El})
                selObjectDNNDictListResolved.append({'channel':'Mu','selObject':MuSelObjResolved2b0Wj,'DNN_Inputs':inputs_array_2b0Wj_Mu,'DNN_Score':DNNScore_2b0Wj_Mu})


                if not self.args.OnlyYield:
                    ChannelDictList.append({'channel':'El','sel':ElSelObjResolved2b0Wj.sel,'lep':ElColl[0],'met':self.corrMET,
                                            'j1':ElSelObjResolved2b0WjJets[0],'j2':ElSelObjResolved2b0WjJets[1],
                                            'j3':None,'j4':None,
                                            'nJet':2,'nbJet':2,'suffix':ElSelObjResolved2b0Wj.selName,
                                            'is_MC':self.is_MC})
                    ChannelDictList.append({'channel':'Mu','sel':MuSelObjResolved2b0Wj.sel,'lep':MuColl[0],'met':self.corrMET,
                                            'j1':MuSelObjResolved2b0WjJets[0],'j2':MuSelObjResolved2b0WjJets[1],
                                            'j3':None,'j4':None,
                                            'nJet':2,'nbJet':2,'suffix':MuSelObjResolved2b0Wj.selName,
                                            'is_MC':self.is_MC})
                    
            if "Res1b2Wj" in jetplot_level or "Resolved" in jetplot_level:
                logger.info ('...... JPA : 1b2Wj Node Selection')
                ElSelObjResolved1b2Wj        = ElResolvedSelObjJetsIdxPerJpaNodeDict.get('1b2Wj')[0]
                ElSelObjResolved1b2WjJets    = ElResolvedSelObjJetsIdxPerJpaNodeDict.get('1b2Wj')[1]
                MuSelObjResolved1b2Wj        = MuResolvedSelObjJetsIdxPerJpaNodeDict.get('1b2Wj')[0]
                MuSelObjResolved1b2WjJets    = MuResolvedSelObjJetsIdxPerJpaNodeDict.get('1b2Wj')[1]
                
                if self.args.onlypost:
                    ElSelObjResolved1b2Wj.record_yields = True
                    MuSelObjResolved1b2Wj.record_yields = True
                    ElSelObjResolved1b2Wj.yieldTitle = 'Resolved1b2Wj Channel $e^{\pm}$'
                    MuSelObjResolved1b2Wj.yieldTitle = 'Resolved1b2Wj Channel $\mu^{\pm}$'

                ElSelObjResolved1b2WjVBFJets = op.sort(op.combine(op.select(self.VBFJets, cleanVBFwithJPA_Resolved(ElSelObjResolved1b2WjJets, 3)), N=2, pred=self.lambda_VBFPair),
                                                       lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))
                MuSelObjResolved1b2WjVBFJets = op.sort(op.combine(op.select(self.VBFJets, cleanVBFwithJPA_Resolved(MuSelObjResolved1b2WjJets, 3)), N=2, pred=self.lambda_VBFPair),
                                                       lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))

                inputsClassic_1b2Wj_El = returnClassicInputs_1b2Wj(self, ElColl[0], ElSelObjResolved1b2WjJets[0], None, ElSelObjResolved1b2WjJets[1], ElSelObjResolved1b2WjJets[2], ElSelObjResolved1b2WjVBFJets, elL2OutList)
                inputsLBN_1b2Wj_El     = returnLBNInputs_1b2Wj(self, ElColl[0], ElSelObjResolved1b2WjJets[0], None, ElSelObjResolved1b2WjJets[1], ElSelObjResolved1b2WjJets[2])
                inputsClassic_1b2Wj_Mu = returnClassicInputs_1b2Wj(self, MuColl[0], MuSelObjResolved1b2WjJets[0], None, MuSelObjResolved1b2WjJets[1], MuSelObjResolved1b2WjJets[2], MuSelObjResolved1b2WjVBFJets, muL2OutList)
                inputsLBN_1b2Wj_Mu     = returnLBNInputs_1b2Wj(self, MuColl[0], MuSelObjResolved1b2WjJets[0], None, MuSelObjResolved1b2WjJets[1], MuSelObjResolved1b2WjJets[2])
                input_names_1b2Wj = [key[0] for key in inputsClassic_1b2Wj_El.keys()] + ['LBN_inputs','eventnr']
                DNN_1b2Wj = op.mvaEvaluator(path_model_resolved,mvaType='Tensorflow',otherArgs=(input_names_1b2Wj, output_name))
                inputs_array_1b2Wj_El = [op.array("double",val) for val in inputStaticCast(inputsClassic_1b2Wj_El,"float")]
                inputs_array_1b2Wj_El.append(op.array("double",*inputStaticCast(inputsLBN_1b2Wj_El,"float")))
                inputs_array_1b2Wj_El.append(op.array("long",*inputStaticCast(inputsEventNr,"long")))
                inputs_array_1b2Wj_Mu = [op.array("double",val) for val in inputStaticCast(inputsClassic_1b2Wj_Mu,"float")]
                inputs_array_1b2Wj_Mu.append(op.array("double",*inputStaticCast(inputsLBN_1b2Wj_Mu,"float")))
                inputs_array_1b2Wj_Mu.append(op.array("long",*inputStaticCast(inputsEventNr,"long")))
                DNNScore_1b2Wj_El = DNN_1b2Wj(*inputs_array_1b2Wj_El)
                DNNScore_1b2Wj_Mu = DNN_1b2Wj(*inputs_array_1b2Wj_Mu)

                selObjectDNNDictListResolved.append({'channel':'El','selObject':ElSelObjResolved1b2Wj,'DNN_Inputs':inputs_array_1b2Wj_El,'DNN_Score':DNNScore_1b2Wj_El})
                selObjectDNNDictListResolved.append({'channel':'Mu','selObject':MuSelObjResolved1b2Wj,'DNN_Inputs':inputs_array_1b2Wj_Mu,'DNN_Score':DNNScore_1b2Wj_Mu})

                if not self.args.OnlyYield:
                    ChannelDictList.append({'channel':'El','sel':ElSelObjResolved1b2Wj.sel,'lep':ElColl[0],'met':self.corrMET,
                                            'j1':ElSelObjResolved1b2WjJets[0],'j2':ElSelObjResolved1b2WjJets[1],
                                            'j3':ElSelObjResolved1b2WjJets[2],'j4':None,
                                            'nJet':3,'nbJet':1,'suffix':ElSelObjResolved1b2Wj.selName,
                                            'is_MC':self.is_MC})
                    ChannelDictList.append({'channel':'Mu','sel':MuSelObjResolved1b2Wj.sel,'lep':MuColl[0],'met':self.corrMET,
                                            'j1':MuSelObjResolved1b2WjJets[0],'j2':MuSelObjResolved1b2WjJets[1],
                                            'j3':MuSelObjResolved1b2WjJets[2],'j4':None,
                                            'nJet':3,'nbJet':1,'suffix':MuSelObjResolved1b2Wj.selName,
                                            'is_MC':self.is_MC})
                    
            if "Res1b1Wj" in jetplot_level or "Resolved" in jetplot_level:
                logger.info ('...... JPA : 1b1Wj Node Selection')
                ElSelObjResolved1b1Wj        = ElResolvedSelObjJetsIdxPerJpaNodeDict.get('1b1Wj')[0]
                ElSelObjResolved1b1WjJets    = ElResolvedSelObjJetsIdxPerJpaNodeDict.get('1b1Wj')[1]
                MuSelObjResolved1b1Wj        = MuResolvedSelObjJetsIdxPerJpaNodeDict.get('1b1Wj')[0]
                MuSelObjResolved1b1WjJets    = MuResolvedSelObjJetsIdxPerJpaNodeDict.get('1b1Wj')[1]
            
                if self.args.onlypost:
                    ElSelObjResolved1b1Wj.record_yields = True
                    MuSelObjResolved1b1Wj.record_yields = True
                    ElSelObjResolved1b1Wj.yieldTitle = 'Resolved1b1Wj Channel $e^{\pm}$'
                    MuSelObjResolved1b1Wj.yieldTitle = 'Resolved1b1Wj Channel $\mu^{\pm}$'

                ElSelObjResolved1b1WjVBFJets = op.sort(op.combine(op.select(self.VBFJets, cleanVBFwithJPA_Resolved(ElSelObjResolved1b1WjJets, 2)), N=2, pred=self.lambda_VBFPair),
                                                       lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))
                MuSelObjResolved1b1WjVBFJets = op.sort(op.combine(op.select(self.VBFJets, cleanVBFwithJPA_Resolved(MuSelObjResolved1b1WjJets, 2)), N=2, pred=self.lambda_VBFPair),
                                                       lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))

                inputsClassic_1b1Wj_El = returnClassicInputs_1b1Wj(self, ElColl[0], ElSelObjResolved1b1WjJets[0], None, ElSelObjResolved1b1WjJets[1], None, ElSelObjResolved1b1WjVBFJets, elL2OutList)
                inputsLBN_1b1Wj_El     = returnLBNInputs_1b1Wj(self, ElColl[0], ElSelObjResolved1b1WjJets[0], None, ElSelObjResolved1b1WjJets[1], None)
                inputsClassic_1b1Wj_Mu = returnClassicInputs_1b1Wj(self, MuColl[0], MuSelObjResolved1b1WjJets[0], None, MuSelObjResolved1b1WjJets[1], None, MuSelObjResolved1b1WjVBFJets, muL2OutList)
                inputsLBN_1b1Wj_Mu     = returnLBNInputs_1b1Wj(self, MuColl[0], MuSelObjResolved1b1WjJets[0], None, MuSelObjResolved1b1WjJets[1], None)
                input_names_1b1Wj = [key[0] for key in inputsClassic_1b1Wj_El.keys()] + ['LBN_inputs','eventnr']
                DNN_1b1Wj = op.mvaEvaluator(path_model_resolved,mvaType='Tensorflow',otherArgs=(input_names_1b1Wj, output_name))
                inputs_array_1b1Wj_El = [op.array("double",val) for val in inputStaticCast(inputsClassic_1b1Wj_El,"float")]
                inputs_array_1b1Wj_El.append(op.array("double",*inputStaticCast(inputsLBN_1b1Wj_El,"float")))
                inputs_array_1b1Wj_El.append(op.array("long",*inputStaticCast(inputsEventNr,"long")))
                inputs_array_1b1Wj_Mu = [op.array("double",val) for val in inputStaticCast(inputsClassic_1b1Wj_Mu,"float")]
                inputs_array_1b1Wj_Mu.append(op.array("double",*inputStaticCast(inputsLBN_1b1Wj_Mu,"float")))
                inputs_array_1b1Wj_Mu.append(op.array("long",*inputStaticCast(inputsEventNr,"long")))
                DNNScore_1b1Wj_El = DNN_1b1Wj(*inputs_array_1b1Wj_El)
                DNNScore_1b1Wj_Mu = DNN_1b1Wj(*inputs_array_1b1Wj_Mu)

                selObjectDNNDictListResolved.append({'channel':'El','selObject':ElSelObjResolved1b1Wj,'DNN_Inputs':inputs_array_1b1Wj_El,'DNN_Score':DNNScore_1b1Wj_El})
                selObjectDNNDictListResolved.append({'channel':'Mu','selObject':MuSelObjResolved1b1Wj,'DNN_Inputs':inputs_array_1b1Wj_Mu,'DNN_Score':DNNScore_1b1Wj_Mu})

                if not self.args.OnlyYield:
                    ChannelDictList.append({'channel':'El','sel':ElSelObjResolved1b1Wj.sel,'lep':ElColl[0],'met':self.corrMET,
                                            'j1':ElSelObjResolved1b1WjJets[0],'j2':ElSelObjResolved1b1WjJets[1],
                                            'j3':None,'j4':None,
                                            'nJet':2,'nbJet':1,'suffix':ElSelObjResolved1b1Wj.selName,
                                            'is_MC':self.is_MC})
                    ChannelDictList.append({'channel':'Mu','sel':MuSelObjResolved1b1Wj.sel,'lep':MuColl[0],'met':self.corrMET,
                                            'j1':MuSelObjResolved1b1WjJets[0],'j2':MuSelObjResolved1b1WjJets[1],
                                            'j3':None,'j4':None,
                                            'nJet':2,'nbJet':1,'suffix':MuSelObjResolved1b1Wj.selName,
                                            'is_MC':self.is_MC})
                    
            if "Res1b0Wj" in jetplot_level or "Resolved" in jetplot_level:
                logger.info('...... JPA : 1b0Wj Node Selection')
                ElSelObjResolved1b0Wj        = ElResolvedSelObjJetsIdxPerJpaNodeDict.get('1b0Wj')[0]
                ElSelObjResolved1b0WjJets    = ElResolvedSelObjJetsIdxPerJpaNodeDict.get('1b0Wj')[1]
                MuSelObjResolved1b0Wj        = MuResolvedSelObjJetsIdxPerJpaNodeDict.get('1b0Wj')[0]
                MuSelObjResolved1b0WjJets    = MuResolvedSelObjJetsIdxPerJpaNodeDict.get('1b0Wj')[1]

                if self.args.onlypost:
                    ElSelObjResolved1b0Wj.record_yields = True
                    MuSelObjResolved1b0Wj.record_yields = True
                    ElSelObjResolved1b0Wj.yieldTitle = 'Resolved1b0Wj Channel $e^{\pm}$'
                    MuSelObjResolved1b0Wj.yieldTitle = 'Resolved1b0Wj Channel $\mu^{\pm}$'

                ElSelObjResolved1b0WjVBFJets = op.sort(op.combine(op.select(self.VBFJets, cleanVBFwithJPA_Resolved(ElSelObjResolved1b0WjJets, 1)), N=2, pred=self.lambda_VBFPair),
                                                       lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))
                MuSelObjResolved1b0WjVBFJets = op.sort(op.combine(op.select(self.VBFJets, cleanVBFwithJPA_Resolved(MuSelObjResolved1b0WjJets, 1)), N=2, pred=self.lambda_VBFPair),
                                                       lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))

                inputsClassic_1b0Wj_El = returnClassicInputs_1b0Wj(self, ElColl[0], ElSelObjResolved1b0WjJets[0], None, None, None, ElSelObjResolved1b0WjVBFJets, elL2OutList)
                inputsLBN_1b0Wj_El     = returnLBNInputs_1b0Wj(self, ElColl[0], ElSelObjResolved1b0WjJets[0], None, None, None)
                inputsClassic_1b0Wj_Mu = returnClassicInputs_1b0Wj(self, MuColl[0], MuSelObjResolved1b0WjJets[0], None, None, None, MuSelObjResolved1b0WjVBFJets, muL2OutList)
                inputsLBN_1b0Wj_Mu     = returnLBNInputs_1b0Wj(self, MuColl[0], MuSelObjResolved1b0WjJets[0], None, None, None)
                input_names_1b0Wj = [key[0] for key in inputsClassic_1b0Wj_El.keys()] + ['LBN_inputs','eventnr']
                DNN_1b0Wj = op.mvaEvaluator(path_model_resolved,mvaType='Tensorflow',otherArgs=(input_names_1b0Wj, output_name))
                inputs_array_1b0Wj_El = [op.array("double",val) for val in inputStaticCast(inputsClassic_1b0Wj_El,"float")]
                inputs_array_1b0Wj_El.append(op.array("double",*inputStaticCast(inputsLBN_1b0Wj_El,"float")))
                inputs_array_1b0Wj_El.append(op.array("long",*inputStaticCast(inputsEventNr,"long")))
                inputs_array_1b0Wj_Mu = [op.array("double",val) for val in inputStaticCast(inputsClassic_1b0Wj_Mu,"float")]
                inputs_array_1b0Wj_Mu.append(op.array("double",*inputStaticCast(inputsLBN_1b0Wj_Mu,"float")))
                inputs_array_1b0Wj_Mu.append(op.array("long",*inputStaticCast(inputsEventNr,"long")))
                DNNScore_1b0Wj_El = DNN_1b0Wj(*inputs_array_1b0Wj_El)
                DNNScore_1b0Wj_Mu = DNN_1b0Wj(*inputs_array_1b0Wj_Mu)

                selObjectDNNDictListResolved.append({'channel':'El','selObject':ElSelObjResolved1b0Wj,'DNN_Inputs':inputs_array_1b0Wj_El,'DNN_Score':DNNScore_1b0Wj_El})
                selObjectDNNDictListResolved.append({'channel':'Mu','selObject':MuSelObjResolved1b0Wj,'DNN_Inputs':inputs_array_1b0Wj_Mu,'DNN_Score':DNNScore_1b0Wj_Mu})
                
                if not self.args.OnlyYield:
                    ChannelDictList.append({'channel':'El','sel':ElSelObjResolved1b0Wj.sel,'lep':ElColl[0],'met':self.corrMET,
                                            'j1':ElSelObjResolved1b0WjJets[0],'j2':None,
                                            'j3':None,'j4':None,
                                            'nJet':1,'nbJet':1,'suffix':ElSelObjResolved1b0Wj.selName,
                                            'is_MC':self.is_MC})
                    ChannelDictList.append({'channel':'Mu','sel':MuSelObjResolved1b0Wj.sel,'lep':MuColl[0],'met':self.corrMET,
                                            'j1':MuSelObjResolved1b0WjJets[0],'j2':None,
                                            'j3':None,'j4':None,
                                            'nJet':1,'nbJet':1,'suffix':MuSelObjResolved1b0Wj.selName,
                                            'is_MC':self.is_MC})


            if "Res0b" in jetplot_level or "Resolved" in jetplot_level:
                logger.info ('...... JPA : 0b Node Selection')
                ElSelObjResolved0b    = ElResolvedSelObjJetsIdxPerJpaNodeDict.get('0b')[0]
                MuSelObjResolved0b    = MuResolvedSelObjJetsIdxPerJpaNodeDict.get('0b')[0]
                print(type(ElSelObjResolved0b))
                if self.args.onlypost:
                    ElSelObjResolved0b.record_yields = True
                    MuSelObjResolved0b.record_yields = True
                    ElSelObjResolved0b.yieldTitle = 'Resolved0b Channel $e^{\pm}$'
                    MuSelObjResolved0b.yieldTitle = 'Resolved0b Channel $\mu^{\pm}$'

                inputsClassic_0b_El = returnClassicInputs_0b(self, ElColl[0], None, None, None, None, None, elL2OutList)
                inputsLBN_0b_El     = returnLBNInputs_0b(self, ElColl[0], None, None, None, None)
                inputsClassic_0b_Mu = returnClassicInputs_0b(self, MuColl[0], None, None, None, None, None, muL2OutList)
                inputsLBN_0b_Mu     = returnLBNInputs_0b(self, MuColl[0], None, None, None, None)
                input_names_0b = [key[0] for key in inputsClassic_0b_El.keys()] + ['LBN_inputs','eventnr']
                DNN_0b = op.mvaEvaluator(path_model_resolved,mvaType='Tensorflow',otherArgs=(input_names_0b, output_name))
                inputs_array_0b_El = [op.array("double",val) for val in inputStaticCast(inputsClassic_0b_El,"float")]
                inputs_array_0b_El.append(op.array("double",*inputStaticCast(inputsLBN_0b_El,"float")))
                inputs_array_0b_El.append(op.array("long",*inputStaticCast(inputsEventNr,"long")))
                inputs_array_0b_Mu = [op.array("double",val) for val in inputStaticCast(inputsClassic_0b_Mu,"float")]
                inputs_array_0b_Mu.append(op.array("double",*inputStaticCast(inputsLBN_0b_Mu,"float")))
                inputs_array_0b_Mu.append(op.array("long",*inputStaticCast(inputsEventNr,"long")))
                DNNScore_0b_El = DNN_0b(*inputs_array_0b_El)
                DNNScore_0b_Mu = DNN_0b(*inputs_array_0b_Mu)

                selObjectDNNDictListResolved.append({'channel':'El','selObject':ElSelObjResolved0b,'DNN_Inputs':inputs_array_0b_El,'DNN_Score':DNNScore_0b_El})
                selObjectDNNDictListResolved.append({'channel':'Mu','selObject':MuSelObjResolved0b,'DNN_Inputs':inputs_array_0b_Mu,'DNN_Score':DNNScore_0b_Mu})


            for channelDict in ChannelDictList:
                # Singlelepton #
                plots.extend(makeSinleptonPlots(**{k:channelDict[k] for k in LeptonKeys}))
                # Number of jets #
                #plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**JetsN))
                #plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**FatJetsN))
                # Ak4 Jets #
                plots.extend(makeAk4JetsPlots(**{k:channelDict[k] for k in JetKeys},HLL=self.HLL))
                # MET #
                plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))
                # High level #
                #plots.extend(makeHighLevelPlotsResolved(**{k:channelDict[k] for k in ResolvedKeys},HLL=self.HLL))

            # ML
            for selObjectDNNDict in selObjectDNNDictListResolved:
                lepton = leptonCont[selObjectDNNDict['channel']]
                inputs = selObjectDNNDict['DNN_Inputs']
                output = selObjectDNNDict['DNN_Score']
                selObj = selObjectDNNDict['selObject']
                selObj_0b = makeExclusiveLooseResolvedJetComboSelection(self, selObj, 0, copy_sel=True)
                selObj_1b = makeExclusiveLooseResolvedJetComboSelection(self, selObj, 1, copy_sel=True)
                selObj_2b = makeExclusiveLooseResolvedJetComboSelection(self, selObj, 2, copy_sel=True)
                selObjNodesDict    = makeDNNOutputNodesSelections(self,selObj,output,suffix='')
                selObjNodesDict_0b = makeDNNOutputNodesSelections(self,selObj_0b,output,suffix='')
                selObjNodesDict_1b = makeDNNOutputNodesSelections(self,selObj_1b,output,suffix='')
                selObjNodesDict_2b = makeDNNOutputNodesSelections(self,selObj_2b,output,suffix='')
                #plots.extend(makeDoubleLeptonMachineLearningInputPlots(selObjectDNNDict['selObject'].sel,selObjectDNNDict['selObject'].selName,selObjectDNNDict['channel'],inputs))
                logger.info('Filling DNN responses')
                plots.extend(makeDoubleLeptonMachineLearningOutputPlots(selObjNodesDict,output,self.nodes,channel=selObjectDNNDict['channel']))
                plots.extend(makeDoubleLeptonMachineLearningOutputPlots(selObjNodesDict_0b,output,self.nodes,channel=selObjectDNNDict['channel']))
                plots.extend(makeDoubleLeptonMachineLearningOutputPlots(selObjNodesDict_1b,output,self.nodes,channel=selObjectDNNDict['channel']))
                plots.extend(makeDoubleLeptonMachineLearningOutputPlots(selObjNodesDict_2b,output,self.nodes,channel=selObjectDNNDict['channel']))

        # ========================== JPA Boosted Categories ========================= #
        if any(item in boosted_args for item in jetsel_level):
            ChannelDictList = []
            selObjectDNNDictListBoosted = []
            # dict = {'key':'Node', 'value' : [refined selObj, [JPAjetIndices]]}    
            elL1OutList, elL2OutList, ElBoostedSelObjJetsIdxPerJpaNodeDict = findJPACategoryBoosted (self, ElSelObjBoosted, ElColl[0], self.muonsPreSel, self.electronsPreSel, 
                                                                                                     self.ak8BJets, self.ak4JetsCleanedFromAk8b, self.ak4BJetsLoose, 
                                                                                                     self.ak4BJets, self.corrMET, boostedModelDict, t.event, self.HLL, 
                                                                                                     BoostedJPANodeList,
                                                                                                     plot_yield=True)
            muL1OutList, muL2OutList, MuBoostedSelObjJetsIdxPerJpaNodeDict = findJPACategoryBoosted (self, MuSelObjBoosted, MuColl[0], self.muonsPreSel, self.electronsPreSel, 
                                                                                                     self.ak8BJets, self.ak4JetsCleanedFromAk8b, self.ak4BJetsLoose, 
                                                                                                     self.ak4BJets, self.corrMET, boostedModelDict, t.event, self.HLL, 
                                                                                                     BoostedJPANodeList, 
                                                                                                     plot_yield=True)

            if "Hbb2Wj" in jetplot_level or "Boosted" in jetplot_level:
                print ('...... JPA : Hbb2Wj Node Selection')
                ElSelObjBoostedHbb2Wj        = ElBoostedSelObjJetsIdxPerJpaNodeDict.get('Hbb2Wj')[0]
                ElSelObjBoostedHbb2WjJets    = ElBoostedSelObjJetsIdxPerJpaNodeDict.get('Hbb2Wj')[1]
                MuSelObjBoostedHbb2Wj        = MuBoostedSelObjJetsIdxPerJpaNodeDict.get('Hbb2Wj')[0]
                MuSelObjBoostedHbb2WjJets    = MuBoostedSelObjJetsIdxPerJpaNodeDict.get('Hbb2Wj')[1]

                # VBF Jets
                ElSelObjBoostedHbb2WjVBFJets = op.sort(op.combine(op.select(self.VBFJets, cleanVBFwithJPA_Boosted(ElSelObjBoostedHbb2WjJets, 2)), N=2, pred=self.lambda_VBFPair), 
                                                       lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))
                MuSelObjBoostedHbb2WjVBFJets = op.sort(op.combine(op.select(self.VBFJets, cleanVBFwithJPA_Boosted(MuSelObjBoostedHbb2WjJets, 2)), N=2, pred=self.lambda_VBFPair), 
                                                       lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))

                inputsClassic_Hbb2Wj_El = returnClassicInputs_Hbb2Wj(self, ElColl[0], self.ak8BJets[0], ElSelObjBoostedHbb2WjJets[0], ElSelObjBoostedHbb2WjJets[1], 
                                                                     ElSelObjBoostedHbb2WjVBFJets, elL2OutList, self.era)
                inputsLBN_Hbb2Wj_El     = returnLBNInputs_Hbb2Wj(self, ElColl[0], self.ak8BJets[0], ElSelObjBoostedHbb2WjJets[0], ElSelObjBoostedHbb2WjJets[1])
                inputsClassic_Hbb2Wj_Mu = returnClassicInputs_Hbb2Wj(self, ElColl[0], self.ak8BJets[0], MuSelObjBoostedHbb2WjJets[0], MuSelObjBoostedHbb2WjJets[1], 
                                                                     MuSelObjBoostedHbb2WjVBFJets, muL2OutList, self.era)
                inputsLBN_Hbb2Wj_Mu     = returnLBNInputs_Hbb2Wj(self, MuColl[0], self.ak8BJets[0], MuSelObjBoostedHbb2WjJets[0], MuSelObjBoostedHbb2WjJets[1])

                input_names_Hbb2Wj      = [key[0] for key in inputsClassic_Hbb2Wj_El.keys()] + ['LBN_inputs','eventnr']
                DNN_Hbb2Wj   = op.mvaEvaluator(path_model_resolved,mvaType='Tensorflow',otherArgs=(input_names_Hbb2Wj, output_name))
                inputs_array_Hbb2Wj_El = [op.array("double",val) for val in inputStaticCast(inputsClassic_Hbb2Wj_El,"float")]
                inputs_array_Hbb2Wj_El.append(op.array("double",*inputStaticCast(inputsLBN_Hbb2Wj_El,"float")))
                inputs_array_Hbb2Wj_El.append(op.array("long",*inputStaticCast(inputsEventNr,"long")))
                inputs_array_Hbb2Wj_Mu = [op.array("double",val) for val in inputStaticCast(inputsClassic_Hbb2Wj_Mu,"float")]
                inputs_array_Hbb2Wj_Mu.append(op.array("double",*inputStaticCast(inputsLBN_Hbb2Wj_Mu,"float")))
                inputs_array_Hbb2Wj_Mu.append(op.array("long",*inputStaticCast(inputsEventNr,"long")))
                DNNScore_Hbb2Wj_El = DNN_Hbb2Wj(*inputs_array_Hbb2Wj_El)
                DNNScore_Hbb2Wj_Mu = DNN_Hbb2Wj(*inputs_array_Hbb2Wj_Mu)

                selObjectDNNDictListBoosted.append({'channel':'El','selObject':ElSelObjBoostedHbb2Wj,'DNN_Inputs':inputs_array_Hbb2Wj_El,'DNN_Score':DNNScore_Hbb2Wj_El})
                selObjectDNNDictListBoosted.append({'channel':'Mu','selObject':MuSelObjBoostedHbb2Wj,'DNN_Inputs':inputs_array_Hbb2Wj_Mu,'DNN_Score':DNNScore_Hbb2Wj_Mu})

                if not self.args.OnlyYield:
                    ChannelDictList.append({'channel':'El','sel':ElSelObjBoostedHbb2Wj.sel,'lep':ElColl[0],'met':self.corrMET,
                                            'j1':self.ak8BJets[0],'j2':ElSelObjBoostedHbb2WjJets[0],'j3':ElSelObjBoostedHbb2WjJets[1],
                                            'has1fat1slim':False,'has1fat2slim':True,'bothAreFat':False,
                                            'suffix':ElSelObjBoostedHbb2Wj.selName,
                                            'is_MC':self.is_MC})
                    ChannelDictList.append({'channel':'Mu','sel':MuSelObjBoostedHbb2Wj.sel,'lep':MuColl[0],'met':self.corrMET,
                                            'j1':self.ak8BJets[0],'j2':MuSelObjBoostedHbb2WjJets[0],'j3':MuSelObjBoostedHbb2WjJets[1],
                                            'has1fat1slim':False,'has1fat2slim':True,'bothAreFat':False,
                                            'suffix':MuSelObjBoostedHbb2Wj.selName,
                                            'is_MC':self.is_MC})
                    
            if "Hbb1Wj" in jetplot_level or "Boosted" in jetplot_level:
                print ('...... JPA : Hbb1Wj Node Selection')
                ElSelObjBoostedHbb1Wj        = ElBoostedSelObjJetsIdxPerJpaNodeDict.get('Hbb1Wj')[0]
                ElSelObjBoostedHbb1WjJets    = ElBoostedSelObjJetsIdxPerJpaNodeDict.get('Hbb1Wj')[1]
                MuSelObjBoostedHbb1Wj        = MuBoostedSelObjJetsIdxPerJpaNodeDict.get('Hbb1Wj')[0]
                MuSelObjBoostedHbb1WjJets    = MuBoostedSelObjJetsIdxPerJpaNodeDict.get('Hbb1Wj')[1]
                
                # VBF Jets
                ElSelObjBoostedHbb1WjVBFJets = op.sort(op.combine(op.select(self.VBFJets, cleanVBFwithJPA_Boosted(ElSelObjBoostedHbb1WjJets, 1)), N=2, pred=self.lambda_VBFPair), 
                                                       lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))
                MuSelObjBoostedHbb1WjVBFJets = op.sort(op.combine(op.select(self.VBFJets, cleanVBFwithJPA_Boosted(MuSelObjBoostedHbb1WjJets, 1)), N=2, pred=self.lambda_VBFPair), 
                                                       lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))

                inputsClassic_Hbb1Wj_El = returnClassicInputs_Hbb1Wj(self, ElColl[0], self.ak8BJets[0], ElSelObjBoostedHbb1WjJets[1], None, ElSelObjBoostedHbb1WjVBFJets, elL2OutList, self.era)
                inputsLBN_Hbb1Wj_El     = returnLBNInputs_Hbb1Wj(self, ElColl[0], self.ak8BJets[0], ElSelObjBoostedHbb2WjJets[1], None)
                inputsClassic_Hbb1Wj_Mu = returnClassicInputs_Hbb1Wj(self, ElColl[0], self.ak8BJets[0], MuSelObjBoostedHbb1WjJets[1], None, MuSelObjBoostedHbb1WjVBFJets, muL2OutList, self.era)
                inputsLBN_Hbb1Wj_Mu     = returnLBNInputs_Hbb1Wj(self, MuColl[0], self.ak8BJets[0], MuSelObjBoostedHbb1WjJets[1], None)

                input_names_Hbb1Wj      = [key[0] for key in inputsClassic_Hbb1Wj_El.keys()] + ['LBN_inputs','eventnr']
                DNN_Hbb1Wj   = op.mvaEvaluator(path_model_resolved,mvaType='Tensorflow',otherArgs=(input_names_Hbb1Wj, output_name))
                inputs_array_Hbb1Wj_El = [op.array("double",val) for val in inputStaticCast(inputsClassic_Hbb1Wj_El,"float")]
                inputs_array_Hbb1Wj_El.append(op.array("double",*inputStaticCast(inputsLBN_Hbb1Wj_El,"float")))
                inputs_array_Hbb1Wj_El.append(op.array("long",*inputStaticCast(inputsEventNr,"long")))
                inputs_array_Hbb1Wj_Mu = [op.array("double",val) for val in inputStaticCast(inputsClassic_Hbb1Wj_Mu,"float")]
                inputs_array_Hbb1Wj_Mu.append(op.array("double",*inputStaticCast(inputsLBN_Hbb1Wj_Mu,"float")))
                inputs_array_Hbb1Wj_Mu.append(op.array("long",*inputStaticCast(inputsEventNr,"long")))
                DNNScore_Hbb1Wj_El = DNN_Hbb1Wj(*inputs_array_Hbb1Wj_El)
                DNNScore_Hbb1Wj_Mu = DNN_Hbb1Wj(*inputs_array_Hbb1Wj_Mu)

                selObjectDNNDictListBoosted.append({'channel':'El','selObject':ElSelObjBoostedHbb1Wj,'DNN_Inputs':inputs_array_Hbb1Wj_El,'DNN_Score':DNNScore_Hbb1Wj_El})
                selObjectDNNDictListBoosted.append({'channel':'Mu','selObject':MuSelObjBoostedHbb1Wj,'DNN_Inputs':inputs_array_Hbb1Wj_Mu,'DNN_Score':DNNScore_Hbb1Wj_Mu})


                if not self.args.OnlyYield:
                    ChannelDictList.append({'channel':'El','sel':ElSelObjBoostedHbb1Wj.sel,'lep':ElColl[0],'met':self.corrMET,
                                            'j1':self.ak8BJets[0],'j2':ElSelObjBoostedHbb1WjJets[0],'j3':None,
                                            'has1fat1slim':True,'has1fat2slim':False,'bothAreFat':False,
                                            'suffix':ElSelObjBoostedHbb1Wj.selName,
                                            'is_MC':self.is_MC})
                    ChannelDictList.append({'channel':'Mu','sel':MuSelObjBoostedHbb1Wj.sel,'lep':MuColl[0],'met':self.corrMET,
                                            'j1':self.ak8BJets[0],'j2':MuSelObjBoostedHbb1WjJets[0],'j3':None,
                                            'has1fat1slim':True,'has1fat2slim':False,'bothAreFat':False,
                                            'suffix':MuSelObjBoostedHbb1Wj.selName,
                                            'is_MC':self.is_MC})
                        
            if "Hbb0Wj" in jetplot_level or "Boosted" in jetplot_level:
                print ('...... JPA : Hbb0Wj Node Selection')
                ElSelObjBoostedHbb0Wj        = ElBoostedSelObjJetsIdxPerJpaNodeDict.get('Hbb1Wj')[0]
                MuSelObjBoostedHbb0Wj        = MuBoostedSelObjJetsIdxPerJpaNodeDict.get('Hbb1Wj')[0]
                
                # VBF Jets
                ElSelObjBoostedHbb0WjVBFJets = op.sort(op.combine(op.select(self.VBFJets, cleanVBFwithJPA_Boosted_rest()), N=2, pred=self.lambda_VBFPair),
                                                       lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))
                MuSelObjBoostedHbb0WjVBFJets = op.sort(op.combine(op.select(self.VBFJets, cleanVBFwithJPA_Boosted_rest()), N=2, pred=self.lambda_VBFPair),
                                                       lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))

                inputsClassic_Hbb0Wj_El = returnClassicInputs_Hbb0Wj(self, ElColl[0], self.ak8BJets[0], None, None, ElSelObjBoostedHbb0WjVBFJets, elL2OutList, self.era)
                inputsLBN_Hbb0Wj_El     = returnLBNInputs_Hbb0Wj(self, ElColl[0], self.ak8BJets[0], None, None)
                inputsClassic_Hbb0Wj_Mu = returnClassicInputs_Hbb0Wj(self, ElColl[0], self.ak8BJets[0], None, None, MuSelObjBoostedHbb0WjVBFJets, muL2OutList, self.era)
                inputsLBN_Hbb0Wj_Mu     = returnLBNInputs_Hbb0Wj(self, MuColl[0], self.ak8BJets[0], None, None)

                input_names_Hbb0Wj      = [key[0] for key in inputsClassic_Hbb0Wj_El.keys()] + ['LBN_inputs','eventnr']
                DNN_Hbb0Wj   = op.mvaEvaluator(path_model_resolved,mvaType='Tensorflow',otherArgs=(input_names_Hbb0Wj, output_name))
                inputs_array_Hbb0Wj_El = [op.array("double",val) for val in inputStaticCast(inputsClassic_Hbb0Wj_El,"float")]
                inputs_array_Hbb0Wj_El.append(op.array("double",*inputStaticCast(inputsLBN_Hbb0Wj_El,"float")))
                inputs_array_Hbb0Wj_El.append(op.array("long",*inputStaticCast(inputsEventNr,"long")))
                inputs_array_Hbb0Wj_Mu = [op.array("double",val) for val in inputStaticCast(inputsClassic_Hbb0Wj_Mu,"float")]
                inputs_array_Hbb0Wj_Mu.append(op.array("double",*inputStaticCast(inputsLBN_Hbb0Wj_Mu,"float")))
                inputs_array_Hbb0Wj_Mu.append(op.array("long",*inputStaticCast(inputsEventNr,"long")))
                DNNScore_Hbb0Wj_El = DNN_Hbb0Wj(*inputs_array_Hbb0Wj_El)
                DNNScore_Hbb0Wj_Mu = DNN_Hbb0Wj(*inputs_array_Hbb0Wj_Mu)

                selObjectDNNDictListBoosted.append({'channel':'El','selObject':ElSelObjBoostedHbb0Wj,'DNN_Inputs':inputs_array_Hbb0Wj_El,'DNN_Score':DNNScore_Hbb0Wj_El})
                selObjectDNNDictListBoosted.append({'channel':'Mu','selObject':MuSelObjBoostedHbb0Wj,'DNN_Inputs':inputs_array_Hbb0Wj_Mu,'DNN_Score':DNNScore_Hbb0Wj_Mu})

            for channelDict in ChannelDictList:
                # Dilepton #
                plots.extend(makeSinleptonPlots(**{k:channelDict[k] for k in LeptonKeys}))
                # Number of jets #
                ##plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**FatJetsN))
                ##plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**SlimJetsN))
                # Ak8 Jets #
                ##plots.extend(makeSingleLeptonAk8JetsPlots(**{k:channelDict[k] for k in FatJetKeys},nMedBJets=self.nMediumBTaggedSubJets, HLL=self.HLL))
                # MET #
                plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))
                # HighLevel #
                ##plots.extend(makeHighLevelPlotsBoosted(**{k:channelDict[k] for k in BoostedKeys}, HLL=self.HLL))
                
            # ML
            for selObjectDNNDict in selObjectDNNDictListBoosted:
                lepton = leptonCont[selObjectDNNDict['channel']]
                inputs = selObjectDNNDict['DNN_Inputs']
                output = selObjectDNNDict['DNN_Score']
                selObj = selObjectDNNDict['selObject']
                selObjNodesDict    = makeDNNOutputNodesSelections(self,selObj,output,suffix='')
                #plots.extend(makeDoubleLeptonMachineLearningInputPlots(selObjectDNNDict['selObject'].sel,selObjectDNNDict['selObject'].selName,selObjectDNNDict['channel'],inputs))
                logger.info('Filling DNN responses Boosted')
                plots.extend(makeDoubleLeptonMachineLearningOutputPlots(selObjNodesDict,output,self.nodes,channel=selObjectDNNDict['channel']))
                      
                
        #----- Add the Yield plots -----#
        plots.append(self.yields)
        #plots.extend(cutFlowPlots)
        return plots

