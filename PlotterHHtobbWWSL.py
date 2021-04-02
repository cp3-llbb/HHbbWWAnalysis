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
        # ----------------------------------- NodeList ------------------------------------- #
        # keep the exact same order of nodes as mentioned in respective xml files
        ResolvedJPANodeList = ['2b2Wj','2b1Wj','1b2Wj','2b0Wj','1b1Wj','1b0Wj','0b']
        BoostedJPANodeList  = ['Hbb2Wj','Hbb1Wj','Hbb0Wj']

        basepath = os.path.join(os.path.abspath(os.path.dirname(__file__)),'MachineLearning','ml-models','JPA_Loose_ttH')
        resolvedModelDict = getResolvedJpaModelDict(basepath, ResolvedJPANodeList, era)        
        boostedModelDict  = getBoostedJpaModelDict(basepath, BoostedJPANodeList, era)
        
        # ---------- LBN+DNN models ----------- #
        #path_model_resolved = os.path.join('/home/ucl/cp3/gsaha/bamboodev/HHbbWWAnalysis/MachineLearning/ml-models/DNN/Resolved','Resolved'+str(era)+'.pb')
        #path_model_boosted  = os.path.join('/home/ucl/cp3/gsaha/bamboodev/HHbbWWAnalysis/MachineLearning/ml-models/DNN/Boosted','Boosted_allEras_x512.pb')
        path_model_resolved_SM  = os.path.join('/home/ucl/cp3/gsaha/bamboodev/HHbbWWAnalysis/MachineLearning/ml-models/DNN/Resolved','Resolved'+str(era)+'SMv2.pb')
        #        path_model_resolved_BSM = os.path.join('/home/ucl/cp3/gsaha/bamboodev/HHbbWWAnalysis/MachineLearning/ml-models/DNN/Resolved','Resolved'+str(era)+'BSMv2.pb')
        #        path_model_resolved_BSM = os.path.join('/home/ucl/cp3/gsaha/bamboodev/HHbbWWAnalysis/MachineLearning/ml-models/DNN/Resolved','Resolved'+str(era)+'BSMv3_w0p1.pb')
        path_model_resolved_BSM = os.path.join('/home/ucl/cp3/gsaha/bamboodev/HHbbWWAnalysis/MachineLearning/ml-models/DNN/Resolved','Resolved'+str(era)+'BSMv3_w1.pb')
        path_model_boosted_SM   = os.path.join('/home/ucl/cp3/gsaha/bamboodev/HHbbWWAnalysis/MachineLearning/ml-models/DNN/Boosted','BoostedSMv2.pb')
        #        path_model_boosted_BSM  = os.path.join('/home/ucl/cp3/gsaha/bamboodev/HHbbWWAnalysis/MachineLearning/ml-models/DNN/Boosted','BoostedBSMv2.pb')
        #        path_model_boosted_BSM  = os.path.join('/home/ucl/cp3/gsaha/bamboodev/HHbbWWAnalysis/MachineLearning/ml-models/DNN/Boosted','BoostedBSMv3_w0p1_256x4.pb')
        path_model_boosted_BSM  = os.path.join('/home/ucl/cp3/gsaha/bamboodev/HHbbWWAnalysis/MachineLearning/ml-models/DNN/Boosted','BoostedBSMv3_w1_256x4.pb')
        logger.info('DNN_Model_Resolved SM : {}'.format(path_model_resolved_SM))
        logger.info('DNN_Model_Resolved BSM : {}'.format(path_model_resolved_BSM))
        logger.info('DNN_Model_Boosted SM : {}'.format(path_model_boosted_SM))
        logger.info('DNN_Model_Boosted BSM : {}'.format(path_model_boosted_BSM))
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
        resolved_args = ["Res2b2Wj","Res2b1Wj","Res1b2Wj","Res2b0Wj","Res1b1Wj","Res1b0Wj","Res0b","Resolved"]
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
            
        # Selections:    
        #---- Lepton selection ----#
        ElColl = [t.Electron[op.switch(op.rng_len(self.electronsTightSel) == 1, 
                                       self.electronsTightSel[0].idx, 
                                       self.electronsFakeSel[0].idx)]]
        MuColl = [t.Muon[op.switch(op.rng_len(self.muonsTightSel) == 1, 
                                   self.muonsTightSel[0].idx, 
                                   self.muonsFakeSel[0].idx)]]

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
            
        #----- Ak4 jets selection -----#
        if "Ak4" in jetsel_level:
            logger.info ("... Processing Ak4Jets Selection for Resolved category : nAk4Jets >= 3")

            ElSelObjResolved = makeResolvedSelection(self,ElSelObj,copy_sel=True,plot_yield=True)
            MuSelObjResolved = makeResolvedSelection(self,MuSelObj,copy_sel=True,plot_yield=True)
        
        #----- Ak8-b jets selection -----#
        if "Ak8" in jetsel_level:
            logger.info ("...... Processing Ak8b jet selection for SemiBoosted & Boosted Category")

            ElSelObjBoosted = makeBoostedSelection(self,ElSelObj,copy_sel=True,plot_yield=True)
            MuSelObjBoosted = makeBoostedSelection(self,MuSelObj,copy_sel=True,plot_yield=True)
            

        self.nodes       = ['Ewk','GGF','H','Top','VBF','WJets']
        inputsEventNr    = returnEventNr(self, t)
        output_name      = "Identity"
        ClassicInputKeys = ['lepton','jpaSelectedJets','L1out','L2out','jpaArg']
        LBNInputKeys     = ['lepton','jet1','jet2','jet3','jet4']

        # ========================== JPA Resolved Categories ========================= #
        if any(item in resolved_args for item in jetsel_level):
            ChannelDictListR = []

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

                if self.args.onlypost:
                    ElSelObjResolved2b2Wj.record_yields = True
                    MuSelObjResolved2b2Wj.record_yields = True
                    ElSelObjResolved2b2Wj.yieldTitle = 'Resolved2b2Wj Channel $e^{\pm}$'
                    MuSelObjResolved2b2Wj.yieldTitle = 'Resolved2b2Wj Channel $\mu^{\pm}$'
                
                if not self.args.OnlyYield:
                    ChannelDictListR.append({'channel':'El','selObj':ElSelObjResolved2b2Wj,'sel':ElSelObjResolved2b2Wj.sel,
                                             'lepton':ElColl[0],'met':self.corrMET,'jpaSelectedJets':ElSelObjResolved2b2WjJets,
                                             'jet1':ElSelObjResolved2b2WjJets[0],'jet2':ElSelObjResolved2b2WjJets[1],
                                             'jet3':ElSelObjResolved2b2WjJets[2],'jet4':ElSelObjResolved2b2WjJets[3],
                                             'nJet':4,'nbJet':2,'suffix':ElSelObjResolved2b2Wj.selName,
                                             'is_MC':self.is_MC,'jpaArg':'Res2b2Wj'})
                    ChannelDictListR.append({'channel':'Mu','selObj':MuSelObjResolved2b2Wj,'sel':MuSelObjResolved2b2Wj.sel,
                                             'lepton':MuColl[0],'met':self.corrMET,'jpaSelectedJets':MuSelObjResolved2b2WjJets,
                                             'jet1':MuSelObjResolved2b2WjJets[0],'jet2':MuSelObjResolved2b2WjJets[1],
                                             'jet3':MuSelObjResolved2b2WjJets[2],'jet4':MuSelObjResolved2b2WjJets[3],
                                             'nJet':4,'nbJet':2,'suffix':MuSelObjResolved2b2Wj.selName,
                                             'is_MC':self.is_MC,'jpaArg':'Res2b2Wj'})
                    
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
                
                if not self.args.OnlyYield:
                    ChannelDictListR.append({'channel':'El','selObj':ElSelObjResolved2b1Wj,'sel':ElSelObjResolved2b1Wj.sel,
                                            'lep':ElColl[0],'met':self.corrMET,'jpaSelectedJets':ElSelObjResolved2b1WjJets,
                                            'jet1':ElSelObjResolved2b1WjJets[0],'jet2':ElSelObjResolved2b1WjJets[1],
                                            'jet3':ElSelObjResolved2b1WjJets[2],'jet4':None,
                                            'nJet':3,'nbJet':2,'suffix':ElSelObjResolved2b1Wj.selName,
                                            'is_MC':self.is_MC, 'jpaArg':'Res2b1Wj'})
                    ChannelDictListR.append({'channel':'Mu','selObj':MuSelObjResolved2b1Wj,'sel':MuSelObjResolved2b1Wj.sel,
                                            'lep':MuColl[0],'met':self.corrMET,'jpaSelectedJets':MuSelObjResolved2b1WjJets,
                                            'jet1':MuSelObjResolved2b1WjJets[0],'jet2':MuSelObjResolved2b1WjJets[1],
                                            'jet3':MuSelObjResolved2b1WjJets[2],'jet4':None,
                                            'nJet':3,'nbJet':2,'suffix':MuSelObjResolved2b1Wj.selName,
                                            'is_MC':self.is_MC, 'jpaArg':'Res2b1Wj'})
                    
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

                if not self.args.OnlyYield:
                    ChannelDictListR.append({'channel':'El','selObj':ElSelObjResolved1b2Wj,'sel':ElSelObjResolved1b2Wj.sel,
                                             'lep':ElColl[0],'met':self.corrMET,'jpaSelectedJets':ElSelObjResolved1b2WjJets,
                                             'jet1':ElSelObjResolved1b2WjJets[0],'jet2':ElSelObjResolved1b2WjJets[1],
                                             'jet3':ElSelObjResolved1b2WjJets[2],'jet4':None,
                                             'nJet':3,'nbJet':1,'suffix':ElSelObjResolved1b2Wj.selName,
                                             'is_MC':self.is_MC, 'jpaArg':'Res1b2Wj'})
                    ChannelDictListR.append({'channel':'Mu','selObj':MuSelObjResolved1b2Wj,'sel':MuSelObjResolved1b2Wj.sel,
                                             'lep':MuColl[0],'met':self.corrMET,'jpaSelectedJets':MuSelObjResolved1b2WjJets,
                                             'jet1':MuSelObjResolved1b2WjJets[0],'jet2':MuSelObjResolved1b2WjJets[1],
                                             'jet3':MuSelObjResolved1b2WjJets[2],'jet4':None,
                                             'nJet':3,'nbJet':1,'suffix':MuSelObjResolved1b2Wj.selName,
                                             'is_MC':self.is_MC, 'jpaArg':'Res1b2Wj'})
                    
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

                if not self.args.OnlyYield:
                    ChannelDictListR.append({'channel':'El','selObj':ElSelObjResolved2b0Wj,'sel':ElSelObjResolved2b0Wj.sel,
                                             'lep':ElColl[0],'met':self.corrMET,'jpaSelectedJets':ElSelObjResolved2b0WjJets,
                                             'jet1':ElSelObjResolved2b0WjJets[0],'jet2':ElSelObjResolved2b0WjJets[1],
                                             'jet3':None,'jet4':None,
                                             'nJet':2,'nbJet':2,'suffix':ElSelObjResolved2b0Wj.selName,
                                             'is_MC':self.is_MC, 'jpaArg':'Res2b0Wj'})
                    ChannelDictListR.append({'channel':'Mu','selObj':MuSelObjResolved2b0Wj,'sel':MuSelObjResolved2b0Wj.sel,
                                             'lep':MuColl[0],'met':self.corrMET,'jpaSelectedJets':MuSelObjResolved2b0WjJets,
                                             'jet1':MuSelObjResolved2b0WjJets[0],'jet2':MuSelObjResolved2b0WjJets[1],
                                             'jet3':None,'jet4':None,
                                             'nJet':2,'nbJet':2,'suffix':MuSelObjResolved2b0Wj.selName,
                                             'is_MC':self.is_MC, 'jpaArg':'Res2b0Wj'})
                    
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

                if not self.args.OnlyYield:
                    ChannelDictListR.append({'channel':'El','selObj':ElSelObjResolved1b1Wj,'sel':ElSelObjResolved1b1Wj.sel,
                                             'lep':ElColl[0],'met':self.corrMET,'jpaSelectedJets':ElSelObjResolved1b1WjJets,
                                             'jet1':ElSelObjResolved1b1WjJets[0],'jet2':ElSelObjResolved1b1WjJets[1],
                                             'jet3':None,'jet4':None,
                                             'nJet':2,'nbJet':1,'suffix':ElSelObjResolved1b1Wj.selName,
                                             'is_MC':self.is_MC, 'jpaArg':'Res1b1Wj'})
                    ChannelDictListR.append({'channel':'Mu','selObj':MuSelObjResolved1b1Wj,'sel':MuSelObjResolved1b1Wj.sel,
                                             'lep':MuColl[0],'met':self.corrMET,'jpaSelectedJets':MuSelObjResolved1b1WjJets,
                                             'jet1':MuSelObjResolved1b1WjJets[0],'jet2':MuSelObjResolved1b1WjJets[1],
                                             'jet3':None,'jet4':None,
                                             'nJet':2,'nbJet':1,'suffix':MuSelObjResolved1b1Wj.selName,
                                             'is_MC':self.is_MC, 'jpaArg':'Res1b1Wj'})
                    
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

                if not self.args.OnlyYield:
                    ChannelDictListR.append({'channel':'El','selObj':ElSelObjResolved1b0Wj,'sel':ElSelObjResolved1b0Wj.sel,
                                             'lep':ElColl[0],'met':self.corrMET,'jpaSelectedJets':ElSelObjResolved1b0WjJets,
                                             'jet1':ElSelObjResolved1b0WjJets[0],'jet2':None,
                                             'jet3':None,'jet4':None,
                                             'nJet':1,'nbJet':1,'suffix':ElSelObjResolved1b0Wj.selName,
                                             'is_MC':self.is_MC, 'jpaArg':'Res1b0Wj'})
                    ChannelDictListR.append({'channel':'Mu','selObj':MuSelObjResolved1b0Wj,'sel':MuSelObjResolved1b0Wj.sel,
                                             'lep':MuColl[0],'met':self.corrMET,'jpaSelectedJets':MuSelObjResolved1b0WjJets,
                                             'jet1':MuSelObjResolved1b0WjJets[0],'jet2':None,
                                             'jet3':None,'jet4':None,
                                             'nJet':1,'nbJet':1,'suffix':MuSelObjResolved1b0Wj.selName,
                                             'is_MC':self.is_MC, 'jpaArg':'Res1b0Wj'})
                    

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

                if not self.args.OnlyYield:
                    ChannelDictListR.append({'channel':'El','selObj':ElSelObjResolved0b,'sel':ElSelObjResolved0b.sel,
                                             'lep':ElColl[0],'met':self.corrMET,'jpaSelectedJets':None,
                                             'jet1':None,'jet2':None,
                                             'jet3':None,'jet4':None,
                                             'nJet':1,'nbJet':1,'suffix':ElSelObjResolved0b.selName,
                                             'is_MC':self.is_MC, 'jpaArg':'Res0b'})
                    ChannelDictListR.append({'channel':'Mu','selObj':MuSelObjResolved0b,'sel':MuSelObjResolved0b.sel,
                                             'lep':MuColl[0],'met':self.corrMET,'jpaSelectedJets':None,
                                             'jet1':None,'jet2':None,
                                             'jet3':None,'jet4':None,
                                             'nJet':1,'nbJet':1,'suffix':MuSelObjResolved1b0Wj.selName,
                                             'is_MC':self.is_MC, 'jpaArg':'Res0b'})
                    
            for channelDict in ChannelDictListR:
                if any(['Res2b2Wj','Res2b1Wj','Res1b2Wj','Res2b0Wj','Res1b1Wj','Res1b0Wj']) in jetplot_level:
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
                    ##plots.extend(makeHighLevelPlotsResolved(**{k:channelDict[k] for k in ResolvedKeys},HLL=self.HLL))
                
                inputsCommon  = returnCommonInputs_Resolved(self)
                inputsClassic = returnClassicInputs_Resolved(self, **{k:channelDict[k] for k in ClassicInputKeys})
                inputsLBN     = returnLBNInputs_Resolved(self, **{k:channelDict[k] for k in LBNInputKeys})
                input_names   = [key[0] for key in inputsClassic.keys()] + ['LBN_inputs','eventnr']
                inputs_array  = [op.array("double",val) for val in inputStaticCast(inputsCommon,"float")]
                inputs_array.append(op.array("double",*inputStaticCast(inputsClassic,"float")))
                inputs_array.append(op.array("double",*inputStaticCast(inputsLBN,"float")))
                inputs_array.append(op.array("long",*inputStaticCast(inputsEventNr,"long")))

                DNN_SM  = op.mvaEvaluator (path_model_resolved_SM, mvaType='Tensorflow',otherArgs=(input_names, output_name)) 
                DNN_BSM = op.mvaEvaluator(path_model_resolved_BSM, mvaType='Tensorflow',otherArgs=(input_names, output_name)) 
                DNN_Score_SM  = DNN_SM (*inputs_array) 
                DNN_Score_BSM = DNN_BSM (*inputs_array)

                selObj = channelDict['selObj']
                selObj_1b = makeExclusiveLooseResolvedJetComboSelection(self, selObj, 1, copy_sel=True)
                selObj_2b = makeExclusiveLooseResolvedJetComboSelection(self, selObj, 2, copy_sel=True)
                
                selObjNodesDict_SM     = makeDNNOutputNodesSelections(self,selObj,   DNN_Score_SM, suffix='_SM_')
                selObjNodesDict_SM_1b  = makeDNNOutputNodesSelections(self,selObj_1b,DNN_Score_SM, suffix='_SM_')
                selObjNodesDict_SM_2b  = makeDNNOutputNodesSelections(self,selObj_2b,DNN_Score_SM, suffix='_SM_')

                selObjNodesDict_BSM    = makeDNNOutputNodesSelections(self,selObj,   DNN_Score_BSM, suffix='_BSM_')
                selObjNodesDict_BSM_1b = makeDNNOutputNodesSelections(self,selObj_1b,DNN_Score_BSM, suffix='_BSM_')
                selObjNodesDict_BSM_2b = makeDNNOutputNodesSelections(self,selObj_2b,DNN_Score_BSM, suffix='_BSM_')

                logger.info('Filling SM DNN responses')
                plots.extend(makeDoubleLeptonMachineLearningOutputPlots(selObjNodesDict_SM,    DNN_Score_SM, self.nodes,channel=channelDict['channel']))
                plots.extend(makeDoubleLeptonMachineLearningOutputPlots(selObjNodesDict_SM_1b, DNN_Score_SM, self.nodes,channel=channelDict['channel']))
                plots.extend(makeDoubleLeptonMachineLearningOutputPlots(selObjNodesDict_SM_2b, DNN_Score_SM, self.nodes,channel=channelDict['channel']))
                logger.info('Filling BSM DNN responses')
                plots.extend(makeDoubleLeptonMachineLearningOutputPlots(selObjNodesDict_BSM,    DNN_Score_BSM, self.nodes,channel=channelDict['channel']))
                plots.extend(makeDoubleLeptonMachineLearningOutputPlots(selObjNodesDict_BSM_1b, DNN_Score_BSM, self.nodes,channel=channelDict['channel']))
                plots.extend(makeDoubleLeptonMachineLearningOutputPlots(selObjNodesDict_BSM_2b, DNN_Score_BSM, self.nodes,channel=channelDict['channel']))

        # ========================== JPA Boosted Categories ========================= #
        if any(item in boosted_args for item in jetsel_level):
            ChannelDictListB = []
            FatJetKeys  = ['channel','sel','jet1','jet2','jet3','jet4','has1fat1slim','has1fat2slim','suffix']

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

                if not self.args.OnlyYield:
                    ChannelDictListB.append({'channel':'El','selObj':ElSelObjBoostedHbb2Wj, 'sel':ElSelObjBoostedHbb2Wj.sel,
                                             'lep':ElColl[0],'met':self.corrMET,'jpaSelectedJets':ElSelObjBoostedHbb2WjJets,
                                             'jet1':self.ak8BJets[0],'jet2':None,'jet3':ElSelObjBoostedHbb2WjJets[0],'jet4':ElSelObjBoostedHbb2WjJets[1],
                                             'has1fat1slim':False,'has1fat2slim':True,'bothAreFat':False,
                                             'suffix':ElSelObjBoostedHbb2Wj.selName,
                                             'is_MC':self.is_MC, 'jpaArg':'Hbb2Wj'})
                    ChannelDictListB.append({'channel':'Mu','selObj':MuSelObjBoostedHbb2Wj,'sel':MuSelObjBoostedHbb2Wj.sel,
                                             'lep':MuColl[0],'met':self.corrMET,'jpaSelectedJets':MuSelObjBoostedHbb2WjJets,
                                             'jet1':self.ak8BJets[0],'jet2':None,'jet3':MuSelObjBoostedHbb2WjJets[0],'jet4':MuSelObjBoostedHbb2WjJets[1],
                                             'has1fat1slim':False,'has1fat2slim':True,'bothAreFat':False,
                                             'suffix':MuSelObjBoostedHbb2Wj.selName,
                                             'is_MC':self.is_MC, 'jpaArg':'Hbb2Wj'})
                    
            if "Hbb1Wj" in jetplot_level or "Boosted" in jetplot_level:
                print ('...... JPA : Hbb1Wj Node Selection')
                ElSelObjBoostedHbb1Wj        = ElBoostedSelObjJetsIdxPerJpaNodeDict.get('Hbb1Wj')[0]
                ElSelObjBoostedHbb1WjJets    = ElBoostedSelObjJetsIdxPerJpaNodeDict.get('Hbb1Wj')[1]
                MuSelObjBoostedHbb1Wj        = MuBoostedSelObjJetsIdxPerJpaNodeDict.get('Hbb1Wj')[0]
                MuSelObjBoostedHbb1WjJets    = MuBoostedSelObjJetsIdxPerJpaNodeDict.get('Hbb1Wj')[1]

                if not self.args.OnlyYield:
                    ChannelDictListB.append({'channel':'El','selObj':ElSelObjBoostedHbb1Wj,'sel':ElSelObjBoostedHbb1Wj.sel,
                                             'lep':ElColl[0],'met':self.corrMET,'jpaSelectedJets':ElSelObjBoostedHbb1WjJets,
                                             'jet1':self.ak8BJets[0],'jet2':None,'jet3':ElSelObjBoostedHbb1WjJets[0],'jet4':None,
                                             'has1fat1slim':True,'has1fat2slim':False,'bothAreFat':False,
                                             'suffix':ElSelObjBoostedHbb1Wj.selName,
                                             'is_MC':self.is_MC, 'jpaArg':'Hbb1Wj'})
                    ChannelDictListB.append({'channel':'Mu','selObj':MuSelObjBoostedHbb1Wj,'sel':MuSelObjBoostedHbb1Wj.sel,
                                             'lep':MuColl[0],'met':self.corrMET,'jpaSelectedJets':MuSelObjBoostedHbb1WjJets,
                                             'jet1':self.ak8BJets[0],'jet2':None,'jet3':MuSelObjBoostedHbb1WjJets[0],'jet4':None,
                                             'has1fat1slim':True,'has1fat2slim':False,'bothAreFat':False,
                                             'suffix':MuSelObjBoostedHbb1Wj.selName,
                                             'is_MC':self.is_MC, 'jpaArg':'Hbb1Wj'})
                    
            if "Hbb0Wj" in jetplot_level or "Boosted" in jetplot_level:
                print ('...... JPA : Hbb0Wj Node Selection')
                ElSelObjBoostedHbb0Wj        = ElBoostedSelObjJetsIdxPerJpaNodeDict.get('Hbb0Wj')[0]
                MuSelObjBoostedHbb0Wj        = MuBoostedSelObjJetsIdxPerJpaNodeDict.get('Hbb0Wj')[0]

                if not self.args.OnlyYield:
                    ChannelDictListB.append({'channel':'El','selObj':ElSelObjBoostedHbb0Wj,'sel':ElSelObjBoostedHbb0Wj.sel,
                                             'lep':ElColl[0],'met':self.corrMET,'jpaSelectedJets':None,
                                             'jet1':self.ak8BJets[0],'jet2':None,'jet3':None,'jet4':None,
                                             'has1fat1slim':False,'has1fat2slim':False,'bothAreFat':False,
                                             'suffix':ElSelObjBoostedHbb0Wj.selName,
                                             'is_MC':self.is_MC, 'jpaArg':'Hbb0Wj'})
                    ChannelDictListB.append({'channel':'Mu','selObj':MuSelObjBoostedHbb0Wj,'sel':MuSelObjBoostedHbb0Wj.sel,
                                             'lep':MuColl[0],'met':self.corrMET,'jpaSelectedJets':None,
                                             'jet1':self.ak8BJets[0],'jet2':None,'jet3':None,'jet4':None,
                                             'has1fat1slim':False,'has1fat2slim':False,'bothAreFat':False,
                                             'suffix':MuSelObjBoostedHbb0Wj.selName,
                                             'is_MC':self.is_MC, 'jpaArg':'Hbb0Wj'})
                    
            for channelDict in ChannelDictListB:
                if any(['Hbb2Wj','Hbb1Wj']) in jetplot_level:
                    # Dilepton #
                    plots.extend(makeSinleptonPlots(**{k:channelDict[k] for k in LeptonKeys}))
                    # Number of jets #
                    ##plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**FatJetsN))
                    ##plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**SlimJetsN))
                    # Ak8 Jets #
                    plots.extend(makeSingleLeptonAk8JetsPlots(**{k:channelDict[k] for k in FatJetKeys},nMedBJets=self.nMediumBTaggedSubJets, HLL=self.HLL))
                    # MET #
                    plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))
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
                
        #----- Add the Yield plots -----#
        plots.append(self.yields)
        #plots.extend(cutFlowPlots)
        return plots

