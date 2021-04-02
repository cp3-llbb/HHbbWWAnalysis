import os
import sys
from copy import copy
from operator import mul
from functools import reduce

from itertools import chain

import logging
logger = logging.getLogger(__name__) 

from bamboo.analysismodules import SkimmerModule
from bamboo import treefunctions as op
from bamboo.analysisutils import makePileupWeight

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)))) # Add scripts in this directory
from BaseHHtobbWW import BaseNanoHHtobbWW
from selectionDef import *
from highlevelLambdas import *
from JPA import *
from bamboo.root import gbl
import ROOT
from DDHelper import DataDrivenFake, DataDrivenDY
from variableMakerSL import *
#from bamboo.root import loadHeader
#loadHeader("jpa.h")

#===============================================================================================#
#                                 SkimmerHHtobbWW                                               #
#===============================================================================================#
class SkimmerNanoHHtobbWWSL(BaseNanoHHtobbWW,SkimmerModule):
    """ Plotter module: HH->bbW(->e/µ nu)W(->e/µ nu) histograms from NanoAOD """
    def __init__(self, args):
        super(SkimmerNanoHHtobbWWSL, self).__init__(args)

    def initialize(self):
        super(SkimmerNanoHHtobbWWSL, self).initialize(True) # avoids doing the pseudo-data for skimmer

    def defineSkimSelection(self, t, noSel, sample=None, sampleCfg=None): 
        noSel = super(SkimmerNanoHHtobbWWSL,self).prepareObjects(t, noSel, sample, sampleCfg, "SL", forSkimmer=True)
        # For the Skimmer, SF must not use defineOnFirstUse -> segmentation fault

        era = sampleCfg['era'] 

        # Initialize varsToKeep dict #
        varsToKeep = dict()  

        #---------------------------------------------------------------------------------------# 
        #                                     Selections                                        #
        #---------------------------------------------------------------------------------------#
        # keep the exact same order of nodes as mentioned in respective JPA model xml files
        ResolvedJPANodeList = ['2b2Wj','2b1Wj','1b2Wj','2b0Wj','1b1Wj','1b0Wj','0b']
        BoostedJPANodeList  = ['Hbb2Wj','Hbb1Wj','Hbb0Wj']

        # JPA Models
        basepath = os.path.join(os.path.abspath(os.path.dirname(__file__)),'MachineLearning','ml-models','JPA_Loose_ttH')
        resolvedModelDict = getResolvedJpaModelDict(basepath, ResolvedJPANodeList, era)
        boostedModelDict  = getBoostedJpaModelDict(basepath, BoostedJPANodeList, era)

        if not self.inclusive_sel:
            #----- Check arguments -----#
            jet_level    = ["Ak4","Ak8","Res2b2Wj","Res2b1Wj","Res2b0Wj","Res1b2Wj","Res1b1Wj","Res1b0Wj","Res0b","Hbb2Wj","Hbb1Wj","Hbb0Wj"]

            if [boolean for (level,boolean) in self.args.__dict__.items() if level in jet_level].count(True) != 1:
                raise RuntimeError("Only one of the jet arguments must be used, check --help")

            if self.args.Channel not in ["El","Mu"]:
                raise RuntimeError("Channel must be either 'El' or 'Mu'")            
            
            #----- Lepton selection -----#
            ElSelObj, MuSelObj = makeSingleLeptonSelection(self,noSel,use_dd=False,fake_selection=self.args.FakeCR)

            if self.args.Channel is None:
                raise RuntimeError("You need to specify --Channel")
            if self.args.Channel == "El":
                selObj = ElSelObj
                lep = self.electronsTightSel[0]
            if self.args.Channel == "Mu":
                selObj = MuSelObj
                lep = self.muonsTightSel[0]
            
            #----- Jet selection -----#
            if any([self.args.__dict__[item] for item in ["Ak4","Res2b2Wj","Res2b1Wj","Res2b0Wj","Res1b2Wj","Res1b1Wj","Res1b0Wj","Res0b"]]):
                makeResolvedSelection(self,selObj)
                if self.args.Channel == "El":
                    print('... Resolved :: El Channel')
                    L1out,L2out,selObjAndJetsPerJpaCatDict = findJPACategoryResolved (self, selObj, lep, 
                                                                                      self.muonsPreSel, self.electronsPreSel, 
                                                                                      self.ak4Jets, self.ak4BJetsLoose, self.ak4BJets, self.corrMET, 
                                                                                      resolvedModelDict, t.event, self.HLL, ResolvedJPANodeList, 
                                                                                      plot_yield=False)
                                        
                if self.args.Channel == "Mu":
                    print('... Resolved :: Mu Channel')
                    L1out,L2out,selObjAndJetsPerJpaCatDict = findJPACategoryResolved (self, selObj, lep, 
                                                                                      self.muonsPreSel, self.electronsPreSel, 
                                                                                      self.ak4Jets, self.ak4BJetsLoose, self.ak4BJets, self.corrMET, 
                                                                                      resolvedModelDict, t.event, self.HLL, ResolvedJPANodeList, 
                                                                                      plot_yield=False)
                    
            if any([self.args.__dict__[item] for item in ["Ak8","Hbb2Wj","Hbb1Wj","Hbb0Wj"]]):
                makeBoostedSelection(self,selObj)
                if self.args.Channel == "El":
                    print('... Boosted :: El Channel')
                    L1out,L2out,selObjAndJetsPerJpaCatDict = findJPACategoryBoosted (self, selObj, lep, 
                                                                                     self.muonsPreSel, self.electronsPreSel, 
                                                                                     self.ak8BJets, self.ak4JetsCleanedFromAk8b, self.ak4BJetsLoose, 
                                                                                     self.ak4BJets, self.corrMET, boostedModelDict, t.event, self.HLL, 
                                                                                     BoostedJPANodeList, plot_yield=False)

                if self.args.Channel == "Mu":
                    print('... Boosted :: Mu Channel')
                    L1out,L2out,selObjAndJetsPerJpaCatDict = findJPACategoryBoosted (self, selObj, lep, 
                                                                                     self.muonsPreSel, self.electronsPreSel, 
                                                                                     self.ak8BJets, self.ak4JetsCleanedFromAk8b, self.ak4BJetsLoose, 
                                                                                     self.ak4BJets, self.corrMET, boostedModelDict, t.event, self.HLL, 
                                                                                     BoostedJPANodeList, plot_yield=False)
                    
            if self.args.Res2b2Wj:
                print("...... 2b2Wj")
                selObj   = selObjAndJetsPerJpaCatDict.get('2b2Wj')[0]
                jpaJets  = selObjAndJetsPerJpaCatDict.get('2b2Wj')[1]
                jpaarg   = "Res2b2Wj"

            if self.args.Res2b1Wj:
                print("...... 2b1Wj")
                selObj  = selObjAndJetsPerJpaCatDict.get('2b1Wj')[0]
                jpaJets = selObjAndJetsPerJpaCatDict.get('2b1Wj')[1]
                jpaarg   = "Res2b1Wj"

            if self.args.Res1b2Wj:
                print("...... 1b2Wj")
                selObj  = selObjAndJetsPerJpaCatDict.get('1b2Wj')[0]
                jpaJets = selObjAndJetsPerJpaCatDict.get('1b2Wj')[1]
                jpaarg   = "Res1b2Wj"

            if self.args.Res2b0Wj:
                print("...... 2b0Wj")
                selObj  = selObjAndJetsPerJpaCatDict.get('2b0Wj')[0]
                jpaJets = selObjAndJetsPerJpaCatDict.get('2b0Wj')[1]
                jpaarg   = "Res2b0Wj"

            if self.args.Res1b1Wj:
                print("...... 1b1Wj")
                selObj  = selObjAndJetsPerJpaCatDict.get('1b1Wj')[0]
                jpaJets = selObjAndJetsPerJpaCatDict.get('1b1Wj')[1]
                jpaarg   = "Res1b1Wj"

            if self.args.Res1b0Wj:
                print("...... 1b0Wj")
                selObj  = selObjAndJetsPerJpaCatDict.get('1b0Wj')[0]
                jpaJets = selObjAndJetsPerJpaCatDict.get('1b0Wj')[1]
                jpaarg   = "Res1b0Wj"

            if self.args.Res0b:
                print("...... 0b")
                selObj  = selObjAndJetsPerJpaCatDict.get('0b')[0]
                jpaJets = None 
                jpaarg   = "Res0b"

            #######################################
            # Hbb2Wj : jet1 jet2   jet3   jet4
            # Hbb1Wj : jet1 jet2   jet3   jet4=0
            # Hbb0Wj : jet1 jet2   jet3=0 jet4=0
            #######################################
            if self.args.Hbb2Wj:
                print("...... Hbb2Wj")
                selObj  = selObjAndJetsPerJpaCatDict.get('Hbb2Wj')[0]
                jpaJets = selObjAndJetsPerJpaCatDict.get('Hbb2Wj')[1]
                jpaarg   = "Hbb2Wj"

            if self.args.Hbb1Wj:
                print("...... Hbb1Wj")
                selObj  = selObjAndJetsPerJpaCatDict.get('Hbb1Wj')[0]
                jpaJets = selObjAndJetsPerJpaCatDict.get('Hbb1Wj')[1]
                jpaarg   = "Hbb1Wj"

            if self.args.Hbb0Wj:
                print("...... Hbb0Wj")
                selObj  = selObjAndJetsPerJpaCatDict.get('Hbb0Wj')[0]
                jpaJets = None
                jpaarg   = "Hbb0Wj"

        else:
            noSel = self.beforeJetselection(noSel)
            
        #---------------------------------------------------------------------------------------# 
        #                                 Synchronization tree                                  #
        #---------------------------------------------------------------------------------------#
        # Removed !!! 
        #---------------------------------------------------------------------------------------# 
        #                                    Selection tree                                     #
        #---------------------------------------------------------------------------------------#
        #----- EVT variables -----#
        varsToKeep["event"]     = None # Already in tree                                               
        varsToKeep["run"]       = None # Already in tree             
        varsToKeep["ls"]        = t.luminosityBlock

        genVarsList = self.HLL.comp_cosThetaSbetBeamAndHiggs(t.GenPart)
        varsToKeep['mHH_gen']        = genVarsList[0]
        varsToKeep['consTheta1_gen'] = genVarsList[1]
        varsToKeep['consTheta2_gen'] = genVarsList[2]

        if any([self.args.__dict__[item] for item in ["Ak4","Res2b2Wj","Res2b1Wj","Res2b0Wj","Res1b2Wj","Res1b1Wj","Res1b0Wj","Res0b"]]):
            commonInputs_Resolved  = returnCommonInputs_Resolved  (self = self)
            classicInputs_Resolved = returnClassicInputs_Resolved (self = self, 
                                                                   lepton = lep, 
                                                                   jpaSelectedJets = jpaJets, 
                                                                   L1out = L1out, L2out = L2out, 
                                                                   jpaArg = jpaarg)
            inputs_resolved        = {**commonInputs_Resolved, **classicInputs_Resolved}

            for (varname,_,_),var in inputs_resolved.items():
                varsToKeep[varname] = var 

        #----- Fatjet variables -----#
        if any([self.args.__dict__[item] for item in ["Ak8","Hbb2Wj","Hbb1Wj","Hbb0Wj"]]):
            commonInputs_Boosted  = returnCommonInputs_Boosted  (self = self)
            classicInputs_Boosted = returnClassicInputs_Boosted (self = self, 
                                                                 lepton = lep, 
                                                                 jpaSelectedJets = jpaJets, 
                                                                 L1out = L1out, L2out = L2out,
                                                                 jpaArg = jpaarg)
            inputs_boosted        = {**commonInputs_Boosted, **classicInputs_Boosted}
            
            for (varname,_,_),var in inputs_boosted.items():
                varsToKeep[varname] = var 

        #----- Additional variables -----#
        varsToKeep["MC_weight"]    = t.genWeight
        varsToKeep['total_weight'] = selObj.sel.weight

        #return leptonSel.sel, varsToKeep
        return selObj.sel, varsToKeep
