import os
import sys
from copy import copy

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
from variableMakerSL_Basic import *
from bamboo.root import gbl
import ROOT

#===============================================================================================#
#                                 SkimmerHHtobbWW                                               #
#===============================================================================================#
class SkimmerNanoHHtobbWWBasicSL(BaseNanoHHtobbWW,SkimmerModule):
    """ Plotter module: HH->bbW(->e/µ nu)W(->e/µ nu) histograms from NanoAOD """
    def __init__(self, args):
        super(SkimmerNanoHHtobbWWBasicSL, self).__init__(args)

    def initialize(self):
        super(SkimmerNanoHHtobbWWBasicSL, self).initialize(True) # avoids doing the pseudo-data for skimmer

    def defineSkimSelection(self, t, noSel, sample=None, sampleCfg=None): 
        noSel = super(SkimmerNanoHHtobbWWBasicSL,self).prepareObjects(t, noSel, sample, sampleCfg, "SL", forSkimmer=True)
        era = sampleCfg['era'] 
        # Initialize varsToKeep dict #
        varsToKeep = dict()  
        if not self.inclusive_sel:
            jet_level    = ["Res2b2Wj","Res1b3Wj","Res2b1Wj","Res1b2Wj","Hbb2Wj","Hbb1Wj"]
            # Only one lepton_level must be in args and Only one jet_level must be in args
            if [boolean for (level,boolean) in self.args.__dict__.items() if level in jet_level].count(True) != 1:
                raise RuntimeError("Only one of the jet arguments must be used, check --help")

            if self.args.Channel not in ["El","Mu"]:
                raise RuntimeError("Channel must be either 'El' or 'Mu'")            
            
            #----- Lepton selection -----#
            # Args are passed within the self #
            #ElSelObj, MuSelObj = makeSingleLeptonSelection(self,noSel,use_dd=False)
            ElSelObj, MuSelObj = makeSingleLeptonSelection(self,noSel,use_dd=False,fake_selection=self.args.FakeCR)
            if self.args.Channel == "El":
                selObj = ElSelObj
                lepton = self.electronsTightSel[0]
                
            if self.args.Channel == "Mu":
                selObj = MuSelObj
                lepton = self.muonsTightSel[0]
            

            #----- Jet selection -----#
            # Since the selections in one line, we can use the non copy option of the selection to modify the selection object internally
            if any([self.args.__dict__[item] for item in ["Res2b2Wj","Res1b3Wj"]]):
                makeCoarseResolvedSelection(self,selObj,4)
                if self.args.Res2b2Wj:
                    makeExclusiveResolvedJetComboSelection(self,selObj,2,4)
                    print('Resolved 2b2Wj')
                    jet1 = self.ak4BJets[0]
                    jet2 = self.ak4BJets[1]
                    # simple basic reco method
                    WjjPairs = chooseWjj(self, self.ak4BJets[0], self.ak4BJets[1], lepton, self.makeJetPairs(self.ak4LightJetsByPt), isResolved=True, isBoosted=False)
                    jet3 = WjjPairs[0][0]
                    jet4 = WjjPairs[0][1]
                    catarg = 'Res2b2Wj'
                    
                if self.args.Res1b3Wj:
                    makeExclusiveResolvedJetComboSelection(self,selObj,1,4)
                    print('Resolved 1b3Wj')
                    jet1 = self.ak4BJets[0]
                    jet2 = self.ak4LightJetsByBtagScore[0]
                    WjjPairs = chooseWjj(self, self.ak4BJets[0], self.ak4LightJetsByBtagScore[0], lepton, self.makeJetPairs(self.remainingJets), 
                                         isResolved=True, isBoosted=False)
                    jet3 = WjjPairs[0][0]
                    jet4 = WjjPairs[0][1]
                    catarg = 'Res1b3Wj'

            if any([self.args.__dict__[item] for item in ["Res2b1Wj","Res1b2Wj"]]):
                makeCoarseResolvedSelection(self,selObj,3)
                if self.args.Res2b1Wj:
                    makeExclusiveResolvedJetComboSelection(self,selObj,2,3)
                    print('Resolved 2b1Wj')
                    jet1  = self.ak4BJets[0]
                    jet2  = self.ak4BJets[1]
                    jet3  = self.ak4LightJetsByPt[0]
                    jet4  = None
                    catarg = 'Res2b1Wj'

                if self.args.Res1b2Wj:
                    makeExclusiveResolvedJetComboSelection(self,selObj,1,3)
                    print('Resolved 1b2Wj')
                    jet1  = self.ak4BJets[0]
                    jet2  = self.ak4LightJetsByBtagScore[0]
                    jet3  = self.remainingJets[0]
                    jet4  = None
                    catarg = 'Res1b2Wj'

            if any([self.args.__dict__[item] for item in ["Hbb2Wj","Hbb1Wj"]]):
                makeBoostedSelection(self,selObj)
                fatjet  = self.ak8BJets[0]
                if self.args.Hbb1Wj:
                    makeSemiBoostedHbbSelection(self,selObj,1)
                    print('Resolved Hbb1Wj')
                    jet3 = self.ak4JetsCleanedFromAk8b[0]
                    jet4 = None
                    catarg = 'Hbb1Wj'

                if self.args.Hbb2Wj:
                    makeSemiBoostedHbbSelection(self,selObj,2)
                    print('Resolved Hbb2Wj')
                    WjjPairs = chooseWjj(self, self.ak8BJets[0], None, lepton, self.makeJetPairs(self.ak4JetsCleanedFromAk8b), isResolved=False, isBoosted=True)
                    jet3 = WjjPairs[0][0]
                    jet4 = WjjPairs[0][1]
                    catarg = 'Hbb2Wj'

        else:
            noSel = self.beforeJetselection(noSel)
            

        #---------------------------------------------------------------------------------------# 
        #                                    Selection tree                                     #
        #---------------------------------------------------------------------------------------#
        #----- EVT variables -----#
        varsToKeep["event"]     = None # Already in tree                                               
        varsToKeep["run"]       = None # Already in tree          
        varsToKeep["ls"]        = t.luminosityBlock

        if any([self.args.__dict__[item] for item in ["Res2b2Wj","Res1b3Wj","Res2b1Wj","Res1b2Wj"]]):
            classicInputs_Resolved = returnClassicInputs_Resolved (self = self,
                                                                   lepton = lepton,
                                                                   jet1=jet1, jet2=jet2, 
                                                                   jet3=jet3, jet4=jet4,
                                                                   category = catarg)
            inputs_resolved        = {**classicInputs_Resolved}
            
            for (varname,_,_),var in inputs_resolved.items():
                varsToKeep[varname] = var

        if any([self.args.__dict__[item] for item in ["Hbb2Wj","Hbb1Wj"]]):
            classicInputs_Boosted = returnClassicInputs_Boosted (self = self,
                                                                 lepton = lepton,
                                                                 jet3=jet3,jet4=jet4,
                                                                 category = catarg)
            inputs_boosted        = {**classicInputs_Boosted}

            for (varname,_,_),var in inputs_boosted.items():
                varsToKeep[varname] = var


        #----- Additional variables -----#
        varsToKeep["MC_weight"]    = t.genWeight
        varsToKeep['total_weight'] = selObj.sel.weight

        #return leptonSel.sel, varsToKeep
        return selObj.sel, varsToKeep
