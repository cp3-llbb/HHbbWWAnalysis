import os
import sys
from operator import mul
from functools import reduce

from bamboo.analysismodules import SkimmerModule
from bamboo import treefunctions as op
from bamboo.analysisutils import makePileupWeight

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)))) # Add scripts in this directory
from BaseHHtobbWW import BaseNanoHHtobbWW
from selectionDef import *

#===============================================================================================#
#                                 SkimmerHHtobbWW                                               #
#===============================================================================================#
class SkimmerHME(BaseNanoHHtobbWW,SkimmerModule):
    """ Plotter module: HH->bbW(->e/µ nu)W(->e/µ nu) histograms from NanoAOD """
    def __init__(self, args):
        super(SkimmerHME, self).__init__(args)

    def initialize(self):
        super(SkimmerHME, self).initialize(True) # avoids doing the pseudo-data for skimmer


    def defineSkimSelection(self, t, noSel, sample=None, sampleCfg=None): 
        noSel = super(SkimmerHME,self).prepareObjects(t, noSel, sample, sampleCfg, "DL", forSkimmer=True)
            # For the Skimmer, SF must not use defineOnFirstUse -> segmentation fault

        era = sampleCfg['era'] 

        # Initialize varsToKeep dict #
        varsToKeep = dict()  

        #---------------------------------------------------------------------------------------# 
        #                                     Selections                                        #
        #---------------------------------------------------------------------------------------#
        #----- Check arguments -----#
        jet_level = ["Ak4","Ak8","Resolved0Btag","Resolved1Btag","Resolved2Btag","Boosted0Btag","Boosted1Btag"] # Only one must be in args

        if [boolean for (level,boolean) in self.args.__dict__.items() if level in jet_level].count(True) != 1:
            raise RuntimeError("Only one of the jet arguments must be used, check --help")
        if self.args.Channel not in ["ElEl","MuMu","ElMu"]:
            raise RuntimeError("Channel must be either 'ElEl', 'MuMu' or 'ElMu'")

        #----- Lepton selection -----#
        # Args are passed within the self #
        ElElSelObj,MuMuSelObj,ElMuSelObj = makeDoubleLeptonSelection(self,noSel,use_dd=False,fake_selection=self.args.FakeCR)

        if self.args.Channel == "ElEl":
            selObj = ElElSelObj
            dilepton = self.ElElFakeSel[0]
        elif self.args.Channel == "MuMu":
            selObj = MuMuSelObj
            dilepton = self.MuMuFakeSel[0]
        elif self.args.Channel == "ElMu":
            selObj = ElMuSelObj
            dilepton = self.ElMuFakeSel[0]
        else:
            raise RuntimeError("You need to specify --Channel")

        #----- HME -----#
        if self.args.analysis != 'res':
            raise RuntimeError("This part of the Skimmer is only for resonant")

        #----- Apply jet corrections -----#
        ElElSelObj.sel = self.beforeJetselection(ElElSelObj.sel,'ElEl')
        MuMuSelObj.sel = self.beforeJetselection(MuMuSelObj.sel,'MuMu')
        ElMuSelObj.sel = self.beforeJetselection(ElMuSelObj.sel,'ElMu')

        #----- Jet selection -----#
        # Since the selections in one line, we can use the non copy option of the selection to modify the selection object internally
        if any([self.args.__dict__[item] for item in ["Ak4","Resolved0Btag","Resolved1Btag","Resolved2Btag"]]):
            makeAtLeastTwoAk4JetSelection(self,selObj,use_dd=False) 
        if any([self.args.__dict__[item] for item in ["Ak8","Boosted0Btag","Boosted1Btag"]]):
           makeAtLeastOneAk8JetSelection(self,selObj,use_dd=False) 
        if self.args.Resolved0Btag:
            makeExclusiveResolvedNoBtagSelection(self,selObj,use_dd=False)
        if self.args.Resolved1Btag:
            makeExclusiveResolvedOneBtagSelection(self,selObj,use_dd=False,dy_selection=self.args.DYCR)
        if self.args.Resolved2Btag:
            makeExclusiveResolvedTwoBtagsSelection(self,selObj,use_dd=False,dy_selection=self.args.DYCR)
        if self.args.Boosted0Btag:
            makeInclusiveBoostedNoBtagSelection(self,selObj,use_dd=False)
        if self.args.Boosted1Btag:
            makeInclusiveBoostedOneBtagSelection(self,selObj,use_dd=False,dy_selection=self.args.DYCR)


        #----- HME ----#
        if self.args.Resolved1Btag or self.args.Resolved2Btag:
            HME,HME_eff = self.computeResolvedHMEAfterLeptonSelections(
                                sel   = selObj.sel,
                                l1    = dilepton[0],
                                l2    = dilepton[1],
                                bjets = self.ak4JetsByBtagScore,
                                met   = self.corrMET)
        elif self.args.Boosted1Btag:
            HME,HME_eff = self.computeBoostedHMEAfterLeptonSelections(
                                sel     = selObj.sel,
                                l1      = dilepton[0],
                                l2      = dilepton[1],
                                fatjets = self.ak8BJets,
                                met     = self.corrMET)
        else:
            raise RuntimeError("Wrong category for resonant HME computations")

        varsToKeep["HME"] = HME
        varsToKeep["HME_eff"] = HME_eff

        #----- Additional variables -----#
        varsToKeep['era']           = op.c_int(int(self.era))
        varsToKeep["event"]         = None # Already in tree
        varsToKeep["run"]           = None # Already in tree 
        varsToKeep["ls"]            = t.luminosityBlock

        return selObj.sel, varsToKeep

    ### PostProcess ###
    def postProcess(self, taskList, config=None, workdir=None, resultsdir=None):
        super(SkimmerHME, self).postProcess(taskList, config, workdir, resultsdir, forSkimmer=True)

