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
from highlevelLambdas import *

#===============================================================================================#
#                                 SkimmerHHtobbWW                                               #
#===============================================================================================#
class SkimmerMEMNanoHHtobbWWDL(BaseNanoHHtobbWW,SkimmerModule):
    """ Plotter module: HH->bbW(->e/µ nu)W(->e/µ nu) histograms from NanoAOD """
    def __init__(self, args):
        super(SkimmerMEMNanoHHtobbWWDL, self).__init__(args)

    def initialize(self):
        super(SkimmerMEMNanoHHtobbWWDL, self).initialize(True) # avoids doing the pseudo-data for skimmer


    def defineSkimSelection(self, t, noSel, sample=None, sampleCfg=None): 
        noSel = super(SkimmerMEMNanoHHtobbWWDL,self).prepareObjects(t, noSel, sample, sampleCfg, "DL", forSkimmer=True)
            # For the Skimmer, SF must not use defineOnFirstUse -> segmentation fault

        era = sampleCfg['era'] 

        # Initialize varsToKeep dict #
        varsToKeep = dict()  

        if self.inclusive_sel:
            raise RuntimeError("Inclusive analysis not possible")

        #---------------------------------------------------------------------------------------# 
        #                                     Selections                                        #
        #---------------------------------------------------------------------------------------#
        #----- Check arguments -----#
        jet_level = ["Resolved0Btag","Resolved1Btag","Resolved2Btag","Boosted0Btag","Boosted1Btag"] # Only one must be in args

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
        if self.args.Channel == "MuMu":
            selObj = MuMuSelObj
            dilepton = self.MuMuFakeSel[0]
        if self.args.Channel == "ElMu":
            selObj = ElMuSelObj
            dilepton = self.ElMuFakeSel[0]

        #----- Jet selection -----#
        # Since the selections in one line, we can use the non copy option of the selection to modify the selection object internally
        if any([self.args.__dict__[item] for item in ["Ak4","Resolved0Btag","Resolved1Btag","Resolved2Btag"]]):
            makeAtLeastTwoAk4JetSelection(self,selObj,use_dd=False) 
        if any([self.args.__dict__[item] for item in ["Ak8","Boosted0Btag","Boosted1Btag"]]):
           makeAtLeastOneAk8JetSelection(self,selObj,use_dd=False) 
        if self.args.Resolved0Btag:
            makeExclusiveResolvedNoBtagSelection(self,selObj,use_dd=False)
        if self.args.Resolved1Btag:
            makeExclusiveResolvedOneBtagSelection(self,selObj,use_dd=False)
        if self.args.Resolved2Btag:
            makeExclusiveResolvedTwoBtagsSelection(self,selObj,use_dd=False)
        if self.args.Boosted0Btag:
            makeInclusiveBoostedNoBtagSelection(self,selObj,use_dd=False)
        if self.args.Boosted1Btag:
            makeInclusiveBoostedOneBtagSelection(self,selObj,use_dd=False)


        #---------------------------------------------------------------------------------------# 
        #                                    Selection tree                                     #
        #---------------------------------------------------------------------------------------#

        #----- MET variables -----#
        MET = self.corrMET

        varsToKeep['met_pt']  = MET.pt
        varsToKeep['met_phi']  = MET.phi
        varsToKeep['met_E'] = MET.p4.E()
        varsToKeep['met_Px'] = MET.p4.Px()
        varsToKeep['met_Py'] = MET.p4.Py()
        varsToKeep['met_Pz'] = MET.p4.Pz()

        #----- Lepton variables -----#
        if self.args.Channel is None:
            raise RuntimeError("You need to specify --Channel")
        if self.args.Channel == "ElEl": dilepton = self.ElElTightSel[0]
        if self.args.Channel == "MuMu": dilepton = self.MuMuTightSel[0]
        if self.args.Channel == "ElMu": dilepton = self.ElMuTightSel[0]

        varsToKeep["is_SR"] = op.static_cast("UInt_t",op.OR(op.rng_len(self.ElElTightSel)>0,
                                                            op.rng_len(self.MuMuTightSel)>0,
                                                            op.rng_len(self.ElMuTightSel)>0))
        varsToKeep['is_ee'] = op.static_cast("UInt_t",op.rng_len(self.ElElTightSel)>0)
        varsToKeep['is_mm'] = op.static_cast("UInt_t",op.rng_len(self.MuMuTightSel)>0)
        varsToKeep['is_em'] = op.static_cast("UInt_t",op.rng_len(self.ElMuTightSel)>0)
        varsToKeep['resolved1b_tag'] = op.static_cast("UInt_t",op.AND(op.rng_len(self.ak4BJets)==1,op.rng_len(self.ak8BJets)==0))
        varsToKeep['resolved2b_tag'] = op.static_cast("UInt_t",op.AND(op.rng_len(self.ak4BJets)>=2,op.rng_len(self.ak8BJets)==0))
        varsToKeep['boosted1b_tag'] = op.static_cast("UInt_t",op.AND(op.rng_len(self.ak8BJets)>0))

        l1 = dilepton[0]
        l2 = dilepton[1]

        varsToKeep['l1_Px']     = l1.p4.Px()
        varsToKeep['l1_Py']     = l1.p4.Py()
        varsToKeep['l1_Pz']     = l1.p4.Pz()
        varsToKeep['l1_E']      = l1.p4.E()
        varsToKeep['l1_pt']     = l1.pt
        varsToKeep['l1_eta']    = l1.eta
        varsToKeep['l1_phi']    = l1.phi
        varsToKeep['l1_pdgId']  = l1.pdgId
        varsToKeep['l1_charge'] = l1.charge

        varsToKeep['l2_Px']     = l2.p4.Px()
        varsToKeep['l2_Py']     = l2.p4.Py()
        varsToKeep['l2_Pz']     = l2.p4.Pz()
        varsToKeep['l2_E']      = l2.p4.E()
        varsToKeep['l2_pt']     = l2.pt
        varsToKeep['l2_eta']    = l2.eta
        varsToKeep['l2_phi']    = l2.phi
        varsToKeep['l2_pdgId']  = l2.pdgId
        varsToKeep['l2_charge'] = l2.charge

        #----- Jet variables -----#
        jets = self.ak4JetsByBtagScore
        for idx in range(1,5):
            varsToKeep[f'j{idx}_Px']  = op.switch(op.rng_len(jets)>=idx, jets[idx-1].p4.Px(), op.c_float(-9999))
            varsToKeep[f'j{idx}_Py']  = op.switch(op.rng_len(jets)>=idx, jets[idx-1].p4.Py(), op.c_float(-9999))
            varsToKeep[f'j{idx}_Pz']  = op.switch(op.rng_len(jets)>=idx, jets[idx-1].p4.Pz(), op.c_float(-9999))
            varsToKeep[f'j{idx}_E']   = op.switch(op.rng_len(jets)>=idx, jets[idx-1].p4.E(), op.c_float(-9999))
            varsToKeep[f'j{idx}_pt']  = op.switch(op.rng_len(jets)>=idx, jets[idx-1].pt, op.c_float(-9999))
            varsToKeep[f'j{idx}_eta'] = op.switch(op.rng_len(jets)>=idx, jets[idx-1].eta, op.c_float(-9999))
            varsToKeep[f'j{idx}_phi'] = op.switch(op.rng_len(jets)>=idx, jets[idx-1].phi, op.c_float(-9999))
            varsToKeep[f'j{idx}_btag']= op.switch(op.rng_len(jets)>=idx, jets[idx-1].btagDeepFlavB, op.c_float(-9999))


        varsToKeep['n_ak4'] = op.static_cast("UInt_t",op.rng_len(self.ak4Jets))
        varsToKeep['n_ak4_btag'] = op.static_cast("UInt_t",op.rng_len(self.ak4BJets))


        #----- Fatjet variables -----#
        fatjets = self.ak8BJets
        subJet1 = fatjets[0].subJet1
        subJet2 = fatjets[0].subJet2

        varsToKeep['fatj_sub1_Px']   = op.switch(op.rng_len(fatjets)>0, subJet1.p4.Px(), op.c_float(-9999))
        varsToKeep['fatj_sub1_Py']   = op.switch(op.rng_len(fatjets)>0, subJet1.p4.Py(), op.c_float(-9999))
        varsToKeep['fatj_sub1_Pz']   = op.switch(op.rng_len(fatjets)>0, subJet1.p4.Pz(), op.c_float(-9999))
        varsToKeep['fatj_sub1_E']    = op.switch(op.rng_len(fatjets)>0, subJet1.p4.E(), op.c_float(-9999))
        varsToKeep['fatj_sub1_pt']   = op.switch(op.rng_len(fatjets)>0, subJet1.pt, op.c_float(-9999))
        varsToKeep['fatj_sub1_eta']  = op.switch(op.rng_len(fatjets)>0, subJet1.eta, op.c_float(-9999))
        varsToKeep['fatj_sub1_phi']  = op.switch(op.rng_len(fatjets)>0, subJet1.phi, op.c_float(-9999))
        varsToKeep['fatj_sub1_btag'] = op.switch(op.rng_len(fatjets)>0, subJet1.btagDeepB, op.c_float(-9999))

        varsToKeep['fatj_sub2_Px']   = op.switch(op.rng_len(fatjets)>0, subJet2.p4.Px(), op.c_float(-9999))
        varsToKeep['fatj_sub2_Py']   = op.switch(op.rng_len(fatjets)>0, subJet2.p4.Py(), op.c_float(-9999))
        varsToKeep['fatj_sub2_Pz']   = op.switch(op.rng_len(fatjets)>0, subJet2.p4.Pz(), op.c_float(-9999))
        varsToKeep['fatj_sub2_E']    = op.switch(op.rng_len(fatjets)>0, subJet2.p4.E(), op.c_float(-9999))
        varsToKeep['fatj_sub2_pt']   = op.switch(op.rng_len(fatjets)>0, subJet2.pt, op.c_float(-9999))
        varsToKeep['fatj_sub2_eta']  = op.switch(op.rng_len(fatjets)>0, subJet2.eta, op.c_float(-9999))
        varsToKeep['fatj_sub2_phi']  = op.switch(op.rng_len(fatjets)>0, subJet2.phi, op.c_float(-9999))
        varsToKeep['fatj_sub2_btag'] = op.switch(op.rng_len(fatjets)>0, subJet2.btagDeepB, op.c_float(-9999))

        varsToKeep['fatj_Px']             = op.switch(op.rng_len(fatjets)>0, fatjets[0].p4.Px(), op.c_float(-9999))
        varsToKeep['fatj_Py']             = op.switch(op.rng_len(fatjets)>0, fatjets[0].p4.Py(), op.c_float(-9999))
        varsToKeep['fatj_Pz']             = op.switch(op.rng_len(fatjets)>0, fatjets[0].p4.Pz(), op.c_float(-9999))
        varsToKeep['fatj_E']              = op.switch(op.rng_len(fatjets)>0, fatjets[0].p4.E(), op.c_float(-9999))
        varsToKeep['fatj_pt']             = op.switch(op.rng_len(fatjets)>0, fatjets[0].pt, op.c_float(-9999))
        varsToKeep['fatj_eta']            = op.switch(op.rng_len(fatjets)>0, fatjets[0].eta, op.c_float(-9999))
        varsToKeep['fatj_phi']            = op.switch(op.rng_len(fatjets)>0, fatjets[0].phi, op.c_float(-9999))
        varsToKeep['fatj_softdropMass']   = op.switch(op.rng_len(fatjets)>0, fatjets[0].msoftdrop, op.c_float(-9999))
        varsToKeep['fatj_btagDeepB']      = op.switch(op.rng_len(fatjets)>0, fatjets[0].btagDeepB, op.c_float(-9999))
        varsToKeep['fatj_btagHbb']        = op.switch(op.rng_len(fatjets)>0, fatjets[0].btagHbb, op.c_float(-9999))

        varsToKeep['n_ak8'] = op.static_cast("UInt_t",op.rng_len(self.ak8Jets))
        varsToKeep['n_ak8_btag'] = op.static_cast("UInt_t",op.rng_len(self.ak8BJets))

        #----- Additional variables -----#
        if self.is_MC:
            varsToKeep["MC_weight"]         = t.genWeight
        varsToKeep['total_weight']      = selObj.sel.weight
        varsToKeep["event"]             = None # Already in tree
        varsToKeep["run"]               = None # Already in tree 
        varsToKeep["ls"]                = t.luminosityBlock

        return selObj.sel, varsToKeep
