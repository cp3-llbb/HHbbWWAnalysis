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
class SkimmerNanoHHtobbWWDL(BaseNanoHHtobbWW,SkimmerModule):
    """ Plotter module: HH->bbW(->e/µ nu)W(->e/µ nu) histograms from NanoAOD """
    def __init__(self, args):
        super(SkimmerNanoHHtobbWWDL, self).__init__(args)

    def initialize(self):
        super(SkimmerNanoHHtobbWWDL, self).initialize(True) # avoids doing the pseudo-data for skimmer


    def defineSkimSelection(self, t, noSel, sample=None, sampleCfg=None): 
        noSel = super(SkimmerNanoHHtobbWWDL,self).prepareObjects(t, noSel, sample, sampleCfg, "DL", forSkimmer=True)
            # For the Skimmer, SF must not use defineOnFirstUse -> segmentation fault

        era = sampleCfg['era'] 

        #self.datadrivenContributions = {} # Avoid all data-driven estimates

        # Initialize varsToKeep dict #
        varsToKeep = dict()  

        #---------------------------------------------------------------------------------------# 
        #                                     Selections                                        #
        #---------------------------------------------------------------------------------------#
        if not self.inclusive_sel:
            #----- Check arguments -----#
            lepton_level = ["Preselected","Fakeable","Tight","FakeExtrapolation"]               # Only one must be in args
            jet_level = ["Ak4","Ak8","Resolved0Btag","Resolved1Btag","Resolved2Btag","Boosted"] # Only one must be in args

            if [boolean for (level,boolean) in self.args.__dict__.items() if level in lepton_level].count(True) != 1:
                raise RuntimeError("Only one of the lepton arguments must be used, check --help")
            if [boolean for (level,boolean) in self.args.__dict__.items() if level in jet_level].count(True) != 1:
                raise RuntimeError("Only one of the jet arguments must be used, check --help")
            if self.args.Channel not in ["ElEl","MuMu","ElMu"]:
                raise RuntimeError("Channel must be either 'ElEl', 'MuMu' or 'ElMu'")

            #----- Lepton selection -----#
            # Args are passed within the self #
            selLeptonDict = makeDoubleLeptonSelection(self,noSel,use_dd=False)
                # makeDoubleLeptonSelection returns dict -> value is list of three selections for 3 channels 
                # [0] -> we take the first and only key and value because restricted to one lepton selection
            selLeptonList = list(selLeptonDict.values())[0]
            if self.args.Channel == "ElEl":
                selObj = selLeptonList[0] # First item of list is ElEl selection
            if self.args.Channel == "MuMu":
                selObj = selLeptonList[1] # Second item of list is MuMu selection
            if self.args.Channel == "ElMu":
                selObj = selLeptonList[2] # Third item of list is ElMu selection

            #----- Jet selection -----#
            # Since the selections in one line, we can use the non copy option of the selection to modify the selection object internally
            if any([self.args.__dict__[item] for item in ["Ak4","Resolved0Btag","Resolved1Btag","Resolved2Btag"]]):
                makeAtLeastTwoAk4JetSelection(self,selObj,use_dd=False) 
            if any([self.args.__dict__[item] for item in ["Ak8","Boosted"]]):
               makeAtLeastOneAk8JetSelection(self,selObj,use_dd=False) 
            if self.args.Resolved0Btag:
                makeExclusiveResolvedNoBtagSelection(self,selObj,use_dd=False)
            if self.args.Resolved1Btag:
                makeExclusiveResolvedOneBtagSelection(self,selObj,use_dd=False)
            if self.args.Resolved2Btag:
                makeExclusiveResolvedTwoBtagsSelection(self,selObj,use_dd=False)
            if self.args.Boosted:
                makeInclusiveBoostedSelection(self,selObj,use_dd=False)


        #---------------------------------------------------------------------------------------# 
        #                                 Synchronization tree                                  #
        #---------------------------------------------------------------------------------------#
        if self.args.Synchronization:
            # Event variables #
            varsToKeep["event"]             = None # Already in tree
            varsToKeep["run"]               = None # Already in tree 
            varsToKeep["ls"]                = t.luminosityBlock
            varsToKeep["n_presel_mu"]       = op.static_cast("UInt_t",op.rng_len(self.muonsPreSel))
            varsToKeep["n_fakeablesel_mu"]  = op.static_cast("UInt_t",op.rng_len(self.muonsFakeSel))
            varsToKeep["n_mvasel_mu"]       = op.static_cast("UInt_t",op.rng_len(self.muonsTightSel))
            varsToKeep["n_presel_ele"]      = op.static_cast("UInt_t",op.rng_len(self.electronsPreSel))
            varsToKeep["n_fakeablesel_ele"] = op.static_cast("UInt_t",op.rng_len(self.electronsFakeSel))
            varsToKeep["n_mvasel_ele"]      = op.static_cast("UInt_t",op.rng_len(self.electronsTightSel))
            varsToKeep["n_presel_ak4Jet"]   = op.static_cast("UInt_t",op.rng_len(self.ak4Jets))    
            varsToKeep["n_presel_ak8Jet"]   = op.static_cast("UInt_t",op.rng_len(self.ak8BJets))    
            varsToKeep["n_medium_ak4BJet"]  = op.static_cast("UInt_t",op.rng_len(self.ak4BJets))    
            varsToKeep["is_SR"]             = op.static_cast("UInt_t",op.OR(op.rng_len(self.ElElDileptonTightSel)>=1,
                                                                            op.rng_len(self.MuMuDileptonTightSel)>=1,
                                                                            op.rng_len(self.ElMuDileptonTightSel)>=1))
            varsToKeep["is_CR"]             = op.static_cast("UInt_t",op.OR(op.rng_len(self.ElElDileptonFakeExtrapolationSel)>=1,
                                                                            op.rng_len(self.MuMuDileptonFakeExtrapolationSel)>=1,
                                                                            op.rng_len(self.ElMuDileptonFakeExtrapolationSel)>=1))
            varsToKeep["is_ee"]             = op.static_cast("UInt_t",op.OR(op.rng_len(self.ElElDileptonTightSel)>=1, op.rng_len(self.ElElDileptonFakeExtrapolationSel)>=1))
            varsToKeep["is_mm"]             = op.static_cast("UInt_t",op.OR(op.rng_len(self.MuMuDileptonTightSel)>=1, op.rng_len(self.MuMuDileptonFakeExtrapolationSel)>=1))
            varsToKeep["is_em"]             = op.static_cast("UInt_t",op.OR(op.rng_len(self.ElMuDileptonTightSel)>=1, op.rng_len(self.ElMuDileptonFakeExtrapolationSel)>=1))
            varsToKeep["is_resolved"]       = op.switch(op.AND(op.rng_len(self.ak4Jets)>=2,op.rng_len(self.ak4BJets)>=1,op.rng_len(self.ak8BJets)==0), op.c_bool(True), op.c_bool(False))
            varsToKeep["is_boosted"]        = op.switch(op.rng_len(self.ak8BJets)>=1, op.c_bool(True), op.c_bool(False))


            # Triggers #
            varsToKeep['n_leadfakeableSel_ele']     = op.static_cast("UInt_t",op.rng_len(self.leadElectronsFakeSel))
            varsToKeep['n_leadfakeableSel_mu']      = op.static_cast("UInt_t",op.rng_len(self.leadMuonsFakeSel))
            varsToKeep["triggers"]                  = self.triggers
            varsToKeep["triggers_SingleElectron"]   = op.OR(*self.triggersPerPrimaryDataset['SingleElectron'])
            varsToKeep["triggers_SingleMuon"]       = op.OR(*self.triggersPerPrimaryDataset['SingleMuon'])
            varsToKeep["triggers_DoubleElectron"]   = op.OR(*self.triggersPerPrimaryDataset['DoubleEGamma'])
            varsToKeep["triggers_DoubleMuon"]       = op.OR(*self.triggersPerPrimaryDataset['DoubleMuon'])
            varsToKeep["triggers_MuonElectron"]     = op.OR(*self.triggersPerPrimaryDataset['MuonEG'])

            # Muons #
            for i in range(1,3): # 2 leading muons
                varsToKeep["mu{}_pt".format(i)]                    = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].pt, op.c_float(-9999., "float"))
                varsToKeep["mu{}_eta".format(i)]                   = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].eta, op.c_float(-9999.))
                varsToKeep["mu{}_phi".format(i)]                   = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].phi, op.c_float(-9999.))
                varsToKeep["mu{}_E".format(i)]                     = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].p4.E(), op.c_float(-9999., "float"))
                varsToKeep["mu{}_charge".format(i)]                = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].charge, op.c_int(-9999.))
                varsToKeep["mu{}_conept".format(i)]                = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muon_conept[self.muonsPreSel[i-1].idx], op.c_float(-9999.))
                varsToKeep["mu{}_miniRelIso".format(i)]            = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].miniPFRelIso_all, op.c_float(-9999.))
                varsToKeep["mu{}_PFRelIso04".format(i)]            = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].pfRelIso04_all, op.c_float(-9999.))
                varsToKeep["mu{}_jetNDauChargedMVASel".format(i)]  = op.c_float(-9999.)
                varsToKeep["mu{}_jetPtRel".format(i)]              = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].jetPtRelv2, op.c_float(-9999.))
                varsToKeep["mu{}_jetRelIso".format(i)]             = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].jetRelIso, op.c_float(-9999.))
                varsToKeep["mu{}_jetDeepJet".format(i)]            = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].jet.btagDeepFlavB, op.c_float(-9999.))
                varsToKeep["mu{}_sip3D".format(i)]                 = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].sip3d, op.c_float(-9999.))
                varsToKeep["mu{}_dxy".format(i)]                   = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].dxy, op.c_float(-9999.))
                varsToKeep["mu{}_dxyAbs".format(i)]                = op.switch(op.rng_len(self.muonsPreSel) >= i, op.abs(self.muonsPreSel[i-1].dxy), op.c_float(-9999.))
                varsToKeep["mu{}_dz".format(i)]                    = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].dz, op.c_float(-9999.))
                varsToKeep["mu{}_segmentCompatibility".format(i)]  = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].segmentComp, op.c_float(-9999.))
                varsToKeep["mu{}_leptonMVA".format(i)]             = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].mvaTTH, op.c_float(-9999.))
                varsToKeep["mu{}_mediumID".format(i)]              = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].mediumId, op.c_float(-9999.,"Bool_t"))
                varsToKeep["mu{}_dpt_div_pt".format(i)]            = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].tunepRelPt, op.c_float(-9999.))  # Not sure
                varsToKeep["mu{}_isfakeablesel".format(i)]         = op.switch(op.rng_len(self.muonsPreSel) >= i, op.switch(self.lambda_muonFakeSel(self.muonsPreSel[i-1]), op.c_int(1), op.c_int(0)), op.c_int(-9999))
                varsToKeep["mu{}_ismvasel".format(i)]              = op.switch(op.rng_len(self.muonsPreSel) >= i, op.switch(op.AND(self.lambda_muonTightSel(self.muonsPreSel[i-1]), self.lambda_muonFakeSel(self.muonsPreSel[i-1])), op.c_int(1), op.c_int(0)), op.c_int(-9999)) # mvasel encompasses fakeablesel
                varsToKeep["mu{}_isGenMatched".format(i)]          = op.switch(op.rng_len(self.muonsPreSel) >= i, op.switch(self.lambda_is_matched(self.muonsPreSel[i-1]), op.c_int(1), op.c_int(0)), op.c_int(-9999))
                varsToKeep["mu{}_genPartFlav".format(i)]           = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].genPartFlav, op.c_int(-9999))
                #varsToKeep["mu{}_FR".format(i)]                    = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonFR(self.muonsPreSel[i-1]), op.c_int(-9999))
                #varsToKeep["mu{}_FRCorr".format(i)]                = op.switch(op.rng_len(self.muonsPreSel) >= i, self.lambda_FF_mu(self.muonsPreSel[i-1]), op.c_int(-9999))
                varsToKeep["mu{}_looseSF".format(i)]               = op.switch(op.rng_len(self.muonsPreSel) >= i, reduce(mul,self.lambda_MuonLooseSF(self.muonsPreSel[i-1])), op.c_int(-9999))
                varsToKeep["mu{}_tightSF".format(i)]               = op.switch(op.rng_len(self.muonsPreSel) >= i, reduce(mul,self.lambda_MuonTightSF(self.muonsPreSel[i-1])), op.c_int(-9999))
            
            # Electrons #
            for i in range(1,3): # 2 leading electrons 
                varsToKeep["ele{}_pt".format(i)]                    = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].pt, op.c_float(-9999.))
                varsToKeep["ele{}_eta".format(i)]                   = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].eta, op.c_float(-9999.))
                varsToKeep["ele{}_phi".format(i)]                   = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].phi, op.c_float(-9999.))
                varsToKeep["ele{}_E".format(i)]                     = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].p4.E(), op.c_float(-9999., "float"))
                varsToKeep["ele{}_charge".format(i)]                = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].charge, op.c_int(-9999.))
                varsToKeep["ele{}_conept".format(i)]                = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electron_conept[self.electronsPreSel[i-1].idx], op.c_float(-9999.))
                varsToKeep["ele{}_miniRelIso".format(i)]            = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].miniPFRelIso_all, op.c_float(-9999.))
                varsToKeep["ele{}_PFRelIso03".format(i)]            = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].pfRelIso03_all, op.c_float(-9999.)) # Iso03, Iso04 not in NanoAOD
                varsToKeep["ele{}_jetNDauChargedMVASel".format(i)]  = op.c_float(-9999.)
                varsToKeep["ele{}_jetPtRel".format(i)]              = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].jetPtRelv2, op.c_float(-9999.))
                varsToKeep["ele{}_jetRelIso".format(i)]             = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].jetRelIso, op.c_float(-9999.))
                varsToKeep["ele{}_jetDeepJet".format(i)]            = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].jet.btagDeepFlavB, op.c_float(-9999.))
                varsToKeep["ele{}_sip3D".format(i)]                 = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].sip3d, op.c_float(-9999.))

                varsToKeep["ele{}_dxy".format(i)]                   = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].dxy, op.c_float(-9999.))
                varsToKeep["ele{}_dxyAbs".format(i)]                = op.switch(op.rng_len(self.electronsPreSel) >= i, op.abs(self.electronsPreSel[i-1].dxy), op.c_float(-9999.))
                varsToKeep["ele{}_dz".format(i)]                    = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].dz, op.c_float(-9999.))
                varsToKeep["ele{}_ntMVAeleID".format(i)]            = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].mvaFall17V2noIso, op.c_float(-9999.))
                varsToKeep["ele{}_leptonMVA".format(i)]             = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].mvaTTH, op.c_float(-9999.))
                varsToKeep["ele{}_passesConversionVeto".format(i)]  = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].convVeto, op.c_float(-9999.,"Bool_t"))
                varsToKeep["ele{}_nMissingHits".format(i)]          = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].lostHits, op.c_float(-9999.,"UChar_t"))
                varsToKeep["ele{}_sigmaEtaEta".format(i)]           = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].sieie, op.c_float(-9999.))
                varsToKeep["ele{}_HoE".format(i)]                   = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].hoe, op.c_float(-9999.))
                varsToKeep["ele{}_OoEminusOoP".format(i)]           = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].eInvMinusPInv, op.c_float(-9999.))
                varsToKeep["ele{}_isfakeablesel".format(i)]         = op.switch(op.rng_len(self.electronsPreSel) >= i, op.switch(self.lambda_electronFakeSel(self.electronsPreSel[i-1]), op.c_int(1), op.c_int(0)), op.c_int(-9999))
                varsToKeep["ele{}_ismvasel".format(i)]              = op.switch(op.rng_len(self.electronsPreSel) >= i, op.switch(op.AND(self.lambda_electronTightSel(self.electronsPreSel[i-1]), self.lambda_electronFakeSel(self.electronsPreSel[i-1])), op.c_int(1), op.c_int(0)), op.c_int(-9999)) # mvasel encompasses fakeablesel
                varsToKeep["ele{}_isGenMatched".format(i)]          = op.switch(op.rng_len(self.electronsPreSel) >= i, op.switch(self.lambda_is_matched(self.electronsPreSel[i-1]), op.c_int(1), op.c_int(0)), op.c_int(-9999))
                varsToKeep["ele{}_genPartFlav".format(i)]           = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].genPartFlav, op.c_int(-9999))
                varsToKeep["ele{}_deltaEtaSC".format(i)]            = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].deltaEtaSC, op.c_int(-9999))
                #varsToKeep["ele{}_FR".format(i)]                    = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronFR(self.electronsPreSel[i-1]), op.c_int(-9999))
                #varsToKeep["ele{}_FF".format(i)]                    = op.switch(op.rng_len(self.electronsPreSel) >= i, self.lambda_FF_el(self.electronsPreSel[i-1]), op.c_int(-9999))
                varsToKeep["ele{}_looseSF".format(i)]               = op.switch(op.rng_len(self.electronsPreSel) >= i, reduce(mul,self.lambda_ElectronLooseSF(self.electronsPreSel[i-1])), op.c_int(-9999))
                varsToKeep["ele{}_tightSF".format(i)]               = op.switch(op.rng_len(self.electronsPreSel) >= i, reduce(mul,self.lambda_ElectronTightSF(self.electronsPreSel[i-1])), op.c_int(-9999))

            # AK4 Jets #
            for i in range(1,5): # 4 leading jets 
                varsToKeep["ak4Jet{}_pt".format(i)]                 = op.switch(op.rng_len(self.ak4Jets) >= i, self.ak4Jets[i-1].pt, op.c_float(-9999.,"float"))
                varsToKeep["ak4Jet{}_eta".format(i)]                = op.switch(op.rng_len(self.ak4Jets) >= i, self.ak4Jets[i-1].eta, op.c_float(-9999.))
                varsToKeep["ak4Jet{}_phi".format(i)]                = op.switch(op.rng_len(self.ak4Jets) >= i, self.ak4Jets[i-1].phi, op.c_float(-9999.))
                varsToKeep["ak4Jet{}_E".format(i)]                  = op.switch(op.rng_len(self.ak4Jets) >= i, self.ak4Jets[i-1].p4.E(), op.c_float(-9999., "float"))
                varsToKeep["ak4Jet{}_CSV".format(i)]                = op.switch(op.rng_len(self.ak4Jets) >= i, self.ak4Jets[i-1].btagDeepFlavB, op.c_float(-9999.))
                varsToKeep["ak4Jet{}_hadronFlavour".format(i)]      = op.switch(op.rng_len(self.ak4Jets) >= i, self.ak4Jets[i-1].hadronFlavour, op.c_float(-9999.))
                varsToKeep["ak4Jet{}_btagSF".format(i)]             = op.switch(op.rng_len(self.ak4Jets) >= i, self.DeepJetDiscReshapingSF(self.ak4Jets[i-1]), op.c_float(-9999.))

            # AK8 Jets #
            for i in range(1,3): # 2 leading fatjets 
                varsToKeep["ak8Jet{}_pt".format(i)]                 = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].pt, op.c_float(-9999.))
                varsToKeep["ak8Jet{}_eta".format(i)]                = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].eta, op.c_float(-9999.))
                varsToKeep["ak8Jet{}_phi".format(i)]                = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].phi, op.c_float(-9999.))
                varsToKeep["ak8Jet{}_E".format(i)]                  = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].p4.E(), op.c_float(-9999., "float"))
                varsToKeep["ak8Jet{}_msoftdrop".format(i)]          = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].msoftdrop, op.c_float(-9999.))
                varsToKeep["ak8Jet{}_tau1".format(i)]               = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].tau1, op.c_float(-9999.))
                varsToKeep["ak8Jet{}_tau2".format(i)]               = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].tau2, op.c_float(-9999.))
                varsToKeep["ak8Jet{}_subjet0_pt".format(i)]         = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].subJet1.pt, op.c_float(-9999.))
                varsToKeep["ak8Jet{}_subjet0_eta".format(i)]        = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].subJet1.eta, op.c_float(-9999.))
                varsToKeep["ak8Jet{}_subjet0_phi".format(i)]        = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].subJet1.phi, op.c_float(-9999.))
                varsToKeep["ak8Jet{}_subjet0_CSV".format(i)]        = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].subJet1.btagDeepB, op.c_float(-9999.))
                varsToKeep["ak8Jet{}_subjet1_pt".format(i)]         = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].subJet2.pt, op.c_float(-9999.))
                varsToKeep["ak8Jet{}_subjet1_eta".format(i)]        = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].subJet2.eta, op.c_float(-9999.))
                varsToKeep["ak8Jet{}_subjet1_phi".format(i)]        = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].subJet2.phi, op.c_float(-9999.))
                varsToKeep["ak8Jet{}_subjet1_CSV".format(i)]        = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].subJet2.btagDeepB, op.c_float(-9999.))

            # MET #
             
            varsToKeep["PFMET"]    = self.corrMET.pt
            varsToKeep["PFMETphi"] = self.corrMET.phi

            # HME #

            # SF #
            electronMuon_cont = op.combine((self.electronsFakeSel, self.muonsFakeSel))
            varsToKeep["trigger_SF"] = op.multiSwitch(
                    (op.AND(op.rng_len(self.electronsTightSel)==1,op.rng_len(self.muonsTightSel)==0) , self.ttH_singleElectron_trigSF(self.electronsTightSel[0])),
                    (op.AND(op.rng_len(self.electronsTightSel)==0,op.rng_len(self.muonsTightSel)==1) , self.ttH_singleMuon_trigSF(self.muonsTightSel[0])),
                    (op.AND(op.rng_len(self.electronsTightSel)>=2,op.rng_len(self.muonsTightSel)==0) , self.lambda_ttH_doubleElectron_trigSF(self.electronsTightSel)),
                    (op.AND(op.rng_len(self.electronsTightSel)==0,op.rng_len(self.muonsTightSel)>=2) , self.lambda_ttH_doubleMuon_trigSF(self.muonsTightSel)),
                    (op.AND(op.rng_len(self.electronsTightSel)>=1,op.rng_len(self.muonsTightSel)>=1) , self.lambda_ttH_electronMuon_trigSF(electronMuon_cont[0])),
                     op.c_float(1.))

            varsToKeep["lepton_IDSF"] = op.rng_product(self.electronsTightSel, lambda el : reduce(mul,self.lambda_ElectronLooseSF(el)+self.lambda_ElectronTightSF(el))) * \
                                        op.rng_product(self.muonsTightSel, lambda mu : reduce(mul,self.lambda_MuonLooseSF(mu)+self.lambda_MuonTightSF(mu))) 

            varsToKeep["lepton_IDSF_recoToLoose"] = op.rng_product(self.electronsTightSel, lambda el : reduce(mul,self.lambda_ElectronLooseSF(el))) * \
                                                    op.rng_product(self.muonsTightSel, lambda mu : reduce(mul,self.lambda_MuonLooseSF(mu)))
            varsToKeep["lepton_IDSF_looseToTight"] = op.rng_product(self.electronsTightSel, lambda el : reduce(mul,self.lambda_ElectronTightSF(el))) * \
                                                     op.rng_product(self.muonsTightSel, lambda mu : reduce(mul,self.lambda_MuonTightSF(mu)))

            # L1 Prefire #
            if era in ["2016","2017"]:
                varsToKeep["L1prefire"] = self.L1Prefiring
            else:
                varsToKeep["L1prefire"] = op.c_float(-9999.)

            # Fake rate #
            if self.args.FakeExtrapolation:
                varsToKeep["fakeRate"] = op.multiSwitch((op.rng_len(self.ElElDileptonFakeExtrapolationSel)>=1,self.ElElFakeFactor(self.ElElDileptonFakeExtrapolationSel[0])),
                                                        (op.rng_len(self.MuMuDileptonFakeExtrapolationSel)>=1,self.MuMuFakeFactor(self.MuMuDileptonFakeExtrapolationSel[0])),
                                                        (op.rng_len(self.ElMuDileptonFakeExtrapolationSel)>=1,self.ElMuFakeFactor(self.ElMuDileptonFakeExtrapolationSel[0])),
                                                        op.c_float(0.))
            else:
                varsToKeep["fakeRate"] = op.c_float(-9999.)

            # Btagging SF #
            varsToKeep["btag_SF"] = self.btagAk4SF
            if "BtagRatioWeight" in self.__dict__.keys():
                varsToKeep["btag_reweighting"] = self.BtagRatioWeight
                varsToKeep["btag_reweighting_SF"] = self.btagAk4SF * self.BtagRatioWeight

            # ttbar PT reweighting #
            if "group" in sampleCfg and sampleCfg["group"] == 'ttbar':
                varsToKeep["topPt_wgt"] = self.ttbar_weight(self.genTop[0],self.genAntitop[0])

           # Event Weight #
            if self.is_MC:
                #varsToKeep["MC_weight"] = op.sign(t.genWeight)
                varsToKeep["MC_weight"] = t.genWeight
                puWeightsFile = os.path.join(os.path.dirname(__file__), "data", "pileup",sample+'_%s.json'%era)
                #puWeightsFile = os.path.join(os.path.dirname(__file__), "data" , "pileup", sampleCfg["pufile"])
                varsToKeep["PU_weight"] = makePileupWeight(puWeightsFile, t.Pileup_nTrueInt, nameHint=f"puweightFromFile{sample}".replace('-','_'))
                varsToKeep["eventWeight"] = noSel.weight if self.inclusive_sel else selObj.sel.weight

           
            if self.inclusive_sel:
                return noSel, varsToKeep
            else:
                return selObj.sel, varsToKeep
                

        #---------------------------------------------------------------------------------------# 
        #                                    Selection tree                                     #
        #---------------------------------------------------------------------------------------#

        #----- MET variables -----#
        MET = self.corrMET

        varsToKeep['MET_pt']  = MET.pt
        varsToKeep['MET_phi']  = MET.phi

        #----- Lepton variables -----#
        if self.args.Channel is None:
            raise RuntimeError("You need to specify --Channel")
        dilepton = None
        if self.args.Preselected:
            if self.args.Channel == "ElEl": dilepton = self.ElElDileptonPreSel[0] 
            if self.args.Channel == "MuMu": dilepton = self.MuMuDileptonPreSel[0]
            if self.args.Channel == "ElMu": dilepton = self.ElMuDileptonPreSel[0]
        if self.args.Fakeable:
            if self.args.Channel == "ElEl": dilepton = self.ElElDileptonFakeSel[0]
            if self.args.Channel == "MuMu": dilepton = self.MuMuDileptonFakeSel[0]
            if self.args.Channel == "ElMu": dilepton = self.ElMuDileptonFakeSel[0]
        if self.args.Tight:
            if self.args.Channel == "ElEl": dilepton = self.ElElDileptonTightSel[0]
            if self.args.Channel == "MuMu": dilepton = self.MuMuDileptonTightSel[0]
            if self.args.Channel == "ElMu": dilepton = self.ElMuDileptonTightSel[0]
        if self.args.FakeExtrapolation:
            if self.args.Channel == "ElEl": dilepton = self.ElElDileptonFakeExtrapolationSel[0]
            if self.args.Channel == "MuMu": dilepton = self.MuMuDileptonFakeExtrapolationSel[0]
            if self.args.Channel == "ElMu": dilepton = self.ElMuDileptonFakeExtrapolationSel[0]
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

        if self.args.Channel == "ElEl":
            varsToKeep['l1_conept'] = self.electron_conept[l1.idx]
            varsToKeep['l2_conept'] = self.electron_conept[l2.idx]
        if self.args.Channel == "MuMu":
            varsToKeep['l1_conept'] = self.muon_conept[l1.idx]
            varsToKeep['l2_conept'] = self.muon_conept[l2.idx]
        if self.args.Channel == "ElMu":
            varsToKeep['l1_conept'] = self.electron_conept[l1.idx]
            varsToKeep['l2_conept'] = self.muon_conept[l2.idx]

        varsToKeep['ll_pt'] = (l1.p4+l2.p4).Pt()
        varsToKeep['ll_DR'] = op.deltaR(l1.p4,l2.p4)
        varsToKeep['ll_DPhi'] = op.deltaPhi(l1.p4,l2.p4) # Might need abs

        varsToKeep['llmet_DPhi'] = self.HLL.DilepMET_deltaPhi(l1,l2,MET)
        varsToKeep['llmet_pt'] = self.HLL.DilepMET_Pt(l1,l2,MET)

        varsToKeep['ll_M'] = op.invariant_mass(l1.p4,l2.p4) 
        varsToKeep['ll_MT'] = self.HLL.MT_ll(l1,l2,MET)

        #----- Jet variables -----#
        if any([self.args.__dict__[item] for item in ["Ak4","Resolved0Btag","Resolved1Btag","Resolved2Btag"]]):
            if self.args.Ak4:
                j1 = self.ak4Jets[0]
                j2 = self.ak4Jets[1]
            if self.args.Resolved0Btag:
                j1 = self.ak4LightJetsByBtagScore[0]
                j2 = self.ak4LightJetsByBtagScore[0]
            if self.args.Resolved1Btag:
                j1 = self.ak4BJets[0]
                j2 = self.ak4LightJetsByBtagScore[0]
            if self.args.Resolved2Btag:
                j1 = self.ak4BJets[0]
                j2 = self.ak4BJets[1]

            varsToKeep['j1_Px']  = j1.p4.Px()
            varsToKeep['j1_Py']  = j1.p4.Py()
            varsToKeep['j1_Pz']  = j1.p4.Pz()
            varsToKeep['j1_E']   = j1.p4.E()
            varsToKeep['j1_pt']  = j1.pt
            varsToKeep['j1_eta'] = j1.eta
            varsToKeep['j1_phi'] = j1.phi
            varsToKeep['j1_btag']= j1.btagDeepFlavB

            varsToKeep['j2_Px']  = j2.p4.Px()
            varsToKeep['j2_Py']  = j2.p4.Py()
            varsToKeep['j2_Pz']  = j2.p4.Pz()
            varsToKeep['j2_E']   = j2.p4.E()
            varsToKeep['j2_pt']  = j2.pt
            varsToKeep['j2_eta'] = j2.eta
            varsToKeep['j2_phi'] = j2.phi
            varsToKeep['j2_btag']= j2.btagDeepFlavB

            varsToKeep['jj_pt'] = (j1.p4+j2.p4).Pt()
            varsToKeep['jj_DR'] = op.deltaR(j1.p4,j2.p4)
            varsToKeep['jj_DPhi'] = op.deltaPhi(j1.p4,j2.p4) # Might need abs
            varsToKeep['jj_M'] = op.invariant_mass(j1.p4,j2.p4) 

            varsToKeep['lljj_M'] = self.HLL.M_lljj(l1,l2,j1,j2)
            varsToKeep['lljj_MT'] = self.HLL.MT_lljj(l1,l2,j1,j2,MET)
            varsToKeep['lj_MinDR'] = self.HLL.MinDR_lj(l1,l2,j1,j2)
            varsToKeep['l1_jetsMinDR'] = self.HLL.JetsMinDR(l1,j1,j2)
            varsToKeep['l2_jetsMinDR'] = self.HLL.JetsMinDR(l2,j1,j2)
            varsToKeep['j1_lepsMinDR'] = self.HLL.LepsMinDR(j1,l1,l2)
            varsToKeep['j2_lepsMinDR'] = self.HLL.LepsMinDR(j2,l1,l2)
            varsToKeep['HT2'] = self.HLL.HT2(l1,l2,j1,j2,MET)
            varsToKeep['HT2R'] = self.HLL.HT2R(l1,l2,j1,j2,MET)
            #varsToKeep['jj_mRegCorr'] = op.invariant_mass(self.HLL.bJetCorrP4(j1), self.HLL.bJetCorrP4(j2)),

            varsToKeep['MET_LD'] = self.HLL.MET_LD_DL(MET,self.ak4Jets,self.electronsFakeSel,self.muonsFakeSel)

            varsToKeep['n_ak4'] = op.static_cast("UInt_t",op.rng_len(self.ak4Jets))
            varsToKeep['n_ak4_btag'] = op.static_cast("UInt_t",op.rng_len(self.ak4BJets))
#
#        #----- Fatjet variables -----#
#        if any([self.args.__dict__[item] for item in ["Ak8","Boosted"]]):
#            if self.args.Ak8:
#                fatjet = self.ak8Jets[0]
#            if self.args.Boosted:
#                fatjet = self.ak8BJets[0]
#
#            j1 = fatjet.subJet1
#            j2 = fatjet.subJet2
#
#            varsToKeep['j1_Px']  = j1.p4.Px()
#            varsToKeep['j1_Py']  = j1.p4.Py()
#            varsToKeep['j1_Pz']  = j1.p4.Pz()
#            varsToKeep['j1_E']   = j1.p4.E()
#            varsToKeep['j1_pt']  = j1.pt
#            varsToKeep['j1_eta'] = j1.eta
#            varsToKeep['j1_phi'] = j1.phi
#            varsToKeep['j1_btag']= j1.btagDeepB
#
#            varsToKeep['j2_Px']  = j2.p4.Px()
#            varsToKeep['j2_Py']  = j2.p4.Py()
#            varsToKeep['j2_Pz']  = j2.p4.Pz()
#            varsToKeep['j2_E']   = j2.p4.E()
#            varsToKeep['j2_pt']  = j2.pt
#            varsToKeep['j2_eta'] = j2.eta
#            varsToKeep['j2_phi'] = j2.phi
#            varsToKeep['j2_btag']= j2.btagDeepB
#
#            varsToKeep['fatjet_Px']  = fatjet.p4.Px()
#            varsToKeep['fatjet_Py']  = fatjet.p4.Py()
#            varsToKeep['fatjet_Pz']  = fatjet.p4.Pz()
#            varsToKeep['fatjet_E']   = fatjet.p4.E()
#            varsToKeep['fatjet_pt']  = fatjet.pt
#            varsToKeep['fatjet_eta'] = fatjet.eta
#            varsToKeep['fatjet_phi'] = fatjet.phi
#            varsToKeep['fatjet_softdropMass'] = fatjet.msoftdrop
#            varsToKeep['fatjet_btagDeepB'] = fatjet.btagDeepB
#            varsToKeep['fatjet_btagHbb'] = fatjet.btagHbb
#
#            varsToKeep['jj_pt'] = (j1.p4+j2.p4).Pt()
#            varsToKeep['jj_DR'] = op.deltaR(j1.p4,j2.p4)
#            varsToKeep['jj_DPhi'] = op.deltaPhi(j1.p4,j2.p4) # Might need abs
#            varsToKeep['jj_M'] = op.invariant_mass(j1.p4,j2.p4) 
#
#            varsToKeep['lljj_M'] = self.HLL.M_lljj(l1,l2,j1,j2)
#            varsToKeep['lljj_MT'] = self.HLL.MT_lljj(l1,l2,j1,j2,MET)
#            varsToKeep['lj_MinDR'] = self.HLL.MinDR_lj(l1,l2,j1,j2)
#            varsToKeep['l1_jetsMinDR'] = self.HLL.JetsMinDR(l1,j1,j2)
#            varsToKeep['l2_jetsMinDR'] = self.HLL.JetsMinDR(l2,j1,j2)
#            varsToKeep['j1_lepsMinDR'] = self.HLL.LepsMinDR(j1,l1,l2)
#            varsToKeep['j2_lepsMinDR'] = self.HLL.LepsMinDR(j2,l1,l2)
#            varsToKeep['HT2'] = self.HLL.HT2(l1,l2,j1,j2,MET)
#            varsToKeep['HT2R'] = self.HLL.HT2R(l1,l2,j1,j2,MET)
#
#            varsToKeep['MET_LD'] = self.HLL.MET_LD_DL(MET,self.ak4Jets,self.electronsFakeSel,self.muonsFakeSel)
#
#            varsToKeep['n_ak8'] = op.static_cast("UInt_t",op.rng_len(self.ak8Jets))
#            varsToKeep['n_ak8_btag'] = op.static_cast("UInt_t",op.rng_len(self.ak8BJets))
#
#        #----- Additional variables -----#
#        varsToKeep["MC_weight"]         = t.genWeight
#        varsToKeep['total_weight']      = selObj.sel.weight
#        varsToKeep["event"]             = None # Already in tree
#        varsToKeep["run"]               = None # Already in tree 
#        varsToKeep["ls"]                = t.luminosityBlock
#
#        # Test #
#        inputs = [0.5]*32
#        input_names = ['input_%d'%i for i in range(32)]
#        for inp, inp_name in zip(inputs,input_names):
#            varsToKeep[inp_name] = op.c_float(inp)
#
#        DNN = op.mvaEvaluator("/home/users/f/b/fbury/bamboodev/HHbbWWAnalysis/MachineLearning/ml-models/models/multi-classification/dnn/01/model/model.pb",mvaType='Tensorflow',otherArgs=(["input_1", "input_2"], "model/output/Softmax"))
#        
#        varsToKeep["output"] = DNN(*inputs)



        return selObj.sel, varsToKeep
