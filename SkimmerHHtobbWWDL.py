import os
import sys

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
class SkimmerNanoHHtobbWWSL(BaseNanoHHtobbWW,SkimmerModule):
    """ Plotter module: HH->bbW(->e/µ nu)W(->e/µ nu) histograms from NanoAOD """
    def __init__(self, args):
        super(SkimmerNanoHHtobbWWSL, self).__init__(args)

    def defineSkimSelection(self, t, noSel, sample=None, sampleCfg=None): 
        noSel = super(SkimmerNanoHHtobbWWSL,self).prepareObjects(t, noSel, sample, sampleCfg, "DL")

        era = sampleCfg['era'] 

        # Initialize varsToKeep dict #
        varsToKeep = dict()  

        #---------------------------------------------------------------------------------------# 
        #                                     Selections                                        #
        #---------------------------------------------------------------------------------------#
        if not self.basic_sync:
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
            selLeptonDict = makeLeptonSelection(self,noSel)
                # makeLeptonSelection returns dict -> value is list of three selections for 3 channels 
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
                makeAk4JetSelection(self,selObj) 
            if any([self.args.__dict__[item] for item in ["Ak8","Boosted"]]):
                makeAk8JetSelection(self,selObj) 
            if self.args.Resolved0Btag:
                makeExclusiveResolvedNoBtagSelection(self,selObj)
            if self.args.Resolved1Btag:
                makeExclusiveResolvedOneBtagSelection(self,selObj)
            if self.args.Resolved2Btag:
                makeExclusiveResolvedTwoBtagsSelection(self,selObj)
            if self.args.Boosted:
                makeInclusiveBoostedSelection(self,selObj)


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

                varsToKeep["muTight{}_pt".format(i)]               = op.switch(op.rng_len(self.muonsTightSel) >= i, self.muonsTightSel[i-1].pt, op.c_float(-9999., "float"))
                varsToKeep["muTight{}_eta".format(i)]              = op.switch(op.rng_len(self.muonsTightSel) >= i, self.muonsTightSel[i-1].eta, op.c_float(-9999.))
                varsToKeep["muTight{}_phi".format(i)]              = op.switch(op.rng_len(self.muonsTightSel) >= i, self.muonsTightSel[i-1].phi, op.c_float(-9999.))
            
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

                varsToKeep["eleTight{}_pt".format(i)]               = op.switch(op.rng_len(self.electronsTightSel) >= i, self.electronsTightSel[i-1].pt, op.c_float(-9999., "float"))
                varsToKeep["eleTight{}_eta".format(i)]              = op.switch(op.rng_len(self.electronsTightSel) >= i, self.electronsTightSel[i-1].eta, op.c_float(-9999.))
                varsToKeep["eleTight{}_phi".format(i)]              = op.switch(op.rng_len(self.electronsTightSel) >= i, self.electronsTightSel[i-1].phi, op.c_float(-9999.))
     
            # AK4 Jets #
            for i in range(1,5): # 4 leading jets 
                varsToKeep["ak4Jet{}_pt".format(i)]                 = op.switch(op.rng_len(self.ak4Jets) >= i, self.ak4Jets[i-1].pt, op.c_float(-9999.,"float"))
                varsToKeep["ak4Jet{}_eta".format(i)]                = op.switch(op.rng_len(self.ak4Jets) >= i, self.ak4Jets[i-1].eta, op.c_float(-9999.))
                varsToKeep["ak4Jet{}_phi".format(i)]                = op.switch(op.rng_len(self.ak4Jets) >= i, self.ak4Jets[i-1].phi, op.c_float(-9999.))
                varsToKeep["ak4Jet{}_E".format(i)]                  = op.switch(op.rng_len(self.ak4Jets) >= i, self.ak4Jets[i-1].p4.E(), op.c_float(-9999., "float"))
                varsToKeep["ak4Jet{}_CSV".format(i)]                = op.switch(op.rng_len(self.ak4Jets) >= i, self.ak4Jets[i-1].btagDeepFlavB, op.c_float(-9999.))
                varsToKeep["ak4Jet{}_hadronFlavour".format(i)]      = op.switch(op.rng_len(self.ak4Jets) >= i, self.ak4Jets[i-1].hadronFlavour, op.c_float(-9999.))
                varsToKeep["ak4Jet{}_btagSF".format(i)]             = op.switch(op.rng_len(self.ak4Jets) >= i, self.DeepJetDiscReshapingSF(self.ak4Jets[i-1]), op.c_float(-9999.))

            # AK4 BJets #
            for i in range(1,3): # 2 leading bjets
                varsToKeep["ak4BJet{}_pt".format(i)]                 = op.switch(op.rng_len(self.ak4BJets) >= i, self.ak4BJets[i-1].pt, op.c_float(-9999.,"float"))
                varsToKeep["ak4BJet{}_eta".format(i)]                = op.switch(op.rng_len(self.ak4BJets) >= i, self.ak4BJets[i-1].eta, op.c_float(-9999.))
                varsToKeep["ak4BJet{}_phi".format(i)]                = op.switch(op.rng_len(self.ak4BJets) >= i, self.ak4BJets[i-1].phi, op.c_float(-9999.))
                varsToKeep["ak4BJet{}_E".format(i)]                  = op.switch(op.rng_len(self.ak4BJets) >= i, self.ak4BJets[i-1].p4.E(), op.c_float(-9999., "float"))
                varsToKeep["ak4BJet{}_CSV".format(i)]                = op.switch(op.rng_len(self.ak4BJets) >= i, self.ak4BJets[i-1].btagDeepFlavB, op.c_float(-9999.))
                varsToKeep["ak4BJet{}_hadronFlavour".format(i)]      = op.switch(op.rng_len(self.ak4BJets) >= i, self.ak4BJets[i-1].hadronFlavour, op.c_float(-9999.))
                varsToKeep["ak4BJet{}_btagSF".format(i)]             = op.switch(op.rng_len(self.ak4BJets) >= i, self.DeepJetDiscReshapingSF(self.ak4BJets[i-1]), op.c_float(-9999.))

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

            # SR sync #
            for channel,dileptons in zip(["ElEl","MuMu","ElMu"],[self.ElElDileptonTightSel,self.MuMuDileptonTightSel,self.ElMuDileptonTightSel]):
                varsToKeep[channel+"_tightLep1_pt"]    = op.switch(op.rng_len(dileptons)>=1 , dileptons[0][0].pt , op.c_float(-9999.))
                varsToKeep[channel+"_tightLep1_eta"]   = op.switch(op.rng_len(dileptons)>=1 , dileptons[0][0].eta, op.c_float(-9999.))
                varsToKeep[channel+"_tightLep1_phi"]   = op.switch(op.rng_len(dileptons)>=1 , dileptons[0][0].phi , op.c_float(-9999.))
                varsToKeep[channel+"_tightLep1_pdgId"] = op.switch(op.rng_len(dileptons)>=1 , dileptons[0][0].pdgId, op.c_float(-9999.))
                varsToKeep[channel+"_tightLep2_pt"]    = op.switch(op.rng_len(dileptons)>=1 , dileptons[0][1].pt , op.c_float(-9999.))
                varsToKeep[channel+"_tightLep2_eta"]   = op.switch(op.rng_len(dileptons)>=1 , dileptons[0][1].eta, op.c_float(-9999.))
                varsToKeep[channel+"_tightLep2_phi"]   = op.switch(op.rng_len(dileptons)>=1 , dileptons[0][1].phi , op.c_float(-9999.))
                varsToKeep[channel+"_tightLep2_pdgId"] = op.switch(op.rng_len(dileptons)>=1 , dileptons[0][1].pdgId, op.c_float(-9999.))

            for channel,dileptons in zip(["ElEl","MuMu","ElMu"],[self.ElElDileptonFakeExtrapolationSel,self.MuMuDileptonFakeExtrapolationSel,self.ElMuDileptonFakeExtrapolationSel]):
                varsToKeep[channel+"_fakeCRLep1_pt"]    = op.switch(op.rng_len(dileptons)>=1 , dileptons[0][0].pt , op.c_float(-9999.))
                varsToKeep[channel+"_fakeCRLep1_eta"]   = op.switch(op.rng_len(dileptons)>=1 , dileptons[0][0].eta, op.c_float(-9999.))
                varsToKeep[channel+"_fakeCRLep1_phi"]   = op.switch(op.rng_len(dileptons)>=1 , dileptons[0][0].phi , op.c_float(-9999.))
                varsToKeep[channel+"_fakeCRLep1_pdgId"] = op.switch(op.rng_len(dileptons)>=1 , dileptons[0][0].pdgId, op.c_float(-9999.))
                varsToKeep[channel+"_fakeCRLep2_pt"]    = op.switch(op.rng_len(dileptons)>=1 , dileptons[0][1].pt , op.c_float(-9999.))
                varsToKeep[channel+"_fakeCRLep2_eta"]   = op.switch(op.rng_len(dileptons)>=1 , dileptons[0][1].eta, op.c_float(-9999.))
                varsToKeep[channel+"_fakeCRLep2_phi"]   = op.switch(op.rng_len(dileptons)>=1 , dileptons[0][1].phi , op.c_float(-9999.))
                varsToKeep[channel+"_fakeCRLep2_pdgId"] = op.switch(op.rng_len(dileptons)>=1 , dileptons[0][1].pdgId, op.c_float(-9999.))


            varsToKeep['len_ElElDileptonFakeSel'] = op.static_cast("UInt_t",op.rng_len(self.ElElDileptonFakeSel))
            varsToKeep['len_MuMuDileptonFakeSel'] = op.static_cast("UInt_t",op.rng_len(self.MuMuDileptonFakeSel))
            varsToKeep['len_ElMuDileptonFakeSel'] = op.static_cast("UInt_t",op.rng_len(self.ElMuDileptonFakeSel))

            varsToKeep['len_lead_electrons'] = op.static_cast("UInt_t",op.rng_len(self.leadElectronsFakeSel))
            varsToKeep['len_lead_muons'] = op.static_cast("UInt_t",op.rng_len(self.leadMuonsFakeSel))
            for i in range(1,3): 
                varsToKeep["lead_ele{}_pt".format(i)]      = op.switch(op.rng_len(self.leadElectronsFakeSel) >= i, self.leadElectronsFakeSel[i-1].pt, op.c_float(-9999.))
                varsToKeep["lead_ele{}_eta".format(i)]     = op.switch(op.rng_len(self.leadElectronsFakeSel) >= i, self.leadElectronsFakeSel[i-1].eta, op.c_float(-9999.))
                varsToKeep["lead_ele{}_phi".format(i)]     = op.switch(op.rng_len(self.leadElectronsFakeSel) >= i, self.leadElectronsFakeSel[i-1].phi, op.c_float(-9999.))
                varsToKeep["lead_ele{}_pdgId".format(i)]   = op.switch(op.rng_len(self.leadElectronsFakeSel) >= i, self.leadElectronsFakeSel[i-1].pdgId, op.c_float(-9999.))
            for i in range(1,3): 
                varsToKeep["lead_mu{}_pt".format(i)]      = op.switch(op.rng_len(self.leadMuonsFakeSel) >= i, self.leadMuonsFakeSel[i-1].pt, op.c_float(-9999.))
                varsToKeep["lead_mu{}_eta".format(i)]     = op.switch(op.rng_len(self.leadMuonsFakeSel) >= i, self.leadMuonsFakeSel[i-1].eta, op.c_float(-9999.))
                varsToKeep["lead_mu{}_phi".format(i)]     = op.switch(op.rng_len(self.leadMuonsFakeSel) >= i, self.leadMuonsFakeSel[i-1].phi, op.c_float(-9999.))
                varsToKeep["lead_mu{}_pdgId".format(i)]   = op.switch(op.rng_len(self.leadMuonsFakeSel) >= i, self.leadMuonsFakeSel[i-1].pdgId, op.c_float(-9999.))


#            lambda_lowMllCut    = lambda dileptons: op.NOT(op.rng_any(dileptons, lambda dilep : op.invariant_mass(dilep[0].p4, dilep[1].p4)<12.))
#            lambda_outZ         = lambda dileptons: op.NOT(op.rng_any(dileptons, lambda dilep : op.in_range(80.,op.invariant_mass(dilep[0].p4, dilep[1].p4),100.)))
#
#            varsToKeep['ElEl_mll12Cut'] = lambda_lowMllCut(self.ElElDileptonPreSel)
#            varsToKeep['MuMu_mll12Cut'] = lambda_lowMllCut(self.MuMuDileptonPreSel)
#            varsToKeep['ElMu_mll12Cut'] = lambda_lowMllCut(self.ElMuDileptonPreSel)
#            varsToKeep['ElEl_mllZVeto'] = lambda_outZ(self.OSElElDileptonPreSel)
#            varsToKeep['MuMu_mllZVeto'] = lambda_outZ(self.OSMuMuDileptonPreSel)
#            varsToKeep['ElMu_mllZVeto'] = lambda_outZ(self.OSElMuDileptonPreSel)
#
#            varsToKeep['ElEl_rng_len_preselected'] = op.static_cast("UInt_t",op.rng_len(self.ElElDileptonPreSel))
#            varsToKeep['ElEl_rng_len_preselected_OS'] = op.static_cast("UInt_t",op.rng_len(self.OSElElDileptonPreSel))
#            varsToKeep['MuMu_rng_len_preselected'] = op.static_cast("UInt_t",op.rng_len(self.MuMuDileptonPreSel))
#            varsToKeep['MuMu_rng_len_preselected_OS'] = op.static_cast("UInt_t",op.rng_len(self.OSMuMuDileptonPreSel))
#            varsToKeep['ElMu_rng_len_preselected'] = op.static_cast("UInt_t",op.rng_len(self.ElMuDileptonPreSel))
#            varsToKeep['ElMu_rng_len_preselected_OS'] = op.static_cast("UInt_t",op.rng_len(self.OSElMuDileptonPreSel))
#
#            for i in range(1,5): 
#                varsToKeep["ElElPreSelPair{}_lep1_pt".format(i)]         = op.switch(op.rng_len(self.ElElDileptonPreSel) >= i, self.ElElDileptonPreSel[i-1][0].pt, op.c_float(-9999.))
#                varsToKeep["ElElPreSelPair{}_lep1_eta".format(i)]        = op.switch(op.rng_len(self.ElElDileptonPreSel) >= i, self.ElElDileptonPreSel[i-1][0].eta, op.c_float(-9999.))
#                varsToKeep["ElElPreSelPair{}_lep1_phi".format(i)]        = op.switch(op.rng_len(self.ElElDileptonPreSel) >= i, self.ElElDileptonPreSel[i-1][0].phi, op.c_float(-9999.))
#                varsToKeep["ElElPreSelPair{}_lep1_pdgId".format(i)]      = op.switch(op.rng_len(self.ElElDileptonPreSel) >= i, self.ElElDileptonPreSel[i-1][0].pdgId, op.c_float(-9999.))
#                varsToKeep["ElElPreSelPair{}_lep2_pt".format(i)]         = op.switch(op.rng_len(self.ElElDileptonPreSel) >= i, self.ElElDileptonPreSel[i-1][1].pt, op.c_float(-9999.))
#                varsToKeep["ElElPreSelPair{}_lep2_eta".format(i)]        = op.switch(op.rng_len(self.ElElDileptonPreSel) >= i, self.ElElDileptonPreSel[i-1][1].eta, op.c_float(-9999.))
#                varsToKeep["ElElPreSelPair{}_lep2_phi".format(i)]        = op.switch(op.rng_len(self.ElElDileptonPreSel) >= i, self.ElElDileptonPreSel[i-1][1].phi, op.c_float(-9999.))
#                varsToKeep["ElElPreSelPair{}_lep2_pdgId".format(i)]      = op.switch(op.rng_len(self.ElElDileptonPreSel) >= i, self.ElElDileptonPreSel[i-1][1].pdgId, op.c_float(-9999.))
#                varsToKeep["ElElPreSelPair{}_invMass".format(i)]         = op.switch(op.rng_len(self.ElElDileptonPreSel) >= i, op.invariant_mass(self.ElElDileptonPreSel[i-1][0].p4,self.ElElDileptonPreSel[i-1][1].p4), op.c_float(-9999.))
#                varsToKeep["ElElOSPreSelPair{}_lep1_pt".format(i)]         = op.switch(op.rng_len(self.OSElElDileptonPreSel) >= i, self.OSElElDileptonPreSel[i-1][0].pt, op.c_float(-9999.))
#                varsToKeep["ElElOSPreSelPair{}_lep1_eta".format(i)]        = op.switch(op.rng_len(self.OSElElDileptonPreSel) >= i, self.OSElElDileptonPreSel[i-1][0].eta, op.c_float(-9999.))
#                varsToKeep["ElElOSPreSelPair{}_lep1_phi".format(i)]        = op.switch(op.rng_len(self.OSElElDileptonPreSel) >= i, self.OSElElDileptonPreSel[i-1][0].phi, op.c_float(-9999.))
#                varsToKeep["ElElOSPreSelPair{}_lep1_pdgId".format(i)]      = op.switch(op.rng_len(self.OSElElDileptonPreSel) >= i, self.OSElElDileptonPreSel[i-1][0].pdgId, op.c_float(-9999.))
#                varsToKeep["ElElOSPreSelPair{}_lep2_pt".format(i)]         = op.switch(op.rng_len(self.OSElElDileptonPreSel) >= i, self.OSElElDileptonPreSel[i-1][1].pt, op.c_float(-9999.))
#                varsToKeep["ElElOSPreSelPair{}_lep2_eta".format(i)]        = op.switch(op.rng_len(self.OSElElDileptonPreSel) >= i, self.OSElElDileptonPreSel[i-1][1].eta, op.c_float(-9999.))
#                varsToKeep["ElElOSPreSelPair{}_lep2_phi".format(i)]        = op.switch(op.rng_len(self.OSElElDileptonPreSel) >= i, self.OSElElDileptonPreSel[i-1][1].phi, op.c_float(-9999.))
#                varsToKeep["ElElOSPreSelPair{}_lep2_pdgId".format(i)]      = op.switch(op.rng_len(self.OSElElDileptonPreSel) >= i, self.OSElElDileptonPreSel[i-1][1].pdgId, op.c_float(-9999.))
#                varsToKeep["ElElOSPreSelPair{}_invMass".format(i)]         = op.switch(op.rng_len(self.OSElElDileptonPreSel) >= i, op.invariant_mass(self.OSElElDileptonPreSel[i-1][0].p4,self.OSElElDileptonPreSel[i-1][1].p4), op.c_float(-9999.))
#            for i in range(1,5): 
#                varsToKeep["MuMuPreSelPair{}_lep1_pt".format(i)]         = op.switch(op.rng_len(self.MuMuDileptonPreSel) >= i, self.MuMuDileptonPreSel[i-1][0].pt, op.c_float(-9999.))
#                varsToKeep["MuMuPreSelPair{}_lep1_eta".format(i)]        = op.switch(op.rng_len(self.MuMuDileptonPreSel) >= i, self.MuMuDileptonPreSel[i-1][0].eta, op.c_float(-9999.))
#                varsToKeep["MuMuPreSelPair{}_lep1_phi".format(i)]        = op.switch(op.rng_len(self.MuMuDileptonPreSel) >= i, self.MuMuDileptonPreSel[i-1][0].phi, op.c_float(-9999.))
#                varsToKeep["MuMuPreSelPair{}_lep1_pdgId".format(i)]      = op.switch(op.rng_len(self.MuMuDileptonPreSel) >= i, self.MuMuDileptonPreSel[i-1][0].pdgId, op.c_float(-9999.))
#                varsToKeep["MuMuPreSelPair{}_lep2_pt".format(i)]         = op.switch(op.rng_len(self.MuMuDileptonPreSel) >= i, self.MuMuDileptonPreSel[i-1][1].pt, op.c_float(-9999.))
#                varsToKeep["MuMuPreSelPair{}_lep2_eta".format(i)]        = op.switch(op.rng_len(self.MuMuDileptonPreSel) >= i, self.MuMuDileptonPreSel[i-1][1].eta, op.c_float(-9999.))
#                varsToKeep["MuMuPreSelPair{}_lep2_phi".format(i)]        = op.switch(op.rng_len(self.MuMuDileptonPreSel) >= i, self.MuMuDileptonPreSel[i-1][1].phi, op.c_float(-9999.))
#                varsToKeep["MuMuPreSelPair{}_lep2_pdgId".format(i)]      = op.switch(op.rng_len(self.MuMuDileptonPreSel) >= i, self.MuMuDileptonPreSel[i-1][1].pdgId, op.c_float(-9999.))
#                varsToKeep["MuMuPreSelPair{}_invMass".format(i)]         = op.switch(op.rng_len(self.MuMuDileptonPreSel) >= i, op.invariant_mass(self.MuMuDileptonPreSel[i-1][0].p4,self.MuMuDileptonPreSel[i-1][1].p4), op.c_float(-9999.))
#                varsToKeep["MuMuOSPreSelPair{}_lep1_pt".format(i)]         = op.switch(op.rng_len(self.OSMuMuDileptonPreSel) >= i, self.OSMuMuDileptonPreSel[i-1][0].pt, op.c_float(-9999.))
#                varsToKeep["MuMuOSPreSelPair{}_lep1_eta".format(i)]        = op.switch(op.rng_len(self.OSMuMuDileptonPreSel) >= i, self.OSMuMuDileptonPreSel[i-1][0].eta, op.c_float(-9999.))
#                varsToKeep["MuMuOSPreSelPair{}_lep1_phi".format(i)]        = op.switch(op.rng_len(self.OSMuMuDileptonPreSel) >= i, self.OSMuMuDileptonPreSel[i-1][0].phi, op.c_float(-9999.))
#                varsToKeep["MuMuOSPreSelPair{}_lep1_pdgId".format(i)]      = op.switch(op.rng_len(self.OSMuMuDileptonPreSel) >= i, self.OSMuMuDileptonPreSel[i-1][0].pdgId, op.c_float(-9999.))
#                varsToKeep["MuMuOSPreSelPair{}_lep2_pt".format(i)]         = op.switch(op.rng_len(self.OSMuMuDileptonPreSel) >= i, self.OSMuMuDileptonPreSel[i-1][1].pt, op.c_float(-9999.))
#                varsToKeep["MuMuOSPreSelPair{}_lep2_eta".format(i)]        = op.switch(op.rng_len(self.OSMuMuDileptonPreSel) >= i, self.OSMuMuDileptonPreSel[i-1][1].eta, op.c_float(-9999.))
#                varsToKeep["MuMuOSPreSelPair{}_lep2_phi".format(i)]        = op.switch(op.rng_len(self.OSMuMuDileptonPreSel) >= i, self.OSMuMuDileptonPreSel[i-1][1].phi, op.c_float(-9999.))
#                varsToKeep["MuMuOSPreSelPair{}_lep2_pdgId".format(i)]      = op.switch(op.rng_len(self.OSMuMuDileptonPreSel) >= i, self.OSMuMuDileptonPreSel[i-1][1].pdgId, op.c_float(-9999.))
#                varsToKeep["MuMuOSPreSelPair{}_invMass".format(i)]         = op.switch(op.rng_len(self.OSMuMuDileptonPreSel) >= i, op.invariant_mass(self.OSMuMuDileptonPreSel[i-1][0].p4,self.OSMuMuDileptonPreSel[i-1][1].p4), op.c_float(-9999.))
#            for i in range(1,5): 
#                varsToKeep["ElMuPreSelPair{}_lep1_pt".format(i)]         = op.switch(op.rng_len(self.ElMuDileptonPreSel) >= i, self.ElMuDileptonPreSel[i-1][0].pt, op.c_float(-9999.))
#                varsToKeep["ElMuPreSelPair{}_lep1_eta".format(i)]        = op.switch(op.rng_len(self.ElMuDileptonPreSel) >= i, self.ElMuDileptonPreSel[i-1][0].eta, op.c_float(-9999.))
#                varsToKeep["ElMuPreSelPair{}_lep1_phi".format(i)]        = op.switch(op.rng_len(self.ElMuDileptonPreSel) >= i, self.ElMuDileptonPreSel[i-1][0].phi, op.c_float(-9999.))
#                varsToKeep["ElMuPreSelPair{}_lep1_pdgId".format(i)]      = op.switch(op.rng_len(self.ElMuDileptonPreSel) >= i, self.ElMuDileptonPreSel[i-1][0].pdgId, op.c_float(-9999.))
#                varsToKeep["ElMuPreSelPair{}_lep2_pt".format(i)]         = op.switch(op.rng_len(self.ElMuDileptonPreSel) >= i, self.ElMuDileptonPreSel[i-1][1].pt, op.c_float(-9999.))
#                varsToKeep["ElMuPreSelPair{}_lep2_eta".format(i)]        = op.switch(op.rng_len(self.ElMuDileptonPreSel) >= i, self.ElMuDileptonPreSel[i-1][1].eta, op.c_float(-9999.))
#                varsToKeep["ElMuPreSelPair{}_lep2_phi".format(i)]        = op.switch(op.rng_len(self.ElMuDileptonPreSel) >= i, self.ElMuDileptonPreSel[i-1][1].phi, op.c_float(-9999.))
#                varsToKeep["ElMuPreSelPair{}_lep2_pdgId".format(i)]      = op.switch(op.rng_len(self.ElMuDileptonPreSel) >= i, self.ElMuDileptonPreSel[i-1][1].pdgId, op.c_float(-9999.))
#                varsToKeep["ElMuPreSelPair{}_invMass".format(i)]         = op.switch(op.rng_len(self.ElMuDileptonPreSel) >= i, op.invariant_mass(self.ElMuDileptonPreSel[i-1][0].p4,self.ElMuDileptonPreSel[i-1][1].p4), op.c_float(-9999.))
#                varsToKeep["ElMuOSPreSelPair{}_lep1_pt".format(i)]         = op.switch(op.rng_len(self.OSElMuDileptonPreSel) >= i, self.OSElMuDileptonPreSel[i-1][0].pt, op.c_float(-9999.))
#                varsToKeep["ElMuOSPreSelPair{}_lep1_eta".format(i)]        = op.switch(op.rng_len(self.OSElMuDileptonPreSel) >= i, self.OSElMuDileptonPreSel[i-1][0].eta, op.c_float(-9999.))
#                varsToKeep["ElMuOSPreSelPair{}_lep1_phi".format(i)]        = op.switch(op.rng_len(self.OSElMuDileptonPreSel) >= i, self.OSElMuDileptonPreSel[i-1][0].phi, op.c_float(-9999.))
#                varsToKeep["ElMuOSPreSelPair{}_lep1_pdgId".format(i)]      = op.switch(op.rng_len(self.OSElMuDileptonPreSel) >= i, self.OSElMuDileptonPreSel[i-1][0].pdgId, op.c_float(-9999.))
#                varsToKeep["ElMuOSPreSelPair{}_lep2_pt".format(i)]         = op.switch(op.rng_len(self.OSElMuDileptonPreSel) >= i, self.OSElMuDileptonPreSel[i-1][1].pt, op.c_float(-9999.))
#                varsToKeep["ElMuOSPreSelPair{}_lep2_eta".format(i)]        = op.switch(op.rng_len(self.OSElMuDileptonPreSel) >= i, self.OSElMuDileptonPreSel[i-1][1].eta, op.c_float(-9999.))
#                varsToKeep["ElMuOSPreSelPair{}_lep2_phi".format(i)]        = op.switch(op.rng_len(self.OSElMuDileptonPreSel) >= i, self.OSElMuDileptonPreSel[i-1][1].phi, op.c_float(-9999.))
#                varsToKeep["ElMuOSPreSelPair{}_lep2_pdgId".format(i)]      = op.switch(op.rng_len(self.OSElMuDileptonPreSel) >= i, self.OSElMuDileptonPreSel[i-1][1].pdgId, op.c_float(-9999.))
#                varsToKeep["ElMuOSPreSelPair{}_invMass".format(i)]         = op.switch(op.rng_len(self.OSElMuDileptonPreSel) >= i, op.invariant_mass(self.OSElMuDileptonPreSel[i-1][0].p4,self.OSElMuDileptonPreSel[i-1][1].p4), op.c_float(-9999.))





            # MET #
             
            varsToKeep["PFMET"]    = self.corrMET.pt
            varsToKeep["PFMETphi"] = self.corrMET.phi

            # HME #

            # SF #
            from operator import mul
            from functools import reduce

#            if era == '2016':
#                electronMuon_cont = op.combine((self.electronsFakeSel, self.muonsFakeSel), pred=self.lambda_leptonOS)
#                varsToKeep["trigger_SF"] = op.multiSwitch(
#                        (op.AND(op.rng_len(self.electronsTightSel)==1,op.rng_len(self.muonsTightSel)==0) , self.ttH_singleElectron_trigSF(self.electronsTightSel[0])),
#                        (op.AND(op.rng_len(self.electronsTightSel)==0,op.rng_len(self.muonsTightSel)==1) , self.ttH_singleMuon_trigSF(self.muonsTightSel[0])),
#                        (op.AND(op.rng_len(self.electronsTightSel)>=2,op.rng_len(self.muonsTightSel)==0) , self.lambda_ttH_doubleElectron_trigSF(self.electronsTightSel)),
#                        (op.AND(op.rng_len(self.electronsTightSel)==0,op.rng_len(self.muonsTightSel)>=2) , self.lambda_ttH_doubleMuon_trigSF(self.muonsTightSel)),
#                        (op.AND(op.rng_len(self.electronsTightSel)>=1,op.rng_len(self.muonsTightSel)>=1) , self.lambda_ttH_electronMuon_trigSF(electronMuon_cont)),
#                         op.c_float(1.))
#            else:
#                raise NotImplementedError

            if era == '2016' or era == '2017':
                varsToKeep["lepton_IDSF"] = op.rng_product(self.electronsTightSel, lambda el : reduce(mul,self.lambda_ElectronLooseSF(el)+self.lambda_ElectronTightSF(el))) * \
                                            op.rng_product(self.muonsTightSel, lambda mu : reduce(mul,self.lambda_MuonLooseSF(mu)+self.lambda_MuonTightSF(mu))) 
            else:
                raise NotImplementedError

            # Btagging SF #
            varsToKeep["btag_SF"] = self.btagSF
            varsToKeep["btag_reweighting"] = self.BtagRatioWeight
            varsToKeep["btag_reweighting_SF"] = self.btagSF * self.BtagRatioWeight

            # ttbar PT reweighting #
            if "group" in sampleCfg and sampleCfg["group"] == 'ttbar':
                varsToKeep["topPt_wgt"] = self.ttbar_weight(self.genTop[0],self.genAntitop[0])

            # Event Weight #
            if self.is_MC:
                #varsToKeep["MC_weight"] = op.sign(t.genWeight)
                varsToKeep["MC_weight"] = t.genWeight
                puWeightsFile = os.path.join(os.path.dirname(__file__), "data" , "pileup", sampleCfg["pufile"])
                varsToKeep["PU_weight"] = makePileupWeight(puWeightsFile, t.Pileup_nTrueInt, nameHint=f"puweightFromFile{sample}".replace('-','_'))
                varsToKeep["eventWeight"] = noSel.weight if self.basic_sync else selObj.sel.weight


            # Triggers #
            varsToKeep['HLT_IsoMu24'] = None
            varsToKeep['HLT_Ele27_WPTight_Gsf'] = None
            varsToKeep['HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL'] = None
            varsToKeep['HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ'] = None
            varsToKeep['HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ'] = None
            varsToKeep['HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ'] = None
            varsToKeep['HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ'] = None
            varsToKeep['HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL'] = None
            varsToKeep['HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL'] = None

            if self.basic_sync:
                return noSel, varsToKeep
            else:
                return selObj.sel, varsToKeep
                

        #---------------------------------------------------------------------------------------# 
        #                                    Selection tree                                     #
        #---------------------------------------------------------------------------------------#

        #----- MET variables -----#
        MET = self.corrMET

        #varsToKeep['MET_pt']  = MET.pt
        #varsToKeep['MET_phi']  = MET.phi

        #----- Lepton variables -----#
        if self.args.Channel is None:
            raise RuntimeError("You need to specify --Channel")
        dilepton = None
        if self.args.Preselected:
            if self.args.Channel == "ElEl": dilepton = self.OSElElDileptonPreSel[0] 
            if self.args.Channel == "MuMu": dilepton = self.OSMuMuDileptonPreSel[0]
            if self.args.Channel == "ElMu": dilepton = self.OSElMuDileptonPreSel[0]
        if self.args.Fakeable:
            if self.args.Channel == "ElEl": dilepton = self.OSElElDileptonFakeSel[0]
            if self.args.Channel == "MuMu": dilepton = self.OSMuMuDileptonFakeSel[0]
            if self.args.Channel == "ElMu": dilepton = self.OSElMuDileptonFakeSel[0]
        if self.args.Tight:
            if self.args.Channel == "ElEl": dilepton = self.OSElElDileptonTightSel[0]
            if self.args.Channel == "MuMu": dilepton = self.OSMuMuDileptonTightSel[0]
            if self.args.Channel == "ElMu": dilepton = self.OSElMuDileptonTightSel[0]
        if self.args.FakeExtrapolation:
            if self.args.Channel == "ElEl": dilepton = self.OSElElDileptonFakeExtrapolationSel[0]
            if self.args.Channel == "MuMu": dilepton = self.OSMuMuDileptonFakeExtrapolationSel[0]
            if self.args.Channel == "ElMu": dilepton = self.OSElMuDileptonFakeExtrapolationSel[0]

        varsToKeep['l1_pt']  = dilepton[0].pt
        varsToKeep['l1_eta'] = dilepton[0].eta
        varsToKeep['l1_phi'] = dilepton[0].phi

        varsToKeep['l2_pt']  = dilepton[1].pt
        varsToKeep['l2_eta'] = dilepton[1].eta
        varsToKeep['l2_phi'] = dilepton[1].phi

        varsToKeep['PT_ll'] = (dilepton[0].p4+dilepton[1].p4).Pt()
        varsToKeep['DR_ll'] = op.deltaR(dilepton[0].p4,dilepton[1].p4)
        varsToKeep['DPhi_ll'] = op.deltaPhi(dilepton[0].p4,dilepton[1].p4) # Might need abs

        varsToKeep['DPhi_llmet'] = DilepMET_deltaPhi(dilepton[0],dilepton[1],MET)
        varsToKeep['PT_llmet'] = DilepMET_Pt(dilepton[0],dilepton[1],MET)

        varsToKeep['M_ll'] = op.invariant_mass(dilepton[0].p4,dilepton[1].p4) 
        varsToKeep['MT_ll'] = MT_ll(dilepton[0],dilepton[1],MET)

        #----- Jet variables -----#
        if any([self.args.__dict__[item] for item in ["Ak4","Resolved0Btag","Resolved1Btag","Resolved2Btag"]]):
            if self.args.Ak4:
                leadjet = self.ak4Jets[0]
                subleadjet = self.ak4Jets[1]
            if self.args.Resolved0Btag:
                leadjet = self.ak4LightJets[0]
                subleadjet = self.ak4LightJets[1]
            if self.args.Resolved1Btag:
                leadjet = self.ak4BJets[0]
                subleadjet = self.ak4LightJets[1]
            if self.args.Resolved2Btag:
                leadjet = self.ak4BJets[0]
                subleadjet = self.ak4BJets[1]

            varsToKeep['j1_pt']  = leadjet.pt
            varsToKeep['j1_eta'] = leadjet.eta
            varsToKeep['j1_phi'] = leadjet.phi

            varsToKeep['j2_pt']  = subleadjet.pt
            varsToKeep['j2_eta'] = subleadjet.eta
            varsToKeep['j2_phi'] = subleadjet.phi

            varsToKeep['PT_jj'] = (leadjet.p4+subleadjet.p4).Pt()
            varsToKeep['DR_jj'] = op.deltaR(leadjet.p4,subleadjet.p4)
            varsToKeep['DPhi_jj'] = op.deltaPhi(leadjet.p4,subleadjet.p4) # Might need abs
            varsToKeep['M_jj'] = op.invariant_mass(leadjet.p4,subleadjet.p4) 

            varsToKeep['M_lljj'] = M_lljj(dilepton[0],dilepton[1],leadjet,subleadjet)
            varsToKeep['MT_lljj'] = MT_lljj(dilepton[0],dilepton[1],leadjet,subleadjet,MET)
            varsToKeep['MinDR_lj'] = MinDR_lj(dilepton[0],dilepton[1],leadjet,subleadjet)
            varsToKeep['HT2'] = HT2(dilepton[0],dilepton[1],leadjet,subleadjet,MET)
            varsToKeep['HT2R'] = HT2R(dilepton[0],dilepton[1],leadjet,subleadjet,MET)

        #----- Fatjet variables -----#
        if any([self.args.__dict__[item] for item in ["Ak8","Boosted"]]):
            if self.args.Ak8:
                fatjet = self.ak8Jets[0]
            if self.args.Boosted:
                fatjet = self.ak8BJets[0]

            varsToKeep['fatjet_pt']  = fatjet.pt
            varsToKeep['fatjet_eta'] = fatjet.eta
            varsToKeep['fatjet_phi'] = fatjet.phi
            varsToKeep['fatjet_softdropMass'] = fatjet.msoftdrop

            varsToKeep['M_lljj'] = op.invariant_mass(dilepton[0],dilepton[1],fatjet.subJet1,fatjet.subJet2)
            varsToKeep['MT_lljj'] = MT_lljj(dilepton[0],dilepton[1],fatjet.subJet1,fatjet.subJet2,MET)
            varsToKeep['MinDR_lj'] = MinDR_lj(dilepton[0],dilepton[1],fatjet.subJet1,fatjet.subJet2)
            varsToKeep['HT2'] = HT2(dilepton[0],dilepton[1],fatjet.subJet1,fatjet.subJet2,MET)
            varsToKeep['HT2R'] = HT2R(dilepton[0],dilepton[1],fatjet.subJet1,fatjet.subJet2,MET)

        #----- Additional variables -----#
        varsToKeep['total_weight'] = selObj.sel.weight

        #return leptonSel.sel, varsToKeep
        return selObj.sel, varsToKeep
