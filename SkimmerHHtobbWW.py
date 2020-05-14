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
class SkimmerNanoHHtobbWW(BaseNanoHHtobbWW,SkimmerModule):
    """ Plotter module: HH->bbW(->e/µ nu)W(->e/µ nu) histograms from NanoAOD """
    def __init__(self, args):
        super(SkimmerNanoHHtobbWW, self).__init__(args)

    def defineSkimSelection(self, t, noSel, sample=None, sampleCfg=None): 
        noSel = self.prepareObjects(t, noSel, sample, sampleCfg)

        era = sampleCfg['era'] 

        # Initialize varsToKeep dict #
        varsToKeep = dict()  

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

            # Muons #
            for i in range(1,3): # 2 leading muons
                varsToKeep["mu{}_pt".format(i)]                    = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].pt, op.c_float(-9999., "float"))
                varsToKeep["mu{}_eta".format(i)]                   = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].eta, op.c_float(-9999.))
                varsToKeep["mu{}_phi".format(i)]                   = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].phi, op.c_float(-9999.))
                varsToKeep["mu{}_E".format(i)]                     = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].p4.E(), op.c_float(-9999., "float"))
                varsToKeep["mu{}_charge".format(i)]                = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].charge, op.c_int(-9999.))
                varsToKeep["mu{}_conept".format(i)]                = op.switch(op.rng_len(self.muonsPreSel) >= i, self.lambda_conept_muon(self.muonsPreSel[i-1]), op.c_float(-9999.))
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
                varsToKeep["mu{}_genPartFlav".format(i)]          = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].genPartFlav, op.c_int(-9999))
            
            # Electrons #
            for i in range(1,3): # 2 leading electrons 
                varsToKeep["ele{}_pt".format(i)]                    = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].pt, op.c_float(-9999.))
                varsToKeep["ele{}_eta".format(i)]                   = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].eta, op.c_float(-9999.))
                varsToKeep["ele{}_phi".format(i)]                   = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].phi, op.c_float(-9999.))
                varsToKeep["ele{}_E".format(i)]                     = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].p4.E(), op.c_float(-9999., "float"))
                varsToKeep["ele{}_charge".format(i)]                = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].charge, op.c_int(-9999.))
                varsToKeep["ele{}_conept".format(i)]                = op.switch(op.rng_len(self.electronsPreSel) >= i, self.lambda_conept_electron(self.electronsPreSel[i-1]), op.c_float(-9999.))
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
                varsToKeep["ele{}_genPartFlav".format(i)]          = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].genPartFlav, op.c_int(-9999))
                varsToKeep["ele{}_deltaEtaSC".format(i)]          = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].deltaEtaSC, op.c_int(-9999))
     
            # AK4 Jets #
            for i in range(1,5): # 4 leading jets 
                varsToKeep["ak4Jet{}_pt".format(i)]                 = op.switch(op.rng_len(self.ak4Jets) >= i, self.ak4Jets[i-1].pt, op.c_float(-9999.,"float"))
                varsToKeep["ak4Jet{}_eta".format(i)]                = op.switch(op.rng_len(self.ak4Jets) >= i, self.ak4Jets[i-1].eta, op.c_float(-9999.))
                varsToKeep["ak4Jet{}_phi".format(i)]                = op.switch(op.rng_len(self.ak4Jets) >= i, self.ak4Jets[i-1].phi, op.c_float(-9999.))
                varsToKeep["ak4Jet{}_E".format(i)]                  = op.switch(op.rng_len(self.ak4Jets) >= i, self.ak4Jets[i-1].p4.E(), op.c_float(-9999., "float"))
                varsToKeep["ak4Jet{}_CSV".format(i)]                = op.switch(op.rng_len(self.ak4Jets) >= i, self.ak4Jets[i-1].btagDeepFlavB, op.c_float(-9999.))

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

            # Event Weight #
            if self.is_MC:
                varsToKeep["MC_weight"] = t.genWeight
                puWeightsFile = os.path.join(os.path.dirname(__file__), "data" , "pileup", sampleCfg["pufile"])
                varsToKeep["PU_weight"] = makePileupWeight(puWeightsFile, t.Pileup_nTrueInt, nameHint=f"puweightFromFile{sample}".replace('-','_'))

            return noSel, varsToKeep

        #---------------------------------------------------------------------------------------# 
        #                                    Selection tree                                     #
        #---------------------------------------------------------------------------------------#
        #----- Check arguments -----#
        lepton_level = ["Preselected","Fakeable","Tight","FakeExtrapolation"]               # Only one must be in args
        jet_level = ["Ak4","Ak8","Resolved0Btag","Resolved1Btag","Resolved2Btag","Boosted"] # Only one must be in args

        if [boolean for (level,boolean) in self.args.__dict__.items() if level in lepton_level].count(True) != 1:
            raise RuntimeError("Only one of the lepton arguments must be used, check --help")
        if [boolean for (level,boolean) in self.args.__dict__.items() if level in jet_level].count(True) != 1:
            raise RuntimeError("Only one of the jet arguments must be used, check --help")
        if self.args.Channel not in ["ElEl","MuMu","ElMu"]:
            raise RuntimeError("Channel mus be either 'ElEl', 'MuMu' or 'ElMu'")

        #----- Lepton selection -----#
        # Args are passed within the self #
        selLeptonDict = makeLeptonSelection(self,noSel)
            # makeLeptonSelection returns dict -> value is list of three selections for 3 channels 
            # [0] -> we take the first and only key and value because restricted to one leptons selection
        selLeptonList = list(selLeptonDict.values())[0]
        if self.args.Channel == "ElEl":
            tupleSel = selLeptonList[0] # First item of list is ElEl selection
        if self.args.Channel == "MuMu":
            tupleSel = selLeptonList[1] # Second item of list is MuMu selection
        if self.args.Channel == "ElMu":
            tupleSel = selLeptonList[2] # Third item of list is ElMu selection

        #----- Jet selection -----#
        # Since the selections in one line, we can use the non copy option of the selection to modify the selection object internally
        if any([self.args.__dict__[item] for item in ["Ak4","Resolved0Btag","Resolved1Btag","Resolved2Btag"]]):
            makeAk4JetSelection(self,tupleSel) 
        if any([self.args.__dict__[item] for item in ["Ak8","Boosted"]]):
            makeAk8JetSelection(self,tupleSel) 
        if self.args.Resolved0Btag:
            makeExclusiveResolvedNoBtagSelection(self,tupleSel)
        if self.args.Resolved1Btag:
            makeExclusiveResolvedOneBtagSelection(self,tupleSel)
        if self.args.Resolved2Btag:
            makeExclusiveResolvedTwoBtagsSelection(self,tupleSel)
        if self.args.Boosted:
            makeInclusiveBoostedSelection(self,tupleSel)

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
            print ( ["Ak8","Boosted"])
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
        varsToKeep['total_weight'] = tupleSel.sel.weight

        #return leptonSel.sel, varsToKeep
        return tupleSel.sel, varsToKeep
