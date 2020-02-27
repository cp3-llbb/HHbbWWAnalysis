import os
import sys

from bamboo import treefunctions as op
from bamboo.analysismodules import SkimmerModule

sys.path.append('/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/') # Add scripts in this directory -- TODO : make cleaner
from BaseHHtobbWW import BaseNanoHHtobbWW

#===============================================================================================#
#                                 PlotterHHtobbWW                                               #
#===============================================================================================#
class SkimmerNanoHHtobbWW(BaseNanoHHtobbWW,SkimmerModule):
    """ Plotter module: HH->bbW(->e/µ nu)W(->e/µ nu) histograms from NanoAOD """
    def __init__(self, args):
        super(SkimmerNanoHHtobbWW, self).__init__(args)
        #self.BaseClass = BaseNanoHHtobbWW(args)

    #def prepareTree(self, tree, sample=None, sampleCfg=None):
    #    #return super(SkimmerNanoHHtobbWW,self).prepareTree(tree, sample, sampleCfg)
    #    return self.BaseClass.prepareTree(tree, sample, sampleCfg)

    def defineSkimSelection(self, t, noSel, sample=None, sampleCfg=None): 
        noSel = self.prepareObjects(t, noSel, sample, sampleCfg)
        #noSel, objDict = self.BaseClass.prepareObjects(t, noSel, sample, sampleCfg)

        # Initialize varsToKeep dict #
        varsToKeep = dict()  

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
        varsToKeep["n_presel_ak4Jet"]   = op.static_cast("UInt_t",op.rng_len(self.ak4JetsPreSel))    # Should I include cleaning ?
        varsToKeep["n_presel_ak8Jet"]   = op.static_cast("UInt_t",op.rng_len(self.ak8JetsPreSel))    # Should I include cleaning ?

        # Muons #
        for i in range(1,3): # 2 leading muons
            varsToKeep["mu{}_pt".format(i)]                    = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].pt, op.c_float(-9999., "float"))
            varsToKeep["mu{}_eta".format(i)]                   = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].eta, op.c_float(-9999.))
            varsToKeep["mu{}_phi".format(i)]                   = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].phi, op.c_float(-9999.))
            varsToKeep["mu{}_E".format(i)]                     = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].p4.E(), op.c_float(-9999., "float"))
            varsToKeep["mu{}_charge".format(i)]                = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].charge, op.c_int(-9999.))
            #varsToKeep["mu{}_conept".format(i)]                = 
            varsToKeep["mu{}_miniRelIso".format(i)]            = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].miniPFRelIso_all, op.c_float(-9999.))
            varsToKeep["mu{}_PFRelIso04".format(i)]            = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].pfRelIso04_all, op.c_float(-9999.))
            #varsToKeep["mu{}_jetNDauChargedMVASel".format(i)]  = 
            varsToKeep["mu{}_jetPtRel".format(i)]              = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].jetPtRelv2, op.c_float(-9999.))
            varsToKeep["mu{}_jetRelIso".format(i)]             = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].jetRelIso, op.c_float(-9999.))
            #varsToKeep["mu{}_jetDeepJet".format(i)]            = 
            varsToKeep["mu{}_sip3D".format(i)]                 = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].sip3d, op.c_float(-9999.))
            varsToKeep["mu{}_dxy".format(i)]                   = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].dxy, op.c_float(-9999.))
            varsToKeep["mu{}_dxyAbs".format(i)]                = op.switch(op.rng_len(self.muonsPreSel) >= i, op.abs(self.muonsPreSel[i-1].dxy), op.c_float(-9999.))
            varsToKeep["mu{}_dz".format(i)]                    = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].dz, op.c_float(-9999.))
            varsToKeep["mu{}_segmentCompatibility".format(i)]  = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].segmentComp, op.c_float(-9999.))
            varsToKeep["mu{}_leptonMVA".format(i)]             = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].mvaTTH, op.c_float(-9999.))
            varsToKeep["mu{}_mediumID".format(i)]              = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].mediumId, op.c_float(-9999.,"Bool_t"))
            varsToKeep["mu{}_dpt_div_pt".format(i)]            = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].tunepRelPt, op.c_float(-9999.))  # Not sure
            varsToKeep["mu{}_isfakeablesel".format(i)]         = op.switch(op.rng_len(self.muonsPreSel) >= i, self.lambda_muonFakeSel(self.muonsPreSel[i-1]), op.c_bool(-9999.))
            varsToKeep["mu{}_ismvasel".format(i)]              = op.switch(op.rng_len(self.muonsPreSel) >= i, self.lambda_muonTightSel(self.muonsPreSel[i-1]), op.c_bool(-9999.))  
            #varsToKeep["mu{}_isGenMatched".format(i)]          = 
        
        # Electrons #
        for i in range(1,3): # 2 leading electrons 
            varsToKeep["ele{}_pt".format(i)]                    = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].pt, op.c_float(-9999.))
            varsToKeep["ele{}_eta".format(i)]                   = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].eta, op.c_float(-9999.))
            varsToKeep["ele{}_phi".format(i)]                   = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].phi, op.c_float(-9999.))
            varsToKeep["ele{}_E".format(i)]                     = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].p4.E(), op.c_float(-9999., "float"))
            varsToKeep["ele{}_charge".format(i)]                = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].charge, op.c_int(-9999.))
            #varsToKeep["ele{}_conept".format(i)]                = 
            varsToKeep["ele{}_miniRelIso".format(i)]            = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].miniPFRelIso_all, op.c_float(-9999.))
            varsToKeep["ele{}_PFRelIso03".format(i)]            = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].pfRelIso03_all, op.c_float(-9999.)) # Iso03, Iso04 not in NanoAOD
            #varsToKeep["ele{}_jetNDauChargedMVASel".format(i)]  = 
            varsToKeep["ele{}_jetPtRel".format(i)]              = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].jetPtRelv2, op.c_float(-9999.))
            varsToKeep["ele{}_jetRelIso".format(i)]             = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].jetRelIso, op.c_float(-9999.))
            #varsToKeep["ele{}_jetDeepJet".format(i)]            = 
            varsToKeep["ele{}_sip3D".format(i)]                 = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].sip3d, op.c_float(-9999.))
            varsToKeep["ele{}_dxy".format(i)]                   = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].dxy, op.c_float(-9999.))
            varsToKeep["ele{}_dxyAbs".format(i)]                = op.switch(op.rng_len(self.electronsPreSel) >= i, op.abs(self.electronsPreSel[i-1].dxy), op.c_float(-9999.))
            varsToKeep["ele{}_dz".format(i)]                    = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].dz, op.c_float(-9999.))
            #varsToKeep["ele{}_ntMVAeleID".format(i)]            = non-triggering POG MVA electron ID discriminator
            varsToKeep["ele{}_leptonMVA".format(i)]             = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].mvaTTH, op.c_float(-9999.))
            varsToKeep["ele{}_passesConversionVeto".format(i)]  = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].convVeto, op.c_float(-9999.,"Bool_t"))
            varsToKeep["ele{}_nMissingHits".format(i)]          = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].lostHits, op.c_float(-9999.,"UChar_t"))
            varsToKeep["ele{}_sigmaEtaEta".format(i)]           = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].sieie, op.c_float(-9999.))
            varsToKeep["ele{}_HoE".format(i)]                   = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].hoe, op.c_float(-9999.))
            varsToKeep["ele{}_OoEminusOoP".format(i)]           = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].eInvMinusPInv, op.c_float(-9999.))
            varsToKeep["ele{}_isfakeablesel".format(i)]         = op.switch(op.rng_len(self.electronsPreSel) >= i, self.lambda_electronFakeSel(self.electronsPreSel[i-1]), op.c_bool(-9999.))
            varsToKeep["ele{}_ismvasel".format(i)]              = op.switch(op.rng_len(self.electronsPreSel) >= i, self.lambda_electronTightSel(self.electronsPreSel[i-1]), op.c_bool(-9999.))  
           #varsToKeep["ele{}_isGenMatched".format(i)]          = 
 
        # AK4 Jets #
        for i in range(1,5): # 4 leading jets 
            varsToKeep["ak4Jet{}_pt".format(i)]                    = op.switch(op.rng_len(self.ak4JetsPreSel) >= i, self.ak4JetsPreSel[i-1].pt, op.c_float(-9999.,"float"))
            varsToKeep["ak4Jet{}_eta".format(i)]                   = op.switch(op.rng_len(self.ak4JetsPreSel) >= i, self.ak4JetsPreSel[i-1].eta, op.c_float(-9999.))
            varsToKeep["ak4Jet{}_phi".format(i)]                   = op.switch(op.rng_len(self.ak4JetsPreSel) >= i, self.ak4JetsPreSel[i-1].phi, op.c_float(-9999.))
            varsToKeep["ak4Jet{}_E".format(i)]                     = op.switch(op.rng_len(self.ak4JetsPreSel) >= i, self.ak4JetsPreSel[i-1].p4.E(), op.c_float(-9999., "float"))
            varsToKeep["ak4Jet{}_CSV".format(i)]                   = op.switch(op.rng_len(self.ak4JetsPreSel) >= i, self.ak4JetsPreSel[i-1].btagDeepFlavB, op.c_float(-9999.))

        # AK8 Jets #
        for i in range(1,5): # 2 leading fatjets 
            varsToKeep["ak8Jet{}_pt".format(i)]                    = op.switch(op.rng_len(self.ak8JetsPreSel) >= i, self.ak8JetsPreSel[i-1].pt, op.c_float(-9999.))
            varsToKeep["ak8Jet{}_eta".format(i)]                   = op.switch(op.rng_len(self.ak8JetsPreSel) >= i, self.ak8JetsPreSel[i-1].eta, op.c_float(-9999.))
            varsToKeep["ak8Jet{}_phi".format(i)]                   = op.switch(op.rng_len(self.ak8JetsPreSel) >= i, self.ak8JetsPreSel[i-1].phi, op.c_float(-9999.))
            varsToKeep["ak8Jet{}_E".format(i)]                     = op.switch(op.rng_len(self.ak8JetsPreSel) >= i, self.ak8JetsPreSel[i-1].p4.E(), op.c_float(-9999., "float"))
            varsToKeep["ak8Jet{}_msoftdrop".format(i)]             = op.switch(op.rng_len(self.ak8JetsPreSel) >= i, self.ak8JetsPreSel[i-1].msoftdrop, op.c_float(-9999.))
            varsToKeep["ak8Jet{}_tau1".format(i)]                  = op.switch(op.rng_len(self.ak8JetsPreSel) >= i, self.ak8JetsPreSel[i-1].tau1, op.c_float(-9999.))
            varsToKeep["ak8Jet{}_tau2".format(i)]                  = op.switch(op.rng_len(self.ak8JetsPreSel) >= i, self.ak8JetsPreSel[i-1].tau2, op.c_float(-9999.))
            varsToKeep["ak8Jet{}_subjet1_pt".format(i)]            = op.switch(op.rng_len(self.ak8JetsPreSel) >= i, self.ak8JetsPreSel[i-1].subJet1.pt, op.c_float(-9999.))
            varsToKeep["ak8Jet{}_subjet1_eta".format(i)]           = op.switch(op.rng_len(self.ak8JetsPreSel) >= i, self.ak8JetsPreSel[i-1].subJet1.eta, op.c_float(-9999.))
            varsToKeep["ak8Jet{}_subjet1_phi".format(i)]           = op.switch(op.rng_len(self.ak8JetsPreSel) >= i, self.ak8JetsPreSel[i-1].subJet1.phi, op.c_float(-9999.))
            varsToKeep["ak8Jet{}_subjet1_CSV".format(i)]           = op.switch(op.rng_len(self.ak8JetsPreSel) >= i, self.ak8JetsPreSel[i-1].subJet1.btagDeepB, op.c_float(-9999.))
            varsToKeep["ak8Jet{}_subjet2_pt".format(i)]            = op.switch(op.rng_len(self.ak8JetsPreSel) >= i, self.ak8JetsPreSel[i-1].subJet2.pt, op.c_float(-9999.))
            varsToKeep["ak8Jet{}_subjet2_eta".format(i)]           = op.switch(op.rng_len(self.ak8JetsPreSel) >= i, self.ak8JetsPreSel[i-1].subJet2.eta, op.c_float(-9999.))
            varsToKeep["ak8Jet{}_subjet2_phi".format(i)]           = op.switch(op.rng_len(self.ak8JetsPreSel) >= i, self.ak8JetsPreSel[i-1].subJet2.phi, op.c_float(-9999.))
            varsToKeep["ak8Jet{}_subjet2_CSV".format(i)]           = op.switch(op.rng_len(self.ak8JetsPreSel) >= i, self.ak8JetsPreSel[i-1].subJet2.btagDeepB, op.c_float(-9999.))

        # MET #
         
        varsToKeep["PFMET"]    = self.corrMET.pt
        varsToKeep["PFMETphi"] = self.corrMET.phi

        # HME #

        # Event Weight #
        varsToKeep["MC_weight"] = t.genWeight
        varsToKeep["PU_weight"] = self.PUWeight

#         varsToKeep["nSelMuons"] = op.static_cast("UInt_t", op.rng_len(muons)) ## TBranch doesn't accept size_t         
#         varsToKeep["selMu_miniPFRelIsoNeu"] = op.map(muons, lambda mu : mu.miniPFRelIso_all - mu.miniPFRelIso_chg)


        return noSel, varsToKeep

