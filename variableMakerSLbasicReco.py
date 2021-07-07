from bamboo import treefunctions as op
from highlevelLambdas import highlevelLambdas
import random 
import sys
import os
from bamboo.root import loadHeader
loadHeader("/home/ucl/cp3/gsaha/bamboodev/HHbbWWAnalysis/jpa.h")

# --------------------- Choose Wjets pair for fullReco categories ---------------------- #
def chooseWjj(self, bjet1, bjet2, lep, lightJetPairs, isResolved=False, isBoosted=False):
    if isResolved : 
        mH1_mH2_diff = lambda dijet: op.abs((dijet[0].p4+dijet[1].p4+lep.p4+self.corrMET.p4).M() - (self.HLL.bJetCorrP4(bjet1) + self.HLL.bJetCorrP4(bjet2)).M())
    elif isBoosted:
        mH1_mH2_diff = lambda dijet: op.abs((dijet[0].p4+dijet[1].p4+lep.p4+self.corrMET.p4).M() - bjet1.p4.M())
    else:
        raise RuntimeError('Please mention isBoosted = True or isResolved=True') 

    findWjj = lambda jetPairs : op.sort(jetPairs, mH1_mH2_diff)
    return findWjj(lightJetPairs)

# ----------------------------------------- VBF ---------------------------------------- #
def findVBFJets_Resolved(self, nAk4JetsSelected, jet1=None, jet2=None, jet3=None, jet4=None):
    if nAk4JetsSelected > 3 : 
        lambda_findVBFjets = lambda j : op.AND(op.deltaR(jet1.p4, j.p4) > 0.8,
                                               op.deltaR(jet2.p4, j.p4) > 0.8,
                                               op.deltaR(jet3.p4, j.p4) > 0.8,
                                               op.deltaR(jet4.p4, j.p4) > 0.8)
    else :
        lambda_findVBFjets = lambda j : op.AND(op.deltaR(jet1.p4, j.p4) > 0.8,
                                               op.deltaR(jet2.p4, j.p4) > 0.8,
                                               op.deltaR(jet3.p4, j.p4) > 0.8)
        
    return op.sort(op.combine(op.select(self.VBFJets, lambda_findVBFjets), N=2, pred=self.lambda_VBFPair), 
                   lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))


def findVBFJets_Boosted(self, nAk4JetsSelected, fatJet, jet3=None, jet4=None):
    if nAk4JetsSelected > 1 : 
        lambda_findVBFjets = lambda j : op.AND(op.deltaR(fatJet.p4, j.p4) > 1.2, 
                                               op.AND(op.deltaR(jet3.p4, j.p4) > 0.8,
                                                      op.deltaR(jet4.p4, j.p4) > 0.8))
    else :
        lambda_findVBFjets = lambda j : op.AND(op.deltaR(fatJet.p4, j.p4) > 1.2, 
                                               op.deltaR(jet3.p4, j.p4) > 0.8)
        
    return op.sort(op.combine(op.select(self.VBFJets, lambda_findVBFjets), N=2, pred=self.lambda_VBFPair), 
                   lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))

# ----------------------------------- Resolved Classic Inputs ------------------------------------- #
def returnClassicInputs_Resolved(self, lepton, jet1=None, jet2=None, jet3=None, jet4=None, category=None):
    if category == 'Res2b2Wj':
        VBFJetPairs = findVBFJets_Resolved(self, 4, jet1, jet2, jet3, jet4)
                      
    elif category == 'Res1b3Wj':
        VBFJetPairs = findVBFJets_Resolved(self, 4, jet1, jet2, jet3, jet4)

    elif category == 'Res2b1Wj':
        VBFJetPairs = findVBFJets_Resolved(self, 3, jet1, jet2, jet3)

    elif category == 'Res1b2Wj':
        VBFJetPairs = findVBFJets_Resolved(self, 3, jet1, jet2, jet3)

    else:
        raise RuntimeError('Resolved JPA category is not mentioned!')

    isFullReco = any([category == 'Res2b2Wj', category == 'Res1b3Wj'])
    isMissReco = any([category == 'Res2b1Wj', category == 'Res1b2Wj'])

    lep_Px = lepton.p4.Px()
    lep_Py = lepton.p4.Py()
    lep_Pz = lepton.p4.Pz()
    lep_E  = lepton.p4.E()
    bj1_Px = self.HLL.bJetCorrP4(jet1).Px()
    bj1_Py = self.HLL.bJetCorrP4(jet1).Py()
    bj1_Pz = self.HLL.bJetCorrP4(jet1).Pz()
    bj1_E  = self.HLL.bJetCorrP4(jet1).E()
    bj2_Px = self.HLL.bJetCorrP4(jet2).Px()
    bj2_Py = self.HLL.bJetCorrP4(jet2).Py()
    bj2_Pz = self.HLL.bJetCorrP4(jet2).Pz()
    bj2_E  = self.HLL.bJetCorrP4(jet2).E()
    wj1_Px = jet3.p4.Px()
    wj1_Py = jet3.p4.Py()
    wj1_Pz = jet3.p4.Pz()
    wj1_E  = jet3.p4.E()
    wj2_Px = jet4.p4.Px() if isFullReco else op.c_float(0.)
    wj2_Py = jet4.p4.Py() if isFullReco else op.c_float(0.)
    wj2_Pz = jet4.p4.Pz() if isFullReco else op.c_float(0.)
    wj2_E  = jet4.p4.E()  if isFullReco else op.c_float(0.)

    # Jet-1 variables :: Available for each category
    bj1_M              = self.HLL.bJetCorrP4(jet1).M()
    bj1_pt             = self.HLL.bJetCorrP4(jet1).Pt()
    bj1_eta            = self.HLL.bJetCorrP4(jet1).Eta()
    bj1_phi            = self.HLL.bJetCorrP4(jet1).Phi()
    bj1_bTagDeepFlavB  = jet1.btagDeepFlavB
    bj1LepDR           = op.deltaR(jet1.p4, lepton.p4)
    bj1LepDEta         = op.abs(jet1.eta - lepton.eta)
    bj1LepDPhi         = op.abs(op.deltaPhi(jet1.p4, lepton.p4))
    bj1MetDPhi         = op.abs(self.HLL.SinglepMet_dPhi(jet1, self.corrMET))

    # Jet-2 variables :: Available for 2b2Wj, 2b1Wj. 2b0Wj
    bj2_M             = self.HLL.bJetCorrP4(jet2).M()
    bj2_pt            = self.HLL.bJetCorrP4(jet2).Pt()
    bj2_eta           = self.HLL.bJetCorrP4(jet2).Eta()
    bj2_phi           = self.HLL.bJetCorrP4(jet2).Phi()
    bj2_bTagDeepFlavB = jet2.btagDeepFlavB
    bj2LepDR          = op.deltaR(jet2.p4, lepton.p4)
    bj2LepDEta        = op.abs(jet2.eta - lepton.eta)
    bj2LepDPhi        = op.abs(op.deltaPhi(jet2.p4, lepton.p4))
    bj2MetDPhi        = op.abs(self.HLL.SinglepMet_dPhi(jet2, self.corrMET))

    # h->bb related
    HbbLepDR          = op.deltaR(jet1.p4+jet2.p4, lepton.p4)
    HbbLepDEta        = op.abs((jet1.p4+jet2.p4).Eta() - lepton.eta)
    HbbLepDPhi        = op.abs(op.deltaPhi(jet1.p4+jet2.p4, lepton.p4))
    HbbMetDPhi        = op.abs(op.deltaPhi(jet1.p4+jet2.p4, self.corrMET.p4))
    cosThetaS_Hbb     = op.extMethod("HHbbWWJPA::cosThetaS", returnType="float")(jet1.p4, jet2.p4)

    # Jet-3 variables :: Available for 2b2Wj, 2b1Wj, 1b2Wj, 1b1Wj
    wj1_M             = jet3.mass
    wj1_pt            = jet3.pt
    wj1_eta           = jet3.eta
    wj1_phi           = jet3.phi
    wj1_bTagDeepFlavB = jet3.btagDeepFlavB
    wj1LepDR          = op.deltaR(jet3.p4, lepton.p4)
    wj1LepDEta        = op.abs(jet3.eta - lepton.eta)
    wj1LepDPhi        = op.abs(op.deltaPhi(jet3.p4, lepton.p4))
    wj1MetDPhi        = op.abs(self.HLL.SinglepMet_dPhi(jet3, self.corrMET))
    
    # Jet-4 variables :: Available for 2b2Wj, 1b2Wj
    wj2_M             = jet4.mass                                                              if isFullReco else op.c_float(0.)
    wj2_pt            = jet4.pt                                                                if isFullReco else op.c_float(0.)
    wj2_eta           = jet4.eta                                                               if isFullReco else op.c_float(0.)
    wj2_phi           = jet4.phi                                                               if isFullReco else op.c_float(0.)
    wj2_bTagDeepFlavB = jet4.btagDeepFlavB                                                     if isFullReco else op.c_float(0.)
    wj2LepDR          = op.deltaR(jet4.p4, lepton.p4)                                          if isFullReco else op.c_float(0.) 
    wj2LepDEta        = op.abs(jet4.eta - lepton.eta)                                          if isFullReco else op.c_float(0.)
    wj2LepDPhi        = op.abs(op.deltaPhi(jet4.p4, lepton.p4))                                if isFullReco else op.c_float(0.)
    wj2MetDPhi        = op.abs(self.HLL.SinglepMet_dPhi(jet4, self.corrMET))                   if isFullReco else op.c_float(0.)

    WjjLepDR          = op.deltaR(jet3.p4+jet4.p4, lepton.p4)                                  if isFullReco else op.c_float(0.)
    WjjLepDEta        = op.abs((jet3.p4+jet4.p4).Eta() - lepton.eta)                           if isFullReco else op.c_float(0.)
    WjjLepDPhi        = op.abs(op.deltaPhi(jet3.p4+jet4.p4, lepton.p4))                        if isFullReco else op.c_float(0.)
    WjjMetDPhi        = op.abs(op.deltaPhi(jet3.p4+jet4.p4, self.corrMET.p4))                  if isFullReco else op.c_float(0.)

    # h->ww related
    HWW_Mass          = self.HLL.HWW_simple(jet3.p4,jet4.p4,lepton.p4,self.corrMET).M()        if isFullReco else op.c_float(0.) # w.r.t. neutrino
    HWW_Simple_Mass   = self.HLL.HWW_met_simple(jet3.p4,jet4.p4,lepton.p4,self.corrMET.p4).M() if isFullReco else op.c_float(0.) # not w.r.t. neutrino
    HWW_dR            = self.HLL.dR_Hww(jet3.p4,jet4.p4,lepton.p4,self.corrMET)                if isFullReco else op.c_float(0.)
    HWW_dEta          = self.HLL.dEta_Hww(jet3.p4,jet4.p4,lepton.p4,self.corrMET)              if isFullReco else op.c_float(0.)
    HWW_dPhi          = self.HLL.dPhi_Hww(jet3.p4,jet4.p4,lepton.p4,self.corrMET)              if isFullReco else op.c_float(0.)
    angleBetWWPlane   = self.HLL.angleWWplane(lepton.p4, self.corrMET, jet3.p4, jet4.p4)       if isFullReco else op.c_float(0.)
    
    cosThetaS_Wjj_simple    = op.extMethod("HHbbWWJPA::cosThetaS", returnType="float")(jet3.p4, jet4.p4)                          if isFullReco else op.c_float(0.)
    cosThetaS_WW_simple_met = op.extMethod("HHbbWWJPA::cosThetaS", returnType="float")(self.HLL.Wjj_simple(jet3.p4,jet4.p4), self.HLL.Wlep_met_simple(lepton.p4, self.corrMET.p4)) if isFullReco else op.c_float(0.)
    cosThetaS_HH_simple_met = op.extMethod("HHbbWWJPA::cosThetaS", returnType="float")(self.HLL.bJetCorrP4(jet1)+self.HLL.bJetCorrP4(jet2), self.HLL.HWW_met_simple(jet3.p4,jet4.p4,lepton.p4, self.corrMET.p4))  if isFullReco else op.c_float(0.)    
    mHH                     = (jet3.p4+jet4.p4+lepton.p4+self.corrMET.p4+self.HLL.bJetCorrP4(jet1)+self.HLL.bJetCorrP4(jet2)).M() if isFullReco else op.c_float(0.)
    angleBetHWPlane         = self.HLL.angleBetPlanes(jet1.p4,jet2.p4,jet3.p4,jet4.p4)                                            if isFullReco else op.c_float(0.)
    
    # neu p4 only for 1b2Wj and 2b2Wj
    neuPx = self.HLL.neuP4(jet3.p4+jet4.p4+lepton.p4, self.corrMET).Px()                       if isFullReco else op.c_float(0.)
    neuPy = self.HLL.neuP4(jet3.p4+jet4.p4+lepton.p4, self.corrMET).Py()                       if isFullReco else op.c_float(0.)
    neuPz = self.HLL.neuP4(jet3.p4+jet4.p4+lepton.p4, self.corrMET).Pz()                       if isFullReco else op.c_float(0.)
    neuE  = self.HLL.neuP4(jet3.p4+jet4.p4+lepton.p4, self.corrMET).E()                        if isFullReco else op.c_float(0.)
    neuPt = self.HLL.neuP4(jet3.p4+jet4.p4+lepton.p4, self.corrMET).Pt()                       if isFullReco else op.c_float(0.)

    # diJet dR & dPhi
    bj1bj2_DR   = op.deltaR(jet1.p4,jet2.p4)  
    bj1bj2_DPhi = op.abs(op.deltaPhi(jet1.p4,jet2.p4)) 
    bj1bj2_DEta = op.abs(jet1.eta - jet2.eta)

    bj1wj1_DR   = op.deltaR(jet1.p4,jet3.p4) 
    bj1wj1_DPhi = op.abs(op.deltaPhi(jet1.p4,jet3.p4))
    bj1wj1_DEta = op.abs(jet1.eta - jet3.eta)

    bj1wj2_DR   = op.deltaR(jet1.p4,jet4.p4)             if isFullReco else op.c_float(0.)
    bj1wj2_DPhi = op.abs(op.deltaPhi(jet1.p4,jet4.p4))   if isFullReco else op.c_float(0.)
    bj1wj2_DEta = op.abs(jet1.eta - jet4.eta)            if isFullReco else op.c_float(0.)

    bj2wj1_DR   = op.deltaR(jet2.p4,jet3.p4)
    bj2wj1_DPhi = op.abs(op.deltaPhi(jet2.p4,jet3.p4))
    bj2wj1_DEta = op.abs(jet2.eta - jet3.eta)

    bj2wj2_DR   = op.deltaR(jet2.p4,jet4.p4)             if isFullReco else op.c_float(0.)
    bj2wj2_DPhi = op.abs(op.deltaPhi(jet2.p4,jet4.p4))   if isFullReco else op.c_float(0.)
    bj2wj2_DEta = op.abs(jet2.eta - jet4.eta)            if isFullReco else op.c_float(0.)

    wj1wj2_DR   = op.deltaR(jet3.p4,jet4.p4)             if isFullReco else op.c_float(0.)
    wj1wj2_DPhi = op.abs(op.deltaPhi(jet3.p4,jet4.p4))   if isFullReco else op.c_float(0.)
    wj1wj2_DEta = op.abs(jet3.eta - jet4.eta)            if isFullReco else op.c_float(0.)

    # VBF variables
    VBFj1Px        = op.switch(op.rng_len(VBFJetPairs) > 0, 
                               op.switch(VBFJetPairs[0][0].pt > VBFJetPairs[0][1].pt,
                                         VBFJetPairs[0][0].p4.Px(), VBFJetPairs[0][1].p4.Px()),
                               op.c_float(0.))
    VBFj2Px        = op.switch(op.rng_len(VBFJetPairs) > 0, 
                               op.switch(VBFJetPairs[0][0].pt > VBFJetPairs[0][1].pt,
                                         VBFJetPairs[0][1].p4.Px(), VBFJetPairs[0][0].p4.Px()),
                               op.c_float(0.))
    VBFj1Py        = op.switch(op.rng_len(VBFJetPairs) > 0, 
                               op.switch(VBFJetPairs[0][0].pt > VBFJetPairs[0][1].pt,
                                         VBFJetPairs[0][0].p4.Py(), VBFJetPairs[0][1].p4.Py()),
                               op.c_float(0.))
    VBFj2Py        = op.switch(op.rng_len(VBFJetPairs) > 0, 
                               op.switch(VBFJetPairs[0][0].pt > VBFJetPairs[0][1].pt,
                                         VBFJetPairs[0][1].p4.Py(), VBFJetPairs[0][0].p4.Py()),
                               op.c_float(0.))
    VBFj1Pz        = op.switch(op.rng_len(VBFJetPairs) > 0, 
                               op.switch(VBFJetPairs[0][0].pt > VBFJetPairs[0][1].pt,
                                         VBFJetPairs[0][0].p4.Pz(), VBFJetPairs[0][1].p4.Pz()),
                               op.c_float(0.))
    VBFj2Pz        = op.switch(op.rng_len(VBFJetPairs) > 0, 
                               op.switch(VBFJetPairs[0][0].pt > VBFJetPairs[0][1].pt,
                                         VBFJetPairs[0][1].p4.Pz(), VBFJetPairs[0][0].p4.Pz()),
                               op.c_float(0.))
    VBFj1E         = op.switch(op.rng_len(VBFJetPairs) > 0, 
                               op.switch(VBFJetPairs[0][0].pt > VBFJetPairs[0][1].pt,
                                         VBFJetPairs[0][0].p4.E(), VBFJetPairs[0][1].p4.E()),
                               op.c_float(0.))
    VBFj2E         = op.switch(op.rng_len(VBFJetPairs) > 0, 
                               op.switch(VBFJetPairs[0][0].pt > VBFJetPairs[0][1].pt,
                                         VBFJetPairs[0][1].p4.E(), VBFJetPairs[0][0].p4.E()),
                               op.c_float(0.))
    VBFj1pt        = op.switch(op.rng_len(VBFJetPairs) > 0, 
                                 op.switch(VBFJetPairs[0][0].pt > VBFJetPairs[0][1].pt,
                                           VBFJetPairs[0][0].pt, VBFJetPairs[0][1].pt),
                                 op.c_float(0.))
    VBFj2pt        = op.switch(op.rng_len(VBFJetPairs) > 0, 
                               op.switch(VBFJetPairs[0][0].pt > VBFJetPairs[0][1].pt,
                                         VBFJetPairs[0][1].pt, VBFJetPairs[0][0].pt),
                               op.c_float(0.))
    VBFj1eta       = op.switch(op.rng_len(VBFJetPairs) > 0, 
                               op.switch(VBFJetPairs[0][0].pt > VBFJetPairs[0][1].pt,
                                         VBFJetPairs[0][0].eta, VBFJetPairs[0][1].eta),
                               op.c_float(0.))
    VBFj2eta       = op.switch(op.rng_len(VBFJetPairs) > 0, 
                               op.switch(VBFJetPairs[0][0].pt > VBFJetPairs[0][1].pt,
                                         VBFJetPairs[0][1].eta, VBFJetPairs[0][0].eta),
                               op.c_float(0.))
    VBFj1j2dEta    = op.switch(op.rng_len(VBFJetPairs) > 0,op.abs(VBFJetPairs[0][0].eta - VBFJetPairs[0][1].eta),
                               op.c_float(0.))
    
    VBFj1j2dPhi    = op.switch(op.rng_len(VBFJetPairs) > 0,op.abs(op.deltaPhi(VBFJetPairs[0][0].p4,VBFJetPairs[0][1].p4)),
                               op.c_float(0.))
    
    VBFj1j2invM    = op.switch(op.rng_len(VBFJetPairs) > 0,op.invariant_mass(VBFJetPairs[0][0].p4,  VBFJetPairs[0][1].p4),
                               op.c_float(0.))
    
    VBF_tag        = op.c_int(op.rng_len(VBFJetPairs) > 0)
    zeppenfeldVar  = op.switch(op.rng_len(VBFJetPairs) > 0, 
                               op.abs(lepton.eta - (VBFJetPairs[0][0].eta + VBFJetPairs[0][1].eta)/2.0)/op.abs(VBFJetPairs[0][0].eta - VBFJetPairs[0][1].eta),
                               op.c_float(0.))

    if isFullReco : 
        pt_taggedJetsAvg = (self.HLL.bJetCorrP4(jet1).Pt()+self.HLL.bJetCorrP4(jet2).Pt()+jet3.pt+jet4.pt)/4.0
        mass_taggedJetsAvg = (self.HLL.bJetCorrP4(jet1).M()+self.HLL.bJetCorrP4(jet2).M()+jet3.mass+jet4.mass)/4.0
        bj1bj2_pt      = (self.HLL.bJetCorrP4(jet1)+self.HLL.bJetCorrP4(jet2)).Pt()
        bj1bj2_M       = op.invariant_mass(self.HLL.bJetCorrP4(jet1),self.HLL.bJetCorrP4(jet2))
        mT_top_3particle  = op.min(self.HLL.mT2(self.HLL.bJetCorrP4(jet1),lepton.p4 ,self.corrMET.p4), 
                                   self.HLL.mT2(self.HLL.bJetCorrP4(jet2), lepton.p4, self.corrMET.p4))
        avgJetDR       = (op.deltaR(jet1.p4,jet2.p4) + op.deltaR(jet1.p4,jet3.p4) + op.deltaR(jet1.p4,jet4.p4) + op.deltaR(jet2.p4,jet3.p4) 
                          + op.deltaR(jet2.p4,jet4.p4) + op.deltaR(jet3.p4,jet4.p4))/6.0
        avgJetDPhi     = (op.abs(op.deltaPhi(jet1.p4,jet2.p4)) + op.abs(op.deltaPhi(jet1.p4,jet3.p4)) + op.abs(op.deltaPhi(jet1.p4,jet4.p4)) 
                          + op.abs(op.deltaPhi(jet2.p4,jet3.p4)) + op.abs(op.deltaPhi(jet2.p4,jet4.p4)) + op.abs(op.deltaPhi(jet3.p4,jet4.p4)))/6.0
        avgJetDEta     = (op.abs(jet1.eta-jet2.eta) + op.abs(jet1.eta-jet3.eta) + op.abs(jet1.eta-jet4.eta) + op.abs(jet2.eta-jet3.eta) 
                          + op.abs(jet2.eta-jet4.eta) + op.abs(jet3.eta-jet4.eta))/6.0
        minJetDR       = self.HLL.MinDiJetDRTight(jet1,jet2,jet3,jet4)
        minJetDPhi     = self.HLL.MinDiJetDPhiTight(jet1,jet2,jet3,jet4)
        minJetDEta     = self.HLL.MinDiJetDEtaTight(jet1,jet2,jet3,jet4)
        maxJetDR       = self.HLL.MaxDiJetDRTight(jet1,jet2,jet3,jet4)
        maxJetDEta     = self.HLL.MaxDiJetDEtaTight(jet1,jet2,jet3,jet4)
        maxJetDPhi     = self.HLL.MaxDiJetDPhiTight(jet1,jet2,jet3,jet4)
        minLepJetDR    = self.HLL.MinDR_lep4j(lepton,jet1,jet2,jet3,jet4)
        minLepJetDPhi  = self.HLL.MinDPhi_lep4j(lepton,jet1,jet2,jet3,jet4)
        minLepJetDEta  = self.HLL.MinDEta_lep4j(lepton,jet1,jet2,jet3,jet4)
        maxLepJetDR    = self.HLL.MaxDR_lep4j(lepton,jet1,jet2,jet3,jet4)
        maxLepJetDPhi  = self.HLL.MaxDPhi_lep4j(lepton,jet1,jet2,jet3,jet4)
        maxLepJetDEta  = self.HLL.MaxDEta_lep4j(lepton,jet1,jet2,jet3,jet4)
        HT2_lepJetMet  = self.HLL.HT2_2b2Wj(lepton,jet1,jet2,jet3,jet4,self.corrMET)
        HT2R_lepJetMet = self.HLL.HT2R_2b2Wj(lepton,jet1,jet2,jet3,jet4,self.corrMET)
        wj1wj2_pt      = (jet3.p4 + jet4.p4).Pt()   
        wj1wj2_M       = op.invariant_mass(jet3.p4, jet4.p4)  
        w1w2_MT        = self.HLL.MT_W1W2_ljj(lepton,jet3,jet4,self.corrMET) # MT of H->WW->lep, jet3, jet4  and Met
        recocat        = op.c_int(1) if category == 'Res2b2Wj' else op.c_int(2)

    if isMissReco : 
        pt_taggedJetsAvg = (self.HLL.bJetCorrP4(jet1).Pt()+self.HLL.bJetCorrP4(jet2).Pt()+jet3.pt)/3.0
        mass_taggedJetsAvg = (self.HLL.bJetCorrP4(jet1).M()+self.HLL.bJetCorrP4(jet2).M()+jet3.mass)/3.0
        bj1bj2_pt      = (self.HLL.bJetCorrP4(jet1)+self.HLL.bJetCorrP4(jet2)).Pt()
        bj1bj2_M       = op.invariant_mass(self.HLL.bJetCorrP4(jet1),self.HLL.bJetCorrP4(jet2))
        mT_top_3particle  = op.min(self.HLL.mT2(self.HLL.bJetCorrP4(jet1),lepton.p4 ,self.corrMET.p4), 
                                   self.HLL.mT2(self.HLL.bJetCorrP4(jet2), lepton.p4, self.corrMET.p4))
        avgJetDR       = (op.deltaR(jet1.p4,jet2.p4) + op.deltaR(jet1.p4,jet3.p4) + op.deltaR(jet2.p4,jet3.p4))/3.0
        avgJetDPhi     = (op.abs(op.deltaPhi(jet1.p4,jet2.p4)) + op.abs(op.deltaPhi(jet1.p4,jet3.p4)) + op.abs(op.deltaPhi(jet2.p4,jet3.p4)))/3.0
        avgJetDEta     = (op.abs(jet1.eta-jet2.eta) + op.abs(jet1.eta-jet3.eta) + op.abs(jet2.eta-jet3.eta))/3.0
        minJetDR       = self.HLL.MinDiJetDRLoose(jet1,jet2,jet3)
        minJetDPhi     = self.HLL.MinDiJetDPhiLoose(jet1,jet2,jet3)
        minJetDEta     = self.HLL.MinDiJetDEtaLoose(jet1,jet2,jet3)
        maxJetDR       = self.HLL.MaxDiJetDRLoose(jet1,jet2,jet3)
        maxJetDEta     = self.HLL.MaxDiJetDEtaLoose(jet1,jet2,jet3)
        maxJetDPhi     = self.HLL.MaxDiJetDPhiLoose(jet1,jet2,jet3)
        minLepJetDR    = self.HLL.MinDR_lep3j(lepton,jet1,jet2,jet3)
        minLepJetDPhi  = self.HLL.MinDPhi_lep3j(lepton,jet1,jet2,jet3)
        minLepJetDEta  = self.HLL.MinDEta_lep3j(lepton,jet1,jet2,jet3)
        maxLepJetDR    = self.HLL.MaxDR_lep3j(lepton,jet1,jet2,jet3)
        maxLepJetDPhi  = self.HLL.MaxDPhi_lep3j(lepton,jet1,jet2,jet3)
        maxLepJetDEta  = self.HLL.MaxDEta_lep3j(lepton,jet1,jet2,jet3)
        HT2_lepJetMet  = self.HLL.HT2_2b1Wj(lepton,jet1,jet2,jet3,self.corrMET)
        HT2R_lepJetMet = self.HLL.HT2R_2b1Wj(lepton,jet1,jet2,jet3,self.corrMET)
        wj1wj2_pt      = jet3.pt   
        wj1wj2_M       = jet3.p4.M()  
        w1w2_MT        = self.HLL.MT_W1W2_lj(lepton,jet3,self.corrMET) # MT of H->WW->lep, jet3  and Met
        recocat        = op.c_int(3) if category == 'Res2b1Wj' else op.c_int(4)
        

    return{('leppdgId',                'leppdgId',                (45,-22.,22.)):   lepton.pdgId,
           ('lepcharge',               'lepcharge',               (3,-1,1.)):       lepton.charge,
           ('lepE',                    'lepE',                    (50,0.,500.)):    lep_E,
           ('lepPx',                   'lepPx',                   (50,-250.,250.)): lep_Px,
           ('lepPy',                   'lepPy',                   (50,-250.,250.)): lep_Py,
           ('lepPz',                   'lepPz',                   (50,-250.,250.)): lep_Pz,
           ('METpx',                   'METpx',                   (50, 0, 300)):    self.corrMET.p4.Px(),
           ('METpy',                   'METpy',                   (50, 0, 300)):    self.corrMET.p4.Py(),
           ('METpz',                   'METpz',                   (50, 0, 300)):    self.corrMET.p4.Pz(),
           ('METenergy',               'METenergy',               (50, 0, 300)):    self.corrMET.p4.E(),
           ('METpt',                   'METpt',                   (50, 0, 300)):    self.corrMET.pt,
           ('METphi',                  'METphi',                  (50, 0, 300)):    self.corrMET.phi,
           ('METLD',                   'MET_LD',                  (40,0,200)):      self.HLL.MET_LD_DL(self.corrMET, self.ak4Jets, 
                                                                                                       self.electronsFakeSel, self.muonsFakeSel),
           ('nfakeableMuons',          'nfakeableMuons',          (5,0,5)):         op.static_cast("UInt_t",op.rng_len(self.muonsFakeSel)),
           ('nfakeableElectrons',      'nfakeableElectrons',      (5,0,5)):         op.static_cast("UInt_t",op.rng_len(self.electronsFakeSel)),
           ('hT',                      'hT',                      (50,0,300)):      op.static_cast("Float_t", self.HLL.HT_SL(self.ak4Jets)),
           ('nAk4Jets',                'nAk4Jets',                (10,0,10)):       op.static_cast("UInt_t",op.rng_len(self.ak4Jets)),
           ('nAk4BJets',               'nAk4BJets',               (10,0,10)):       op.static_cast("UInt_t",op.rng_len(self.ak4BJets)),
           ('nAk4BJetsLoose',          'nAk4BJetsLoose',          (10,0,10)):       op.static_cast("UInt_t",op.rng_len(self.ak4BJetsLoose)),
           ('nAk8Jets',                'nAk8Jets',                (5,0,5)):         op.static_cast("UInt_t",op.rng_len(self.ak8Jets)),
           ('nAk8BJets',               'nAk8BJets',               (5,0,5)):         op.static_cast("UInt_t",op.rng_len(self.ak8BJets)),
           ('nAk4JetsCleandWithAk8b',  'nAk4JetsCleandWithAk8b',  (5,0,5)):         op.static_cast("UInt_t",op.rng_len(self.ak4JetsCleanedFromAk8b)),
           ('bj1E',                    'bj1E',                    (50,0.,500.)):    bj1_E,
           ('bj1Px',                   'bj1Px',                   (50,-250.,250.)): bj1_Px,
           ('bj1Py',                   'bj1Py',                   (50,-250.,250.)): bj1_Py,
           ('bj1Pz',                   'bj1Pz',                   (50,-250.,250.)): bj1_Pz,
           ('bj2E',                    'bj2E',                    (50,0.,500.)):    bj2_E,
           ('bj2Px',                   'bj2Px',                   (50,-250.,250.)): bj2_Px,
           ('bj2Py',                   'bj2Py',                   (50,-250.,250.)): bj2_Py,
           ('bj2Pz',                   'bj2Pz',                   (50,-250.,250.)): bj2_Pz,
           ('wj1E',                    'wj1E',                    (50,0.,500.)):    wj1_E,
           ('wj1Px',                   'wj1Px',                   (50,-250.,250.)): wj1_Px,
           ('wj1Py',                   'wj1Py',                   (50,-250.,250.)): wj1_Py,
           ('wj1Pz',                   'wj1Pz',                   (50,-250.,250.)): wj1_Pz,
           ('wj2E',                    'wj2E',                    (50,0.,500.)):    wj2_E,
           ('wj2Px',                   'wj2Px',                   (50,-250.,250.)): wj2_Px,
           ('wj2Py',                   'wj2Py',                   (50,-250.,250.)): wj2_Py,
           ('wj2Pz',                   'wj2Pz',                   (50,-250.,250.)): wj2_Pz,
           ('lepmetDPhi',              'lepmet_DPhi',             (20,0,3.2)):      op.abs(self.HLL.SinglepMet_dPhi(lepton, self.corrMET)),
           ('lepmetPt',                'lepmet_pt',               (50,0,300)):      self.HLL.SinglepMet_Pt(lepton, self.corrMET),
           ('lepMT',                   'lep_MT',                  (40,0,200)):      self.HLL.MT(lepton, self.corrMET),
           ('lepConept',               'lep_conept',              (40,0,200)):      self.HLL.lambdaConePt(lepton),
           ('minDRLepAllJets',         'minDR_lep_allJets',       (5,0,5)):         self.HLL.MinDR_part1_partCont(lepton, self.ak4Jets),
           ('minDEtaLepAllJets',       'minDEta_lep_allJets',     (5,0,5)):         self.HLL.MinDEta_part1_partCont(lepton, self.ak4Jets),
           ('bj1M',                    'bj1_M',                   (50,0.,500.)):    bj1_M,
           ('bj1Pt',                   'bj1_Pt',                  (50,-250.,250.)): bj1_pt,
           ('bj1Eta',                  'bj1_Eta',                 (50,0.,5.)):      bj1_eta,
           ('bj1Phi',                  'bj1Phi',                  (50,0.,5.)):      bj1_phi,
           ('bj1bTagDeepFlavB',        'bj1_bTagDeepFlavB',       (50,0.,1.)):      bj1_bTagDeepFlavB,
           ('bj1LepDR',                'bj1LepDR',                (50,0.,5.)):      bj1LepDR,
           ('bj1LepDEta',              'bj1LepDEta',              (50,0.,5.)):      bj1LepDEta,
           ('bj1LepDPhi',              'bj1LepDPhi',              (50,0.,5.)):      bj1LepDPhi,
           ('bj1MetDPhi',              'bj1MetDPhi',              (50,0.,5.)):      bj1MetDPhi,
           ('bj2M',                    'bj2_M',                   (50,0.,500.)):    bj2_M,
           ('bj2Pt',                   'bj2_Pt',                  (50,-250.,250.)): bj2_pt,
           ('bj2Eta',                  'bj2_Eta',                 (50,0.,5.)):      bj2_eta,
           ('bj2Phi',                  'bj2Phi',                  (50,0.,5.)):      bj2_phi,
           ('bj2bTagDeepFlavB',        'bj2_bTagDeepFlavB',       (50,0.,1.)):      bj2_bTagDeepFlavB,
           ('bj2LepDR',                'bj2LepDR',                (50,0.,5.)):      bj2LepDR,
           ('bj2LepDEta',              'bj2LepDEta',              (50,0.,5.)):      bj2LepDEta,
           ('bj2LepDPhi',              'bj2LepDPhi',              (50,0.,5.)):      bj2LepDPhi,
           ('bj2MetDPhi',              'bj2MetDPhi',              (50,0.,5.)):      bj2MetDPhi,
           ('HbbLepDR',                'HbbLepDR',                (50,0.,5.)):      HbbLepDR,
           ('HbbLepDEta',              'HbbLepDEta',              (50,0.,5.)):      HbbLepDEta,
           ('HbbLepDPhi',              'HbbLepDPhi',              (50,0.,5.)):      HbbLepDPhi,
           ('HbbMetDPhi',              'HbbMetDPhi',              (50,0.,5.)):      HbbMetDPhi,
           ('cosThetaSHbb',            'cosThetaS_Hbb',           (50,0.,1.)):      cosThetaS_Hbb,
           ('wj1M',                    'wj1_M',                   (50,0.,500.)):    wj1_M,
           ('wj1Pt',                   'wj1_Pt',                  (50,-250.,250.)): wj1_pt,
           ('wj1Eta',                  'wj1_Eta',                 (50,0.,5.)):      wj1_eta,
           ('wj1Phi',                  'wj1Phi',                  (50,0.,5.)):      wj1_phi,
           ('wj1bTagDeepFlavB',        'wj1_bTagDeepFlavB',       (50,0.,1.)):      wj1_bTagDeepFlavB,
           ('wj1LepDR',                'wj1LepDR',                (50,0.,5.)):      wj1LepDR,
           ('wj1LepDEta',              'wj1LepDEta',              (50,0.,5.)):      wj1LepDEta,
           ('wj1LepDPhi',              'wj1LepDPhi',              (50,0.,5.)):      wj1LepDPhi,
           ('wj1MetDPhi',              'wj1MetDPhi',              (50,0.,5.)):      wj1MetDPhi,
           ('wj2M',                    'wj2_M',                   (50,0.,500.)):    wj2_M,
           ('wj2Pt',                   'wj2_Pt',                  (50,-250.,250.)): wj2_pt,
           ('wj2Eta',                  'wj2_Eta',                 (50,0.,5.)):      wj2_eta,
           ('wj2Phi',                  'wj2Phi',                  (50,0.,5.)):      wj2_phi,
           ('wj2bTagDeepFlavB',        'wj2_bTagDeepFlavB',       (50,0.,1.)):      wj2_bTagDeepFlavB,
           ('wj2LepDR',                'wj2LepDR',                (50,0.,5.)):      wj2LepDR,
           ('wj2LepDEta',              'wj2LepDEta',              (50,0.,5.)):      wj2LepDEta,
           ('wj2LepDPhi',              'wj2LepDPhi',              (50,0.,5.)):      wj2LepDPhi,
           ('wj2MetDPhi',              'w21MetDPhi',              (50,0.,5.)):      wj2MetDPhi,
           ('WjjLepDR',                'wLepDR',                  (50,0.,5.)):      WjjLepDR,
           ('WjjLepDEta',              'wLepDEta',                (50,0.,5.)):      WjjLepDEta,
           ('WjjLepDPhi',              'wLepDPhi',                (50,0.,5.)):      WjjLepDPhi,
           ('WjjMetDPhi',              'wMetDPhi',                (50,0.,5.)):      WjjMetDPhi,
           ('HWWMass',                 'HWW_Mass',                (50,0.,300.)):    HWW_Mass,
           ('HWWSimpleMass',           'HWW_Simple_Mass',         (50,0.,300.)):    HWW_Simple_Mass,
           ('HWWdR',                   'HWW_dR',                  (50,0.,5.)):      HWW_dR,
           ('HWWdPhi',                 'HWW_dPhi',                (50,0.,5.)):      HWW_dPhi,
           ('HWWdEta',                 'HWW_dEta',                (50,0.,5.)):      HWW_dEta,
           ('angleBetWWPlane',         'angleBetWWPlane',         (50,0.,5.)):      angleBetWWPlane,           
           ('cosThetaSWjjSimple',      'cosThetaS_Wjj_simple',    (20,0.,1.)):      cosThetaS_Wjj_simple,
           ('cosThetaSWWSimpleMet',    'cosThetaS_WW_simple_met', (20,0.,1.)):      cosThetaS_WW_simple_met,
           ('cosThetaSHHSimpleMet',    'cosThetaS_HH_simple_met', (20,0.,1.)):      cosThetaS_HH_simple_met,
           ('mHH',                     'mHH',                     (50,0.,300.)):    mHH,
           ('angleBetHWPlane',         'angleBetHWPlane',         (50,0.,5.)):      angleBetHWPlane,
           ('neuPx',                   'neuPx',                   (50,-100.,100.)): neuPx,
           ('neuPy',                   'neuPy',                   (50,-100.,100.)): neuPy,
           ('neuPz',                   'neuPz',                   (50,-100.,100.)): neuPz,
           ('neuPt',                   'neuPt',                   (50,-100.,100.)): neuPt,
           ('neuE',                    'neuE',                    (50,-100.,100.)): neuE,
           ('bj1bj2DR',                'bj1bj2_DR',               (50,0.,5)):       bj1bj2_DR,
           ('bj1bj2DPhi',              'bj1bj2_DPhi',             (50,0.,5)):       bj1bj2_DPhi,
           ('bj1bj2DEta',              'bj1bj2_DEta',             (50,0.,5)):       bj1bj2_DEta,
           ('bj1wj1DR',                'bj1wj1_DR',               (50,0.,5)):       bj1wj1_DR,
           ('bj1wj1DPhi',              'bj1wj1_DPhi',             (50,0.,5)):       bj1wj1_DPhi,
           ('bj1wj1DEta',              'bj1wj1_DEta',             (50,0.,5)):       bj1wj1_DEta,
           ('bj1wj2DR',                'bj1wj2_DR',               (50,0.,5)):       bj1wj2_DR,
           ('bj1wj2DPhi',              'bj1wj2_DPhi',             (50,0.,5)):       bj1wj2_DPhi,
           ('bj1wj2DEta',              'bj1wj2_DEta',             (50,0.,5)):       bj1wj2_DEta,
           ('bj2wj1DR',                'bj2wj1_DR',               (50,0.,5)):       bj2wj1_DR,
           ('bj2wj1DPhi',              'bj2wj1_DPhi',             (50,0.,5)):       bj2wj1_DPhi,
           ('bj2wj1DEta',              'bj2wj1_DEta',             (50,0.,5)):       bj2wj1_DEta,
           ('bj2wj2DR',                'bj2wj2_DR',               (50,0.,5)):       bj2wj2_DR,
           ('bj2wj2DPhi',              'bj2wj2_DPhi',             (50,0.,5)):       bj2wj2_DPhi,
           ('bj2wj2DEta',              'bj2wj2_DEta',             (50,0.,5)):       bj2wj2_DEta,
           ('wj1wj2DR',                'wj1wj2_DR',               (50,0.,5)):       wj1wj2_DR,
           ('wj1wj2DPhi',              'wj1wj2_DPhi',             (50,0.,5)):       wj1wj2_DPhi,
           ('wj1wj2DEta',              'wj1wj2_DEta',             (50,0.,5)):       wj1wj2_DEta,
           ('VBFj1Px',                 'VBFj1Px',                 (50,0.,200)):     VBFj1Px,           
           ('VBFj1Py',                 'VBFj1Py',                 (50,0.,200)):     VBFj1Py,           
           ('VBFj1Pz',                 'VBFj1Pz',                 (50,0.,200)):     VBFj1Pz,           
           ('VBFj1E',                  'VBFj1E',                  (50,0.,200)):     VBFj1E,           
           ('VBFj1pt',                 'VBFj1pt',                 (50,0.,200)):     VBFj1pt,           
           ('VBFj1eta',                'VBFj1eta',                (50,0.,200)):     VBFj1eta,           
           ('VBFj2Px',                 'VBFj2Px',                 (50,0.,200)):     VBFj2Px,           
           ('VBFj2Py',                 'VBFj2Py',                 (50,0.,200)):     VBFj2Py,           
           ('VBFj2Pz',                 'VBFj2Pz',                 (50,0.,200)):     VBFj2Pz,           
           ('VBFj2E',                  'VBFj2E',                  (50,0.,200)):     VBFj2E,           
           ('VBFj2pt',                 'VBFj2pt',                 (50,0.,200)):     VBFj2pt,           
           ('VBFj2eta',                'VBFj2eta',                (50,0.,200)):     VBFj2eta,
           ('VBFj1j2dEta',             'VBFj1j2dEta',             (50,0.,200)):     VBFj1j2dEta,
           ('VBFj1j2dPhi',             'VBFj1j2dPhi',             (50,0.,200)):     VBFj1j2dPhi,
           ('VBFj1j2invM',             'VBFj1j2invM',             (50,0.,200)):     VBFj1j2invM,           
           ('VBFtag',                  'VBF_tag',                 (2, 0, 2)):       VBF_tag,
           ('zeppenfeldVar',           'zeppenfeldVar',           (50,0.,200)):     zeppenfeldVar,           
           ('ptTaggedJetsAvg',         'pt_taggedJetsAvg',        (50,0.,200)):     pt_taggedJetsAvg,
           ('massTaggedJetsAvg',       'mass_taggedJetsAvg',      (50,0.,200)):     mass_taggedJetsAvg,
           ('bj1bj2Pt',                'bj1bj2_pt',               (50,0.,200)):     bj1bj2_pt,
           ('bj1bj2M',                 'bj1bj2_M',                (50,0.,200)):     bj1bj2_M,
           ('mTtop3particle',          'mT_top_3particle',        (50,0.,200)):     mT_top_3particle,
           ('avgJetDR',                'avgJetDR',                (50,0.,200)):     avgJetDR,
           ('minJetDR',                'minJetDR',                (50,0.,200)):     minJetDR,
           ('maxJetDR',                'maxJetDR',                (50,0.,200)):     maxJetDR,
           ('minLepJetDR',             'minLepJetDR',             (50,0.,200)):     minLepJetDR,
           ('maxLepJetDR',             'maxLepJetDR',             (50,0.,200)):     maxLepJetDR,
           ('avgJetDEta',              'avgJetDEta',              (50,0.,200)):     avgJetDEta,
           ('minJetDEta',              'minJetDEta',              (50,0.,200)):     minJetDEta,
           ('maxJetDEta',              'maxJetDEta',              (50,0.,200)):     maxJetDEta,
           ('minLepJetDEta',           'minLepJetDEta',           (50,0.,200)):     minLepJetDEta,
           ('maxLepJetDEta',           'maxLepJetDEta',           (50,0.,200)):     maxLepJetDEta,
           ('avgJetDPhi',              'avgJetDPhi',              (50,0.,200)):     avgJetDPhi,
           ('minJetDPhi',              'minJetDPhi',              (50,0.,200)):     minJetDPhi,
           ('maxJetDPhi',              'maxJetDPhi',              (50,0.,200)):     maxJetDPhi,
           ('minLepJetDPhi',           'minLepJetDPhi',           (50,0.,200)):     minLepJetDPhi,
           ('maxLepJetDPhi',           'maxLepJetDPhi',           (50,0.,200)):     maxLepJetDPhi,
           ('HT2LepJetMet',            'HT2_lepJetMet',           (50,0.,200)):     HT2_lepJetMet,
           ('HT2RLepJetMet',           'HT2R_lepJetMet',          (50,0.,200)):     HT2R_lepJetMet,
           ('wj1wj2Pt',                'wj1wj2_pt',               (50,0.,200)):     wj1wj2_pt,
           ('wj1wj2M',                 'wj1wj2_M',                (50,0.,200)):     wj1wj2_M,
           ('w1w2MT',                  'w1w2_MT',                 (50,0.,200)):     w1w2_MT,
           ('recocat',                 'recocat',                 (8,0.,8)):        recocat 
    }

# ----------------------------------- Resolved LBN Inputs ------------------------------------- #
def returnLBNInputs_Resolved(self, lepton, jet1=None, jet2=None, jet3=None, jet4=None, category=None) :
    isFullReco = any([category=='Res2b2Wj', category=='Res1b3Wj'])

    lep_Px = lepton.p4.Px()
    lep_Py = lepton.p4.Py()
    lep_Pz = lepton.p4.Pz()
    lep_E  = lepton.p4.E()
    bj1_Px = self.HLL.bJetCorrP4(jet1).Px()
    bj1_Py = self.HLL.bJetCorrP4(jet1).Py()
    bj1_Pz = self.HLL.bJetCorrP4(jet1).Pz()
    bj1_E  = self.HLL.bJetCorrP4(jet1).E() 
    bj2_Px = self.HLL.bJetCorrP4(jet2).Px()
    bj2_Py = self.HLL.bJetCorrP4(jet2).Py()
    bj2_Pz = self.HLL.bJetCorrP4(jet2).Pz()
    bj2_E  = self.HLL.bJetCorrP4(jet2).E() 
    wj1_Px = jet3.p4.Px()
    wj1_Py = jet3.p4.Py()
    wj1_Pz = jet3.p4.Pz()
    wj1_E  = jet3.p4.E() 
    wj2_Px = jet4.p4.Px() if isFullReco else op.c_float(0.)
    wj2_Py = jet4.p4.Py() if isFullReco else op.c_float(0.)
    wj2_Pz = jet4.p4.Pz() if isFullReco else op.c_float(0.)
    wj2_E  = jet4.p4.E()  if isFullReco else op.c_float(0.)
    return {('lepE',   'lepE',   (50,0.,500.)):    lep_E,
            ('lepPx',  'lepPx',  (50,-250.,250.)): lep_Px,
            ('lepPy',  'lepPy',  (50,-250.,250.)): lep_Py,
            ('lepPz',  'lepPz',  (50,-250.,250.)): lep_Pz,
            ('bj1E',   'bj1E',   (50,0.,500.)):    bj1_E,
            ('bj1Px',  'bj1Px',  (50,-250.,250.)): bj1_Px,
            ('bj1Py',  'bj1Py',  (50,-250.,250.)): bj1_Py,
            ('bj1Pz',  'bj1Pz',  (50,-250.,250.)): bj1_Pz,
            ('bj2E',   'bj2E',   (50,0.,500.)):    bj2_E,
            ('bj2Px',  'bj2Px',  (50,-250.,250.)): bj2_Px,
            ('bj2Py',  'bj2Py',  (50,-250.,250.)): bj2_Py,
            ('bj2Pz',  'bj2Pz',  (50,-250.,250.)): bj2_Pz,
            ('wj1E',   'wj1E',   (50,0.,500.)):    wj1_E,
            ('wj1Px',  'wj1Px',  (50,-250.,250.)): wj1_Px,
            ('wj1Py',  'wj1Py',  (50,-250.,250.)): wj1_Py,
            ('wj1Pz',  'wj1Pz',  (50,-250.,250.)): wj1_Pz,
            ('wj2E',   'wj2E',   (50,0.,500.)):    wj2_E,
            ('wj2Px',  'wj2Px',  (50,-250.,250.)): wj2_Px,
            ('wj2Py',  'wj2Py',  (50,-250.,250.)): wj2_Py,
            ('wj2Pz',  'wj2Pz',  (50,-250.,250.)): wj2_Pz
    }
    

# ----------------------------------- Boosted Classic Inputs ------------------------------------- #
def returnClassicInputs_Boosted(self, lepton, jet3=None, jet4=None, category=None):
    isFullReco = False
    isMissReco = False

    fjet = self.ak8BJets[0]

    if category == 'Hbb2Wj':
        isFullReco = True
        VBFJetPairs = findVBFJets_Boosted(self, 2, fjet, jet3, jet4)

    elif category == 'Hbb1Wj':
        isMissReco = True
        VBFJetPairs = findVBFJets_Boosted(self, 1, fjet, jet3, jet4)

    else:
        raise RuntimeError('Boosted JPA category is not mentioned!')
    
    fatbj_E  = fjet.p4.E()
    fatbj_Px = fjet.p4.Px()
    fatbj_Py = fjet.p4.Py()
    fatbj_Pz = fjet.p4.Pz()
    wj1_Px   = jet3.p4.Px()
    wj1_Py   = jet3.p4.Py()
    wj1_Pz   = jet3.p4.Pz()
    wj1_E    = jet3.p4.E() 
    wj2_Px   = jet4.p4.Px() if isFullReco else op.c_float(0.)
    wj2_Py   = jet4.p4.Py() if isFullReco else op.c_float(0.)
    wj2_Pz   = jet4.p4.Pz() if isFullReco else op.c_float(0.)
    wj2_E    = jet4.p4.E()  if isFullReco else op.c_float(0.)
    
    # FatJet Variables
    fatbj_pt             = fjet.pt
    fatbj_eta            = fjet.eta
    fatbj_phi            = fjet.phi
    fatbj_sub1pt         = op.switch(fjet.subJet1.pt > fjet.subJet2.pt, fjet.subJet1.pt, fjet.subJet2.pt)
    fatbj_sub1eta        = op.switch(fjet.subJet1.pt > fjet.subJet2.pt, fjet.subJet1.eta, fjet.subJet2.eta)
    fatbj_sub2pt         = op.switch(fjet.subJet1.pt > fjet.subJet2.pt, fjet.subJet2.pt, fjet.subJet1.pt)
    fatbj_sub2eta        = op.switch(fjet.subJet1.pt > fjet.subJet2.pt, fjet.subJet2.eta, fjet.subJet1.eta)
    fatbj_softdropMass   = fjet.msoftdrop
    fatbj_tau1           = fjet.tau1 
    fatbj_tau2           = fjet.tau2 
    fatbj_tau3           = fjet.tau3 
    fatbj_tau4           = fjet.tau4 
    fatbj_tau43          = op.static_cast("Float_t",(fjet.tau4/fjet.tau3))
    fatbj_tau32          = op.static_cast("Float_t",(fjet.tau3/fjet.tau2))
    fatbj_tau21          = op.static_cast("Float_t",(fjet.tau2/fjet.tau1))
    fatbj_btagDDBvL      = fjet.btagDDBvL
    fatbj_btagDDBvL_noMD = fjet.btagDDBvL_noMD
    fatbj_btagDDCvB      = fjet.btagDDCvB
    fatbj_btagDDCvB_noMD = fjet.btagDDCvB_noMD
    fatbj_btagDDCvL      = fjet.btagDDCvL
    fatbj_btagDDCvL_noMD = fjet.btagDDCvL_noMD
    fatbj_btagDeepB      = fjet.btagDeepB
    
    cosThetaS_Hbb     = op.extMethod("HHbbWWJPA::cosThetaS", returnType="float")(fjet.subJet1.p4, fjet.subJet2.p4)
    mT_top_3particle  = op.min(self.HLL.mT2(fjet.subJet1.p4,lepton.p4,self.corrMET.p4), self.HLL.mT2(fjet.subJet2.p4,lepton.p4,self.corrMET.p4))
    # fatjet-lepton
    fatbj_lepDR       = op.deltaR(fjet.p4, lepton.p4)
    fatbjSub1_lepDR   = op.switch(fjet.subJet1.pt > fjet.subJet2.pt, op.deltaR(fjet.subJet1.p4, lepton.p4), op.deltaR(fjet.subJet2.p4, lepton.p4))
    fatbjSub2_lepDR   = op.switch(fjet.subJet1.pt > fjet.subJet2.pt, op.deltaR(fjet.subJet2.p4, lepton.p4), op.deltaR(fjet.subJet1.p4, lepton.p4))
    fatbj_lepDPhi     = op.abs(op.deltaPhi(fjet.p4, lepton.p4))
    fatbjSub1_lepDPhi = op.switch(fjet.subJet1.pt > fjet.subJet2.pt, op.abs(op.deltaPhi(fjet.subJet1.p4, lepton.p4)), op.abs(op.deltaPhi(fjet.subJet2.p4, lepton.p4)))
    fatbjSub2_lepDPhi = op.switch(fjet.subJet1.pt > fjet.subJet2.pt, op.abs(op.deltaPhi(fjet.subJet2.p4, lepton.p4)), op.abs(op.deltaPhi(fjet.subJet1.p4, lepton.p4)))
    minSubJetLepDR    = op.min(op.deltaR(fjet.subJet1.p4, lepton.p4), op.deltaR(fjet.subJet2.p4, lepton.p4))
    fatbj_lepDEta     = op.abs(fjet.eta - lepton.eta)
    fatbj_metDPhi     = op.abs(op.deltaPhi(fjet.p4, self.corrMET.p4))    

    # Wj1 variables
    wj1_M             = jet3.p4.M() 
    wj1_pt            = jet3.pt     
    wj1_eta           = jet3.eta    
    wj1_phi           = jet3.phi    
    wj1_bTagDeepFlavB = jet3.btagDeepFlavB

    # Wj2 variables
    wj2_M             = jet4.p4.M()                  if isFullReco else op.c_float(0.)
    wj2_pt            = jet4.pt                      if isFullReco else op.c_float(0.)
    wj2_eta           = jet4.eta                     if isFullReco else op.c_float(0.)
    wj2_phi           = jet4.phi                     if isFullReco else op.c_float(0.)
    wj2_bTagDeepFlavB = jet4.btagDeepFlavB           if isFullReco else op.c_float(0.)
    
    # fatjet-Wj1
    fatbj_Wj1DR       = op.deltaR(fjet.p4, jet3.p4)                   
    fatbj_Wj1DPhi     = op.abs(op.deltaPhi(fjet.p4, jet3.p4))         
    fatbj_Wj1DEta     = op.abs(fjet.eta - jet3.eta)                   
    fatbjSub1_Wj1DR   = op.switch(fjet.subJet1.pt > fjet.subJet2.pt, 
                                  op.deltaR(fjet.subJet1.p4, jet3.p4), 
                                  op.deltaR(fjet.subJet2.p4, jet3.p4))
    fatbjSub1_Wj1DPhi = op.switch(fjet.subJet1.pt > fjet.subJet2.pt, 
                                    op.abs(op.deltaPhi(fjet.subJet1.p4, jet3.p4)), 
                                    op.abs(op.deltaPhi(fjet.subJet2.p4, jet3.p4)))
    fatbjSub2_Wj1DR   = op.switch(fjet.subJet1.pt > fjet.subJet2.pt,
                                    op.deltaR(fjet.subJet2.p4, jet3.p4),
                                    op.deltaR(fjet.subJet1.p4, jet3.p4))  
    fatbjSub2_Wj1DPhi = op.switch(fjet.subJet1.pt > fjet.subJet2.pt,
                                  op.abs(op.deltaPhi(fjet.subJet2.p4, jet3.p4)),
                                  op.abs(op.deltaPhi(fjet.subJet1.p4, jet3.p4))) 
    # Wjs-lepton
    wj1_lepDR         = op.deltaR(jet3.p4, lepton.p4)                     
    wj1_lepDEta       = op.abs(jet3.eta - lepton.eta)                     
    wj1_lepDPhi       = op.abs(op.deltaPhi(jet3.p4, lepton.p4))           
    wj1_metDPhi       = op.abs(op.deltaPhi(jet3.p4, self.corrMET.p4))     

    wj2_lepDR         = op.deltaR(jet4.p4, lepton.p4)                     if isFullReco else op.c_float(0.)
    wj2_lepDEta       = op.abs(jet4.eta - lepton.eta)                     if isFullReco else op.c_float(0.)
    wj2_lepDPhi       = op.abs(op.deltaPhi(jet4.p4, lepton.p4))           if isFullReco else op.c_float(0.)
    wj2_metDPhi       = op.abs(op.deltaPhi(jet4.p4, self.corrMET.p4))     if isFullReco else op.c_float(0.)
    
    # fatjet-Wj2
    fatbj_wj2DR       = op.deltaR(fjet.p4, jet4.p4)                       if isFullReco else op.c_float(0.) 
    fatbj_wj2DPhi     = op.abs(op.deltaPhi(fjet.p4, jet4.p4))             if isFullReco else op.c_float(0.)
    fatbj_wj2DEta     = op.abs(fjet.eta - jet4.eta)                       if isFullReco else op.c_float(0.)
    fatbjSub1_wj2DR   = op.deltaR(fjet.subJet1.p4, jet4.p4)               if isFullReco else op.c_float(0.)
    fatbjSub1_wj2DPhi = op.abs(op.deltaPhi(fjet.subJet1.p4, jet4.p4))     if isFullReco else op.c_float(0.)
    fatbjSub2_wj2DR   = op.deltaR(fjet.subJet2.p4, jet4.p4)               if isFullReco else op.c_float(0.)
    fatbjSub2_wj2DPhi = op.abs(op.deltaPhi(fjet.subJet2.p4, jet4.p4))     if isFullReco else op.c_float(0.)
    # Wj1-Wj2
    wj1wj2DR             = op.deltaR(jet3.p4, jet4.p4)                                                                       if isFullReco else op.c_float(0.)
    wj1wj2DPhi           = op.abs(op.deltaPhi(jet3.p4, jet4.p4))                                                             if isFullReco else op.c_float(0.)
    wj1wj2DEta           = op.abs(jet3.eta - jet4.eta)                                                                       if isFullReco else op.c_float(0.)
    WWplaneAngle         = self.HLL.angleBetPlanes(lepton.p4,self.HLL.neuP4(jet3.p4+jet4.p4+lepton.p4, self.corrMET),
                                                   jet3.p4, jet4.p4)                                                         if isFullReco else op.c_float(0.) 
    WWplaneAngle_withMET = self.HLL.angleWWplane(lepton.p4, self.corrMET, jet3.p4, jet4.p4)                                  if isFullReco else op.c_float(0.)
    HWplaneAngle         = self.HLL.angleBetPlanes(fjet.subJet1.p4,fjet.subJet2.p4,jet3.p4,jet4.p4)                          if isFullReco else op.c_float(0.)
    
    HWW_Mass          = op.extMethod("HHbbWWJPA::HWW_SimpleMassWithRecoNeutrino", returnType="float")(jet3.p4, jet4.p4, 
                                                                                                      lepton.p4, self.corrMET.p4)  if isFullReco else op.c_float(0.)
    HWW_Simple_Mass   = self.HLL.HWW_met_simple(jet3.p4,jet4.p4,lepton.p4,self.corrMET.p4).M()                               if isFullReco else op.c_float(0.)
    HWW_dR            = self.HLL.dR_Hww(jet3.p4,jet4.p4,lepton.p4,self.corrMET)                                              if isFullReco else op.c_float(0.)
    HWW_dEta          = self.HLL.dEta_Hww(jet3.p4,jet4.p4,lepton.p4,self.corrMET)                                            if isFullReco else op.c_float(0.)
    HWW_dPhi          = self.HLL.dPhi_Hww(jet3.p4,jet4.p4,lepton.p4,self.corrMET)                                            if isFullReco else op.c_float(0.)

    WjjLepDR          = op.deltaR((jet3.p4 + jet4.p4), lepton.p4)                  if isFullReco else op.c_float(0.)
    WjjLepDEta        = op.abs((jet3.p4 + jet4.p4).Eta() - lepton.eta)             if isFullReco else op.c_float(0.)
    WjjLepDPhi        = op.abs(op.deltaPhi((jet3.p4 + jet4.p4), lepton.p4))        if isFullReco else op.c_float(0.)
    WjjMetDPhi        = op.abs(op.deltaPhi((jet3.p4 + jet4.p4), self.corrMET.p4))  if isFullReco else op.c_float(0.)

    cosThetaS_Wjj_simple    = op.extMethod("HHbbWWJPA::cosThetaS", returnType="float")(jet3.p4, jet4.p4)                     if isFullReco else op.c_float(0.)
    cosThetaS_WW_simple_met = op.extMethod("HHbbWWJPA::cosThetaS", returnType="float")(self.HLL.Wjj_simple(jet3.p4,jet4.p4), 
                                                                                       self.HLL.Wlep_met_simple(
                                                                                           lepton.p4, self.corrMET.p4))      if isFullReco else op.c_float(0.)
    cosThetaS_HH_simple_met = op.extMethod("HHbbWWJPA::cosThetaS", returnType="float")(fjet.subJet1.p4 + fjet.subJet2.p4, 
                                                                                       self.HLL.HWW_met_simple(jet3.p4,jet4.p4,
                                                                                                               lepton.p4,
                                                                                                               self.corrMET.p4)) if isFullReco else op.c_float(0.)
    
    
    VBFj1Px        = op.switch(op.rng_len(VBFJetPairs) > 0, 
                                 op.switch(VBFJetPairs[0][0].pt > VBFJetPairs[0][1].pt,
                                           VBFJetPairs[0][0].p4.Px(), VBFJetPairs[0][1].p4.Px()), op.c_float(0.)) 
    VBFj2Px        = op.switch(op.rng_len(VBFJetPairs) > 0, 
                                 op.switch(VBFJetPairs[0][0].pt > VBFJetPairs[0][1].pt,
                                           VBFJetPairs[0][1].p4.Px(), VBFJetPairs[0][0].p4.Px()), op.c_float(0.)) 
    VBFj1Py        = op.switch(op.rng_len(VBFJetPairs) > 0, 
                                 op.switch(VBFJetPairs[0][0].pt > VBFJetPairs[0][1].pt,
                                           VBFJetPairs[0][0].p4.Py(), VBFJetPairs[0][1].p4.Py()), op.c_float(0.)) 
    VBFj2Py        = op.switch(op.rng_len(VBFJetPairs) > 0, 
                                 op.switch(VBFJetPairs[0][0].pt > VBFJetPairs[0][1].pt,
                                           VBFJetPairs[0][1].p4.Py(), VBFJetPairs[0][0].p4.Py()), op.c_float(0.)) 
    VBFj1Pz        = op.switch(op.rng_len(VBFJetPairs) > 0, 
                                 op.switch(VBFJetPairs[0][0].pt > VBFJetPairs[0][1].pt,
                                           VBFJetPairs[0][0].p4.Pz(), VBFJetPairs[0][1].p4.Pz()), op.c_float(0.)) 
    VBFj2Pz        = op.switch(op.rng_len(VBFJetPairs) > 0, 
                                 op.switch(VBFJetPairs[0][0].pt > VBFJetPairs[0][1].pt,
                                           VBFJetPairs[0][1].p4.Pz(), VBFJetPairs[0][0].p4.Pz()), op.c_float(0.)) 
    VBFj1E         = op.switch(op.rng_len(VBFJetPairs) > 0, 
                                 op.switch(VBFJetPairs[0][0].pt > VBFJetPairs[0][1].pt,
                                           VBFJetPairs[0][0].p4.E(), VBFJetPairs[0][1].p4.E()), op.c_float(0.))   
    VBFj2E         = op.switch(op.rng_len(VBFJetPairs) > 0, 
                                 op.switch(VBFJetPairs[0][0].pt > VBFJetPairs[0][1].pt, 
                                           VBFJetPairs[0][1].p4.E(), VBFJetPairs[0][0].p4.E()), op.c_float(0.))   
    VBFj1pt        = op.switch(op.rng_len(VBFJetPairs) > 0, 
                                 op.switch(VBFJetPairs[0][0].pt > VBFJetPairs[0][1].pt,
                                           VBFJetPairs[0][0].pt, VBFJetPairs[0][1].pt), op.c_float(0.))    
    VBFj2pt        = op.switch(op.rng_len(VBFJetPairs) > 0, 
                                 op.switch(VBFJetPairs[0][0].pt > VBFJetPairs[0][1].pt,
                                           VBFJetPairs[0][1].pt, VBFJetPairs[0][0].pt), op.c_float(0.))    
    VBFj1eta       = op.switch(op.rng_len(VBFJetPairs) > 0, 
                                 op.switch(VBFJetPairs[0][0].pt > VBFJetPairs[0][1].pt,
                                           VBFJetPairs[0][0].eta, VBFJetPairs[0][1].eta), op.c_float(0.))  
    VBFj2eta       = op.switch(op.rng_len(VBFJetPairs) > 0, 
                                 op.switch(VBFJetPairs[0][0].pt > VBFJetPairs[0][1].pt,
                                           VBFJetPairs[0][1].eta, VBFJetPairs[0][0].eta), op.c_float(0.))  
    
    VBFj1j2dEta    = op.switch(op.rng_len(VBFJetPairs) > 0,op.abs(VBFJetPairs[0][0].eta - VBFJetPairs[0][1].eta), op.c_float(0.))
    VBFj1j2dPhi    = op.switch(op.rng_len(VBFJetPairs) > 0,op.abs(op.deltaPhi(VBFJetPairs[0][0].p4,VBFJetPairs[0][1].p4)), op.c_float(0.))
    VBFj1j2invM    = op.switch(op.rng_len(VBFJetPairs) > 0,op.invariant_mass(VBFJetPairs[0][0].p4,  VBFJetPairs[0][1].p4), op.c_float(0.))
    VBF_tag        = op.c_int(op.rng_len(VBFJetPairs)>0)
    zeppenfeldVar  = op.switch(op.rng_len(VBFJetPairs) > 0, op.abs(lepton.eta - (VBFJetPairs[0][0].eta + VBFJetPairs[0][1].eta)/2.0)/op.abs(VBFJetPairs[0][0].eta - VBFJetPairs[0][1].eta),op.c_float(0.))
    
    if isFullReco:
        wj1wj2_pt   = (jet3.p4+jet4.p4).Pt()
        wj1wj2invM  = op.invariant_mass(jet3.p4,jet4.p4)    
        #jetMinDR    = self.HLL.MinDiJetDRTight(fjet.subJet1,fjet.subJet2,jet3,jet4)
        jetMinDR    = op.min(op.deltaR(jet3.p4, jet4.p4), op.min(op.deltaR(fjet.p4, jet3.p4), op.deltaR(fjet.p4, jet4.p4)))
        jetMaxDR    = op.max(op.deltaR(jet3.p4, jet4.p4), op.max(op.deltaR(fjet.p4, jet3.p4), op.deltaR(fjet.p4, jet4.p4)))
        jetLepMinDR = self.HLL.MinDR_lep4j(lepton, fjet.subJet1, fjet.subJet2, jet3, jet4) 
        jetLepMaxDR = self.HLL.MaxDR_lep4j(lepton, fjet.subJet1, fjet.subJet2, jet3, jet4) 
        HT2         = self.HLL.HT2_2b2Wj(lepton, fjet.subJet1, fjet.subJet2,jet3,jet4,self.corrMET)
        HT2R        = self.HLL.HT2R_2b2Wj(lepton, fjet.subJet1, fjet.subJet2, jet3, jet4,self.corrMET)
        MT_W1W2     = self.HLL.MT_W1W2_ljj(lepton,jet3,jet4,self.corrMET)
        HT_taggedJets = fjet.pt+jet3.pt+jet4.pt
        M_taggedJetsAvg = (fjet.msoftdrop+jet3.mass+jet4.mass)/3.0
        E_taggedJetsAvg = (fjet.p4.E()+jet3.p4.E()+jet4.p4.E())/3.0
        recocat      = op.c_int(1)                
        
    if isMissReco:
        wj1wj2_pt   = jet3.p4.Pt()
        wj1wj2invM  = jet3.p4.M()    
        #jetMinDR    = self.HLL.MinDiJetDRLoose(fjet.subJet1,fjet.subJet2,jet3)
        jetMinDR    = op.deltaR(fjet.p4, jet3.p4)
        jetMaxDR    = op.deltaR(fjet.p4, jet3.p4)
        jetLepMinDR = self.HLL.MinDR_lep3j(lepton,fjet.subJet1,fjet.subJet2,jet3) 
        jetLepMaxDR = self.HLL.MaxDR_lep3j(lepton,fjet.subJet1,fjet.subJet2,jet3) 
        HT2         = self.HLL.HT2_2b1Wj(lepton,fjet.subJet1,fjet.subJet2,jet3,self.corrMET)
        HT2R        = self.HLL.HT2R_2b1Wj(lepton,fjet.subJet1,fjet.subJet2,jet3,self.corrMET)
        MT_W1W2     = self.HLL.MT_W1W2_lj(lepton,jet3,self.corrMET)
        HT_taggedJets = fjet.pt+jet3.pt
        M_taggedJetsAvg = (fjet.msoftdrop+jet3.mass)/2.0
        E_taggedJetsAvg = (fjet.p4.E()+jet3.p4.E())/2.0
        recocat      = op.c_int(2)                
        
    return{('leppdgId',                'leppdgId',                (45,-22.,22.)):   lepton.pdgId,
           ('lepcharge',               'lepcharge',               (3,-1,1.)):       lepton.charge,
           ('lepE',                    'lep_E',                   (50,0.,500.)):    lepton.p4.E(),
           ('lepPx',                   'lep_Px',                  (50,-250.,250.)): lepton.p4.Px(),
           ('lepPy',                   'lep_Py',                  (50,-250.,250.)): lepton.p4.Py(),
           ('lepPz',                   'lep_Pz',                  (50,-250.,250.)): lepton.p4.Pz(),
           ('METpx',                   'METpx',                   (50, 0, 300)):    self.corrMET.p4.Px(),
           ('METpy',                   'METpy',                   (50, 0, 300)):    self.corrMET.p4.Py(),
           ('METpz',                   'METpz',                   (50, 0, 300)):    self.corrMET.p4.Pz(),
           ('METenergy',               'METenergy',               (50, 0, 300)):    self.corrMET.p4.E(),
           ('METpt',                   'METpt',                   (50, 0, 300)):    self.corrMET.pt,
           ('METphi',                  'METphi',                  (50, 0, 300)):    self.corrMET.phi,
           ('METLD',                   'MET_LD',                  (40,0,200)):      self.HLL.MET_LD_DL(self.corrMET, self.ak4Jets, self.electronsFakeSel, 
                                                                                                       self.muonsFakeSel),
           ('nfakeableMuons',          'nfakeableMuons',          (5,0,5)):         op.static_cast("UInt_t",op.rng_len(self.muonsFakeSel)),
           ('nfakeableElectrons',      'nfakeableElectrons',      (5,0,5)):         op.static_cast("UInt_t",op.rng_len(self.electronsFakeSel)),
           ('hT',                      'hT',                      (50,0,300)):      op.static_cast("Float_t", self.HLL.HT_SL(self.ak4Jets)),
           ('hTwithCleanAk4',          'hT_withCleanAk4',         (50,0,300)):      op.static_cast("Float_t", self.HLL.HT_SL(self.ak4JetsCleanedFromAk8b)),
           ('hTwithAk8',               'hT_withAk8',              (50,0,300)):      op.static_cast("Float_t", (self.ak8BJets[0].pt+
                                                                                                               self.HLL.HT_SL(self.ak4JetsCleanedFromAk8b))),
           ('nAk4Jets',                'nAk4Jets',                (10,0,10)):       op.static_cast("UInt_t",op.rng_len(self.ak4Jets)),
           ('nAk4BJets',               'nAk4BJets',               (10,0,10)):       op.static_cast("UInt_t",op.rng_len(self.ak4BJets)),
           ('nAk4BJetsLoose',          'nAk4BJetsLoose',          (10,0,10)):       op.static_cast("UInt_t",op.rng_len(self.ak4BJetsLoose)),
           ('nAk8Jets',                'nAk8Jets',                (5,0,5)):         op.static_cast("UInt_t",op.rng_len(self.ak8Jets)),
           ('nAk8BJets',               'nAk8BJets',               (5,0,5)):         op.static_cast("UInt_t",op.rng_len(self.ak8BJets)),
           ('nAk4JetsCleandWithAk8b',  'nAk4JetsCleandWithAk8b',  (5,0,5)):         op.static_cast("UInt_t",op.rng_len(self.ak4JetsCleanedFromAk8b)),
           ('fatbjE',                  'fatbj_E',                 (50,0.,500.)):    fjet.p4.E(),
           ('fatbjPx',                 'bj1_Px',                  (50,-250.,250.)): fjet.p4.Px(),
           ('fatbjPy',                 'bj1_Py',                  (50,-250.,250.)): fjet.p4.Py(),
           ('fatbjPz',                 'bj1_Pz',                  (50,-250.,250.)): fjet.p4.Pz(),
           ('wj1E',                    'wj1_E',                   (50,0.,500.)):    wj1_E,
           ('wj1Px',                   'wj1_Px',                  (50,-250.,250.)): wj1_Px,
           ('wj1Py',                   'wj1_Py',                  (50,-250.,250.)): wj1_Py,
           ('wj1Pz',                   'wj1_Pz',                  (50,-250.,250.)): wj1_Pz,
           ('wj2E',                    'wj2_E',                   (50,0.,500.)):    wj2_E,
           ('wj2Px',                   'wj2_Px',                  (50,-250.,250.)): wj2_Px,
           ('wj2Py',                   'wj2_Py',                  (50,-250.,250.)): wj2_Py,
           ('wj2Pz',                   'wj2_Pz',                  (50,-250.,250.)): wj2_Pz,
           ('lepmetDPhi',              'lepmet_DPhi',             (20,0,3.2)):      op.abs(self.HLL.SinglepMet_dPhi(lepton, self.corrMET)),
           ('lepmetPt',                'lepmet_pt',               (50,0,300)):      self.HLL.SinglepMet_Pt(lepton, self.corrMET),
           ('lepMT',                   'lep_MT',                  (40,0,200)):      self.HLL.MT(lepton, self.corrMET),
           ('lepConept',               'lep_conept',              (40,0,200)):      self.HLL.lambdaConePt(lepton),
           ('minDRLepAllJets',         'minDR_lep_allJets',       (5,0,5)):         self.HLL.MinDR_part1_partCont(lepton, self.ak4JetsCleanedFromAk8b),
           ('minDEtaLepAllJets',       'minDEta_lep_allJets',     (5,0,5)):         self.HLL.MinDEta_part1_partCont(lepton, self.ak4JetsCleanedFromAk8b),
           ('fatbjPt',                 'fatbj_pt',                (50,0.,500.)):    fatbj_pt,
           ('fatbjEta',                'fatbj_eta',               (50,0.,500.)):    fatbj_eta,
           ('fatbjPhi',                'fatbj_phi',               (50,0.,500.)):    fatbj_phi,
           ('fatbjSub1pt',             'fatbj_sub1pt',            (50,0.,500.)):    fatbj_sub1pt,
           ('fatbjSub1eta',            'fatbj_sub1eta',           (50,0.,500.)):    fatbj_sub1eta,
           ('fatbjSub2pt',             'fatbj_sub2pt',            (50,0.,500.)):    fatbj_sub2pt,
           ('fatbjSub2eta',            'fatbj_sub2eta',           (50,0.,500.)):    fatbj_sub2eta,
           ('fatbjSoftdropMass',       'fatbj_softdropMass',      (50,0.,500.)):    fatbj_softdropMass,
           ('fatbjTau1',               'fatbj_tau1',              (50,0.,500.)):    fatbj_tau1,
           ('fatbjTau2',               'fatbj_tau2',              (50,0.,500.)):    fatbj_tau2,
           ('fatbjTau3',               'fatbj_tau3',              (50,0.,500.)):    fatbj_tau3,
           ('fatbjTau4',               'fatbj_tau4',              (50,0.,500.)):    fatbj_tau4,
           ('fatbjTau43',              'fatbj_tau43',             (50,0.,500.)):    fatbj_tau43,
           ('fatbjTau32',              'fatbj_tau32',             (50,0.,500.)):    fatbj_tau32,
           ('fatbjTau21',              'fatbj_tau21',             (50,0.,500.)):    fatbj_tau21,
           ('fatbjBtagDDBvL',          'fatbj_btagDDBvL',         (50,0.,500.)):    fatbj_btagDDBvL,
           ('fatbjBtagDDBvL_noMD',     'fatbj_btagDDBvL_noMD',    (50,0.,500.)):    fatbj_btagDDBvL_noMD,
           ('fatbjBtagDDCvB',          'fatbj_btagDDCvB',         (50,0.,500.)):    fatbj_btagDDCvB,
           ('fatbjBtagDDCvB_noMD',     'fatbj_btagDDCvB_noMD',    (50,0.,500.)):    fatbj_btagDDCvB_noMD,
           ('fatbjBtagDDCvL',          'fatbj_btagDDCvL',         (50,0.,500.)):    fatbj_btagDDCvL,
           ('fatbjBtagDDCvL_noMD',     'fatbj_btagDDCvL_noMD',    (50,0.,500.)):    fatbj_btagDDCvL_noMD,
           ('fatbjBtagDeepB',          'fatbj_btagDeepB',         (50,0.,500.)):    fatbj_btagDeepB,
           ('fatbjLepDR',              'fatbj_lepDR',             (50,0.,500.)):    fatbj_lepDR,                     
           ('fatbjSub1LepDR',          'fatbjSub1_lepDR',         (50,0.,500.)):    fatbjSub1_lepDR,                     
           ('fatbjSub2LepDR',          'fatbjSub2_lepDR',         (50,0.,500.)):    fatbjSub2_lepDR,                     
           ('fatbjLepDPhi',            'fatbj_lepDPhi',           (50,0.,500.)):    fatbj_lepDPhi,
           ('fatbjLepDEta',            'fatbj_lepDEta',           (50,0.,500.)):    fatbj_lepDEta,
           ('fatbjMetDPhi',            'fatbj_metDPhi',           (50,0.,500.)):    fatbj_metDPhi,
           ('minSubJetLepDR',          'minSubJetLepDR',          (50,0.,500.)):    minSubJetLepDR,
           ('cosThetaSHbb',            'cosThetaS_Hbb',           (50,0.,500.)):    cosThetaS_Hbb,
           ('mTtop3particle',          'mT_top_3particle',        (50,0.,500.)):    mT_top_3particle,
           ('wj1M',                    'wj1M',                    (50,0.,500.)):    wj1_M,
           ('wj1Pt',                   'wj1_pt',                  (50,0.,500.)):    wj1_pt,           
           ('wj1Eta',                  'wj1_eta',                 (50,0.,500.)):    wj1_eta,           
           ('wj1Phi',                  'wj1_phi',                 (50,0.,500.)):    wj1_phi,
           ('wj1BTagDeepFlavB',        'wj1_bTagDeepFlavB',       (50,0.,500.)):    wj1_bTagDeepFlavB,
           ('wj1LepDR',                'wj1_lepDR',               (50,0.,500.)):    wj1_lepDR,
           ('wj1LepDEta',              'wj1_lepDEta',             (50,0.,500.)):    wj1_lepDEta,
           ('wj1LepDPhi',              'wj1_lepDPhi',             (50,0.,500.)):    wj1_lepDPhi,
           ('wj1MetDPhi',              'wj1_MetDPhi',             (50,0.,500.)):    wj1_metDPhi,
           ('wj2M',                    'wj2M',                    (50,0.,500.)):    wj2_M,
           ('wj2Pt',                   'wj2_pt',                  (50,0.,500.)):    wj2_pt,           
           ('wj2Eta',                  'wj2_eta',                 (50,0.,500.)):    wj2_eta,           
           ('wj2Phi',                  'wj2_phi',                 (50,0.,500.)):    wj2_phi,
           ('wj2bTagDeepFlavB',        'wj2_bTagDeepFlavB',       (50,0.,500.)):    wj2_bTagDeepFlavB,
           ('wj2LepDR',                'wj2_lepDR',               (50,0.,500.)):    wj2_lepDR,
           ('wj2LepDEta',              'wj2_lepDEta',             (50,0.,500.)):    wj2_lepDEta,
           ('wj2LepDPhi',              'wj2_lepDPhi',             (50,0.,500.)):    wj2_lepDPhi,
           ('wj2MetDPhi',              'wj2_MetDPhi',             (50,0.,500.)):    wj2_metDPhi,
           ('fatbjWj1DR',              'fatbj_Wj1DR',             (50,0.,500.)):    fatbj_Wj1DR,
           ('fatbjWj1DPhi',            'fatbj_Wj1DPhi',           (50,0.,500.)):    fatbj_Wj1DPhi,
           ('fatbjWj1DEta',            'fatbj_Wj1DEta',           (50,0.,500.)):    fatbj_Wj1DEta,
           ('fatbjSub1Wj1DR',          'fatbjSub1_Wj1DR',         (50,0.,500.)):    fatbjSub1_Wj1DR,
           ('fatbjSub1Wj1DPhi',        'fatbjSub1_Wj1DPhi',       (50,0.,500.)):    fatbjSub1_Wj1DPhi,           
           ('fatbjSub2Wj1DR',          'fatbjSub2_Wj1DR',         (50,0.,500.)):    fatbjSub2_Wj1DR,
           ('fatbjSub2Wj1DPhi',        'fatbjSub2_Wj1DPhi',       (50,0.,500.)):    fatbjSub2_Wj1DPhi,
           ('fatbjWj2DR',              'fatbj_wj2DR',             (50,0.,500.)):    fatbj_wj2DR,
           ('fatbjWj2DPhi',            'fatbj_wj2DPhi',           (50,0.,500.)):    fatbj_wj2DPhi,
           ('fatbjWj2DEta',            'fatbj_wj2DEta',           (50,0.,500.)):    fatbj_wj2DEta,
           ('fatbjSub1Wj2DR',          'fatbjSub1_wj2DR',         (50,0.,500.)):    fatbjSub1_wj2DR,
           ('fatbjSub1Wj2DPhi',        'fatbjSub1_wj2DPhi',       (50,0.,500.)):    fatbjSub1_wj2DPhi,
           ('fatbjSub2Wj2DR',          'fatbjSub2_wj2DR',         (50,0.,500.)):    fatbjSub2_wj2DR,
           ('fatbjSub2Wj2DPhi',        'fatbjSub2_wj2DPhi',       (50,0.,500.)):    fatbjSub2_wj2DPhi,
           ('wj1wj2Pt',                'wj1wj2_pt',               (50,0.,500.)):    wj1wj2_pt,
           ('wj1wj2DR',                'wj1wj2DR',                (50,0.,500.)):    wj1wj2DR,
           ('wj1wj2DPhi',              'wj1wj2DPhi',              (50,0.,500.)):    wj1wj2DPhi,
           ('wj1wj2DEta',              'wj1wj2DEta',              (50,0.,500.)):    wj1wj2DEta,
           ('wj1wj2invM',              'wj1wj2invM',              (50,0.,500.)):    wj1wj2invM,
           ('WWplaneAngle',            'WWplaneAngle',            (50,0.,500.)):    WWplaneAngle,
           ('WWplaneAngleWithMET',     'WWplaneAngle_withMET',    (50,0.,500.)):    WWplaneAngle_withMET,
           ('HWplaneAngle',            'HWplaneAngle',            (50,0.,500.)):    HWplaneAngle,
           ('HWWMass',                 'HWW_Mass',                (50,0.,500.)):    HWW_Mass,
           ('HWWSimpleMass',           'HWW_Simple_Mass',         (50,0.,500.)):    HWW_Simple_Mass,
           ('HWWdR',                   'HWW_dR',                  (50,0.,500.)):    HWW_dR,
           ('HWWdPhi',                 'HWW_dPhi',                (50,0.,500.)):    HWW_dPhi,
           ('HWWdEta',                 'HWW_dEta',                (50,0.,500.)):    HWW_dEta,
           ('WjjLepDR',                'WjjLepDR',                (50,0,5)):        WjjLepDR,
           ('WjjLepDEta',              'WjjLepDEta',              (50,0,5)):        WjjLepDEta,
           ('WjjLepDPhi',              'WjjLepDPhi',              (50,0,5)):        WjjLepDPhi,
           ('WjjMetDPhi',              'WjjMetDPhi',              (50,0,5)):        WjjMetDPhi,
           ('cosThetaSWjjSimple',      'cosThetaS_Wjj_simple',    (50,0.,1.)):      cosThetaS_Wjj_simple,
           ('cosThetaSWWSimpleMet',    'cosThetaS_WW_simple_met', (50,0.,1.)):      cosThetaS_WW_simple_met,
           ('cosThetaSHHSimpleMet',    'cosThetaS_HH_simple_met', (50,0.,1.)):      cosThetaS_HH_simple_met,
           ('VBFj1Px',                 'VBFj1Px',                 (50,0.,200)):     VBFj1Px,           
           ('VBFj1Py',                 'VBFj1Py',                 (50,0.,200)):     VBFj1Py,           
           ('VBFj1Pz',                 'VBFj1Pz',                 (50,0.,200)):     VBFj1Pz,           
           ('VBFj1E',                  'VBFj1E',                  (50,0.,200)):     VBFj1E,           
           ('VBFj1pt',                 'VBFj1pt',                 (50,0.,200)):     VBFj1pt,           
           ('VBFj1eta',                'VBFj1eta',                (50,0.,200)):     VBFj1eta,           
           ('VBFj2Px',                 'VBFj2Px',                 (50,0.,200)):     VBFj2Px,           
           ('VBFj2Py',                 'VBFj2Py',                 (50,0.,200)):     VBFj2Py,           
           ('VBFj2Pz',                 'VBFj2Pz',                 (50,0.,200)):     VBFj2Pz,           
           ('VBFj2E',                  'VBFj2E',                  (50,0.,200)):     VBFj2E,           
           ('VBFj2pt',                 'VBFj2pt',                 (50,0.,200)):     VBFj2pt,           
           ('VBFj2eta',                'VBFj2eta',                (50,0.,200)):     VBFj2eta,
           ('VBFj1j2dEta',             'VBFj1j2dEta',             (50,0.,200)):     VBFj1j2dEta,
           ('VBFj1j2dPhi',             'VBFj1j2dPhi',             (50,0.,200)):     VBFj1j2dPhi,
           ('VBFj1j2invM',             'VBFj1j2invM',             (50,0.,200)):     VBFj1j2invM,           
           ('VBFtag',                  'VBF_tag',                 (2, 0, 2)):       VBF_tag,
           ('zeppenfeldVar',           'zeppenfeldVar',           (50,0.,200)):     zeppenfeldVar,           
           ('jetMinDR',                'jetMinDR',                (50,0.,5)):       jetMinDR,           
           ('jetMaxDR',                'jetMaxDR',                (50,0.,5)):       jetMaxDR,           
           ('jetLepMinDR',             'jetLepMinDR',             (50,0.,5)):       jetLepMinDR,           
           ('jetLepMaxDR',             'jetLepMaxDR',             (50,0.,5)):       jetLepMaxDR,           
           ('HT2',                     'HT2',                     (50,0.,300)):     HT2,           
           ('HT2R',                    'HT2R',                    (50,0.,300)):     HT2R,
           ('MTW1W2',                  'MT_W1W2',                 (50,0.,300)):     MT_W1W2,
           ('HTtaggedJets',            'HT_taggedJets',           (50,0.,300)):     HT_taggedJets,
           ('MtaggedJetsAvg',          'M_taggedJetsAvg',         (50,0.,300)):     M_taggedJetsAvg,
           ('EtaggedJetsAvg',          'E_taggedJetsAvg',         (50,0.,300)):     E_taggedJetsAvg,
           ('recocat',                 'JPAcat',                  (3,0.,3)):        recocat           
    }

# ----------------------------------- Boosted LBN Inputs ------------------------------------- #
def returnLBNInputs_Boosted(self,lepton,jet3=None,jet4=None,category=None):
    fjet = self.ak8BJets[0]

    wj1_Px = jet3.p4.Px() 
    wj1_Py = jet3.p4.Py() 
    wj1_Pz = jet3.p4.Pz() 
    wj1_E  = jet3.p4.E()  
    wj2_Px = jet4.p4.Px() if category=='Hbb2Wj' else op.c_float(0.)
    wj2_Py = jet4.p4.Py() if category=='Hbb2Wj' else op.c_float(0.)
    wj2_Pz = jet4.p4.Pz() if category=='Hbb2Wj' else op.c_float(0.)
    wj2_E  = jet4.p4.E()  if category=='Hbb2Wj' else op.c_float(0.)
    return{('lepE',   'lep_E',   (50,0.,500.)):    lepton.p4.E(),
           ('lepPx',  'lep_Px',  (50,-250.,250.)): lepton.p4.Px(),
           ('lepPy',  'lep_Py',  (50,-250.,250.)): lepton.p4.Py(),
           ('lepPz',  'lep_Pz',  (50,-250.,250.)): lepton.p4.Pz(),
           ('fatbjE', 'fatbj_E', (50,0.,500.)):    fjet.p4.E(),
           ('fatbjPx','bj1_Px',  (50,-250.,250.)): fjet.p4.Px(),
           ('fatbjPy','bj1_Py',  (50,-250.,250.)): fjet.p4.Py(),
           ('fatbjPz','bj1_Pz',  (50,-250.,250.)): fjet.p4.Pz(),
           ('wj1E',   'wj1_E',   (50,0.,500.)):    wj1_E,
           ('wj1Px',  'wj1_Px',  (50,-250.,250.)): wj1_Px,
           ('wj1Py',  'wj1_Py',  (50,-250.,250.)): wj1_Py,
           ('wj1Pz',  'wj1_Pz',  (50,-250.,250.)): wj1_Pz,
           ('wj2E',   'wj2_E',   (50,0.,500.)):    wj2_E,
           ('wj2Px',  'wj2_Px',  (50,-250.,250.)): wj2_Px,
           ('wj2Py',  'wj2_Py',  (50,-250.,250.)): wj2_Py,
           ('wj2Pz',  'wj2_Pz',  (50,-250.,250.)): wj2_Pz
    }

### --------------------------------------------------------------------------------------------------- ###

def returnEventNr(self,t):
    return {('evennr','Event number',(100,0.,1e6)): t.event}
    
def inputStaticCast(inputDict,cast='float'):
    return [op.static_cast(cast,v) for v in inputDict.values()]
