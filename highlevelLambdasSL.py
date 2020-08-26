from bamboo import treefunctions as op

# 4-Momentum association #
lep1j_p4 = lambda lep,j1 : lep.p4+j1.p4
lep2j_p4 = lambda lep,j1,j2 : lep.p4+j1.p4+j2.p4
lep3j_p4 = lambda lep,j1,j2,j3 : lep.p4+j1.p4+j2.p4+j3.p4
lep4j_p4 = lambda lep,j1,j2,j3,j4 : lep.p4+j1.p4+j2.p4+j3.p4+j4.p4

# bReg corr 4 momenta of ak4-bTagged jet #
bJetCorrP4 = lambda j : op._to.Construct("ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >", (j.pt*j.bRegCorr, j.eta, j.phi, j.mass)).result

# SingleLep-Met variables
SinglepMet_dPhi = lambda lep, met : lep.p4.Phi()-met.phi

# Jet-Met DeltaPhi


# Dilep-Met variables #
SinglepMet_Pt = lambda lep,met : op.sqrt(op.pow(met.pt*op.cos(met.phi)+lep.p4.Px(),2)+op.pow(met.pt*op.sin(met.phi)+lep.p4.Py(),2))

# Transverse mass #
MT = lambda lep,met : op.sqrt(2*lep.p4.Pt()*met.pt*(1-op.cos(lep.p4.Phi()-met.phi)))
MT_W1W2_ljj = lambda lep,j1,j2,met : op.sqrt(2*lep2j_p4(lep,j1,j2).Pt()*met.pt*(1-op.cos(lep2j_p4(lep,j1,j2).Phi()-met.phi)))
MT_W1W2_lj  = lambda lep,j1,met : op.sqrt(2*lep1j_p4(lep,j1).Pt()*met.pt*(1-op.cos(lep1j_p4(lep,j1).Phi()-met.phi)))

# dilep + dijet #
#M_lljj = lambda l1,l2,j1,j2 : op.invariant_mass(ll_p4(l1,l2),jj_p4(j1,j2))
MinDR_lep3j = lambda lep,j1,j2,j3 : op.min(op.min(op.deltaR(lep.p4,j1.p4),op.deltaR(lep.p4,j2.p4)),op.deltaR(lep.p4,j3.p4))

# Higgs related variables #
HT2_l3jmet  = lambda l,j1,j2,j3,met : op.sqrt(op.pow(met.pt*op.cos(met.phi)+l.p4.Px(),2)+op.pow(met.pt*op.sin(met.phi)+l.p4.Py(),2)) + op.abs((j1.p4+j2.p4+j3.p4).Pt())
HT2R_l3jmet = lambda l,j1,j2,j3,met : HT2_l3jmet(met,l,j1,j2,j3)/(met.pt+l.p4.Pt()+j1.p4.Pt()+j2.p4.Pt()+j3.p4.Pt())
HT2_l4jmet  = lambda l,j1,j2,j3,j4,met : op.sqrt(op.pow(met.pt*op.cos(met.phi)+l.p4.Px(),2)+op.pow(met.pt*op.sin(met.phi)+l.p4.Py(),2)) + op.abs((j1.p4+j2.p4+j3.p4+j4.p4).Pt())
HT2R_l4jmet = lambda l,j1,j2,j3,j4,met : HT2_l4jmet(met,l,j1,j2,j3,j4)/(met.pt+l.p4.Pt()+j1.p4.Pt()+j2.p4.Pt()+j3.p4.Pt()+j4.p4.Pt())


#min j1j2DR
MinDiJetDRLoose = lambda j1,j2,j3: op.min(op.min(op.deltaR(j1.p4,j2.p4), op.deltaR(j2.p4,j3.p4)), op.deltaR(j1.p4,j3.p4))





# ------------------------------------ lambdas for BDT variables ------------------------------------ #
mindr_lep1_jet = lambda lep, jets : op.deltaR(lep.p4, op.sort(jets, lambda j : op.deltaR(lep.p4, j.p4))[0].p4)
HT = lambda jets : op.rng_sum(jets, lambda j : j.p4.Pt())

# mT2
ET  = lambda lep : op.sqrt(op.pow(lep.p4.M(),2) + op.pow(lep.p4.Pt(),2))
mT2 = lambda jet, lep, met : (op.pow(jet.p4.M(),2) + op.pow(lep.p4.M(),2) + op.pow(met.p4.M(),2) + 
                              2*(ET(lep)*ET(jet) - (lep.p4.Px()*jet.p4.Px() + lep.p4.Py()*jet.p4.Py())) +
                              2*(ET(lep)*ET(met) - (lep.p4.Px()*met.p4.Px() + lep.p4.Py()*met.p4.Py())) +
                              2*(ET(jet)*ET(met) - (jet.p4.Px()*met.p4.Px() + jet.p4.Py()*met.p4.Py())))

# pZ component of met
# https://github.com/HEP-KBFI/hh-bbww/blob/f4ab60f81a920268a3f2187b97a58ec449b26883/src/comp_metP4_B2G_18_008.cc
# some necessary constants (visP4 = lepP4 + Wjj_simple)
# - - - - - used to compute neuP4 - - - - - #
_a = lambda visP4, met, mH : (op.pow(mH, 2) - op.pow(visP4.M(), 2) + 2.*visP4.Px()*met.p4.Px() + 2.*visP4.Py()*met.p4.Py())
_A = lambda visP4 : 4.0 * op.pow(visP4.E(), 2) - op.pow(visP4.Pz(), 2)
_B = lambda visP4, met, mH : -4.0 * _a(visP4, met, mH) * visP4.Pz()
_C = lambda visP4, met, mH : 4.0*op.pow(visP4.E(), 2)*(op.pow(met.p4.Px(), 2) + op.pow(met.p4.Py(), 2)) - op.pow(_a(visP4, met, mH), 2)
_D = lambda visP4, met, mH : (op.pow(_B(visP4, met, mH), 2) - 4.0*_A(visP4)*_C(visP4, met, mH))
_pos  = lambda visP4, met, mH : (-_B(visP4, met, mH) + op.sqrt(_D(visP4, met, mH)))/(2.*_A(visP4))
_neg  = lambda visP4, met, mH : (-_B(visP4, met, mH) - op.sqrt(_D(visP4, met, mH)))/(2.*_A(visP4))
neuPz = lambda visP4, met, mH : (op.switch(_D(visP4, met, mH) < 0., - _B(visP4, met, mH)/(2.*_A(visP4)), 
                                           op.switch(op.abs(_pos(visP4, met, mH)) < op.abs(_neg(visP4, met, mH)), 
                                                     _pos(visP4, met, mH), _neg(visP4, met, mH))))
# - - - - - - - - - - - - - - - - - - - - - #
neuP4 = lambda visP4, met, mH : op._to.Construct("ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", 
                                                 (met.p4.Px(), 
                                                  met.p4.Py(), 
                                                  neuPz(visP4, met, mH), 
                                                  op.sqrt(op.pow(met.p4.Px(),2)+op.pow(met.p4.Py(),2)+op.pow(neuPz(visP4, met, mH),2)))).result

# P4 of W1 (l,neu)
Wlep_simple = lambda j1P4,j2P4,lepP4,met,mH : lepP4 + neuP4(j1P4+j2P4+lepP4, met, mH)
# P4 of W2 (j,j)
Wjj_simple  = lambda j1P4,j2P4 : j1P4 + j2P4 
# P4 of HWW (W1 + W2)
HWW_simple  = lambda j1P4,j2P4,lepP4,met,mH : Wjj_simple(j1P4,j2P4) + Wlep_simple(lepP4,neuP4(j1P4+j2P4+lepP4, met, mH))
# dR_HWW
dR_Hww   = lambda j1P4,j2P4,lepP4,met,mH : op.deltaR(Wjj_simple(j1P4,j2P4), Wlep_simple(j1P4,j2P4,lepP4,met,mH))
# P4 of lep + met
Wlep_met_simple = lambda lepP4, metP4 : lepP4 + metP4
# SimpleP4 of HWW (W1 + W2)
HWW_met_simple  = lambda j1P4,j2P4,lepP4,metP4 : Wjj_simple(j1P4, j2P4) + Wlep_met_simple(lepP4,metP4)
# Total P4
HHP4_simple_met = lambda HbbRegP4, j1P4, j2P4, lepP4,metP4 : HbbRegP4 + Wjj_simple(j1P4, j2P4) + Wlep_met_simple(lepP4,metP4)

# CosThetaS calculation
#comp_cosThetaS = lambda ob1p4, ob2p4 : op.abs(ob1p4.Boost(-(ob1p4+ob2p4).BoostVector()).CosTheta())
motherPx = lambda ob1p4,ob2p4 : (ob1p4.Px()+ob2p4.Px())
motherPy = lambda ob1p4,ob2p4 : (ob1p4.Py()+ob2p4.Py())
motherPz = lambda ob1p4,ob2p4 : (ob1p4.Pz()+ob2p4.Pz())
motherE  = lambda ob1p4,ob2p4 : (ob1p4.E()+ob2p4.E())
BoostP4 = lambda ob1p4,ob2p4 : op._to.Construct("ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", (motherPx(ob1p4,ob2p4), motherPy(ob1p4,ob2p4), motherPz(ob1p4,ob2p4), motherE(ob1p4,ob2p4))).result
comp_cosThetaS = lambda ob1p4,ob2p4 : op.abs(op.cos(op.deltaR(BoostP4(ob1p4,ob2p4), ob1p4)))

# MET_LD
# Equation 3 (page 33) of AN-2019/111 v13
# Similar to MET, but more robust against pileup
jetSumPx = lambda jets : op.rng_sum(jets, lambda j : j.p4.Px())
jetSumPy = lambda jets : op.rng_sum(jets, lambda j : j.p4.Py())
lepSumPx = lambda leps : op.rng_sum(leps, lambda l : l.p4.Px())
lepSumPy = lambda leps : op.rng_sum(leps, lambda l : l.p4.Py())
MET_LD = lambda met, jets, leps : 0.6*met.pt + 0.4*op.sqrt(op.pow(jetSumPx(jets)+lepSumPx(leps),2)+op.pow(jetSumPy(jets)+lepSumPy(leps),2)) 



# Lepton conePt
# from Base
lambda_conept_electron = lambda lep : op.multiSwitch((op.AND(op.abs(lep.pdgId)!=11 , op.abs(lep.pdgId)!=13) , op.static_cast("Float_t",lep.pt)),
                                                     # if (abs(lep.pdgId)!=11 and abs(lep.pdgId)!=13): return lep.pt : anything that is not muon or electron
                                                     (op.AND(op.abs(lep.pdgId)==11 , lep.mvaTTH > 0.80) , op.static_cast("Float_t",lep.pt)),
                                                     # if electron, check above MVA 
                                                     op.static_cast("Float_t",0.9*lep.pt*(1.+lep.jetRelIso)))
                                                     # else: return 0.90 * lep.pt / jetPtRatio where jetPtRatio = 1./(Electron_jetRelIso + 1.)
# from Base
lambda_conept_muon = lambda lep : op.multiSwitch((op.AND(op.abs(lep.pdgId)!=11 , op.abs(lep.pdgId)!=13) , op.static_cast("Float_t",lep.pt)),
                                                 # if (abs(lep.pdgId)!=11 and abs(lep.pdgId)!=13): return lep.pt : anything that is not muon or electron
                                                 (op.AND(op.abs(lep.pdgId)==13 , lep.mediumId ,lep.mvaTTH > 0.85) , op.static_cast("Float_t",lep.pt)),
                                                 # if muon, check that passes medium and above MVA
                                                 op.static_cast("Float_t",0.9*lep.pt*(1.+lep.jetRelIso)))
                                                 # else: return 0.90 * lep.pt / lep.jetPtRatiov2

lepConePt = lambda lep: op.switch(op.abs(lep.pdgId) == 11, lambda_conept_electron(lep), lambda_conept_muon(lep))
