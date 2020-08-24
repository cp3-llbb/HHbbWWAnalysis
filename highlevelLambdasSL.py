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
