import os
from bamboo import treefunctions as op
from bamboo.root import loadHeader 
loadHeader(os.path.join(os.path.dirname(os.path.abspath(__file__)),"jpa.h"))

class highlevelLambdas:
    def __init__(self,HHself):
        # All the attributes of the BaseHH are contained in HHself object
        # All the lambdas will be saved in the highlevelLambdas object to avoid confusions of all the attributes of HH base object

        # conept #
        self.conept = lambda lep : op.switch(op.abs(lep.pdgId)==11,
                                             HHself.electron_conept[lep.idx],
                                             HHself.muon_conept[lep.idx])

        self.electron_conept = lambda ele : HHself.electron_conept[ele.idx]
        self.muon_conept     = lambda mu  : HHself.muon_conept[mu.idx]

        # 4-Momentum association #
        self.ll_p4 = lambda l1,l2 : l1.p4+l2.p4
        self.lljj_p4 = lambda l1,l2,j1,j2 : l1.p4+l2.p4+j1.p4+j2.p4
        self.lep1j_p4 = lambda lep,j1 : lep.p4+j1.p4
        self.lep2j_p4 = lambda lep,j1,j2 : lep.p4+j1.p4+j2.p4
        self.lep3j_p4 = lambda lep,j1,j2,j3 : lep.p4+j1.p4+j2.p4+j3.p4
        self.lep4j_p4 = lambda lep,j1,j2,j3,j4 : lep.p4+j1.p4+j2.p4+j3.p4+j4.p4
        
        # bReg corr 4 momenta of ak4-bTagged jet #
        self.bJetCorrP4 = lambda j : op._to.Construct("ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >", (j.pt*j.bRegCorr, j.eta, j.phi, j.mass*j.bRegCorr)).result
        
        # Dilep-Met variables #
        self.DilepMET_deltaPhi = lambda l1,l2,met : self.ll_p4(l1,l2).Phi()-met.phi
        self.DilepMET_Pt = lambda l1,l2,met : op.sqrt(op.pow(met.pt*op.cos(met.phi)+self.ll_p4(l1,l2).Px(),2)+op.pow(met.pt*op.sin(met.phi)+self.ll_p4(l1,l2).Py(),2))
        # SingleLep-Met variables
        #self.SinglepMet_Pt = lambda lep,met : op.sqrt(op.pow(met.pt*op.cos(met.phi)+lep.p4.Px(),2)+op.pow(met.pt*op.sin(met.phi)+lep.p4.Py(),2))
        self.SinglepMet_Pt = lambda lep,met : (lep.p4 + met.p4).Pt()
        #self.SinglepMet_dPhi = lambda lep, met : lep.p4.Phi()-met.phi
        self.SinglepMet_dPhi = lambda lep, met : op.deltaPhi(lep.p4,met.p4)
        
        # Transverse mass #
        self.MT_ll = lambda l1,l2,met : op.sqrt(2*self.ll_p4(l1,l2).Pt()*met.pt*(1-op.cos(self.ll_p4(l1,l2).Phi()-met.phi)))
        self.MT_lljj = lambda l1,l2,j1,j2,met : op.sqrt(2*self.lljj_p4(l1,l2,j1,j2).Pt()*met.pt*(1-op.cos(self.lljj_p4(l1,l2,j1,j2).Phi()-met.phi)))
        self.MT = lambda lep,met : op.sqrt(2*lep.p4.Pt()*met.pt*(1-op.cos(lep.p4.Phi()-met.phi)))
        self.MT_W1W2_ljj = lambda lep,j1,j2,met : op.sqrt(2*self.lep2j_p4(lep,j1,j2).Pt()*met.pt*(1-op.cos(self.lep2j_p4(lep,j1,j2).Phi()-met.phi)))
        self.MT_W1W2_lj  = lambda lep,j1,met : op.sqrt(2*self.lep1j_p4(lep,j1).Pt()*met.pt*(1-op.cos(self.lep1j_p4(lep,j1).Phi()-met.phi)))
        # TODO : clean different versions (eg MT)

        # dilep + dijet #
        self.M_lljj = lambda l1,l2,j1,j2 : op.invariant_mass(self.lljj_p4(l1,l2,j1,j2))
        self.M_HH = lambda l1,l2,j1,j2,met : op.invariant_mass(l1.p4,l2.p4,j1.p4,j2.p4,met.p4)
        self.MinDR_lj = lambda l1,l2,j1,j2 : op.min(op.min(op.deltaR(l1.p4,j1.p4),op.deltaR(l1.p4,j2.p4)),
                                               op.min(op.deltaR(l2.p4,j1.p4),op.deltaR(l2.p4,j2.p4)))
        self.MinDR_part1_partCont   = lambda part1,partCont : op.rng_min(partCont, lambda part2 : op.deltaR(part1.p4, part2.p4))
        self.MinDEta_part1_partCont = lambda part1,partCont : op.rng_min(partCont, lambda part2 : op.abs(part1.eta - part2.eta))
        self.MinDPhi_part1_partCont = lambda part1,partCont : op.rng_min(partCont, lambda part2 : op.abs(op.deltaPhi(part1.p4, part2.p4)))

        self.MinDR_part1_dipart = lambda part1,dipart: op.min(*(op.deltaR(part1.p4, dipart[i2].p4) for i2 in range(2)))

        self.JetsMinDR = lambda l,j1,j2 : op.min(op.deltaR(l.p4,j1.p4),op.deltaR(l.p4,j2.p4))
        self.LepsMinDR = lambda j,l1,l2 : op.min(op.deltaR(j.p4,l1.p4),op.deltaR(j.p4,l2.p4))
        
        self.MinDR_lep3j = lambda lep,j1,j2,j3 : op.min(op.min(op.deltaR(lep.p4,j1.p4),op.deltaR(lep.p4,j2.p4)),op.deltaR(lep.p4,j3.p4))
        

        self.MinDR_lj = lambda l1,l2,j1,j2 : op.min(op.min(op.deltaR(l1.p4,j1.p4),op.deltaR(l1.p4,j2.p4)),
                                               op.min(op.deltaR(l2.p4,j1.p4),op.deltaR(l2.p4,j2.p4)))
        
        self.MinDR_lep2j = lambda lep,j1,j2 : op.min(op.deltaR(lep.p4,j1.p4),op.deltaR(lep.p4,j2.p4))
        self.MinDR_lep3j = lambda lep,j1,j2,j3 : op.min(op.min(op.deltaR(lep.p4,j1.p4),op.deltaR(lep.p4,j2.p4)),op.deltaR(lep.p4,j3.p4))
        self.MinDR_lep4j = lambda lep,j1,j2,j3,j4 : op.min(op.min(op.min(op.deltaR(lep.p4,j1.p4),op.deltaR(lep.p4,j2.p4)),op.deltaR(lep.p4,j3.p4)),op.deltaR(lep.p4,j4.p4))
        self.MinDPhi_lep2j = lambda lep,j1,j2 : op.min(op.abs(op.deltaPhi(lep.p4,j1.p4)), op.abs(op.deltaPhi(lep.p4,j2.p4)))
        self.MinDPhi_lep3j = lambda lep,j1,j2,j3 : op.min(op.min(op.abs(op.deltaPhi(lep.p4,j1.p4)), op.abs(op.deltaPhi(lep.p4,j2.p4))), op.abs(op.deltaPhi(lep.p4,j3.p4)))
        self.MinDPhi_lep4j = lambda lep,j1,j2,j3,j4 : op.min(op.min(op.min(op.abs(op.deltaPhi(lep.p4,j1.p4)), op.abs(op.deltaPhi(lep.p4,j2.p4))), op.abs(op.deltaPhi(lep.p4,j3.p4))), op.abs(op.deltaPhi(lep.p4,j4.p4)))
        self.MinDEta_lep2j = lambda lep,j1,j2 : op.min(op.abs(lep.eta-j1.eta),op.abs(lep.eta-j2.eta))
        self.MinDEta_lep3j = lambda lep,j1,j2,j3 : op.min(op.min(op.abs(lep.eta-j1.eta),op.abs(lep.eta-j2.eta)),op.abs(lep.eta-j3.eta))
        self.MinDEta_lep4j = lambda lep,j1,j2,j3,j4 : op.min(op.min(op.min(op.abs(lep.eta-j1.eta),op.abs(lep.eta-j2.eta)),op.abs(lep.eta-j3.eta)),op.abs(lep.eta-j4.eta))
        
        
        self.MaxDR_lep2j = lambda lep,j1,j2 : op.max(op.deltaR(lep.p4,j1.p4),op.deltaR(lep.p4,j2.p4))
        self.MaxDR_lep3j = lambda lep,j1,j2,j3 : op.max(op.max(op.deltaR(lep.p4,j1.p4),op.deltaR(lep.p4,j2.p4)),op.deltaR(lep.p4,j3.p4))
        self.MaxDR_lep4j = lambda lep,j1,j2,j3,j4 : op.max(op.max(op.max(op.deltaR(lep.p4,j1.p4),op.deltaR(lep.p4,j2.p4)),op.deltaR(lep.p4,j3.p4)),op.deltaR(lep.p4,j4.p4))
        self.MaxDPhi_lep2j = lambda lep,j1,j2 : op.max(op.abs(op.deltaPhi(lep.p4,j1.p4)), op.abs(op.deltaPhi(lep.p4,j2.p4)))
        self.MaxDPhi_lep3j = lambda lep,j1,j2,j3 : op.max(op.max(op.abs(op.deltaPhi(lep.p4,j1.p4)), op.abs(op.deltaPhi(lep.p4,j2.p4))), op.abs(op.deltaPhi(lep.p4,j3.p4)))
        self.MaxDPhi_lep4j = lambda lep,j1,j2,j3,j4 : op.max(op.max(op.max(op.abs(op.deltaPhi(lep.p4,j1.p4)), op.abs(op.deltaPhi(lep.p4,j2.p4))), op.abs(op.deltaPhi(lep.p4,j3.p4))), op.abs(op.deltaPhi(lep.p4,j4.p4)))
        self.MaxDEta_lep2j = lambda lep,j1,j2 : op.max(op.abs(lep.eta-j1.eta),op.abs(lep.eta-j2.eta))
        self.MaxDEta_lep3j = lambda lep,j1,j2,j3 : op.max(op.max(op.abs(lep.eta-j1.eta),op.abs(lep.eta-j2.eta)),op.abs(lep.eta-j3.eta))
        self.MaxDEta_lep4j = lambda lep,j1,j2,j3,j4 : op.max(op.max(op.max(op.abs(lep.eta-j1.eta),op.abs(lep.eta-j2.eta)),op.abs(lep.eta-j3.eta)),op.abs(lep.eta-j4.eta))

        # Higgs related variables #
        #self.HT2 = lambda l1,l2,j1,j2,met : op.sqrt(op.pow(met.pt*op.cos(met.phi)+l1.p4.Px()+l2.p4.Px(),2)+op.pow(met.pt*op.sin(met.phi)+l1.p4.Py()+l2.p4.Py(),2)) + op.abs((j1.p4+j2.p4).Pt())
        #self.HT2R = lambda l1,l2,j1,j2,met : self.HT2(met,l1,l2,j1,j2)/(met.pt+l1.p4.Pt()+l2.p4.Pt()+j1.p4.Pt()+j2.p4.Pt())
        #self.HT2_l1jmet  = lambda l,j1,met : op.sqrt(op.pow(met.pt*op.cos(met.phi)+l.p4.Px(),2)+op.pow(met.pt*op.sin(met.phi)+l.p4.Py(),2)) + op.abs(j1.p4.Pt())
        #self.HT2R_l1jmet = lambda l,j1,met : self.HT2_l1jmet(l,j1,met)/(met.pt+l.p4.Pt()+j1.p4.Pt())
        #self.HT2_l2jmet  = lambda l,j1,j2,met : op.sqrt(op.pow(met.pt*op.cos(met.phi)+l.p4.Px(),2)+op.pow(met.pt*op.sin(met.phi)+l.p4.Py(),2)) + op.abs((j1.p4+j2.p4).Pt())
        #self.HT2R_l2jmet = lambda l,j1,j2,met : self.HT2_l2jmet(l,j1,j2,met)/(met.pt+l.p4.Pt()+j1.p4.Pt()+j2.p4.Pt())
        #self.HT2_l3jmet  = lambda l,j1,j2,j3,met : op.sqrt(op.pow(met.pt*op.cos(met.phi)+l.p4.Px(),2)+op.pow(met.pt*op.sin(met.phi)+l.p4.Py(),2)) + op.abs((j1.p4+j2.p4+j3.p4).Pt())

        self.HT_SL = lambda jets : op.rng_sum(jets, lambda j : j.pt)
        # 0b
        self.HT2_0b  = lambda l,met : op.abs((met.p4+l.p4).Pt())
        self.HT2R_0b = lambda l,met : self.HT2_0b(l,met)/(met.pt+l.pt)

        # 1b0Wj
        self.HT2_1b0Wj  = lambda l,j1,met : op.abs((met.p4+l.p4).Pt()) + j1.pt
        self.HT2R_1b0Wj = lambda l,j1,met : self.HT2_1b0Wj(l,j1,met)/(met.pt+l.pt+j1.pt)

        # 1b1Wj
        self.HT2_1b1Wj  = lambda l,j1,j3,met : op.abs((met.p4+l.p4+j3.p4).Pt()) + j1.pt
        self.HT2R_1b1Wj = lambda l,j1,j3,met : self.HT2_1b1Wj(l,j1,j3,met)/(met.pt+l.pt+j1.pt+j3.pt)

        #2b0Wj
        self.HT2_2b0Wj  = lambda l,j1,j2,met : op.abs((met.p4+l.p4).Pt()) + op.abs((j1.p4+j2.p4).Pt())
        self.HT2R_2b0Wj = lambda l,j1,j2,met : self.HT2_2b0Wj(l,j1,j2,met)/(met.pt+l.pt+j1.pt+j2.pt)

        # 1b2Wj
        self.HT2_1b2Wj  = lambda l,j1,j3,j4,met : op.abs((met.p4+l.p4+j3.p4+j4.p4).Pt()) + j1.pt
        self.HT2R_1b2Wj = lambda l,j1,j3,j4,met : self.HT2_1b2Wj(l,j1,j3,j4,met)/(met.pt+l.pt+j1.pt+j3.pt+j4.pt)

        # 2b1Wj
        self.HT2_2b1Wj  = lambda l,j1,j2,j3,met : op.abs((met.p4+l.p4+j3.p4).Pt()) + op.abs((j1.p4+j2.p4).Pt())
        self.HT2R_2b1Wj = lambda l,j1,j2,j3,met : self.HT2_2b1Wj(l,j1,j2,j3,met)/(met.pt+l.pt+j1.pt+j2.pt+j3.pt)

        #self.HT2_l4jmet  = lambda l,j1,j2,j3,j4,met : op.sqrt(op.pow(met.pt*op.cos(met.phi)+l.p4.Px(),2)+op.pow(met.pt*op.sin(met.phi)+l.p4.Py(),2)) + op.abs((j1.p4+j2.p4+j3.p4+j4.p4).Pt())
        # 2b2Wj
        self.HT2_2b2Wj  = lambda l,j1,j2,j3,j4,met : op.abs((met.p4+l.p4+j3.p4+j4.p4).Pt()) + op.abs((j1.p4+j2.p4).Pt())
        self.HT2R_2b2Wj = lambda l,j1,j2,j3,j4,met : self.HT2_2b2Wj(l,j1,j2,j3,j4,met)/(met.pt+l.pt+j1.pt+j2.pt+j3.pt+j4.pt)
        
        #min j1j2DR
        self.MinDiJetDRLoose   = lambda j1,j2,j3: op.min(op.min(op.deltaR(j1.p4,j2.p4), op.deltaR(j2.p4,j3.p4)), op.deltaR(j1.p4,j3.p4))
        self.MinDiJetDRTight   = lambda j1,j2,j3,j4: op.min(op.min(op.min(self.MinDiJetDRLoose(j1,j2,j3), op.deltaR(j1.p4,j4.p4)), op.deltaR(j2.p4,j4.p4)),op.deltaR(j3.p4,j4.p4))
        self.MinDiJetDEtaLoose = lambda j1,j2,j3: op.min(op.min(op.abs(j1.eta-j2.eta), op.abs(j2.eta-j3.eta)), op.abs(j1.eta-j3.eta))
        self.MinDiJetDEtaTight = lambda j1,j2,j3,j4: op.min(op.min(op.min(self.MinDiJetDEtaLoose(j1,j2,j3), op.abs(j1.eta-j4.eta)), op.abs(j2.eta-j4.eta)),op.abs(j3.eta-j4.eta))
        self.MinDiJetDPhiLoose = lambda j1,j2,j3: op.min(op.min(op.abs(op.deltaPhi(j1.p4,j2.p4)), op.abs(op.deltaPhi(j2.p4,j3.p4))), op.abs(op.deltaPhi(j1.p4,j3.p4)))
        self.MinDiJetDPhiTight = lambda j1,j2,j3,j4: op.min(op.min(op.min(self.MinDiJetDPhiLoose(j1,j2,j3), op.abs(op.deltaPhi(j1.p4,j4.p4))), op.abs(op.deltaPhi(j2.p4,j4.p4))), op.abs(op.deltaPhi(j3.p4,j4.p4)))

        self.MaxDiJetDRLoose   = lambda j1,j2,j3: op.max(op.max(op.deltaR(j1.p4,j2.p4), op.deltaR(j2.p4,j3.p4)), op.deltaR(j1.p4,j3.p4))
        self.MaxDiJetDRTight   = lambda j1,j2,j3,j4: op.max(op.max(op.max(self.MaxDiJetDRLoose(j1,j2,j3), op.deltaR(j1.p4,j4.p4)), op.deltaR(j2.p4,j4.p4)),op.deltaR(j3.p4,j4.p4))
        self.MaxDiJetDEtaLoose = lambda j1,j2,j3: op.max(op.max(op.abs(j1.eta-j2.eta), op.abs(j2.eta-j3.eta)), op.abs(j1.eta-j3.eta))
        self.MaxDiJetDEtaTight = lambda j1,j2,j3,j4: op.max(op.max(op.max(self.MaxDiJetDEtaLoose(j1,j2,j3), op.abs(j1.eta-j4.eta)), op.abs(j2.eta-j4.eta)),op.abs(j3.eta-j4.eta))
        self.MaxDiJetDPhiLoose = lambda j1,j2,j3: op.max(op.max(op.abs(op.deltaPhi(j1.p4,j2.p4)), op.abs(op.deltaPhi(j2.p4,j3.p4))), op.abs(op.deltaPhi(j1.p4,j3.p4)))
        self.MaxDiJetDPhiTight = lambda j1,j2,j3,j4: op.max(op.max(op.max(self.MaxDiJetDPhiLoose(j1,j2,j3), op.abs(op.deltaPhi(j1.p4,j4.p4))), op.abs(op.deltaPhi(j2.p4,j4.p4))), op.abs(op.deltaPhi(j3.p4,j4.p4)))

        
        # ------------------------------------ lambdas for BDT variables ------------------------------------ #
        # min jet-lep DR
        self.mindr_lep1_jet = lambda lep, jets : op.deltaR(lep.p4, op.sort(jets, lambda j : op.deltaR(lep.p4, j.p4))[0].p4)
        
        # HT
        self.HTfull = lambda fleps, j1p4,j2p4,j3p4,j4p4 : j1p4.Pt()+j2p4.Pt()+j3p4.Pt()+j4p4.Pt()+op.rng_sum(fleps, lambda l : l.p4.Pt())
        self.HTmiss = lambda fleps, j1p4,j2p4,j3p4 : j1p4.Pt()+j2p4.Pt()+j3p4.Pt()+op.rng_sum(fleps, lambda l : l.p4.Pt())
        
        # mT2
        ET  = lambda lepp4 : op.sqrt(op.pow(lepp4.M(),2) + op.pow(lepp4.Pt(),2))
        self.mT2 = lambda jetp4, lepp4, metp4 : (op.pow(jetp4.M(),2) + op.pow(lepp4.M(),2) + op.pow(metp4.M(),2) + 
                                                 2*(ET(lepp4)*ET(jetp4) - (lepp4.Px()*jetp4.Px() + lepp4.Py()*jetp4.Py())) +
                                                 2*(ET(lepp4)*ET(metp4) - (lepp4.Px()*metp4.Px() + lepp4.Py()*metp4.Py())) +
                                                 2*(ET(jetp4)*ET(metp4) - (jetp4.Px()*metp4.Px() + jetp4.Py()*metp4.Py())))
        
        # pZ component of met
        # https://github.com/HEP-KBFI/hh-bbww/blob/f4ab60f81a920268a3f2187b97a58ec449b26883/src/comp_metP4_B2G_18_008.cc
        # some necessary constants (visP4 = lepP4 + Wjj_simple)
        # - - - - - used to compute neuP4 - - - - - #
        ax = lambda visP4, met : 125.18*125.18 - op.pow(visP4.M(), 2) + 2.*visP4.Px()*met.p4.Px() + 2.*visP4.Py()*met.p4.Py()
        A  = lambda visP4 : 4.0 * op.pow(visP4.E(), 2) - op.pow(visP4.Pz(), 2)
        B  = lambda visP4, met : -4.0 * ax(visP4, met) * visP4.Pz()
        C  = lambda visP4, met : 4.0*op.pow(visP4.E(), 2)*(op.pow(met.p4.Px(), 2) + op.pow(met.p4.Py(), 2)) - op.pow(ax(visP4, met), 2)
        D  = lambda visP4, met : (op.pow(B(visP4, met), 2) - 4.0*A(visP4)*C(visP4, met))
        pos   = lambda visP4, met : (-B(visP4, met) + op.sqrt(D(visP4, met)))/(2.*A(visP4))
        neg   = lambda visP4, met : (-B(visP4, met) - op.sqrt(D(visP4, met)))/(2.*A(visP4))
        neuPz = lambda visP4, met : (op.switch(D(visP4, met) < 0., -B(visP4, met)/(2*A(visP4)), op.switch(op.abs(pos(visP4, met)) < op.abs(neg(visP4, met)), pos(visP4, met), neg(visP4, met))))
        # - - - - - - - - - - - - - - - - - - - - - #
        self.neuP4 = lambda visP4, met : op._to.Construct("ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", 
                                                          (met.p4.Px(), 
                                                           met.p4.Py(), 
                                                           neuPz(visP4, met), 
                                                           op.sqrt(op.pow(met.p4.Px(),2)+op.pow(met.p4.Py(),2)+op.pow(neuPz(visP4, met),2)))).result
        
        # P4 of W1 (l,neu)
        self.Wlep_simple = lambda wj1P4,wj2P4,lepP4,met : lepP4 + self.neuP4(wj1P4+wj2P4+lepP4, met)
        # P4 of W2 (j,j)
        self.Wjj_simple  = lambda j1P4,j2P4 : j1P4 + j2P4 
        # DR_HadW_bJet
        self.dR_HadW_bjet = lambda bP4,j1P4,j2P4 : op.deltaR(self.Wjj_simple(j1P4,j2P4), bP4)
        # P4 of HWW (W1 + W2)
        self.HWW_simple  = lambda wj1P4,wj2P4,lepP4,met : self.Wjj_simple(wj1P4,wj2P4) + self.Wlep_simple(wj1P4,wj2P4,lepP4,met)
        # dR_HWW
        self.dR_Hww   = lambda j1P4,j2P4,lepP4,met : op.deltaR(self.Wjj_simple(j1P4,j2P4), self.Wlep_simple(j1P4,j2P4,lepP4,met))
        self.dEta_Hww = lambda j1P4,j2P4,lepP4,met : op.abs(self.Wjj_simple(j1P4,j2P4).Eta() - self.Wlep_simple(j1P4,j2P4,lepP4,met).Eta())
        self.dPhi_Hww = lambda j1P4,j2P4,lepP4,met : op.abs(op.deltaPhi(self.Wjj_simple(j1P4,j2P4), self.Wlep_simple(j1P4,j2P4,lepP4,met)))
        # P4 of lep + met
        self.Wlep_met_simple = lambda lepP4, metP4 : lepP4 + metP4
        # SimpleP4 of HWW (W1 + W2)
        self.HWW_met_simple  = lambda j1P4,j2P4,lepP4,metP4 : self.Wjj_simple(j1P4, j2P4) + self.Wlep_met_simple(lepP4,metP4)
        # Total P4
        self.HHP4_simple_met = lambda HbbRegP4,j1P4,j2P4,lepP4,metP4 : HbbRegP4 + self.Wjj_simple(j1P4, j2P4) + self.Wlep_met_simple(lepP4,metP4)
        
        # CosThetaS calculation
        #comp_cosThetaS = lambda ob1p4, ob2p4 : op.abs(ob1p4.Boost(-(ob1p4+ob2p4).BoostVector()).CosTheta())
        motherPx = lambda ob1p4,ob2p4 : (ob1p4.Px()+ob2p4.Px())
        motherPy = lambda ob1p4,ob2p4 : (ob1p4.Py()+ob2p4.Py())
        motherPz = lambda ob1p4,ob2p4 : (ob1p4.Pz()+ob2p4.Pz())
        motherE  = lambda ob1p4,ob2p4 : (ob1p4.E()+ob2p4.E())
        betaX    = lambda ob1p4,ob2p4 : motherPx(ob1p4,ob2p4)/motherE(ob1p4,ob2p4)
        betaY    = lambda ob1p4,ob2p4 : motherPy(ob1p4,ob2p4)/motherE(ob1p4,ob2p4)
        betaZ    = lambda ob1p4,ob2p4 : motherPz(ob1p4,ob2p4)/motherE(ob1p4,ob2p4)
        beta2    = lambda ob1p4,ob2p4 : op.pow(betaX(ob1p4,ob2p4),2)+op.pow(betaY(ob1p4,ob2p4),2)+op.pow(betaZ(ob1p4,ob2p4),2)
        gamma    = lambda ob1p4,ob2p4 : 1.0/op.sqrt(1-beta2(ob1p4,ob2p4))
        betap    = lambda ob1p4,ob2p4 : betaX(ob1p4,ob2p4)*motherPx(ob1p4,ob2p4) + betaY(ob1p4,ob2p4)*motherPy(ob1p4,ob2p4) + betaZ(ob1p4,ob2p4)*motherPz(ob1p4,ob2p4)
        gamma2   = lambda ob1p4,ob2p4 : op.switch(beta2(ob1p4,ob2p4) > 0, (gamma(ob1p4,ob2p4)-1)/beta2(ob1p4,ob2p4) , op.c_float(0.0))
        boostPx  = lambda ob1p4,ob2p4 : ob1p4.Px() + gamma2(ob1p4,ob2p4)*betap(ob1p4,ob2p4)*betaX(ob1p4,ob2p4) + gamma(ob1p4,ob2p4)*betaX(ob1p4,ob2p4)*ob1p4.E()
        boostPy  = lambda ob1p4,ob2p4 : ob1p4.Px() + gamma2(ob1p4,ob2p4)*betap(ob1p4,ob2p4)*betaY(ob1p4,ob2p4) + gamma(ob1p4,ob2p4)*betaY(ob1p4,ob2p4)*ob1p4.E()
        boostPz  = lambda ob1p4,ob2p4 : ob1p4.Pz() + gamma2(ob1p4,ob2p4)*betap(ob1p4,ob2p4)*betaZ(ob1p4,ob2p4) + gamma(ob1p4,ob2p4)*betaZ(ob1p4,ob2p4)*ob1p4.E()
        boostP   = lambda ob1p4,ob2p4 : op.sqrt(op.pow(boostPx(ob1p4,ob2p4),2) + op.pow(boostPy(ob1p4,ob2p4),2) + op.pow(boostPz(ob1p4,ob2p4),2))
        self.comp_cosThetaS = lambda ob1p4,ob2p4 : op.abs(boostPz(ob1p4,ob2p4)/boostP(ob1p4,ob2p4))

        #BoostP3  = lambda ob1p4,ob2p4 : op._to.Construct("ROOT::Math::TVector<ROOT::Math::PxPyPz3D<float>>",(-motherPx(ob1p4,ob2p4), -motherPy(ob1p4,ob2p4), -motherPz(ob1p4,ob2p4))).result
        #boost = lambda ob1p4,ob2p4 : op.construct("ROOT::Math::Boost", (-motherPx(ob1p4,ob2p4)/motherE(ob1p4,ob2p4), 
        #                                                                -motherPy(ob1p4,ob2p4)/motherE(ob1p4,ob2p4), 
        #                                                                -motherPz(ob1p4,ob2p4)/motherE(ob1p4,ob2p4)))
        #self.comp_cosThetaS = lambda ob1p4,ob2p4 : op.abs(boost(ob1p4,ob2p4)(ob1p4).CosTheta()) 

        #p4_boosted = lambda ob1p4,ob2p4 : op.extMethod("ROOT::Math::Boost{-motherPx(ob1p4,ob2p4)/motherE(ob1p4,ob2p4), -motherPy(ob1p4,ob2p4)/motherE(ob1p4,ob2p4), -motherPz(ob1p4,ob2p4)/motherE(ob1p4,ob2p4)}", returnType=(ob1p4+ob2p4)._typeName)(ob1p4+ob2p4)
        #self.comp_cosThetaS = lambda ob1p4,ob2p4 : op.deltaR(ob1p4, p4_boosted(ob1p4,ob2p4))  

        #boost = lambda ob1p4, ob2p4: op.construct("ROOT::Math::Boost", (-betaX(ob1p4, ob2p4), -betaY(ob1p4, ob2p4), -betaZ(ob1p4, ob2p4)))
        #boostP4 = lambda ob1p4,ob2p4 : boost(ob1p4,ob2p4)(ob1p4)
        #self.comp_cosThetaS = lambda ob1p4,ob2p4 : op.abs(boostP4(ob1p4,ob2p4).Pz()/op.sqrt(op.pow(boostP4(ob1p4,ob2p4).Px(),2) + op.pow(boostP4(ob1p4,ob2p4).Py(),2) + op.pow(boostP4(ob1p4,ob2p4).Pz(),2)))


        # MET_LD
        # Equation 3 (page 33) of AN-2019/111 v13
        # Similar to MET, but more robust against pileup
        jetSumPx = lambda jets : op.rng_sum(jets, lambda j : j.p4.Px())
        jetSumPy = lambda jets : op.rng_sum(jets, lambda j : j.p4.Py())
        #lepSumPx = lambda leps : op.rng_sum(leps, lambda l : l.p4.Px())
        #lepSumPy = lambda leps : op.rng_sum(leps, lambda l : l.p4.Py())
        lepSumPx = lambda mus, els : op.rng_sum(mus, lambda l : l.p4.Px()) + op.rng_sum(els, lambda l : l.p4.Px())
        lepSumPy = lambda mus, els : op.rng_sum(mus, lambda l : l.p4.Py()) + op.rng_sum(els, lambda l : l.p4.Py())
        self.MET_LD = lambda met, jets, mus, els : 0.6*met.pt + 0.4*op.sqrt(op.pow(jetSumPx(jets)+lepSumPx(mus,els),2)+op.pow(jetSumPy(jets)+lepSumPy(mus,els),2)) 
        
        empty_p4 = op.construct("ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >",([op.c_float(0.),op.c_float(0.),op.c_float(0.),op.c_float(0.)]))
        self.MET_LD_DL = lambda met, jets, electrons, muons : 0.6 * met.pt +\
                    0.4* (op.rng_sum(jets, (lambda j : j.p4), start=empty_p4) + op.rng_sum(electrons, (lambda e : e.p4), start=empty_p4) + op.rng_sum(muons, (lambda m : m.p4), start=empty_p4)).Pt()


        # conept
        self.lambdaConePt = lambda lep : op.switch(op.abs(lep.pdgId) == 13, HHself.muon_conept[lep.idx], HHself.electron_conept[lep.idx])
       

        # angle between 2 planes
        aDotB = lambda a, b : a.Px()*b.Px() + a.Py()*b.Py() + a.Pz()*b.Pz()
        aMagB = lambda a, b : (op.sqrt(op.pow(a.Px(),2)+op.pow(a.Py(),2)+op.pow(a.Pz(),2)))*(op.sqrt(op.pow(b.Px(),2)+op.pow(b.Py(),2)+op.pow(b.Pz(),2)))
        self.angleWWplane = lambda lp4, met, j3p4, j4p4 : op.acos(aDotB(j3p4+j4p4, self.neuP4(j3p4+j4p4+lp4, met)+lp4)/aMagB(j3p4+j4p4, self.neuP4(j3p4+j4p4+lp4, met)+lp4))
        #self.angleWWplane = lambda lp4, met, j3p4, j4p4 : ((j3p4+j4p4).Vect().Unit()).Angle((self.neuP4(j3p4+j4p4+lp4, met)+lp4).Vect().Unit())
        self.angleBetPlanes = lambda j1p4,j2p4,j3p4,j4p4 : op.acos(op.c_float(aDotB(j1p4+j2p4, j3p4+j4p4)/aMagB(j1p4+j2p4, j3p4+j4p4)))

        self.empty_p4 = op.construct("ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >",([op.c_float(0.),op.c_float(0.),op.c_float(0.),op.c_float(0.)]))
        self.MET_LD_DL = lambda met, jets, electrons, muons : 0.6 * met.pt +\
                    0.4* (op.rng_sum(jets, (lambda j : j.p4), start=self.empty_p4) + op.rng_sum(electrons, (lambda e : e.p4), start=self.empty_p4) + op.rng_sum(muons, (lambda m : m.p4), start=self.empty_p4)).Pt()

        self.isBoosted  = op.AND(op.rng_len(HHself.ak8BJets) >= 1, op.rng_len(HHself.ak4JetsCleanedFromAk8b) >= 1)
        #self.isBoosted  = op.rng_len(HHself.ak8BJets) >= 1
        self.isResolved = op.AND(op.rng_len(HHself.ak4Jets)  >= 3,op.rng_len(HHself.ak4BJets) >= 1,op.rng_len(HHself.ak8BJets) == 0)
        self.has1Wj     = op.rng_len(HHself.probableWJets) == 1
        self.has2Wj     = op.rng_len(HHself.wJetsPairs) >= 1
        self.isFullReco = op.AND(op.rng_len(HHself.bJetsByScore) >= 2, op.rng_len(HHself.wJetsPairs) >= 1)
        self.isMissReco = op.AND(op.rng_len(HHself.bJetsByScore) >= 2, op.rng_len(HHself.probableWJets) == 1)

        self.comp_m_hh_bregcorr = lambda bjets, wjets, lep, met : (op.rng_sum(bjets, (lambda bj : self.bJetCorrP4(bj)), start=empty_p4) + 
                                                                   op.rng_sum(wjets, (lambda wj : self.bJetCorrP4(wj)), start=empty_p4) + 
                                                                   met.p4 + 
                                                                   lep.p4).M()
        self.comp_pt_hh         = lambda bjets, wjets, lep, met : (op.rng_sum(bjets, (lambda bj : bj.p4), start=empty_p4) + 
                                                                   op.rng_sum(wjets, (lambda wj : wj.p4), start=empty_p4) + 
                                                                   met.p4 + 
                                                                   lep.p4).Pt()
        self.comp_dphi_hbb_hww  = lambda bjets, wjets, lep, met : op.deltaPhi((op.rng_sum(wjets, (lambda wj : wj.p4), start=empty_p4) + met.p4 + lep.p4),
                                                                              op.rng_sum(bjets, (lambda bj : bj.p4), start=empty_p4))
        self.comp_dphi_hbb_hwwvis = lambda bjets, wjets, lep : op.deltaPhi((op.rng_sum(wjets, (lambda wj : wj.p4), start=empty_p4) + lep.p4),
                                                                           op.rng_sum(bjets, (lambda bj : bj.p4), start=empty_p4))

    def getCorrBp4(self,j):
        return  op.construct("ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >", (j.pt*j.bRegCorr, j.eta, j.phi, j.mass*j.bRegCorr)) 

    # check florian's snippet
    def comp_cosThetaSbetBeamAndHiggs(self,genColl):
        genh = op.select(genColl,lambda g : op.AND(g.pdgId==25, g.statusFlags & ( 0x1 << 13)))
        HH_p4 = genh[0].p4 + genh[1].p4
        cm = HH_p4.BoostToCM()
        boosted_h1 = op.extMethod("ROOT::Math::VectorUtil::boost", returnType=genh[0].p4._typeName)(genh[0].p4,cm)
        boosted_h2 = op.extMethod("ROOT::Math::VectorUtil::boost", returnType=genh[1].p4._typeName)(genh[1].p4,cm)
        mHH = op.switch(op.rng_len(genh)==2, op.invariant_mass(genh[0].p4,genh[1].p4) , op.c_float(-9999))
        cosTheta1 = op.switch(op.rng_len(genh)==2, op.abs(boosted_h1.Pz()/boosted_h1.P()) , op.c_float(-9999))
        cosTheta2 = op.switch(op.rng_len(genh)==2, op.abs(boosted_h1.Pz()/boosted_h2.P()) , op.c_float(-9999))
        return [mHH, cosTheta1, cosTheta2]

    def comp_smin(self,lep,met,jets,bjets,wjets):
        self.blank_p4 = op.construct("ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >",([op.c_float(0.),op.c_float(0.),op.c_float(0.),op.c_float(0.)]))
        #vis = op.multiSwitch((op.rng_len(jets) >= 4, bjets[0].p4+bjets[1].p4+wjets[0].p4+wjets[1].p4+lep.p4),
        #                     (op.rng_len(jets) == 3, bjets[0].p4+bjets[1].p4+wjets[0].p4+lep.p4),
        #                     op.construct("ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >", (op.c_float(0.),
        #                                                                                                 op.c_float(0.),
        #                                                                                                 op.c_float(0.),
        #                                                                                                 op.c_float(0.))))
        vis = op.rng_sum(bjets, (lambda bj : bj.p4), start=self.blank_p4) + op.rng_sum(wjets, (lambda wj : wj.p4), start=self.blank_p4) + lep.p4

        vis_pt = vis.Pt()
        vis_m  = vis.M()
        vis_et = op.sqrt(op.pow(vis_m, 2) + op.pow(vis_pt, 2))
        met_et = met.p4.E()
        return op.sqrt(op.pow(vis_m, 2) + 2 * (vis_et * met_et - (vis.Px() * met.p4.Px() + vis.Py() * met.p4.Py())))
        


    
