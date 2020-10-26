from bamboo import treefunctions as op

class highlevelLambdas:
    def __init__(self,HHself):
        # All the attributes of the BaseHH are contained in HHself object
        # All the lambdas will be saved in the highlevelLambdas object to avoid confusions of all the attributes of HH base object

        # conept #
        self.conept = lambda lep : op.switch(op.abs(lep.pdgId)==11,
                                             HHself.electron_conept[lep.idx],
                                             HHself.muon_conept[lep.idx])

        # 4-Momentum association #
        self.ll_p4 = lambda l1,l2 : l1.p4+l2.p4
        self.lljj_p4 = lambda l1,l2,j1,j2 : l1.p4+l2.p4+j1.p4+j2.p4
        self.lep1j_p4 = lambda lep,j1 : lep.p4+j1.p4
        self.lep2j_p4 = lambda lep,j1,j2 : lep.p4+j1.p4+j2.p4
        self.lep3j_p4 = lambda lep,j1,j2,j3 : lep.p4+j1.p4+j2.p4+j3.p4
        self.lep4j_p4 = lambda lep,j1,j2,j3,j4 : lep.p4+j1.p4+j2.p4+j3.p4+j4.p4
        

        # bReg corr 4 momenta of ak4-bTagged jet #
        # https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsWG/BJetRegression
        self.bJetCorrP4 = lambda j : op._to.Construct("ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >", (j.pt*j.bRegCorr, j.eta, j.phi, j.mass)).result
        
        # Dilep-Met variables #
        self.DilepMET_deltaPhi = lambda l1,l2,met : self.ll_p4(l1,l2).Phi()-met.phi
        self.DilepMET_Pt = lambda l1,l2,met : op.sqrt(op.pow(met.pt*op.cos(met.phi)+self.ll_p4(l1,l2).Px(),2)+op.pow(met.pt*op.sin(met.phi)+self.ll_p4(l1,l2).Py(),2))
        # SingleLep-Met variables
        self.SinglepMet_Pt = lambda lep,met : op.sqrt(op.pow(met.pt*op.cos(met.phi)+lep.p4.Px(),2)+op.pow(met.pt*op.sin(met.phi)+lep.p4.Py(),2))
        self.SinglepMet_dPhi = lambda lep, met : lep.p4.Phi()-met.phi
        
        # Transverse mass #
        self.MT_ll = lambda l1,l2,met : op.sqrt(2*self.ll_p4(l1,l2).Pt()*met.pt*(1-op.cos(self.ll_p4(l1,l2).Phi()-met.phi)))
        self.MT_lljj = lambda l1,l2,j1,j2,met : op.sqrt(2*self.lljj_p4(l1,l2,j1,j2).Pt()*met.pt*(1-op.cos(self.lljj_p4(l1,l2,j1,j2).Phi()-met.phi)))
        self.MT = lambda lep,met : op.sqrt(2*lep.pt*met.pt*(1-op.cos(lep.p4.Phi()-met.phi)))
        self.MT_W1W2_ljj = lambda lep,j1,j2,met : op.sqrt(2*self.lep2j_p4(lep,j1,j2).Pt()*met.pt*(1-op.cos(self.lep2j_p4(lep,j1,j2).Phi()-met.phi)))
        self.MT_W1W2_lj  = lambda lep,j1,met : op.sqrt(2*self.lep1j_p4(lep,j1).Pt()*met.pt*(1-op.cos(self.lep1j_p4(lep,j1).Phi()-met.phi)))
            # TODO : clean different versions (eg MT)
        
        # dilep + dijet #
        self.M_lljj = lambda l1,l2,j1,j2 : op.invariant_mass(self.lljj_p4(l1,l2,j1,j2))
        self.MinDR_lj = lambda l1,l2,j1,j2 : op.min(op.min(op.deltaR(l1.p4,j1.p4),op.deltaR(l1.p4,j2.p4)),
                                               op.min(op.deltaR(l2.p4,j1.p4),op.deltaR(l2.p4,j2.p4)))
        self.MinDR_part1_partCont = lambda part1,partCont : op.rng_min(partCont, lambda part2 : op.deltaR(part1.p4, part2.p4))
        self.MinDR_part1_dipart = lambda part1,dipart: op.min(*(op.deltaR(part1.p4, dipart[i2].p4) for i2 in range(2)))

        self.JetsMinDR = lambda l,j1,j2 : op.min(op.deltaR(l.p4,j1.p4),op.deltaR(l.p4,j2.p4))
        self.LepsMinDR = lambda j,l1,l2 : op.min(op.deltaR(j.p4,l1.p4),op.deltaR(j.p4,l2.p4))
        
        self.MinDR_lep3j = lambda lep,j1,j2,j3 : op.min(op.min(op.deltaR(lep.p4,j1.p4),op.deltaR(lep.p4,j2.p4)),op.deltaR(lep.p4,j3.p4))
        
        # Higgs related variables #
        self.HT2 = lambda l1,l2,j1,j2,met : op.sqrt(op.pow(met.pt*op.cos(met.phi)+l1.p4.Px()+l2.p4.Px(),2)+op.pow(met.pt*op.sin(met.phi)+l1.p4.Py()+l2.p4.Py(),2)) + op.abs((j1.p4+j2.p4).Pt())
        self.HT2R = lambda l1,l2,j1,j2,met : self.HT2(met,l1,l2,j1,j2)/(met.pt+l1.pt+l2.pt+j1.pt+j2.pt)
        self.HT2_l3jmet  = lambda l,j1,j2,j3,met : op.sqrt(op.pow(met.pt*op.cos(met.phi)+l.p4.Px(),2)+op.pow(met.pt*op.sin(met.phi)+l.p4.Py(),2)) + op.abs((j1.p4+j2.p4+j3.p4).Pt())
        self.HT2R_l3jmet = lambda l,j1,j2,j3,met : self.HT2_l3jmet(met,l,j1,j2,j3)/(met.pt+l.pt+j1.pt+j2.pt+j3.pt)
        self.HT2_l4jmet  = lambda l,j1,j2,j3,j4,met : op.sqrt(op.pow(met.pt*op.cos(met.phi)+l.p4.Px(),2)+op.pow(met.pt*op.sin(met.phi)+l.p4.Py(),2)) + op.abs((j1.p4+j2.p4+j3.p4+j4.p4).Pt())
        self.HT2R_l4jmet = lambda l,j1,j2,j3,j4,met : self.HT2_l4jmet(met,l,j1,j2,j3,j4)/(met.pt+l.pt+j1.pt+j2.pt+j3.pt+j4.pt)
        
        #min j1j2DR
        self.MinDiJetDRLoose = lambda j1,j2,j3: op.min(op.min(op.deltaR(j1.p4,j2.p4), op.deltaR(j2.p4,j3.p4)), op.deltaR(j1.p4,j3.p4))
        
        # ------------------------------------ lambdas for BDT variables ------------------------------------ #
        self.mindr_lep1_jet = lambda lep, jets : op.deltaR(lep.p4, op.sort(jets, lambda j : op.deltaR(lep.p4, j.p4))[0].p4)
        self.HT = lambda jets : op.rng_sum(jets, lambda j : j.pt)
        
        # mT2
        self.ET  = lambda lep : op.sqrt(op.pow(lep.p4.M(),2) + op.pow(lep.pt,2))
        self.mT2 = lambda jet, lep, met : (op.pow(jet.p4.M(),2) + op.pow(lep.p4.M(),2) + op.pow(met.p4.M(),2) + 
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
        self.Wlep_simple = lambda j1P4,j2P4,lepP4,met,mH : lepP4 + neuP4(j1P4+j2P4+lepP4, met, mH)
        # P4 of W2 (j,j)
        self.Wjj_simple  = lambda j1P4,j2P4 : j1P4 + j2P4 
        # P4 of HWW (W1 + W2)
        self.HWW_simple  = lambda j1P4,j2P4,lepP4,met,mH : Wjj_simple(j1P4,j2P4) + Wlep_simple(lepP4,neuP4(j1P4+j2P4+lepP4, met, mH))
        # dR_HWW
        self.dR_Hww   = lambda j1P4,j2P4,lepP4,met,mH : op.deltaR(Wjj_simple(j1P4,j2P4), Wlep_simple(j1P4,j2P4,lepP4,met,mH))
        # P4 of lep + met
        self.Wlep_met_simple = lambda lepP4, metP4 : lepP4 + metP4
        # SimpleP4 of HWW (W1 + W2)
        self.HWW_met_simple  = lambda j1P4,j2P4,lepP4,metP4 : Wjj_simple(j1P4, j2P4) + Wlep_met_simple(lepP4,metP4)
        # Total P4
        self.HHP4_simple_met = lambda HbbRegP4, j1P4, j2P4, lepP4,metP4 : HbbRegP4 + Wjj_simple(j1P4, j2P4) + Wlep_met_simple(lepP4,metP4)
        
        # CosThetaS calculation
        #comp_cosThetaS = lambda ob1p4, ob2p4 : op.abs(ob1p4.Boost(-(ob1p4+ob2p4).BoostVector()).CosTheta())
        motherPx = lambda ob1p4,ob2p4 : (ob1p4.Px()+ob2p4.Px())
        motherPy = lambda ob1p4,ob2p4 : (ob1p4.Py()+ob2p4.Py())
        motherPz = lambda ob1p4,ob2p4 : (ob1p4.Pz()+ob2p4.Pz())
        motherE  = lambda ob1p4,ob2p4 : (ob1p4.E()+ob2p4.E())
        BoostP4 = lambda ob1p4,ob2p4 : op._to.Construct("ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >", (motherPx(ob1p4,ob2p4), motherPy(ob1p4,ob2p4), motherPz(ob1p4,ob2p4), motherE(ob1p4,ob2p4))).result
        self.comp_cosThetaS = lambda ob1p4,ob2p4 : op.abs(op.cos(op.deltaR(BoostP4(ob1p4,ob2p4), ob1p4)))
        
        # MET_LD
        # Equation 3 (page 33) of AN-2019/111 v13
        # Similar to MET, but more robust against pileup
        jetSumPx = lambda jets : op.rng_sum(jets, lambda j : j.p4.Px())
        jetSumPy = lambda jets : op.rng_sum(jets, lambda j : j.p4.Py())
        lepSumPx = lambda leps : op.rng_sum(leps, lambda l : l.p4.Px())
        lepSumPy = lambda leps : op.rng_sum(leps, lambda l : l.p4.Py())
        self.MET_LD = lambda met, jets, leps : 0.6*met.pt + 0.4*op.sqrt(op.pow(jetSumPx(jets)+lepSumPx(leps),2)+op.pow(jetSumPy(jets)+lepSumPy(leps),2)) 

        SumPx = lambda objs : op.rng_sum(objs, lambda obj : obj.p4.Px())
        SumPy = lambda objs : op.rng_sum(objs, lambda obj : obj.p4.Py())

        #self.MET_LD_DL = lambda met, jets, l1 ,l2 : 0.6*met.pt + 0.4*op.sqrt(op.pow(SumPx(jets)+l1.p4.Px()+l2.p4.Px(),2),
        #                                                                    +op.pow(SumPy(jets)+l1.p4.Py()+l2.p4.Py(),2))
        self.MET_LD_DL = lambda met, jets, l1 ,l2 : 0.6*met.pt + 0.4*op.rng_sum(jets, (lambda j : j.p4), start=(l1.p4+l2.p4)).Pt()
                                                                    
    def getCorrBp4(self,j):
        return  op.construct("ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >", (j.pt*j.bRegCorr, j.eta, j.phi, j.mass*j.bRegCorr)) 
   
