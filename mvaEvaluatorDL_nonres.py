from bamboo import treefunctions as op

def returnLeptonsMVAInputs(self,l1,l2,channel):
    if channel == "ElEl":
        cone_l1 = self.getElectronConeP4(l1) 
        cone_l2 = self.getElectronConeP4(l2)
    elif channel == "MuMu":
        cone_l1 = self.getMuonConeP4(l1)
        cone_l2 = self.getMuonConeP4(l2)
    elif channel == "ElMu":
        cone_l1 = self.getElectronConeP4(l1)
        cone_l2 = self.getMuonConeP4(l2)
    else:
        raise RuntimeError("Wrong channel")

    return {('l1_E',        'Lead lepton E [GeV]',          (50,0.,500.))        : op.switch(cone_l1.Pt() >= cone_l2.Pt() , cone_l1.E(),  cone_l2.E()),
            ('l1_Px',       'Lead lepton P_x [GeV]',        (40,-200.,200.))     : op.switch(cone_l1.Pt() >= cone_l2.Pt() , cone_l1.Px(), cone_l2.Px()),
            ('l1_Py',       'Lead lepton P_y [GeV]',        (40,-200.,200.))     : op.switch(cone_l1.Pt() >= cone_l2.Pt() , cone_l1.Py(), cone_l2.Py()),
            ('l1_Pz',       'Lead lepton P_z [GeV]',        (40,-200.,200.))     : op.switch(cone_l1.Pt() >= cone_l2.Pt() , cone_l1.Pz(), cone_l2.Pz()),
            ('l1_pdgId',    'Lead lepton pdg ID',           (45,-22.,22.) )      : op.switch(cone_l1.Pt() >= cone_l2.Pt() , l1.pdgId,     l2.pdgId),
            ('l1_charge',   'Lead lepton charge',           (2,0.,2.)   )        : op.switch(cone_l1.Pt() >= cone_l2.Pt() , l1.charge,    l2.charge),
            ('l2_E',        'Sublead lepton E [GeV]',       (50,0.,500.))        : op.switch(cone_l1.Pt() >= cone_l2.Pt() , cone_l2.E(),  cone_l1.E()),
            ('l2_Px',       'Sublead lepton P_x [GeV]',     (40,-200.,200.))     : op.switch(cone_l1.Pt() >= cone_l2.Pt() , cone_l2.Px(), cone_l1.Px()),
            ('l2_Py',       'Sublead lepton P_y [GeV]',     (40,-200.,200.))     : op.switch(cone_l1.Pt() >= cone_l2.Pt() , cone_l2.Py(), cone_l1.Py()),
            ('l2_Pz',       'Sublead lepton P_z [GeV]',     (40,-200.,200.))     : op.switch(cone_l1.Pt() >= cone_l2.Pt() , cone_l2.Pz(), cone_l1.Pz()),
            ('l2_pdgId',    'Sublead lepton pdg ID',        (45,-22.,22.)  )     : op.switch(cone_l1.Pt() >= cone_l2.Pt() , l2.pdgId,     l1.pdgId),
            ('l2_charge',   'Sublead lepton charge',        (2,0.,2.)  )         : op.switch(cone_l1.Pt() >= cone_l2.Pt() , l2.charge,    l1.charge)}


def returnJetsMVAInputs(self,jets):
    return {('j1_E',        'Lead jet E [GeV]',             (50,0.,500.))        : op.switch(op.rng_len(jets)>0,jets[0].p4.E(),op.c_float(0.)),
            ('j1_Px',       'Lead jet P_x [GeV]',           (40,-200.,200.))     : op.switch(op.rng_len(jets)>0,jets[0].p4.Px(),op.c_float(0.)),
            ('j1_Py',       'Lead jet P_y [GeV]',           (40,-200.,200.))     : op.switch(op.rng_len(jets)>0,jets[0].p4.Py(),op.c_float(0.)),
            ('j1_Pz',       'Lead jet P_z [GeV]',           (40,-200.,200.))     : op.switch(op.rng_len(jets)>0,jets[0].p4.Pz(),op.c_float(0.)),
            ('j1_btag',     'Lead jet btag score',          (50,0.,1.)  )        : op.switch(op.rng_len(jets)>0,jets[0].btagDeepFlavB,op.c_float(0.)),
            ('j2_E',        'Sublead jet E [GeV]',          (50,0.,500.))        : op.switch(op.rng_len(jets)>1,jets[1].p4.E(),op.c_float(0.)),
            ('j2_Px',       'Sublead jet P_x [GeV]',        (40,-200.,200.))     : op.switch(op.rng_len(jets)>1,jets[1].p4.Px(),op.c_float(0.)),
            ('j2_Py',       'Sublead jet P_y [GeV]',        (40,-200.,200.))     : op.switch(op.rng_len(jets)>1,jets[1].p4.Py(),op.c_float(0.)),
            ('j2_Pz',       'Sublead jet P_z [GeV]',        (40,-200.,200.))     : op.switch(op.rng_len(jets)>1,jets[1].p4.Pz(),op.c_float(0.)),
            ('j2_btag',     'Sublead jet btag score',       (50,0.,1.)  )        : op.switch(op.rng_len(jets)>1,jets[1].btagDeepFlavB,op.c_float(0.)),
            ('j3_E',        'Subsublead jet E [GeV]',       (50,0.,500.))        : op.switch(op.rng_len(jets)>2,jets[2].p4.E(),op.c_float(0.)),
            ('j3_Px',       'Subsublead jet P_x [GeV]',     (40,-200.,200.))     : op.switch(op.rng_len(jets)>2,jets[2].p4.Px(),op.c_float(0.)),
            ('j3_Py',       'Subsublead jet P_y [GeV]',     (40,-200.,200.))     : op.switch(op.rng_len(jets)>2,jets[2].p4.Py(),op.c_float(0.)),
            ('j3_Pz',       'Subsublead jet P_z [GeV]',     (40,-200.,200.))     : op.switch(op.rng_len(jets)>2,jets[2].p4.Pz(),op.c_float(0.)),
            ('j3_btag',     'Subsublead jet btag score',    (50,0.,1.)  )        : op.switch(op.rng_len(jets)>2,jets[2].btagDeepFlavB,op.c_float(0.)),
            ('j4_E',        'Subsubsublead jet E [GeV]',    (50,0.,500.))        : op.switch(op.rng_len(jets)>3,jets[3].p4.E(),op.c_float(0.)),
            ('j4_Px',       'Subsubsublead jet P_x [GeV]',  (40,-200.,200.))     : op.switch(op.rng_len(jets)>3,jets[3].p4.Px(),op.c_float(0.)),
            ('j4_Py',       'Subsubsublead jet P_y [GeV]',  (40,-200.,200.))     : op.switch(op.rng_len(jets)>3,jets[3].p4.Py(),op.c_float(0.)),
            ('j4_Pz',       'Subsubsublead jet P_z [GeV]',  (40,-200.,200.))     : op.switch(op.rng_len(jets)>3,jets[3].p4.Pz(),op.c_float(0.)),
            ('j4_btag',     'Subsubsublead jet btag score', (50,0.,1.)  )        : op.switch(op.rng_len(jets)>3,jets[3].btagDeepFlavB,op.c_float(0.))}


def returnFatjetMVAInputs(self,fatjets):
    return {('fatjet_E',        'Fatjet E [GeV]',               (50,0.,500.))        : op.switch(op.rng_len(fatjets)>0, fatjets[0].p4.E(), op.c_float(0.)),
            ('fatjet_Px',       'Fatjet P_x [GeV]',             (40,-200.,200.))     : op.switch(op.rng_len(fatjets)>0, fatjets[0].p4.Px(), op.c_float(0.)),
            ('fatjet_Py',       'Fatjet P_y [GeV]',             (40,-200.,200.))     : op.switch(op.rng_len(fatjets)>0, fatjets[0].p4.Py(), op.c_float(0.)),
            ('fatjet_Pz',       'Fatjet P_z [GeV]',             (40,-200.,200.))     : op.switch(op.rng_len(fatjets)>0, fatjets[0].p4.Pz(), op.c_float(0.)),
            ('fatjet_tau1',     'Fatjet #tau_1',                (50,0.,1.))          : op.switch(op.rng_len(fatjets)>0, fatjets[0].tau1, op.c_float(0.)),
            ('fatjet_tau2',     'Fatjet #tau_2',                (50,0.,1.))          : op.switch(op.rng_len(fatjets)>0, fatjets[0].tau2, op.c_float(0.)), 
            ('fatjet_tau3',     'Fatjet #tau_3',                (50,0.,1.))          : op.switch(op.rng_len(fatjets)>0, fatjets[0].tau3, op.c_float(0.)),
            ('fatjet_tau4',     'Fatjet #tau_4',                (50,0.,1.))          : op.switch(op.rng_len(fatjets)>0, fatjets[0].tau4, op.c_float(0.)), 
            ('fatjet_softdrop', 'Fatjet softdrop mass [GeV]',   (50,0.,1000.))       : op.switch(op.rng_len(fatjets)>0, fatjets[0].msoftdrop, op.c_float(0.))}


def returnMETMVAInputs(self,met):
    return {('met_E',       'MET Energy',                   (50,0.,500.))        : met.p4.E(),
            ('met_Px',      'MET P_x',                      (40,-200.,200.))     : met.p4.Px(),
            ('met_Py',      'MET P_y',                      (40,-200.,200.))     : met.p4.Py(),
            ('met_Pz',      'MET P_z',                      (40,-200.,200.))     : met.p4.Pz()}

def returnHighLevelMVAInputs(self,l1,l2,met,jets,bjets,electrons,muons,channel):
    if channel == "ElEl":
        cone_l1 = self.getElectronConeP4(l1) 
        cone_l2 = self.getElectronConeP4(l2)
    elif channel == "MuMu":
        cone_l1 = self.getMuonConeP4(l1)
        cone_l2 = self.getMuonConeP4(l2)
    elif channel == "ElMu":
        cone_l1 = self.getElectronConeP4(l1)
        cone_l2 = self.getMuonConeP4(l2)
    else:
        raise RuntimeError("Wrong channel")

    dijets = op.combine(jets, N=2)

    import bamboo.treeoperations as _to
    def rng_min(rng, fun=(lambda x : x), typeName="float"):
        return op._to.Reduce.fromRngFun(rng, op.c_float(float("+inf"), typeName), ( lambda fn : (
                lambda res, elm : op.extMethod("std::min", returnType="Float_t")(res, fn(elm))
                        ) )(fun) )

    if self.args.Boosted0Btag or self.args.Boosted1Btag:
        VBFJetPairs = self.VBFJetPairsBoosted
    elif self.args.Resolved0Btag or self.args.Resolved1Btag or self.args.Resolved2Btag:
        VBFJetPairs = self.VBFJetPairsResolved
    else:
        raise RuntimeError("Wrong selection to be used by the DNN")

    return {
            ('m_bb_bregcorr',          'Di-bjet invariant mass (regcorr) [GeV]',  (100,0.,1000.))   : op.multiSwitch((op.rng_len(bjets) == 0, op.c_float(0.)),
                                                                                                                     (op.rng_len(bjets) == 1, self.HLL.getCorrBp4(bjets[0]).M()),
                                                                                                                     op.invariant_mass(self.HLL.getCorrBp4(bjets[0]),self.HLL.getCorrBp4(bjets[1]))),
            ('ht',                     'HT(jets) [GeV]',                          (100,0.,1000.))   : op.rng_sum(jets, lambda j : j.pt),
            ('min_dr_jets_lep1',       'Min(#Delta R(lead lepton,jets))',         (25,0.,5.))       : op.switch(op.rng_len(jets) > 0,
                                                                                                                op.switch(cone_l1.Pt() >= cone_l2.Pt(),
                                                                                                                          self.HLL.MinDR_part1_partCont(cone_l1,jets),
                                                                                                                          self.HLL.MinDR_part1_partCont(cone_l2,jets)),
                                                                                                                op.c_float(0.)),
            ('min_dr_jets_lep2',       'Min(#Delta R(sublead lepton,jets))',      (25,0.,5.))       : op.switch(op.rng_len(jets) > 0,
                                                                                                                op.switch(cone_l1.Pt() >= cone_l2.Pt(),
                                                                                                                          self.HLL.MinDR_part1_partCont(cone_l2,jets),
                                                                                                                          self.HLL.MinDR_part1_partCont(cone_l1,jets)),
                                                                                                                op.c_float(0.)),
            ('m_ll',                   'Dilepton invariant mass [GeV]',           (100,0.,1000.))   : op.invariant_mass(cone_l1,cone_l2),
            ('dr_ll',                  'Dilepton #Delta R',                       (25,0.,5.))       : op.deltaR(cone_l1,cone_l2),
            ('min_dr_jet',             'Min(#Delta R(jets))',                     (25,0.,5.))       : op.switch(op.rng_len(dijets) > 0,
                                                                                                                op.rng_min(dijets,lambda dijet : op.deltaR(dijet[0].p4,dijet[1].p4)),
                                                                                                                op.c_float(0.)),
            ('min_dhi_jet',            'Min(#Delta #Phi(jets))',                  (16,0.,3.2))      : op.switch(op.rng_len(dijets) > 0,
                                                                                                                rng_min(dijets,lambda dijet : op.abs(op.deltaPhi(dijet[0].p4,dijet[1].p4)),typeName='double'),
                                                                                                                op.c_float(0.)),
            ('m_hh_simplemet_bregcorr','M_{HH} (simple MET) (regcorr) [GeV]',     (100,0.,1000.))   : op.invariant_mass(op.rng_sum(bjets,
                                                                                                                                   lambda bjet : self.HLL.getCorrBp4(bjet),
                                                                                                                                   start=self.HLL.empty_p4),
                                                                                                                        cone_l1,
                                                                                                                        cone_l2,
                                                                                                                        met.p4),
            ('met_ld',                 'MET_{LD}',                                (100,0.,1000.))   : self.HLL.MET_LD_DL(met,jets,electrons,muons),
            ('dr_bb',                  'Di-bjet #Delta R',                        (25,0.,5.))       : op.switch(op.rng_len(bjets)>=2,
                                                                                                                op.deltaR(bjets[0].p4,bjets[1].p4),
                                                                                                                op.c_float(0.)),
            ('min_dr_leps_b1',         'Min(#Delta R(lead bjet,dilepton))',       (25,0.,5.))       : op.switch(op.rng_len(bjets)>=1,
                                                                                                                self.HLL.MinDR_part1_dipart(bjets[0].p4,[cone_l1,cone_l2]),
                                                                                                                op.c_float(0.)),
            ('min_dr_leps_b2',         'Min(#Delta R(sublead bjet,dilepton))',    (25,0.,5.))       : op.switch(op.rng_len(bjets)>=2,
                                                                                                                self.HLL.MinDR_part1_dipart(bjets[1].p4,[cone_l1,cone_l2]),
                                                                                                                op.c_float(0.)),
            ('lep1_conept',            'Lead lepton cone-P_T [GeV]',              (40,0.,200.))     : op.switch(cone_l1.Pt() >= cone_l2.Pt() , cone_l1.Pt() , cone_l2.Pt()),
            ('lep2_conept',            'Sublead lepton cone-P_T [GeV]',           (40,0.,200.))     : op.switch(cone_l1.Pt() >= cone_l2.Pt() , cone_l2.Pt() , cone_l1.Pt()),
            ('mww_simplemet',          'M_{WW} (simple MET) [GeV]',               (100,0.,1000.))   : op.invariant_mass(cone_l1,cone_l2,met.p4),
            ('vbf_tag',                'VBF tag',                                 (2,0.,2.))        : op.c_int(op.rng_len(VBFJetPairs)>0),
            ('boosted_tag',            'Boosted tag',                             (2,0.,2.))        : op.c_int(op.OR(op.rng_len(self.ak8BJets) > 0,           # Boosted 1B
                                                                                                                     op.AND(op.rng_len(self.ak8BJets) == 0,   # Boosted 0B
                                                                                                                            op.rng_len(self.ak8Jets) > 0,
                                                                                                                            op.rng_len(self.ak4BJets) == 0))),
            ('dphi_met_dilep',         'Dilepton-MET #Delta #Phi',                (32,-3.2,3.2))    : op.abs(op.deltaPhi(met.p4,(cone_l1+cone_l2))),
            ('dphi_met_dibjet',        'Dibjet-MET #Delta #Phi',                  (32,-3.2,3.2))    : op.multiSwitch((op.rng_len(bjets) == 0, op.c_float(0.)),
                                                                                                                     (op.rng_len(bjets) == 1, op.abs(op.deltaPhi(met.p4,bjets[0].p4))),
                                                                                                                     op.abs(op.deltaPhi(met.p4,(bjets[0].p4+bjets[1].p4)))),
            ('dr_dilep_dijet',         'Dilepton-dijet #Delta R',                 (25,0.,5.))       : op.multiSwitch((op.rng_len(jets) == 0, op.c_float(0.)),
                                                                                                                     (op.rng_len(jets) == 1, op.deltaR((cone_l1+cone_l2),jets[0].p4)),
                                                                                                                     op.deltaR((cone_l1+cone_l2),(jets[0].p4+jets[1].p4))),
            ('dr_dilep_dibjet',        'Dilepton-dibjet #Delta R',                (25,0.,5.))       : op.multiSwitch((op.rng_len(bjets) == 0, op.c_float(0.)),
                                                                                                                     (op.rng_len(bjets) == 1, op.deltaR((cone_l1+cone_l2),bjets[0].p4)),
                                                                                                                     op.deltaR((cone_l1+cone_l2),(bjets[0].p4+bjets[1].p4))),
            ('vbf_pair_mass',          'VBF pair M_{jj}',                         (100,0.,1000.))   : op.switch(op.rng_len(VBFJetPairs)>0,
                                                                                                                op.invariant_mass(VBFJetPairs[0][0].p4,VBFJetPairs[0][1].p4),
                                                                                                                op.c_float(0.)),
            ('vbf_pairs_absdeltaeta',  'VBF pair #Delta#eta',                     (25,0.,5.))       : op.switch(op.rng_len(VBFJetPairs)>0,
                                                                                                                op.abs(VBFJetPairs[0][0].eta-VBFJetPairs[0][1].eta),
                                                                                                                op.c_float(0.)),
            ('sphericity',             'None',                                    (1,0.,1.))        : op.c_float(0.),
            ('sphericity_T',           'None',                                    (1,0.,1.))        : op.c_float(0.),
            ('aplanarity',             'None',                                    (1,0.,1.))        : op.c_float(0.),
            ('eventshape_C',           'None',                                    (1,0.,1.))        : op.c_float(0.),
            ('eventshape_D',           'None',                                    (1,0.,1.))        : op.c_float(0.),
            ('eventshape_Y',           'None',                                    (1,0.,1.))        : op.c_float(0.),
            ('foxwolfram1',            'None',                                    (1,0.,1.))        : op.c_float(0.),
            ('foxwolfram2',            'None',                                    (1,0.,1.))        : op.c_float(0.),
            ('foxwolfram3',            'None',                                    (1,0.,1.))        : op.c_float(0.),
            ('foxwolfram4',            'None',                                    (1,0.,1.))        : op.c_float(0.),
            ('foxwolfram5',            'None',                                    (1,0.,1.))        : op.c_float(0.),
            ('centrality',             'None',                                    (1,0.,1.))        : op.c_float(0.),
            ('centrality_jets',        'None',                                    (1,0.,1.))        : op.c_float(0.),
            ('eigenvalue1',            'None',                                    (1,0.,1.))        : op.c_float(0.),
            ('eigenvalue2',            'None',                                    (1,0.,1.))        : op.c_float(0.),
            ('eigenvalue3',            'None',                                    (1,0.,1.))        : op.c_float(0.),
            }

def returnParamMVAInputs(self):
    if self.args.era is None:
        return {('year',     'Year',         (3,2016.,2019.))    : op.c_int(int(self.era))}
    else:
        return {('year',     'Year',         (3,2016.,2019.))    : op.c_int(int(self.args.era))}

def returnEventNrMVAInputs(self,t):
    return {('eventnr',  'Event number', (int(1e6+1),0.,1e6))     : t.event}

def inputStaticCast(inputDict,cast='float'):
    return [op.static_cast(cast,v) for v in inputDict.values()]
    
