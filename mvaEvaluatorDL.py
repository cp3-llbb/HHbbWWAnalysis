from bamboo import treefunctions as op

def returnLeptonsMVAInputs(self,l1,l2,channel):
    if channel == "ElEl":
        l1conept = lambda l1 : self.electron_conept[l1.idx]
        l2conept = lambda l2 : self.electron_conept[l2.idx]
    if channel == "MuMu":
        l1conept = lambda l1 : self.muon_conept[l1.idx]
        l2conept = lambda l2 : self.muon_conept[l2.idx]
    elif channel == "ElMu":
        l1conept = lambda l1 : self.electron_conept[l1.idx]
        l2conept = lambda l2 : self.muon_conept[l2.idx]
    return {('l1_E',        'Lead lepton E [GeV]',          (50,0.,500.))        : op.switch(l1conept(l1) >= l2conept(l2) , l1.p4.E(), l2.p4.E()),
            ('l1_Px',       'Lead lepton P_x [GeV]',        (40,-200.,200.))     : op.switch(l1conept(l1) >= l2conept(l2) , l1.p4.Px(), l2.p4.Px()),
            ('l1_Py',       'Lead lepton P_y [GeV]',        (40,-200.,200.))     : op.switch(l1conept(l1) >= l2conept(l2) , l1.p4.Py(), l2.p4.Py()),
            ('l1_Pz',       'Lead lepton P_z [GeV]',        (40,-200.,200.))     : op.switch(l1conept(l1) >= l2conept(l2) , l1.p4.Pz(), l2.p4.Pz()),
            ('l1_charge',   'Lead lepton charge',           (2,0.,2.)   )        : op.switch(l1conept(l1) >= l2conept(l2) , l1.charge, l2.charge),
            ('l1_pdgId',    'Lead lepton pdg ID',           (45,-22.,22.) )      : op.switch(l1conept(l1) >= l2conept(l2) , l1.pdgId, l2.pdgId),
            ('l2_E',        'Sublead lepton E [GeV]',       (50,0.,500.))        : op.switch(l1conept(l1) >= l2conept(l2) , l2.p4.E(), l1.p4.E()),
            ('l2_Px',       'Sublead lepton P_x [GeV]',     (40,-200.,200.))     : op.switch(l1conept(l1) >= l2conept(l2) , l2.p4.Px(), l1.p4.Px()),
            ('l2_Py',       'Sublead lepton P_y [GeV]',     (40,-200.,200.))     : op.switch(l1conept(l1) >= l2conept(l2) , l2.p4.Py(), l1.p4.Py()),
            ('l2_Pz',       'Sublead lepton P_z [GeV]',     (40,-200.,200.))     : op.switch(l1conept(l1) >= l2conept(l2) , l2.p4.Pz(), l1.p4.Pz()),
            ('l2_charge',   'Sublead lepton charge',        (2,0.,2.)  )         : op.switch(l1conept(l1) >= l2conept(l2) , l2.charge, l1.charge), 
            ('l2_pdgId',    'Sublead lepton pdg ID',        (45,-22.,22.)  )     : op.switch(l1conept(l1) >= l2conept(l2) , l2.pdgId, l1.pdgId)}

def returnJetsMVAInputs(self,jets):
    return {('j1_E',        'Lead jet E [GeV]',             (50,0.,500.))        : jets[0].p4.E(),
            ('j1_Px',       'Lead jet P_x [GeV]',           (40,-200.,200.))     : jets[0].p4.Px(),
            ('j1_Py',       'Lead jet P_y [GeV]',           (40,-200.,200.))     : jets[0].p4.Py(),
            ('j1_Pz',       'Lead jet P_z [GeV]',           (40,-200.,200.))     : jets[0].p4.Pz(),
            ('j1_btag',     'Lead jet btag score',          (50,0.,1.)  )        : jets[0].btagDeepFlavB,
            ('j2_E',        'Sublead jet E [GeV]',          (50,0.,500.))        : jets[1].p4.E(),
            ('j2_Px',       'Sublead jet P_x [GeV]',        (40,-200.,200.))     : jets[1].p4.Px(),
            ('j2_Py',       'Sublead jet P_y [GeV]',        (40,-200.,200.))     : jets[1].p4.Py(),
            ('j2_Pz',       'Sublead jet P_z [GeV]',        (40,-200.,200.))     : jets[1].p4.Pz(),
            ('j2_btag',     'Sublead jet btag score',       (50,0.,1.)  )        : jets[1].btagDeepFlavB,
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
            ('met_Py',      'MET P_y',                      (40,-200.,200.))     : met.p4.Px(),
            ('met_Pz',      'MET P_z',                      (40,-200.,200.))     : met.p4.Pz()}

def returnHighLevelMVAInputs(self,l1,l2,b1,b2,met,jets,electrons,muons,channel):
    if channel == "ElEl":
        l1conept = lambda l1 : self.electron_conept[l1.idx]
        l2conept = lambda l2 : self.electron_conept[l2.idx]
    if channel == "MuMu":
        l1conept = lambda l1 : self.muon_conept[l1.idx]
        l2conept = lambda l2 : self.muon_conept[l2.idx]
    elif channel == "ElMu":
        l1conept = lambda l1 : self.electron_conept[l1.idx]
        l2conept = lambda l2 : self.muon_conept[l2.idx]
    dijets = op.combine(jets, N=2, samePred=(lambda j1,j2 : j1.idx != j2.idx))

    import bamboo.treeoperations as _to
    def rng_min(rng, fun=(lambda x : x), typeName="float"):
        return op._to.Reduce.fromRngFun(rng, op.c_float(float("+inf"), typeName), ( lambda fn : (
                lambda res, elm : op.extMethod("std::min", returnType="Float_t")(res, fn(elm))
                        ) )(fun) )
    if hasattr(b1,'bRegCorr') and hasattr(b2,'bRegCorr'):
        bRegCorr = self.HLL.getCorrBp4
    else:
        def bRegCorr(obj):
            return obj.p4

    if self.args.Boosted0Btag or self.args.Boosted1Btag:
        VBF_tag = op.c_int(op.rng_len(self.VBFJetPairsBoosted)>0)
    if self.args.Resolved0Btag or self.args.Resolved1Btag or self.args.Resolved2Btag:
        VBF_tag = op.c_int(op.rng_len(self.VBFJetPairsResolved)>0)
    else:
        raise RuntimeError("Wrong selection to be used by the DNN")


    return {('m_bb',                   'Di-bjet invariant mass [GeV]',            (100,0.,1000.))   : op.invariant_mass(b1.p4,b2.p4),
            ('m_bb_bregcorr',          'Di-bjet invariant mass (regcorr) [GeV]',  (100,0.,1000.))   : op.invariant_mass(bRegCorr(b1),bRegCorr(b2)),
            ('ht',                     'HT(jets) [GeV]',                          (100,0.,1000.))   : op.rng_sum(jets, lambda j : j.pt),
            ('min_dr_jets_lep1',       'Min(#Delta R(lead lepton,jets))',         (25,0.,5.))       : op.switch(l1conept(l1) >= l2conept(l2),
                                                                                                                self.HLL.MinDR_part1_partCont(l1,jets),
                                                                                                                self.HLL.MinDR_part1_partCont(l2,jets)),
            ('min_dr_jets_lep2',       'Min(#Delta R(sublead lepton,jets))',      (25,0.,5.))       : op.switch(l1conept(l1) >= l2conept(l2),
                                                                                                                self.HLL.MinDR_part1_partCont(l2,jets),
                                                                                                                self.HLL.MinDR_part1_partCont(l1,jets)),
            ('m_ll',                   'Dilepton invariant mass [GeV]',           (100,0.,1000.))   : op.invariant_mass(l1.p4,l2.p4),
            ('dr_ll',                  'Dilepton #Delta R',                       (25,0.,5.))       : op.deltaR(l1.p4,l2.p4),
            ('dphi_ll',                'Dilepton #Delta #Phi',                    (32,-3.2,3.2))    : op.abs(op.deltaPhi(l1.p4,l2.p4)),
            ('min_dr_jet',             'Min(#Delta R(jets))',                     (25,0.,5.))       : op.rng_min(dijets,lambda dijet : op.deltaR(dijet[0].p4,dijet[1].p4)),
            ('min_dphi_jet',           'Min(#Delta #Phi(jets))',                  (16,0.,3.2))      : rng_min(dijets,lambda dijet : op.abs(op.deltaPhi(dijet[0].p4,dijet[1].p4)),typeName='double'),
            ('m_hh_simplemet',         'M_{HH} (simple MET) [GeV]',               (100,0.,1000.))   : op.invariant_mass(b1.p4,b2.p4,l1.p4,l2.p4,met.p4),
            ('m_hh_simplemet_bregcorr','M_{HH} (simple MET) (regcorr) [GeV]',     (100,0.,1000.))   : op.invariant_mass(bRegCorr(b1),bRegCorr(b2),l1.p4,l2.p4,met.p4),
            ('met_ld',                 'MET_{LD}',                                (100,0.,1000.))   : self.HLL.MET_LD_DL(met,jets,electrons,muons),
            ('dr_bb',                  'Di-bjet #Delta R',                        (25,0.,5.))       : op.deltaR(b1.p4,b2.p4),
            ('dphi_bb',                'Di-bjet #Delta #Phi',                     (32,-3.2,3.2))    : op.abs(op.deltaPhi(b1.p4,b2.p4)), 
            ('min_dr_leps_b1',         'Min(#Delta R(lead bjet,dilepton))',       (25,0.,5.))       : self.HLL.MinDR_part1_dipart(b1,[l1,l2]),
            ('min_dr_leps_b2',         'Min(#Delta R(sublead bjet,dilepton))',    (25,0.,5.))       : self.HLL.MinDR_part1_dipart(b2,[l1,l2]),
            ('lep1_conept',            'Lead lepton cone-P_T [GeV]',              (40,0.,200.))     : op.switch(l1conept(l1) >= l2conept(l2) , l1conept(l1) , l2conept(l2)),
            ('lep2_conept',            'Sublead lepton cone-P_T [GeV]',           (40,0.,200.))     : op.switch(l1conept(l1) >= l2conept(l2) , l2conept(l2) , l1conept(l1)),
            ('mww_simplemet',          'M_{WW} (simple MET) [GeV]',               (100,0.,1000.))   : op.invariant_mass(l1.p4,l2.p4,met.p4),
            #('max_m_jj',               'Max(M_{jj}) [GeV]',                       (100,0.,1000.))   : op.invariant_mass(jets[0].p4,jets[1].p4),
            ('n_btag',                 'N_b',                                     (6,0.,5.))        : op.static_cast("UInt_t",op.rng_len(self.ak4BJets)),
            ('VBF_tag',                'VBF tag',                                 (2,0.,2.))        : VBF_tag}

def returnParamMVAInputs(self):
    return {('year',     'Year',         (3,2016.,2019.))    : op.c_int(int(self.era))}

def returnEventNrMVAInputs(self,t):
    return {('eventnr',  'Event number', (int(1e6+1),0.,1e6))     : op.c_int(t.event)}

