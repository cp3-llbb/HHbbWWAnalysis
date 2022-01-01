from bamboo import treefunctions as op
from highlevelLambdas import highlevelLambdas

def returnLeptonsMVAInputs(self,lep):
    return {('lep_E',        'Lead lepton E [GeV]',          (50,0.,500.))        : lep.p4.E(),
            ('lep_Px',       'Lead lepton P_x [GeV]',        (40,-200.,200.))     : lep.p4.Px(),
            ('lep_Py',       'Lead lepton P_y [GeV]',        (40,-200.,200.))     : lep.p4.Py(),
            ('lep_Pz',       'Lead lepton P_z [GeV]',        (40,-200.,200.))     : lep.p4.Pz(),
            ('lep_pdgId',    'Lead lepton pdg ID',           (45,-22.,22.) )      : lep.pdgId,
            ('lep_charge',   'Lead lepton charge',           (2,0.,2.)   )        : lep.charge}

def returnJetsMVAInputs(self, bjets, jets):
    return {('j1_E',        'Lead jet E [GeV]',             (50,0.,500.))        : op.switch(op.rng_len(bjets)>0,bjets[0].p4.E(),op.c_float(0.)),
            ('j1_Px',       'Lead jet P_x [GeV]',           (40,-200.,200.))     : op.switch(op.rng_len(bjets)>0,bjets[0].p4.Px(),op.c_float(0.)),
            ('j1_Py',       'Lead jet P_y [GeV]',           (40,-200.,200.))     : op.switch(op.rng_len(bjets)>0,bjets[0].p4.Py(),op.c_float(0.)),
            ('j1_Pz',       'Lead jet P_z [GeV]',           (40,-200.,200.))     : op.switch(op.rng_len(bjets)>0,bjets[0].p4.Pz(),op.c_float(0.)),
            ('j1_btag',     'Lead jet btag score',          (50,0.,1.)  )        : op.switch(op.rng_len(bjets)>0,bjets[0].btagDeepFlavB,op.c_float(0.)),
            ('j2_E',        'Sublead jet E [GeV]',          (50,0.,500.))        : op.switch(op.rng_len(bjets)>1,bjets[1].p4.E(),op.c_float(0.)),
            ('j2_Px',       'Sublead jet P_x [GeV]',        (40,-200.,200.))     : op.switch(op.rng_len(bjets)>1,bjets[1].p4.Px(),op.c_float(0.)),
            ('j2_Py',       'Sublead jet P_y [GeV]',        (40,-200.,200.))     : op.switch(op.rng_len(bjets)>1,bjets[1].p4.Py(),op.c_float(0.)),
            ('j2_Pz',       'Sublead jet P_z [GeV]',        (40,-200.,200.))     : op.switch(op.rng_len(bjets)>1,bjets[1].p4.Pz(),op.c_float(0.)),
            ('j2_btag',     'Sublead jet btag score',       (50,0.,1.)  )        : op.switch(op.rng_len(bjets)>1,bjets[1].btagDeepFlavB,op.c_float(0.)),
            ('j3_E',        'Subsublead jet E [GeV]',       (50,0.,500.))        : op.switch(op.rng_len(jets)>0,jets[0].p4.E(),op.c_float(0.)),
            ('j3_Px',       'Subsublead jet P_x [GeV]',     (40,-200.,200.))     : op.switch(op.rng_len(jets)>0,jets[0].p4.Px(),op.c_float(0.)),
            ('j3_Py',       'Subsublead jet P_y [GeV]',     (40,-200.,200.))     : op.switch(op.rng_len(jets)>0,jets[0].p4.Py(),op.c_float(0.)),
            ('j3_Pz',       'Subsublead jet P_z [GeV]',     (40,-200.,200.))     : op.switch(op.rng_len(jets)>0,jets[0].p4.Pz(),op.c_float(0.)),
            ('j3_btag',     'Subsublead jet btag score',    (50,0.,1.)  )        : op.switch(op.rng_len(jets)>0,jets[0].btagDeepFlavB,op.c_float(0.)),
            ('j4_E',        'Subsubsublead jet E [GeV]',    (50,0.,500.))        : op.switch(op.rng_len(jets)>1,jets[1].p4.E(),op.c_float(0.)),
            ('j4_Px',       'Subsubsublead jet P_x [GeV]',  (40,-200.,200.))     : op.switch(op.rng_len(jets)>1,jets[1].p4.Px(),op.c_float(0.)),
            ('j4_Py',       'Subsubsublead jet P_y [GeV]',  (40,-200.,200.))     : op.switch(op.rng_len(jets)>1,jets[1].p4.Py(),op.c_float(0.)),
            ('j4_Pz',       'Subsubsublead jet P_z [GeV]',  (40,-200.,200.))     : op.switch(op.rng_len(jets)>1,jets[1].p4.Pz(),op.c_float(0.)),
            ('j4_btag',     'Subsubsublead jet btag score', (50,0.,1.)  )        : op.switch(op.rng_len(jets)>1,jets[1].btagDeepFlavB,op.c_float(0.)),
            ('j5_E',        'Subsublead jet E [GeV]',       (50,0.,500.))        : op.switch(op.rng_len(jets)>2,jets[2].p4.E(),op.c_float(0.)),
            ('j5_Px',       'Subsublead jet P_x [GeV]',     (40,-200.,200.))     : op.switch(op.rng_len(jets)>2,jets[2].p4.Px(),op.c_float(0.)),
            ('j5_Py',       'Subsublead jet P_y [GeV]',     (40,-200.,200.))     : op.switch(op.rng_len(jets)>2,jets[2].p4.Py(),op.c_float(0.)),
            ('j5_Pz',       'Subsublead jet P_z [GeV]',     (40,-200.,200.))     : op.switch(op.rng_len(jets)>2,jets[2].p4.Pz(),op.c_float(0.)),
            ('j5_btag',     'Subsublead jet btag score',    (50,0.,1.)  )        : op.switch(op.rng_len(jets)>2,jets[2].btagDeepFlavB,op.c_float(0.)),
            ('j6_E',        'Subsublead jet E [GeV]',       (50,0.,500.))        : op.switch(op.rng_len(jets)>3,jets[3].p4.E(),op.c_float(0.)),
            ('j6_Px',       'Subsublead jet P_x [GeV]',     (40,-200.,200.))     : op.switch(op.rng_len(jets)>3,jets[3].p4.Px(),op.c_float(0.)),
            ('j6_Py',       'Subsublead jet P_y [GeV]',     (40,-200.,200.))     : op.switch(op.rng_len(jets)>3,jets[3].p4.Py(),op.c_float(0.)),
            ('j6_Pz',       'Subsublead jet P_z [GeV]',     (40,-200.,200.))     : op.switch(op.rng_len(jets)>3,jets[3].p4.Pz(),op.c_float(0.)),
            ('j6_btag',     'Subsublead jet btag score',    (50,0.,1.)  )        : op.switch(op.rng_len(jets)>3,jets[3].btagDeepFlavB,op.c_float(0.))}

# this is the leading fat bJet or normal fat jet?
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

def returnHighLevelMVAInputs(self, lep, bjets, wjets, VBFJetPairs, channel):
    if channel == 'El':
        lepconept  = self.electron_conept[lep.idx]
    elif channel == 'Mu':
        lepconept  = self.muon_conept[lep.idx]
    else:
        raise RuntimeError('Please mention the correct channel name!')

    import bamboo.treeoperations as _to
    def rng_min(rng, fun=(lambda x : x), typeName="float"):
        return op._to.Reduce.fromRngFun(rng, op.c_float(float("+inf"), typeName), ( lambda fn : (
            lambda res, elm : op.extMethod("std::min", returnType="Float_t")(res, fn(elm))
        ) )(fun) )

    
    return{
        ('m_hh_bregcorr',    'm_hh_bregcorr',     (50,0,400))  : self.HLL.comp_m_hh_bregcorr(bjets, wjets, lep, self.corrMET),
        ('pt_hh',            'pt_hh',             (50,0,400))  : self.HLL.comp_pt_hh(bjets, wjets, lep, self.corrMET),
        ('m_hbb_bregcorr',   'm_hbb_bregcorr',    (25,0,200))  : op.multiSwitch((op.rng_len(bjets) == 0, op.c_float(0.)),
                                                                                (op.rng_len(bjets) == 1, self.HLL.getCorrBp4(bjets[0]).M()),
                                                                                op.invariant_mass(self.HLL.getCorrBp4(bjets[0]),self.HLL.getCorrBp4(bjets[1]))),
        ('pt_hbb',           'pt_hbb',            (25,0,200))  : op.multiSwitch((op.rng_len(bjets) == 0, op.c_float(0.)),
                                                                                (op.rng_len(bjets) == 1, bjets[0].pt),
                                                                                (bjets[0].p4 + bjets[1].p4).Pt()),
        ('m_hww',            'm_hww',             (25,0,200))  : op.multiSwitch((op.rng_len(wjets) >= 2, op.invariant_mass(wjets[0].p4 ,wjets[1].p4, self.corrMET.p4, lep.p4)),
                                                                                (op.rng_len(wjets) == 1, op.invariant_mass(wjets[0].p4 ,self.corrMET.p4, lep.p4)),
                                                                                op.invariant_mass(self.corrMET.p4, lep.p4)),
        ('pt_hww',           'pt_hww',            (25,0,200))  : op.multiSwitch((op.rng_len(wjets) >= 2, (wjets[0].p4 + wjets[1].p4 + self.corrMET.p4 + lep.p4).Pt()),
                                                                                (op.rng_len(wjets) == 1, (wjets[0].p4 + self.corrMET.p4 + lep.p4).Pt()),
                                                                                (self.corrMET.p4 + lep.p4).Pt()),
        ('m_wjj',            'm_wjj',             (25,0,200))  : op.multiSwitch((op.rng_len(wjets) >= 2, op.invariant_mass(wjets[0].p4 ,wjets[1].p4)),
                                                                                (op.rng_len(wjets) == 1, wjets[0].mass),
                                                                                op.c_float(0.)),
        ('pt_wjj',           'pt_wjj',            (25,0,200))  : op.multiSwitch((op.rng_len(wjets) >= 2, (wjets[0].p4 + wjets[1].p4).Pt()),
                                                                                (op.rng_len(wjets) == 1, wjets[0].pt),
                                                                                op.c_float(0.)),
        ('m_wlep',           'm_wlep',            (25,0,200))  : op.invariant_mass(self.corrMET.p4, lep.p4),    
        ('pt_wlep',          'pt_wlep',           (25,0,200))  : (self.corrMET.p4 + lep.p4).Pt(),
        ('min_dr_lepbjets',  'min_dr_lepbjets',   (25,0,5))    : op.switch(op.rng_len(bjets) > 0, self.HLL.MinDR_part1_partCont(lep, bjets), op.c_float(0.)),
        ('dphi_hbb_hww',     'dphi_hbb_hww',      (22,0,3.2))  : self.HLL.comp_dphi_hbb_hww(bjets, wjets, lep, self.corrMET),
        ('dphi_hbb_hwwvis',  'dphi_hbb_hwwvis',   (22,0,3.2))  : self.HLL.comp_dphi_hbb_hwwvis(bjets, wjets, lep),
        ('dphi_met_lep',     'dphi_met_lep',      (22,0,3.2))  : op.abs(op.deltaPhi(self.corrMET.p4, lep.p4)),
        ('dphi_met_hbb',     'dphi_met_hbb',      (22,0,3.2))  : op.multiSwitch((op.rng_len(bjets) >= 2, op.abs(op.deltaPhi(self.corrMET.p4,(bjets[0].p4 + bjets[1].p4)))),
                                                                                (op.rng_len(bjets) == 1, op.abs(op.deltaPhi(self.corrMET.p4, bjets[0].p4))),
                                                                                op.c_float(0.)),
        ('dphi_met_wjj',     'dphi_met_wjj',      (22,0,3.2))  : op.multiSwitch((op.rng_len(wjets) >= 2, op.abs(op.deltaPhi(self.corrMET.p4,(wjets[0].p4 + wjets[1].p4)))),
                                                                                (op.rng_len(wjets) == 1, op.abs(op.deltaPhi(self.corrMET.p4,wjets[0].p4))),
                                                                                op.c_float(0.)),
        ('dr_lep_hbb',       'dr_lep_hbb',        (25,0,5))    : op.multiSwitch((op.rng_len(bjets) >= 2, op.deltaR(lep.p4, (bjets[0].p4+bjets[1].p4))), 
                                                                                (op.rng_len(bjets) == 1, op.deltaR(lep.p4, bjets[0].p4)),
                                                                                op.c_float(0.)),
        ('dr_lep_wjj',       'dr_lep_wjj',        (25,0,5))    : op.multiSwitch((op.rng_len(wjets) >= 2, op.deltaR(lep.p4, (wjets[0].p4+wjets[1].p4))),
                                                                                (op.rng_len(wjets) == 1, op.deltaR(lep.p4, wjets[0].p4)),
                                                                                op.c_float(0.)),
        ('min_dr_wjets',     'min_dr_wjets',      (25,0,5))    : op.switch(op.rng_len(wjets) >= 2, op.deltaR(wjets[0].p4, wjets[1].p4), op.c_float(0.)),
        ('min_dhi_wjets',    'min_dphi_wjets',    (22,0,3.2))  : op.switch(op.rng_len(wjets) >= 2, op.abs(op.deltaPhi(wjets[0].p4,wjets[1].p4)),op.c_float(0.)),
        ('min_dr_bjets',     'min_dr_bjets',      (25,0,5))    : op.switch(op.rng_len(bjets) >= 2, op.deltaR(bjets[0].p4,bjets[1].p4), op.c_float(0.)),
        ('min_dphi_bjets',   'min_dphi_bjets',    (22,0,3.2))  : op.switch(op.rng_len(bjets) >= 2, op.abs(op.deltaPhi(bjets[0].p4,bjets[1].p4)),op.c_float(0.)),
        ('ht',               'ht',                (50,0,300))  : op.rng_sum(self.ak4Jets, lambda j : j.pt),
        ('smin',             'smin',              (50,0,300))  : self.HLL.comp_smin(lep,self.corrMET,self.ak4Jets,bjets,wjets), 
        ('vbf_pair_mass',    'vbf_pair_mass',     (50,0,300))  : op.switch(op.rng_len(VBFJetPairs) > 0,op.invariant_mass(VBFJetPairs[0][0].p4,  VBFJetPairs[0][1].p4),
                                                                           op.c_float(0.)),
        ('vbf_pairs_absdeltaeta', 'vbf_pairs_absdeltaeta', (22,0,3.2)) : op.switch(op.rng_len(VBFJetPairs) > 0,op.abs(VBFJetPairs[0][0].eta - VBFJetPairs[0][1].eta),
                                                                                   op.c_float(0.)),
        ('lep_conept',       'lep_conept',        (25,0,200))  : lepconept,
        ('VBF_tag',          'vbf_tag',           (2,0,2))     : op.c_int(op.rng_len(VBFJetPairs) > 0),
        ('boosted_tag',      'boosted_tag',       (2,0,2))     : op.c_int(op.rng_len(self.ak8BJets) > 0),
        ('n_btag',           'n_btag',            (5,0,5))     : op.c_float(op.rng_len(self.ak4BJets))
    }
    

def returnParamMVAInputs(self):
    return {('year',     'Year',         (3,2016.,2019.))    : op.c_int(int(self.era))}

def returnEventNrMVAInputs(self,t):
    return {('eventnr',  'Event number', (int(1e6+1),0.,1e6))     : t.event}

def inputStaticCast(inputDict,cast='float'):
    return [op.static_cast(cast,v) for v in inputDict.values()]

'''
def returnEventNr(self,t):
    return {('evennr','Event number',(100,0.,1e6)): t.event}
    
def inputStaticCast(inputDict,cast='float'):
    return [op.static_cast(cast,v) for v in inputDict.values()]
'''
