from bamboo import treefunctions as op
from highlevelLambdas import highlevelLambdas
import random 
import sys
import os

# Get JPA Model Dict
# -- key : node
# -- val : model-path
def getResolvedJpaModelDict(basepath, nodelist, era) :
    resolvedModelDict = dict()
    resolvedModelDict[nodelist[0]] = [os.path.join(basepath, 'bb1l_jpa_4jet_resolved_even_all_mass_all_node_'+str(era)+'.xml'), 
                                      os.path.join(basepath, 'bb1l_jpa_4jet_resolved_odd_all_mass_all_node_'+str(era)+'.xml')]
    resolvedModelDict[nodelist[1]] = [os.path.join(basepath, 'bb1l_jpa_missingWJet_resolved_even_all_mass_all_node_'+str(era)+'.xml'), 
                                      os.path.join(basepath, 'bb1l_jpa_missingWJet_resolved_odd_all_mass_all_node_'+str(era)+'.xml')]
    resolvedModelDict[nodelist[2]] = [os.path.join(basepath, 'bb1l_jpa_missingBJet_even_all_mass_all_node_'+str(era)+'.xml'), 
                                      os.path.join(basepath, 'bb1l_jpa_missingBJet_odd_all_mass_all_node_'+str(era)+'.xml')]
    resolvedModelDict[nodelist[3]] = [os.path.join(basepath, 'bb1l_jpa_missingAllWJet_even_all_mass_all_node_'+str(era)+'.xml'), 
                                      os.path.join(basepath, 'bb1l_jpa_missingAllWJet_odd_all_mass_all_node_'+str(era)+'.xml')]
    resolvedModelDict[nodelist[4]] = [os.path.join(basepath, 'bb1l_jpa_missingBJet_missingWJet_even_all_mass_all_node_'+str(era)+'.xml'), 
                                      os.path.join(basepath, 'bb1l_jpa_missingBJet_missingWJet_odd_all_mass_all_node_'+str(era)+'.xml')]
    resolvedModelDict[nodelist[5]] = [os.path.join(basepath, 'bb1l_jpa_missingBJet_missingAllWJet_even_all_mass_all_node_'+str(era)+'.xml'), 
                                      os.path.join(basepath, 'bb1l_jpa_missingBJet_missingAllWJet_odd_all_mass_all_node_'+str(era)+'.xml')]
    resolvedModelDict['evCat']     = [os.path.join(basepath, 'bb1l_jpa_evtCat_resolved_even_all_mass_all_node_'+str(era)+'.xml'), 
                                      os.path.join(basepath, 'bb1l_jpa_evtCat_resolved_odd_all_mass_all_node_'+str(era)+'.xml')]
    return resolvedModelDict

def getBoostedJpaModelDict(basepath, nodelist, era):
    boostedModelDict  = dict()
    boostedModelDict[nodelist[0]] = [os.path.join(basepath, 'bb1l_jpa_4jet_boosted_even_all_mass_all_node_'+str(era)+'.xml'), 
                                     os.path.join(basepath, 'bb1l_jpa_4jet_boosted_odd_all_mass_all_node_'+str(era)+'.xml')]
    boostedModelDict[nodelist[1]] = [os.path.join(basepath, 'bb1l_jpa_missingWJet_boosted_even_all_mass_all_node_'+str(era)+'.xml'), 
                                     os.path.join(basepath, 'bb1l_jpa_missingWJet_boosted_odd_all_mass_all_node_'+str(era)+'.xml')]
    boostedModelDict['evCat']     = [os.path.join(basepath, 'bb1l_jpa_evtCat_boosted_even_all_mass_all_node_'+str(era)+'.xml'), 
                                     os.path.join(basepath, 'bb1l_jpa_evtCat_boosted_odd_all_mass_all_node_'+str(era)+'.xml')]
    return boostedModelDict

# Select VBF jets 
# To use in both Skimmer and Plotter
def cleanVBFwithJPA_Resolved(jpaJets, nJpaJets):
    return lambda j : op.AND(*(op.deltaR(jpaJets[i].p4, j.p4) > 0.8 for i in range(nJpaJets)))

def cleanVBFwithJPA_Boosted(fatjets, jpaJets, nJpaJets):
    return lambda j : op.AND(op.rng_len(fatjets) >= 1, op.AND(op.AND(*(op.deltaR(jpaJets[i].p4, j.p4) > 0.8 for i in range(nJpaJets))), 
                                                              op.deltaR(fatjets[0].p4, j.p4) > 1.2))
def cleanVBFwithJPA_Boosted_rest(fatjets):
    return lambda j : op.AND(op.rng_len(fatjets) >= 1, op.deltaR(fatjets[0].p4, j.p4) > 1.2)


# Format : Dict with key = (variable name, title for plot, binning for plot)
#                    val = bamboo variable
### --------------------------------------------------------------------------------------------------------- ###
###                                                  Resolved                                                 ###
### --------------------------------------------------------------------------------------------------------- ###
def returnCommonInputs_Resolved(self):
    return{('METpx',                   'METpx',                   (50, 0, 300)):    self.corrMET.p4.Px(),
           ('METpy',                   'METpy',                   (50, 0, 300)):    self.corrMET.p4.Py(),
           ('METpz',                   'METpz',                   (50, 0, 300)):    self.corrMET.p4.Pz(),
           ('METenergy',               'METenergy',               (50, 0, 300)):    self.corrMET.p4.E(),
           ('METpt',                   'METpt',                   (50, 0, 300)):    self.corrMET.pt,
           ('METphi',                  'METphi',                  (50, 0, 300)):    self.corrMET.phi,
           ('MET_LD',                  'MET_LD',                  (40,0,200)):      self.HLL.MET_LD_DL(self.corrMET, self.ak4Jets, self.electronsFakeSel, self.muonsFakeSel),
           ('nfakeableMuons',          'nfakeableMuons',          (5,0,5)):         op.static_cast("UInt_t",op.rng_len(self.muonsFakeSel)),
           ('nfakeableElectrons',      'nfakeableElectrons',      (5,0,5)):         op.static_cast("UInt_t",op.rng_len(self.electronsFakeSel)),
           ('hT',                      'hT',                      (50,0,300)):      op.static_cast("Float_t", self.HLL.HT_SL(self.ak4Jets)),
           ('nAk4Jets',                'nAk4Jets',                (10,0,10)):       op.static_cast("UInt_t",op.rng_len(self.ak4Jets)),
           ('nAk4BJets',               'nAk4BJets',               (10,0,10)):       op.static_cast("UInt_t",op.rng_len(self.ak4BJets)),
           ('nAk4BJetsLoose',          'nAk4BJetsLoose',          (10,0,10)):       op.static_cast("UInt_t",op.rng_len(self.ak4BJetsLoose)),
           ('nAk8Jets',                'nAk8Jets',                (5,0,5)):         op.static_cast("UInt_t",op.rng_len(self.ak8Jets)),
           ('nAk8BJets',               'nAk8BJets',               (5,0,5)):         op.static_cast("UInt_t",op.rng_len(self.ak8BJets)),
           ('nAk4JetsCleandWithAk8b',  'nAk4JetsCleandWithAk8b',  (5,0,5)):         op.static_cast("UInt_t",op.rng_len(self.ak4JetsCleanedFromAk8b))
       }

def returnClassicInputs_Resolved(self, lepton, jpaSelectedJets, L1out, L2out, jpaArg):
    ##########################################################################################################
    # Res2b2Wj : jet1   jet2     jet3     jet4
    # Res2b1Wj : jet1   jet2     jet3     jet4=0
    # Res2b0Wj : jet1   jet2     jet3=0   jet4=0
    # Res1b2Wj : jet1   jet2=0   jet3     jet4
    # Res1b1Wj : jet1   jet2=0   jet3     jet4=0
    # Res1b0Wj : jet1   jet2=0   jet3=0   jet4=0
    # Res0b    : all jet = 0 because, this is not a JPA category. Only lepton and MET variables will be saved.
    ##########################################################################################################
    Res2b2Wj = False
    Res2b1Wj = False
    Res1b2Wj = False
    Res2b0Wj = False
    Res1b1Wj = False
    Res1b0Wj = False
    Res0b    = False
    if jpaArg == 'Res2b2Wj':
        Res2b2Wj = True
        jet1 = jpaSelectedJets[0]
        jet2 = jpaSelectedJets[1]
        jet3 = jpaSelectedJets[2]
        jet4 = jpaSelectedJets[3]
        VBFJetPairsJPA = op.sort(op.combine(op.select(self.VBFJets, cleanVBFwithJPA_Resolved(jpaSelectedJets, 4)), N=2, pred=self.lambda_VBFPair), 
                                 lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))
    elif jpaArg == 'Res2b1Wj':
        Res2b1Wj = True
        jet1 = jpaSelectedJets[0]
        jet2 = jpaSelectedJets[1]
        jet3 = jpaSelectedJets[2]
        jet4 = None
        VBFJetPairsJPA = op.sort(op.combine(op.select(self.VBFJets, cleanVBFwithJPA_Resolved(jpaSelectedJets, 3)), N=2, pred=self.lambda_VBFPair), 
                                 lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))
    elif jpaArg == 'Res1b2Wj':
        Res1b2Wj = True
        jet1 = jpaSelectedJets[0]
        jet2 = None
        jet3 = jpaSelectedJets[1]
        jet4 = jpaSelectedJets[2]
        VBFJetPairsJPA = op.sort(op.combine(op.select(self.VBFJets, cleanVBFwithJPA_Resolved(jpaSelectedJets, 3)), N=2, pred=self.lambda_VBFPair), 
                                 lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))
    elif jpaArg == 'Res2b0Wj':
        Res2b0Wj = True
        jet1 = jpaSelectedJets[0]
        jet2 = jpaSelectedJets[1]
        jet3 = None
        jet4 = None
        VBFJetPairsJPA = op.sort(op.combine(op.select(self.VBFJets, cleanVBFwithJPA_Resolved(jpaSelectedJets, 2)), N=2, pred=self.lambda_VBFPair), 
                                 lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))
    elif jpaArg == 'Res1b1Wj':
        Res1b1Wj = True
        jet1 = jpaSelectedJets[0]
        jet2 = None
        jet3 = jpaSelectedJets[1]
        jet4 = None
        VBFJetPairsJPA = op.sort(op.combine(op.select(self.VBFJets, cleanVBFwithJPA_Resolved(jpaSelectedJets, 2)), N=2, pred=self.lambda_VBFPair), 
                                 lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))
    elif jpaArg == 'Res1b0Wj':
        Res1b0Wj = True
        jet1 = jpaSelectedJets[0]
        jet2 = None
        jet3 = None
        jet4 = None
        VBFJetPairsJPA = op.sort(op.combine(op.select(self.VBFJets, cleanVBFwithJPA_Resolved(jpaSelectedJets, 1)), N=2, pred=self.lambda_VBFPair), 
                                 lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))
    elif jpaArg == 'Res0b':
        Res0b = True
        jet1 = None
        jet2 = None
        jet3 = None
        jet4 = None        
        VBFJetPairsJPA = None
    else:
        raise RuntimeError('Resolved JPA category is not mentioned!')

    # 4 bools to make this easy
    has2bj = any([Res2b2Wj, Res2b1Wj, Res2b0Wj])
    has1bj = any([Res1b2Wj, Res1b1Wj, Res1b0Wj])
    has2wj = any([Res2b2Wj, Res1b2Wj])
    has1wj = any([Res2b1Wj, Res1b1Wj])
    # Jet-1 variables :: Available for each category
    bj1_Px             = self.HLL.bJetCorrP4(jet1).Px()                       if not Res0b else op.c_float(0.)
    bj1_Py             = self.HLL.bJetCorrP4(jet1).Py()                       if not Res0b else op.c_float(0.)
    bj1_Pz             = self.HLL.bJetCorrP4(jet1).Pz()                       if not Res0b else op.c_float(0.)
    bj1_E              = self.HLL.bJetCorrP4(jet1).E()                        if not Res0b else op.c_float(0.)
    bj1_M              = self.HLL.bJetCorrP4(jet1).M()                        if not Res0b else op.c_float(0.)
    bj1_pt             = self.HLL.bJetCorrP4(jet1).Pt()                       if not Res0b else op.c_float(0.)
    bj1_eta            = self.HLL.bJetCorrP4(jet1).Eta()                      if not Res0b else op.c_float(0.)
    bj1_phi            = self.HLL.bJetCorrP4(jet1).Phi()                      if not Res0b else op.c_float(0.)
    bj1_bTagDeepFlavB  = jet1.btagDeepFlavB                                   if not Res0b else op.c_float(0.)
    bj1LepDR           = op.deltaR(jet1.p4, lepton.p4)                        if not Res0b else op.c_float(0.)
    bj1LepDEta         = op.abs(jet1.eta - lepton.eta)                        if not Res0b else op.c_float(0.)
    bj1LepDPhi         = op.abs(op.deltaPhi(jet1.p4, lepton.p4))              if not Res0b else op.c_float(0.)  
    bj1MetDPhi         = op.abs(self.HLL.SinglepMet_dPhi(jet1, self.corrMET)) if not Res0b else op.c_float(0.)

    # Jet-2 variables :: Available for 2b2Wj, 2b1Wj. 2b0Wj
    bj2_Px            = self.HLL.bJetCorrP4(jet2).Px()                        if has2bj else op.c_float(0.)
    bj2_Py            = self.HLL.bJetCorrP4(jet2).Py()                        if has2bj else op.c_float(0.)
    bj2_Pz            = self.HLL.bJetCorrP4(jet2).Pz()                        if has2bj else op.c_float(0.)
    bj2_E             = self.HLL.bJetCorrP4(jet2).E()                         if has2bj else op.c_float(0.)
    bj2_M             = self.HLL.bJetCorrP4(jet2).M()                         if has2bj else op.c_float(0.)
    bj2_pt            = self.HLL.bJetCorrP4(jet2).Pt()                        if has2bj else op.c_float(0.)
    bj2_eta           = self.HLL.bJetCorrP4(jet2).Eta()                       if has2bj else op.c_float(0.)
    bj2_phi           = self.HLL.bJetCorrP4(jet2).Phi()                       if has2bj else op.c_float(0.)
    bj2_bTagDeepFlavB = jet2.btagDeepFlavB                                    if has2bj else op.c_float(0.)
    bj2LepDR          = op.deltaR(jet2.p4, lepton.p4)                         if has2bj else op.c_float(0.)
    bj2LepDEta        = op.abs(jet2.eta - lepton.eta)                         if has2bj else op.c_float(0.)
    bj2LepDPhi        = op.abs(op.deltaPhi(jet2.p4, lepton.p4))               if has2bj else op.c_float(0.)
    bj2MetDPhi        = op.abs(self.HLL.SinglepMet_dPhi(jet2, self.corrMET))  if has2bj else op.c_float(0.)

    # h->bb related
    HbbLepDR          = op.deltaR(jet1.p4+jet2.p4, lepton.p4)                                        if has2bj else op.c_float(0.) 
    HbbLepDEta        = op.abs((jet1.p4+jet2.p4).Eta() - lepton.eta)                                 if has2bj else op.c_float(0.) 
    HbbLepDPhi        = op.abs(op.deltaPhi(jet1.p4+jet2.p4, lepton.p4))                              if has2bj else op.c_float(0.) 
    HbbMetDPhi        = op.abs(op.deltaPhi(jet1.p4+jet2.p4, self.corrMET.p4))                        if has2bj else op.c_float(0.) 
    cosThetaS_Hbb     = op.extMethod("HHbbWWJPA::cosThetaS", returnType="float")(jet1.p4, jet2.p4)   if has2bj else op.c_float(0.)

    # Jet-3 variables :: Available for 2b2Wj, 2b1Wj, 1b2Wj, 1b1Wj
    wj1_Px            = jet3.p4.Px()                                  if has1wj or has2wj else op.c_float(0.)
    wj1_Py            = jet3.p4.Py()                                  if has1wj or has2wj else op.c_float(0.)
    wj1_Pz            = jet3.p4.Pz()                                  if has1wj or has2wj else op.c_float(0.)
    wj1_E             = jet3.p4.E()                                   if has1wj or has2wj else op.c_float(0.)
    wj1_M             = jet3.mass                                     if has1wj or has2wj else op.c_float(0.)
    wj1_pt            = jet3.pt                                       if has1wj or has2wj else op.c_float(0.)
    wj1_eta           = jet3.eta                                      if has1wj or has2wj else op.c_float(0.)
    wj1_phi           = jet3.phi                                      if has1wj or has2wj else op.c_float(0.)
    wj1_bTagDeepFlavB = jet3.btagDeepFlavB                            if has1wj or has2wj else op.c_float(0.)
    wj1LepDR          = op.deltaR(jet3.p4, lepton.p4)                 if has1wj or has2wj else op.c_float(0.) 
    wj1LepDEta        = op.abs(jet3.eta - lepton.eta)                 if has1wj or has2wj else op.c_float(0.)
    wj1LepDPhi        = op.abs(op.deltaPhi(jet3.p4, lepton.p4))       if has1wj or has2wj else op.c_float(0.)
    wj1MetDPhi        = op.abs(self.HLL.SinglepMet_dPhi(jet3, self.corrMET))   if has1wj or has2wj else op.c_float(0.)
    
    # Jet-4 variables :: Available for 2b2Wj, 1b2Wj
    wj2_Px            = jet4.p4.Px()                                  if has2wj else op.c_float(0.)
    wj2_Py            = jet4.p4.Py()                                  if has2wj else op.c_float(0.)
    wj2_Pz            = jet4.p4.Pz()                                  if has2wj else op.c_float(0.)
    wj2_E             = jet4.p4.E()                                   if has2wj else op.c_float(0.)
    wj2_M             = jet4.mass                                     if has2wj else op.c_float(0.)
    wj2_pt            = jet4.pt                                       if has2wj else op.c_float(0.)
    wj2_eta           = jet4.eta                                      if has2wj else op.c_float(0.)
    wj2_phi           = jet4.phi                                      if has2wj else op.c_float(0.)
    wj2_bTagDeepFlavB = jet4.btagDeepFlavB                            if has2wj else op.c_float(0.)
    wj2LepDR          = op.deltaR(jet4.p4, lepton.p4)                 if has2wj else op.c_float(0.) 
    wj2LepDEta        = op.abs(jet4.eta - lepton.eta)                 if has2wj else op.c_float(0.)
    wj2LepDPhi        = op.abs(op.deltaPhi(jet4.p4, lepton.p4))       if has2wj else op.c_float(0.)
    wj2MetDPhi        = op.abs(self.HLL.SinglepMet_dPhi(jet4, self.corrMET))   if has2wj else op.c_float(0.)

    WjjLepDR          = op.deltaR(jet3.p4+jet4.p4, lepton.p4)                         if has2wj else op.c_float(0.) ###
    WjjLepDEta        = op.abs((jet3.p4+jet4.p4).Eta() - lepton.eta)                   if has2wj else op.c_float(0.) ###
    WjjLepDPhi        = op.abs(op.deltaPhi(jet3.p4+jet4.p4, lepton.p4))               if has2wj else op.c_float(0.) ###
    WjjMetDPhi        = op.abs(op.deltaPhi(jet3.p4+jet4.p4, self.corrMET.p4))         if has2wj else op.c_float(0.) ###

    # h->ww related
    HWW_Mass          = self.HLL.HWW_simple(jet3.p4,jet4.p4,lepton.p4,self.corrMET).M()        if has2wj else op.c_float(0.) # w.r.t. neutrino
    HWW_Simple_Mass   = self.HLL.HWW_met_simple(jet3.p4,jet4.p4,lepton.p4,self.corrMET.p4).M() if has2wj else op.c_float(0.) # not w.r.t. neutrino
    HWW_dR            = self.HLL.dR_Hww(jet3.p4,jet4.p4,lepton.p4,self.corrMET)                if has2wj else op.c_float(0.)
    HWW_dEta          = self.HLL.dEta_Hww(jet3.p4,jet4.p4,lepton.p4,self.corrMET)              if has2wj else op.c_float(0.)
    HWW_dPhi          = self.HLL.dPhi_Hww(jet3.p4,jet4.p4,lepton.p4,self.corrMET)              if has2wj else op.c_float(0.)
    angleBetWWPlane   = self.HLL.angleWWplane(lepton.p4, self.corrMET, jet3.p4, jet4.p4)       if has2wj else op.c_float(0.)
    
    cosThetaS_Wjj_simple    = op.extMethod("HHbbWWJPA::cosThetaS", returnType="float")(jet3.p4, jet4.p4)                                                   if has2wj else op.c_float(0.)
    cosThetaS_WW_simple_met = op.extMethod("HHbbWWJPA::cosThetaS", returnType="float")(self.HLL.Wjj_simple(jet3.p4,jet4.p4), 
                                                                                       self.HLL.Wlep_met_simple(lepton.p4, self.corrMET.p4))               if has2wj else op.c_float(0.)
    cosThetaS_HH_simple_met = op.extMethod("HHbbWWJPA::cosThetaS", returnType="float")(self.HLL.bJetCorrP4(jet1)+self.HLL.bJetCorrP4(jet2), 
                                                                                       self.HLL.HWW_met_simple(jet3.p4,jet4.p4,lepton.p4,self.corrMET.p4)) if Res2b2Wj else op.c_float(0.)    
    mHH                     = (jet3.p4+jet4.p4+lepton.p4+self.corrMET.p4+self.HLL.bJetCorrP4(jet1)+self.HLL.bJetCorrP4(jet2)).M()                          if Res2b2Wj else op.c_float(0.)
    angleBetHWPlane         = self.HLL.angleBetPlanes(jet1.p4,jet2.p4,jet3.p4,jet4.p4)                                                                     if Res2b2Wj else op.c_float(0.)
    
    # neu p4 only for 1b2Wj and 2b2Wj
    neuPx = self.HLL.neuP4(jet3.p4+jet4.p4+lepton.p4, self.corrMET).Px()  if has2wj else op.c_float(0.)
    neuPy = self.HLL.neuP4(jet3.p4+jet4.p4+lepton.p4, self.corrMET).Py()  if has2wj else op.c_float(0.)
    neuPz = self.HLL.neuP4(jet3.p4+jet4.p4+lepton.p4, self.corrMET).Pz()  if has2wj else op.c_float(0.)
    neuE  = self.HLL.neuP4(jet3.p4+jet4.p4+lepton.p4, self.corrMET).E()   if has2wj else op.c_float(0.)
    neuPt = self.HLL.neuP4(jet3.p4+jet4.p4+lepton.p4, self.corrMET).Pt()  if has2wj else op.c_float(0.)

    
    # diJet dR & dPhi
    bj1bj2_DR   = op.deltaR(jet1.p4,jet2.p4)             if has2bj else op.c_float(0.)
    bj1bj2_DPhi = op.abs(op.deltaPhi(jet1.p4,jet2.p4))   if has2bj else op.c_float(0.)
    bj1bj2_DEta = op.abs(jet1.eta - jet2.eta)            if has2bj else op.c_float(0.)

    bj1wj1_DR   = op.deltaR(jet1.p4,jet3.p4)             if has1wj or has2wj else op.c_float(0.)
    bj1wj1_DPhi = op.abs(op.deltaPhi(jet1.p4,jet3.p4))   if has1wj or has2wj else op.c_float(0.)
    bj1wj1_DEta = op.abs(jet1.eta - jet3.eta)            if has1wj or has2wj else op.c_float(0.)

    bj1wj2_DR   = op.deltaR(jet1.p4,jet4.p4)             if any([Res2b2Wj, Res1b2Wj]) else op.c_float(0.)
    bj1wj2_DPhi = op.abs(op.deltaPhi(jet1.p4,jet4.p4))   if any([Res2b2Wj, Res1b2Wj]) else op.c_float(0.)
    bj1wj2_DEta = op.abs(jet1.eta - jet4.eta)            if any([Res2b2Wj, Res1b2Wj]) else op.c_float(0.)

    bj2wj1_DR   = op.deltaR(jet2.p4,jet3.p4)             if any([Res2b2Wj, Res2b1Wj]) else op.c_float(0.)
    bj2wj1_DPhi = op.abs(op.deltaPhi(jet2.p4,jet3.p4))   if any([Res2b2Wj, Res2b1Wj]) else op.c_float(0.)
    bj2wj1_DEta = op.abs(jet2.eta - jet3.eta)            if any([Res2b2Wj, Res2b1Wj]) else op.c_float(0.)

    bj2wj2_DR   = op.deltaR(jet2.p4,jet4.p4)             if Res2b2Wj else op.c_float(0.)
    bj2wj2_DPhi = op.abs(op.deltaPhi(jet2.p4,jet4.p4))   if Res2b2Wj else op.c_float(0.)
    bj2wj2_DEta = op.abs(jet2.eta - jet4.eta)            if Res2b2Wj else op.c_float(0.)

    wj1wj2_DR   = op.deltaR(jet3.p4,jet4.p4)             if has2wj else op.c_float(0.)
    wj1wj2_DPhi = op.abs(op.deltaPhi(jet3.p4,jet4.p4))   if has2wj else op.c_float(0.)
    wj1wj2_DEta = op.abs(jet3.eta - jet4.eta)            if has2wj else op.c_float(0.)
    
    # VBF variables
    VBFj1Px        = op.switch(op.rng_len(VBFJetPairsJPA) > 0, 
                               op.switch(VBFJetPairsJPA[0][0].p4.Px() > VBFJetPairsJPA[0][1].p4.Px(),
                                         VBFJetPairsJPA[0][0].p4.Px(), VBFJetPairsJPA[0][1].p4.Px()),
                               op.c_float(0.)) if not Res0b else op.c_float(0.)
    VBFj2Px        = op.switch(op.rng_len(VBFJetPairsJPA) > 0, 
                               op.switch(VBFJetPairsJPA[0][0].p4.Px() > VBFJetPairsJPA[0][1].p4.Px(),
                                         VBFJetPairsJPA[0][1].p4.Px(), VBFJetPairsJPA[0][0].p4.Px()),
                               op.c_float(0.)) if not Res0b else op.c_float(0.)
    VBFj1Py        = op.switch(op.rng_len(VBFJetPairsJPA) > 0, 
                               op.switch(VBFJetPairsJPA[0][0].p4.Py() > VBFJetPairsJPA[0][1].p4.Py(),
                                         VBFJetPairsJPA[0][0].p4.Py(), VBFJetPairsJPA[0][1].p4.Py()),
                               op.c_float(0.)) if not Res0b else op.c_float(0.)
    VBFj2Py        = op.switch(op.rng_len(VBFJetPairsJPA) > 0, 
                               op.switch(VBFJetPairsJPA[0][0].p4.Py() > VBFJetPairsJPA[0][1].p4.Py(),
                                         VBFJetPairsJPA[0][1].p4.Py(), VBFJetPairsJPA[0][0].p4.Py()),
                               op.c_float(0.)) if not Res0b else op.c_float(0.)
    VBFj1Pz        = op.switch(op.rng_len(VBFJetPairsJPA) > 0, 
                               op.switch(VBFJetPairsJPA[0][0].p4.Pz() > VBFJetPairsJPA[0][1].p4.Pz(),
                                         VBFJetPairsJPA[0][0].p4.Pz(), VBFJetPairsJPA[0][1].p4.Pz()),
                               op.c_float(0.)) if not Res0b else op.c_float(0.)
    VBFj2Pz        = op.switch(op.rng_len(VBFJetPairsJPA) > 0, 
                               op.switch(VBFJetPairsJPA[0][0].p4.Pz() > VBFJetPairsJPA[0][1].p4.Pz(),
                                         VBFJetPairsJPA[0][1].p4.Pz(), VBFJetPairsJPA[0][0].p4.Pz()),
                               op.c_float(0.)) if not Res0b else op.c_float(0.)
    VBFj1E         = op.switch(op.rng_len(VBFJetPairsJPA) > 0, 
                               op.switch(VBFJetPairsJPA[0][0].p4.E() > VBFJetPairsJPA[0][1].p4.E(),
                                         VBFJetPairsJPA[0][0].p4.E(), VBFJetPairsJPA[0][1].p4.E()),
                               op.c_float(0.)) if not Res0b else op.c_float(0.)
    VBFj2E         = op.switch(op.rng_len(VBFJetPairsJPA) > 0, 
                               op.switch(VBFJetPairsJPA[0][0].p4.E() > VBFJetPairsJPA[0][1].p4.E(),
                                         VBFJetPairsJPA[0][1].p4.E(), VBFJetPairsJPA[0][0].p4.E()),
                               op.c_float(0.)) if not Res0b else op.c_float(0.)
    VBFj1pt        = op.switch(op.rng_len(VBFJetPairsJPA) > 0, 
                                 op.switch(VBFJetPairsJPA[0][0].pt > VBFJetPairsJPA[0][1].pt,
                                           VBFJetPairsJPA[0][0].pt, VBFJetPairsJPA[0][1].pt),
                                 op.c_float(0.)) if not Res0b else op.c_float(0.)
    VBFj2pt        = op.switch(op.rng_len(VBFJetPairsJPA) > 0, 
                               op.switch(VBFJetPairsJPA[0][0].pt > VBFJetPairsJPA[0][1].pt,
                                         VBFJetPairsJPA[0][1].pt, VBFJetPairsJPA[0][0].pt),
                               op.c_float(0.)) if not Res0b else op.c_float(0.)
    VBFj1eta       = op.switch(op.rng_len(VBFJetPairsJPA) > 0, 
                               op.switch(VBFJetPairsJPA[0][0].pt > VBFJetPairsJPA[0][1].pt,
                                         VBFJetPairsJPA[0][0].eta, VBFJetPairsJPA[0][1].eta),
                               op.c_float(0.)) if not Res0b else op.c_float(0.)
    VBFj2eta       = op.switch(op.rng_len(VBFJetPairsJPA) > 0, 
                               op.switch(VBFJetPairsJPA[0][0].pt > VBFJetPairsJPA[0][1].pt,
                                         VBFJetPairsJPA[0][1].eta, VBFJetPairsJPA[0][0].eta),
                               op.c_float(0.)) if not Res0b else op.c_float(0.)
    VBFj1j2dEta    = op.switch(op.rng_len(VBFJetPairsJPA) > 0,op.abs(VBFJetPairsJPA[0][0].eta - VBFJetPairsJPA[0][1].eta),
                               op.c_float(0.)) if not Res0b else op.c_float(0.)
    
    VBFj1j2dPhi    = op.switch(op.rng_len(VBFJetPairsJPA) > 0,op.abs(op.deltaPhi(VBFJetPairsJPA[0][0].p4,VBFJetPairsJPA[0][1].p4)),
                               op.c_float(0.)) if not Res0b else op.c_float(0.)
    
    VBFj1j2invM    = op.switch(op.rng_len(VBFJetPairsJPA) > 0,op.invariant_mass(VBFJetPairsJPA[0][0].p4,  VBFJetPairsJPA[0][1].p4),
                               op.c_float(0.)) if not Res0b else op.c_float(0.)
    
    VBF_tag        = op.c_int(op.rng_len(VBFJetPairsJPA) > 0) if not Res0b else op.c_int(0)
    zeppenfeldVar  = op.switch(op.rng_len(VBFJetPairsJPA) > 0, 
                               op.abs(lepton.eta - (VBFJetPairsJPA[0][0].eta + VBFJetPairsJPA[0][1].eta)/2.0)/op.abs(VBFJetPairsJPA[0][0].eta - VBFJetPairsJPA[0][1].eta),
                               op.c_float(0.)) if not Res0b else op.c_float(0.)
    
    if Res2b2Wj : 
        pt_taggedJetsAvg = (self.HLL.bJetCorrP4(jet1).Pt()+self.HLL.bJetCorrP4(jet2).Pt()+jet3.pt+jet4.pt)/4.0
        mass_taggedJetsAvg = (self.HLL.bJetCorrP4(jet1).M()+self.HLL.bJetCorrP4(jet2).M()+jet3.mass+jet4.mass)/4.0
        bj1bj2_pt      = (self.HLL.bJetCorrP4(jet1)+self.HLL.bJetCorrP4(jet2)).Pt()
        bj1bj2_M       = op.invariant_mass(self.HLL.bJetCorrP4(jet1),self.HLL.bJetCorrP4(jet2))
        mT_top_3particle  = op.min(self.HLL.mT2(self.HLL.bJetCorrP4(jet1),lepton.p4 ,self.corrMET.p4), 
                                   self.HLL.mT2(self.HLL.bJetCorrP4(jet2), lepton.p4, self.corrMET.p4))
        avgJetDR       = (op.deltaR(jet1.p4,jet2.p4) + op.deltaR(jet1.p4,jet3.p4) + op.deltaR(jet1.p4,jet4.p4) + op.deltaR(jet2.p4,jet3.p4) + op.deltaR(jet2.p4,jet4.p4) + op.deltaR(jet3.p4,jet4.p4))/6.0
        avgJetDPhi     = (op.abs(op.deltaPhi(jet1.p4,jet2.p4)) + op.abs(op.deltaPhi(jet1.p4,jet3.p4)) + op.abs(op.deltaPhi(jet1.p4,jet4.p4)) + op.abs(op.deltaPhi(jet2.p4,jet3.p4)) + op.abs(op.deltaPhi(jet2.p4,jet4.p4)) + op.abs(op.deltaPhi(jet3.p4,jet4.p4)))/6.0
        avgJetDEta     = (op.abs(jet1.eta-jet2.eta) + op.abs(jet1.eta-jet3.eta) + op.abs(jet1.eta-jet4.eta) + op.abs(jet2.eta-jet3.eta) + op.abs(jet2.eta-jet4.eta) + op.abs(jet3.eta-jet4.eta))/6.0
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
        JPAcat         = op.c_int(1)
        
    if Res2b1Wj : 
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
        JPAcat         = op.c_int(2)
        
    if Res1b2Wj : 
        pt_taggedJetsAvg = (self.HLL.bJetCorrP4(jet1).Pt()+jet3.pt+jet4.pt)/3.0
        mass_taggedJetsAvg = (self.HLL.bJetCorrP4(jet1).M()+jet3.mass+jet4.mass)/3.0
        bj1bj2_pt      = self.HLL.bJetCorrP4(jet1).Pt()
        bj1bj2_M       = self.HLL.bJetCorrP4(jet1).M()
        mT_top_3particle  = self.HLL.mT2(self.HLL.bJetCorrP4(jet1),lepton.p4 ,self.corrMET.p4)
        avgJetDR       = (op.deltaR(jet1.p4,jet3.p4) + op.deltaR(jet1.p4,jet4.p4) + op.deltaR(jet3.p4,jet4.p4))/3.0
        avgJetDPhi     = (op.abs(op.deltaPhi(jet1.p4,jet3.p4)) + op.abs(op.deltaPhi(jet1.p4,jet4.p4)) + op.abs(op.deltaPhi(jet3.p4,jet4.p4)))/3.0
        avgJetDEta     = (op.abs(jet1.eta-jet3.eta) + op.abs(jet1.eta-jet4.eta) + op.abs(jet3.eta-jet4.eta))/3.0
        minJetDR       = self.HLL.MinDiJetDRLoose(jet1,jet3,jet4)
        minJetDPhi     = self.HLL.MinDiJetDPhiLoose(jet1,jet3,jet4)
        minJetDEta     = self.HLL.MinDiJetDEtaLoose(jet1,jet3,jet4)
        maxJetDR       = self.HLL.MaxDiJetDRLoose(jet1,jet3,jet4)
        maxJetDPhi     = self.HLL.MaxDiJetDPhiLoose(jet1,jet3,jet4)
        maxJetDEta     = self.HLL.MaxDiJetDEtaLoose(jet1,jet3,jet4)
        minLepJetDR    = self.HLL.MinDR_lep3j(lepton,jet1,jet3,jet4)
        minLepJetDPhi  = self.HLL.MinDPhi_lep3j(lepton,jet1,jet3,jet4)
        minLepJetDEta  = self.HLL.MinDEta_lep3j(lepton,jet1,jet3,jet4)
        maxLepJetDR    = self.HLL.MaxDR_lep3j(lepton,jet1,jet3,jet4)
        maxLepJetDPhi  = self.HLL.MaxDPhi_lep3j(lepton,jet1,jet3,jet4)
        maxLepJetDEta  = self.HLL.MaxDEta_lep3j(lepton,jet1,jet3,jet4)
        HT2_lepJetMet  = self.HLL.HT2_1b2Wj(lepton,jet1,jet3,jet4,self.corrMET)
        HT2R_lepJetMet = self.HLL.HT2R_1b2Wj(lepton,jet1,jet3,jet4,self.corrMET)
        wj1wj2_pt      = (jet3.p4 + jet4.p4).Pt()
        wj1wj2_M       = op.invariant_mass(jet3.p4, jet4.p4)
        w1w2_MT        = self.HLL.MT_W1W2_ljj(lepton,jet3,jet4,self.corrMET) # MT of H->WW->lep, jet3, jet4  and Met
        JPAcat         = op.c_int(3)
        
    if Res2b0Wj : 
        pt_taggedJetsAvg = (self.HLL.bJetCorrP4(jet1).Pt()+self.HLL.bJetCorrP4(jet2).Pt())/2.0
        mass_taggedJetsAvg = (self.HLL.bJetCorrP4(jet1).M()+self.HLL.bJetCorrP4(jet2).M())/2.0
        bj1bj2_pt      = (self.HLL.bJetCorrP4(jet1)+self.HLL.bJetCorrP4(jet2)).Pt()
        bj1bj2_M       = op.invariant_mass(self.HLL.bJetCorrP4(jet1),self.HLL.bJetCorrP4(jet2))
        mT_top_3particle  = op.min(self.HLL.mT2(self.HLL.bJetCorrP4(jet1),lepton.p4 ,self.corrMET.p4), 
                                   self.HLL.mT2(self.HLL.bJetCorrP4(jet2), lepton.p4, self.corrMET.p4))
        avgJetDR       = op.deltaR(jet1.p4,jet2.p4)
        avgJetDPhi     = op.abs(op.deltaPhi(jet1.p4,jet2.p4))
        avgJetDEta     = op.abs(jet1.eta - jet2.eta)
        minJetDR       = op.deltaR(jet1.p4,jet2.p4)
        minJetDPhi     = op.abs(op.deltaPhi(jet1.p4,jet2.p4))
        minJetDEta     = op.abs(jet1.eta - jet2.eta)
        maxJetDR       = op.deltaR(jet1.p4,jet2.p4)
        maxJetDPhi       = op.abs(op.deltaPhi(jet1.p4,jet2.p4))
        maxJetDEta       = op.abs(jet1.eta-jet2.eta)
        minLepJetDR    = self.HLL.MinDR_lep2j(lepton,jet1,jet2)
        minLepJetDPhi  = self.HLL.MinDPhi_lep2j(lepton,jet1,jet2)
        minLepJetDEta  = self.HLL.MinDEta_lep2j(lepton,jet1,jet2)
        maxLepJetDR    = self.HLL.MaxDR_lep2j(lepton,jet1,jet2)
        maxLepJetDPhi  = self.HLL.MaxDPhi_lep2j(lepton,jet1,jet2)
        maxLepJetDEta  = self.HLL.MaxDEta_lep2j(lepton,jet1,jet2)
        HT2_lepJetMet  = self.HLL.HT2_2b0Wj(lepton,jet1,jet2,self.corrMET)
        HT2R_lepJetMet = self.HLL.HT2R_2b0Wj(lepton,jet1,jet2,self.corrMET)
        wj1wj2_pt      = op.c_float(0.)
        wj1wj2_M       = op.c_float(0.) 
        w1w2_MT        = op.c_float(0.) # MT of H->WW->lep, jet3  and Met
        JPAcat         = op.c_int(4)
        
    if Res1b1Wj : 
        pt_taggedJetsAvg = (self.HLL.bJetCorrP4(jet1).Pt()+jet3.pt)/2.0
        mass_taggedJetsAvg = (self.HLL.bJetCorrP4(jet1).M()+jet3.mass)/2.0
        bj1bj2_pt      = self.HLL.bJetCorrP4(jet1).Pt()
        bj1bj2_M       = self.HLL.bJetCorrP4(jet1).M()
        mT_top_3particle  = self.HLL.mT2(self.HLL.bJetCorrP4(jet1),lepton.p4 ,self.corrMET.p4)
        avgJetDR       = op.deltaR(jet1.p4,jet3.p4)
        avgJetDPhi     = op.abs(op.deltaPhi(jet1.p4,jet3.p4))
        avgJetDEta     = op.abs(jet1.eta - jet3.eta)
        minJetDR       = op.deltaR(jet1.p4,jet3.p4)
        minJetDPhi     = op.abs(op.deltaPhi(jet1.p4,jet3.p4))
        minJetDEta     = op.abs(jet1.eta-jet3.eta)
        maxJetDR       = op.deltaR(jet1.p4,jet3.p4)
        maxJetDPhi     = op.abs(op.deltaPhi(jet1.p4,jet3.p4))
        maxJetDEta     = op.abs(jet1.eta-jet3.eta)
        minLepJetDR    = self.HLL.MinDR_lep2j(lepton,jet1,jet3)
        minLepJetDPhi  = self.HLL.MinDPhi_lep2j(lepton,jet1,jet3)
        minLepJetDEta  = self.HLL.MinDEta_lep2j(lepton,jet1,jet3)
        maxLepJetDR    = self.HLL.MaxDR_lep2j(lepton,jet1,jet3)
        maxLepJetDPhi  = self.HLL.MaxDPhi_lep2j(lepton,jet1,jet3)
        maxLepJetDEta  = self.HLL.MaxDEta_lep2j(lepton,jet1,jet3)
        HT2_lepJetMet  = self.HLL.HT2_1b1Wj(lepton,jet1,jet3,self.corrMET)
        HT2R_lepJetMet = self.HLL.HT2R_1b1Wj(lepton,jet1,jet3,self.corrMET)
        wj1wj2_pt      = jet3.pt   
        wj1wj2_M       = jet3.p4.M()  
        w1w2_MT        = self.HLL.MT_W1W2_lj(lepton,jet3,self.corrMET) # MT of H->WW->lep, jet3  and Met
        JPAcat         = op.c_int(5)
        
    if Res1b0Wj : 
        pt_taggedJetsAvg = self.HLL.bJetCorrP4(jet1).Pt()
        mass_taggedJetsAvg = self.HLL.bJetCorrP4(jet1).M()
        bj1bj2_pt      = self.HLL.bJetCorrP4(jet1).Pt()
        bj1bj2_M       = self.HLL.bJetCorrP4(jet1).M()
        mT_top_3particle  = self.HLL.mT2(self.HLL.bJetCorrP4(jet1),lepton.p4 ,self.corrMET.p4)
        avgJetDR       = op.c_float(0.)
        avgJetDPhi     = op.c_float(0.)
        avgJetDEta     = op.c_float(0.)
        minJetDR       = op.c_float(0.)
        minJetDPhi     = op.c_float(0.)
        minJetDEta     = op.c_float(0.)
        maxJetDR       = op.c_float(0.)
        maxJetDPhi     = op.c_float(0.)
        maxJetDEta     = op.c_float(0.)
        minLepJetDR    = op.deltaR(lepton.p4, jet1.p4)
        minLepJetDPhi  = op.abs(op.deltaPhi(lepton.p4, jet1.p4))
        minLepJetDEta  = op.abs(lepton.eta- jet1.eta)
        maxLepJetDR    = op.deltaR(lepton.p4, jet1.p4)
        maxLepJetDPhi  = op.abs(op.deltaPhi(lepton.p4, jet1.p4))
        maxLepJetDEta  = op.abs(lepton.eta- jet1.eta)
        HT2_lepJetMet  = self.HLL.HT2_1b0Wj(lepton,jet1,self.corrMET)
        HT2R_lepJetMet = self.HLL.HT2R_1b0Wj(lepton,jet1,self.corrMET)
        wj1wj2_pt      = op.c_float(0.)
        wj1wj2_M       = op.c_float(0.)
        w1w2_MT        = op.c_float(0.)
        JPAcat         = op.c_int(6)
        
    if Res0b : 
        pt_taggedJetsAvg = op.c_float(0.)
        mass_taggedJetsAvg = op.c_float(0.)
        bj1bj2_pt      = op.c_float(0.)
        bj1bj2_M       = op.c_float(0.)
        mT_top_3particle = op.c_float(0.)
        avgJetDR       = op.c_float(0.)
        avgJetDPhi     = op.c_float(0.)
        avgJetDEta     = op.c_float(0.)
        minJetDR       = op.c_float(0.)
        minJetDPhi     = op.c_float(0.)
        minJetDEta     = op.c_float(0.)
        maxJetDR       = op.c_float(0.)
        maxJetDPhi     = op.c_float(0.)
        maxJetDEta     = op.c_float(0.)
        minLepJetDR    = op.c_float(0.)
        minLepJetDPhi  = op.c_float(0.)
        minLepJetDEta  = op.c_float(0.)
        maxLepJetDR    = op.c_float(0.)
        maxLepJetDPhi  = op.c_float(0.)
        maxLepJetDEta  = op.c_float(0.)
        HT2_lepJetMet  = self.HLL.HT2_0b(lepton,self.corrMET)
        HT2R_lepJetMet = self.HLL.HT2R_0b(lepton,self.corrMET)
        wj1wj2_pt      = op.c_float(0.)   
        wj1wj2_M       = op.c_float(0.)
        w1w2_MT        = op.c_float(0.)
        JPAcat         = op.c_int(7)
        
    return{('L1_2b2Wj',                'L1_2b2Wj',                (25,0,1)):        L1out[0],
           ('L1_2b1Wj',                'L1_2b1Wj',                (25,0,1)):        L1out[1],
           ('L1_1b2Wj',                'L1_1b2Wj',                (25,0,1)):        L1out[2],
           ('L1_2b0Wj',                'L1_2b0Wj',                (25,0,1)):        L1out[3],
           ('L1_1b1Wj',                'L1_1b1Wj',                (25,0,1)):        L1out[4],
           ('L1_1b0Wj',                'L1_1b0Wj',                (25,0,1)):        L1out[5],
           ('L2_2b2Wj',                'L2_2b2Wj',                (25,0,1)):        L2out[0],
           ('L2_2b1Wj',                'L2_2b1Wj',                (25,0,1)):        L2out[1],
           ('L2_1b2Wj',                'L2_1b2Wj',                (25,0,1)):        L2out[2],
           ('L2_2b0Wj',                'L2_2b0Wj',                (25,0,1)):        L2out[3],
           ('L2_1b1Wj',                'L2_1b1Wj',                (25,0,1)):        L2out[4],
           ('L2_1b0Wj',                'L2_1b0Wj',                (25,0,1)):        L2out[5],
           ('L2_0b',                   'L2_0b',                   (25,0,1)):        L2out[6],
           ('leppdgId',                'leppdgId',                (45,-22.,22.)):   lepton.pdgId,
           ('lepcharge',               'lepcharge',               (3,-1,1.)):       lepton.charge,
           ('lepmet_DPhi',             'lepmet_DPhi',             (20,0,3.2)):      op.abs(self.HLL.SinglepMet_dPhi(lepton, self.corrMET)),
           ('lepmet_pt',               'lepmet_pt',               (50,0,300)):      self.HLL.SinglepMet_Pt(lepton, self.corrMET),
           ('lep_MT',                  'lep_MT',                  (40,0,200)):      self.HLL.MT(lepton, self.corrMET),
           ('lep_conept',              'lep_conept',              (40,0,200)):      self.HLL.lambdaConePt(lepton),
           ('minDR_lep_allJets',       'minDR_lep_allJets',       (5,0,5)):         self.HLL.MinDR_part1_partCont(lepton, self.ak4Jets),
           ('minDEta_lep_allJets',     'minDEta_lep_allJets',     (5,0,5)):         self.HLL.MinDEta_part1_partCont(lepton, self.ak4Jets),
           #('minDPhi_lep_allJets',     'minDPhi_lep_allJets',     (5,0,5)):         self.HLL.MinDPhi_part1_partCont(lepton, self.ak4Jets),
           ('bj1Px',                   'bj1_Px',                  (50,-250.,250.)): bj1_Px,
           ('bj1Py',                   'bj1_Py',                  (50,-250.,250.)): bj1_Py,
           ('bj1Pz',                   'bj1_Pz',                  (50,-250.,250.)): bj1_Pz,
           ('bj1E',                    'bj1_E',                   (50,0.,500.)):    bj1_E,
           ('bj1M',                    'bj1_M',                   (50,0.,500.)):    bj1_M,
           ('bj1Pt',                   'bj1_Pt',                  (50,-250.,250.)): bj1_pt,
           ('bj1Eta',                  'bj1_Eta',                 (50,0.,5.)):      bj1_eta,
           ('bj1Phi',                  'bj1Phi',                  (50,0.,5.)):      bj1_phi,
           ('bj1bTagDeepFlavB',        'bj1_bTagDeepFlavB',       (50,0.,1.)):      bj1_bTagDeepFlavB,
           ('bj1LepDR',                'bj1LepDR',                (50,0.,5.)):      bj1LepDR,
           ('bj1LepDEta',              'bj1LepDEta',              (50,0.,5.)):      bj1LepDEta,
           ('bj1LepDPhi',              'bj1LepDPhi',              (50,0.,5.)):      bj1LepDPhi,
           ('bj1MetDPhi',              'bj1MetDPhi',              (50,0.,5.)):      bj1MetDPhi,
           ('bj2Px',                   'bj2_Px',                  (50,-250.,250.)): bj2_Px,
           ('bj2Py',                   'bj2_Py',                  (50,-250.,250.)): bj2_Py,
           ('bj2Pz',                   'bj2_Pz',                  (50,-250.,250.)): bj2_Pz,
           ('bj2E',                    'bj2_E',                   (50,0.,500.)):    bj2_E,
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
           ('cosThetaS_Hbb',           'cosThetaS_Hbb',           (50,0.,1.)):      cosThetaS_Hbb,
           ('wj1Px',                   'wj1_Px',                  (50,-250.,250.)): wj1_Px,
           ('wj1Py',                   'wj1_Py',                  (50,-250.,250.)): wj1_Py,
           ('wj1Pz',                   'wj1_Pz',                  (50,-250.,250.)): wj1_Pz,
           ('wj1E',                    'wj1_E',                   (50,0.,500.)):    wj1_E,
           ('wj1M',                    'wj1_M',                   (50,0.,500.)):    wj1_M,
           ('wj1Pt',                   'wj1_Pt',                  (50,-250.,250.)): wj1_pt,
           ('wj1Eta',                  'wj1_Eta',                 (50,0.,5.)):      wj1_eta,
           ('wj1Phi',                  'wj1Phi',                  (50,0.,5.)):      wj1_phi,
           ('wj1bTagDeepFlavB',        'wj1_bTagDeepFlavB',       (50,0.,1.)):      wj1_bTagDeepFlavB,
           ('wj1LepDR',                'wj1LepDR',                (50,0.,5.)):      wj1LepDR,
           ('wj1LepDEta',              'wj1LepDEta',              (50,0.,5.)):      wj1LepDEta,
           ('wj1LepDPhi',              'wj1LepDPhi',              (50,0.,5.)):      wj1LepDPhi,
           ('wj1MetDPhi',              'wj1MetDPhi',              (50,0.,5.)):      wj1MetDPhi,
           ('wj2Px',                   'wj2_Px',                  (50,-250.,250.)): wj2_Px,
           ('wj2Py',                   'wj2_Py',                  (50,-250.,250.)): wj2_Py,
           ('wj2Pz',                   'wj2_Pz',                  (50,-250.,250.)): wj2_Pz,
           ('wj2E',                    'wj2_E',                   (50,0.,500.)):    wj2_E,
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
           ('HWW_Mass',                'HWW_Mass',                (50,0.,300.)):    HWW_Mass,
           ('HWW_Simple_Mass',         'HWW_Simple_Mass',         (50,0.,300.)):    HWW_Simple_Mass,
           ('HWW_dR',                  'HWW_dR',                  (50,0.,5.)):      HWW_dR,
           ('HWW_dPhi',                'HWW_dPhi',                (50,0.,5.)):      HWW_dPhi,
           ('HWW_dEta',                'HWW_dEta',                (50,0.,5.)):      HWW_dEta,
           ('angleBetWWPlane',         'angleBetWWPlane',         (50,0.,5.)):      angleBetWWPlane,           
           ('cosThetaS_Wjj_simple',    'cosThetaS_Wjj_simple',    (20,0.,1.)):      cosThetaS_Wjj_simple,
           ('cosThetaS_WW_simple_met', 'cosThetaS_WW_simple_met', (20,0.,1.)):      cosThetaS_WW_simple_met,
           ('cosThetaS_HH_simple_met', 'cosThetaS_HH_simple_met', (20,0.,1.)):      cosThetaS_HH_simple_met,
           ('mHH',                     'mHH',                     (50,0.,300.)):    mHH,
           ('angleBetHWPlane',         'angleBetHWPlane',         (50,0.,5.)):      angleBetHWPlane,
           ('neuPx',                   'neuPx',                   (50,-100.,100.)): neuPx,
           ('neuPy',                   'neuPy',                   (50,-100.,100.)): neuPy,
           ('neuPz',                   'neuPz',                   (50,-100.,100.)): neuPz,
           ('neuPt',                   'neuPt',                   (50,-100.,100.)): neuPt,
           ('neuE',                    'neuE',                    (50,-100.,100.)): neuE,
           ('bj1bj2_DR',               'bj1bj2_DR',               (50,0.,5)):       bj1bj2_DR,
           ('bj1bj2_DPhi',             'bj1bj2_DPhi',             (50,0.,5)):       bj1bj2_DPhi,
           ('bj1bj2_DEta',             'bj1bj2_DEta',             (50,0.,5)):       bj1bj2_DEta,
           ('bj1wj1_DR',               'bj1wj1_DR',               (50,0.,5)):       bj1wj1_DR,
           ('bj1wj1_DPhi',             'bj1wj1_DPhi',             (50,0.,5)):       bj1wj1_DPhi,
           ('bj1wj1_DEta',             'bj1wj1_DEta',             (50,0.,5)):       bj1wj1_DEta,
           ('bj1wj2_DR',               'bj1wj2_DR',               (50,0.,5)):       bj1wj2_DR,
           ('bj1wj2_DPhi',             'bj1wj2_DPhi',             (50,0.,5)):       bj1wj2_DPhi,
           ('bj1wj2_DEta',             'bj1wj2_DEta',             (50,0.,5)):       bj1wj2_DEta,
           ('bj2wj1_DR',               'bj2wj1_DR',               (50,0.,5)):       bj2wj1_DR,
           ('bj2wj1_DPhi',             'bj2wj1_DPhi',             (50,0.,5)):       bj2wj1_DPhi,
           ('bj2wj1_DEta',             'bj2wj1_DEta',             (50,0.,5)):       bj2wj1_DEta,
           ('bj2wj2_DR',               'bj2wj2_DR',               (50,0.,5)):       bj2wj2_DR,
           ('bj2wj2_DPhi',             'bj2wj2_DPhi',             (50,0.,5)):       bj2wj2_DPhi,
           ('bj2wj2_DEta',             'bj2wj2_DEta',             (50,0.,5)):       bj2wj2_DEta,
           ('wj1wj2_DR',               'wj1wj2_DR',               (50,0.,5)):       wj1wj2_DR,
           ('wj1wj2_DPhi',             'wj1wj2_DPhi',             (50,0.,5)):       wj1wj2_DPhi,
           ('wj1wj2_DEta',             'wj1wj2_DEta',             (50,0.,5)):       wj1wj2_DEta,
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
           ('VBF_tag',                 'VBF_tag',                 (2, 0, 2)):       VBF_tag,
           ('zeppenfeldVar',           'zeppenfeldVar',           (50,0.,200)):     zeppenfeldVar,           
           ('pt_taggedJetsAvg',        'pt_taggedJetsAvg',        (50,0.,200)):     pt_taggedJetsAvg,
           ('mass_taggedJetsAvg',      'mass_taggedJetsAvg',      (50,0.,200)):     mass_taggedJetsAvg,
           ('bj1bj2_pt',               'bj1bj2_pt',               (50,0.,200)):     bj1bj2_pt,
           ('bj1bj2_M',                'bj1bj2_M',                (50,0.,200)):     bj1bj2_M,
           ('mT_top_3particle',        'mT_top_3particle',        (50,0.,200)):     mT_top_3particle,
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
           ('HT2_lepJetMet',           'HT2_lepJetMet',           (50,0.,200)):     HT2_lepJetMet,
           ('HT2R_lepJetMet',          'HT2R_lepJetMet',          (50,0.,200)):     HT2R_lepJetMet,
           ('wj1wj2_pt',               'wj1wj2_pt',               (50,0.,200)):     wj1wj2_pt,
           ('wj1wj2_M',                'wj1wj2_M',                (50,0.,200)):     wj1wj2_M,
           ('w1w2_MT',                 'w1w2_MT',                 (50,0.,200)):     w1w2_MT,
           ('JPAcat',                  'JPAcat',                  (8,0.,8)):        JPAcat 
    }


def returnLBNInputs_Resolved(self,lepton,jet1=None,jet2=None,jet3=None,jet4=None) :
    lep_Px = lepton.p4.Px();
    lep_Py = lepton.p4.Py();
    lep_Pz = lepton.p4.Pz();
    lep_E  = lepton.p4.E();
    bj1_Px = op.c_float(0.) if jet1 == None else jet1.p4.Px()
    bj1_Py = op.c_float(0.) if jet1 == None else jet1.p4.Py()
    bj1_Pz = op.c_float(0.) if jet1 == None else jet1.p4.Pz()
    bj1_E  = op.c_float(0.) if jet1 == None else jet1.p4.E()
    bj2_Px = op.c_float(0.) if jet2 == None else jet2.p4.Px()
    bj2_Py = op.c_float(0.) if jet2 == None else jet2.p4.Py()
    bj2_Pz = op.c_float(0.) if jet2 == None else jet2.p4.Pz()
    bj2_E  = op.c_float(0.) if jet2 == None else jet2.p4.E()
    wj1_Px = op.c_float(0.) if jet3 == None else jet3.p4.Px()
    wj1_Py = op.c_float(0.) if jet3 == None else jet3.p4.Py()
    wj1_Pz = op.c_float(0.) if jet3 == None else jet3.p4.Pz()
    wj1_E  = op.c_float(0.) if jet3 == None else jet3.p4.E()
    wj2_Px = op.c_float(0.) if jet4 == None else jet4.p4.Px()
    wj2_Py = op.c_float(0.) if jet4 == None else jet4.p4.Py()
    wj2_Pz = op.c_float(0.) if jet4 == None else jet4.p4.Pz()
    wj2_E  = op.c_float(0.) if jet4 == None else jet4.p4.E()
    return {('lep_E',   'lep_E',   (50,0.,500.)):    lep_E,
            ('lep_Px',  'lep_Px',  (50,-250.,250.)): lep_Px,
            ('lep_Py',  'lep_Py',  (50,-250.,250.)): lep_Py,
            ('lep_Pz',  'lep_Pz',  (50,-250.,250.)): lep_Pz,
            ('bj1_E',   'bj1_E',   (50,0.,500.)):    bj1_E,
            ('bj1_Px',  'bj1_Px',  (50,-250.,250.)): bj1_Px,
            ('bj1_Py',  'bj1_Py',  (50,-250.,250.)): bj1_Py,
            ('bj1_Pz',  'bj1_Pz',  (50,-250.,250.)): bj1_Pz,
            ('bj2_E',   'bj2_E',   (50,0.,500.)):    bj2_E,
            ('bj2_Px',  'bj2_Px',  (50,-250.,250.)): bj2_Px,
            ('bj2_Py',  'bj2_Py',  (50,-250.,250.)): bj2_Py,
            ('bj2_Pz',  'bj2_Pz',  (50,-250.,250.)): bj2_Pz,
            ('wj1_E',   'wj1_E',   (50,0.,500.)):    wj1_E,
            ('wj1_Px',  'wj1_Px',  (50,-250.,250.)): wj1_Px,
            ('wj1_Py',  'wj1_Py',  (50,-250.,250.)): wj1_Py,
            ('wj1_Pz',  'wj1_Pz',  (50,-250.,250.)): wj1_Pz,
            ('wj2_E',   'wj2_E',   (50,0.,500.)):    wj2_E,
            ('wj2_Px',  'wj2_Px',  (50,-250.,250.)): wj2_Px,
            ('wj2_Py',  'wj2_Py',  (50,-250.,250.)): wj2_Py,
            ('wj2_Pz',  'wj2_Pz',  (50,-250.,250.)): wj2_Pz
    }
    

### --------------------------------------------------------------------------------------------------------- ###
###                                                    Boosted                                                ###
### --------------------------------------------------------------------------------------------------------- ###
def returnCommonInputs_Boosted(self):
    return{('METpx',                   'METpx',                   (50, 0, 300)):    self.corrMET.p4.Px(),
           ('METpy',                   'METpy',                   (50, 0, 300)):    self.corrMET.p4.Py(),
           ('METpz',                   'METpz',                   (50, 0, 300)):    self.corrMET.p4.Pz(),
           ('METenergy',               'METenergy',               (50, 0, 300)):    self.corrMET.p4.E(),
           ('METpt',                   'METpt',                   (50, 0, 300)):    self.corrMET.pt,
           ('METphi',                  'METphi',                  (50, 0, 300)):    self.corrMET.phi,
           ('MET_LD',                  'MET_LD',                  (40,0,200)):      self.HLL.MET_LD_DL(self.corrMET, self.ak4Jets, self.electronsFakeSel, self.muonsFakeSel),
           ('nfakeableMuons',          'nfakeableMuons',          (5,0,5)):         op.static_cast("UInt_t",op.rng_len(self.muonsFakeSel)),
           ('nfakeableElectrons',      'nfakeableElectrons',      (5,0,5)):         op.static_cast("UInt_t",op.rng_len(self.electronsFakeSel)),
           ('hT',                      'hT',                      (50,0,300)):      op.static_cast("Float_t", self.HLL.HT_SL(self.ak4Jets)),
           ('hT_withCleanAk4',         'hT_withCleanAk4',         (50,0,300)):      op.static_cast("Float_t", self.HLL.HT_SL(self.ak4JetsCleanedFromAk8b)),
           ('hT_withAk8',              'hT_withAk8',              (50,0,300)):      op.static_cast("Float_t", (self.ak8BJets[0].pt+self.HLL.HT_SL(self.ak4JetsCleanedFromAk8b))),
           ('nAk4Jets',                'nAk4Jets',                (10,0,10)):       op.static_cast("UInt_t",op.rng_len(self.ak4Jets)),
           ('nAk4BJets',               'nAk4BJets',               (10,0,10)):       op.static_cast("UInt_t",op.rng_len(self.ak4BJets)),
           ('nAk4BJetsLoose',          'nAk4BJetsLoose',          (10,0,10)):       op.static_cast("UInt_t",op.rng_len(self.ak4BJetsLoose)),
           ('nAk8Jets',                'nAk8Jets',                (5,0,5)):         op.static_cast("UInt_t",op.rng_len(self.ak8Jets)),
           ('nAk8BJets',               'nAk8BJets',               (5,0,5)):         op.static_cast("UInt_t",op.rng_len(self.ak8BJets)),
           ('nAk4JetsCleandWithAk8b',  'nAk4JetsCleandWithAk8b',  (5,0,5)):         op.static_cast("UInt_t",op.rng_len(self.ak4JetsCleanedFromAk8b))
       }


def returnClassicInputs_Boosted(self, lepton, jpaSelectedJets, L1out, L2out, jpaArg):
    #######################################
    # Hbb2Wj : jet1 jet2   jet3   jet4
    # Hbb1Wj : jet1 jet2   jet3   jet4=0
    # Hbb0Wj : jet1 jet2   jet3=0 jet4=0
    #######################################
    Hbb2Wj = False
    Hbb1Wj = False
    Hbb0Wj = False

    fjet = self.ak8BJets[0]

    if jpaArg == 'Hbb2Wj':
        Hbb2Wj = True
        jet3 = jpaSelectedJets[0]
        jet4 = jpaSelectedJets[1]
        VBFJetPairsJPA = op.sort(op.combine(op.select(self.VBFJets, cleanVBFwithJPA_Boosted(self.ak8BJets, jpaSelectedJets, 2)), N=2, pred=self.lambda_VBFPair), 
                                 lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))
    elif jpaArg == 'Hbb1Wj':
        Hbb1Wj = True
        jet3 = jpaSelectedJets[0]
        jet4 = None
        VBFJetPairsJPA = op.sort(op.combine(op.select(self.VBFJets, cleanVBFwithJPA_Boosted(self.ak8BJets, jpaSelectedJets, 1)), N=2, pred=self.lambda_VBFPair), 
                                 lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))
    elif jpaArg == 'Hbb0Wj':
        Hbb0Wj = True
        jet3 = None
        jet4 = None
        VBFJetPairsJPA = op.sort(op.combine(op.select(self.VBFJets, cleanVBFwithJPA_Boosted_rest(self.ak8BJets)), N=2, pred=self.lambda_VBFPair),
                                 lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))
    else:
        raise RuntimeError('Boosted JPA category is not mentioned!')

    # FatJet Variables
    fatbj_Px             = fjet.p4.Px()
    fatbj_Py             = fjet.p4.Py()
    fatbj_Pz             = fjet.p4.Pz()
    fatbj_E              = fjet.p4.E()
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
    
    # Wj1 variables
    wj1_Px            = jet3.p4.Px()                 if not Hbb0Wj else op.c_float(0.)
    wj1_Py            = jet3.p4.Py()                 if not Hbb0Wj else op.c_float(0.)
    wj1_Pz            = jet3.p4.Pz()                 if not Hbb0Wj else op.c_float(0.)
    wj1_E             = jet3.p4.E()                  if not Hbb0Wj else op.c_float(0.)
    wj1_pt            = jet3.pt                      if not Hbb0Wj else op.c_float(0.)
    wj1_eta           = jet3.eta                     if not Hbb0Wj else op.c_float(0.)
    wj1_bTagDeepFlavB = jet3.btagDeepFlavB           if not Hbb0Wj else op.c_float(0.)

    # Wj2 variables
    wj2_Px            = jet4.p4.Px()                 if Hbb2Wj else op.c_float(0.)
    wj2_Py            = jet4.p4.Py()                 if Hbb2Wj else op.c_float(0.)
    wj2_Pz            = jet4.p4.Pz()                 if Hbb2Wj else op.c_float(0.)
    wj2_E             = jet4.p4.E()                  if Hbb2Wj else op.c_float(0.)
    wj2_pt            = jet4.pt                      if Hbb2Wj else op.c_float(0.)
    wj2_eta           = jet4.eta                     if Hbb2Wj else op.c_float(0.)
    wj2_bTagDeepFlavB = jet4.btagDeepFlavB           if Hbb2Wj else op.c_float(0.)
    
    # fatjet-Wj1
    fatbj_Wj1DR       = op.deltaR(fjet.p4, jet3.p4)                                if not Hbb0Wj else op.c_float(0.)
    fatbj_Wj1DPhi     = op.abs(op.deltaPhi(fjet.p4, jet3.p4))                      if not Hbb0Wj else op.c_float(0.)
    fatbjSub1_Wj1DR   = op.switch(fjet.subJet1.pt > fjet.subJet2.pt, 
                                  op.deltaR(fjet.subJet1.p4, jet3.p4), 
                                  op.deltaR(fjet.subJet2.p4, jet3.p4))             if not Hbb0Wj else op.c_float(0.)
    fatbjSub1_Wj1DPhi = op.switch(fjet.subJet1.pt > fjet.subJet2.pt, 
                                    op.abs(op.deltaPhi(fjet.subJet1.p4, jet3.p4)), 
                                    op.abs(op.deltaPhi(fjet.subJet2.p4, jet3.p4))) if not Hbb0Wj else op.c_float(0.)
    fatbjSub2_Wj1DR   = op.switch(fjet.subJet1.pt > fjet.subJet2.pt,
                                    op.deltaR(fjet.subJet2.p4, jet3.p4),
                                    op.deltaR(fjet.subJet1.p4, jet3.p4))           if not Hbb0Wj else op.c_float(0.)
    fatbjSub2_Wj1DPhi = op.switch(fjet.subJet1.pt > fjet.subJet2.pt,
                                  op.abs(op.deltaPhi(fjet.subJet2.p4, jet3.p4)),
                                  op.abs(op.deltaPhi(fjet.subJet1.p4, jet3.p4)))   if not Hbb0Wj else op.c_float(0.)
    # Wjs-lepton
    wj1_lepDR         = op.deltaR(jet3.p4, lepton.p4)                     if not Hbb0Wj else op.c_float(0.)
    wj1_lepDPhi       = op.abs(op.deltaPhi(jet3.p4, lepton.p4))           if not Hbb0Wj else op.c_float(0.)
    wj2_lepDR         = op.deltaR(jet4.p4, lepton.p4)                     if Hbb2Wj else op.c_float(0.)
    wj2_lepDPhi       = op.abs(op.deltaPhi(jet4.p4, lepton.p4))           if Hbb2Wj else op.c_float(0.)
    # fatjet-Wj2
    fatbj_wj2DR       = op.deltaR(fjet.p4, jet4.p4)                       if Hbb2Wj else op.c_float(0.) 
    fatbj_wj2DPhi     = op.abs(op.deltaPhi(fjet.p4, jet4.p4))             if Hbb2Wj else op.c_float(0.)
    fatbjSub1_wj2DR   = op.deltaR(fjet.subJet1.p4, jet4.p4)               if Hbb2Wj else op.c_float(0.)
    fatbjSub1_wj2DPhi = op.abs(op.deltaPhi(fjet.subJet1.p4, jet4.p4))     if Hbb2Wj else op.c_float(0.)
    fatbjSub2_wj2DR   = op.deltaR(fjet.subJet2.p4, jet4.p4)               if Hbb2Wj else op.c_float(0.)
    fatbjSub2_wj2DPhi = op.abs(op.deltaPhi(fjet.subJet2.p4, jet4.p4))     if Hbb2Wj else op.c_float(0.)
    # Wj1-Wj2
    wj1wj2DR             = op.deltaR(jet3.p4, jet4.p4)                                                                                               if Hbb2Wj else op.c_float(0.)
    wj1wj2DPhi           = op.abs(op.deltaPhi(jet3.p4, jet4.p4))                                                                                     if Hbb2Wj else op.c_float(0.)
    WWplaneAngle         = self.HLL.angleBetPlanes(lepton.p4,self.HLL.neuP4(jet3.p4+jet4.p4+lepton.p4, self.corrMET),jet3.p4, jet4.p4)               if Hbb2Wj else op.c_float(0.) 
    WWplaneAngle_withMET = self.HLL.angleWWplane(lepton.p4, self.corrMET, jet3.p4, jet4.p4)                                                          if Hbb2Wj else op.c_float(0.)
    HWplaneAngle         = self.HLL.angleBetPlanes(fjet.subJet1.p4,fjet.subJet2.p4,jet3.p4,jet4.p4)                                                  if Hbb2Wj else op.c_float(0.)
    
    HWW_Mass          = op.extMethod("HHbbWWJPA::HWW_SimpleMassWithRecoNeutrino", returnType="float")(jet3.p4, jet4.p4, lepton.p4, self.corrMET.p4)  if Hbb2Wj else op.c_float(0.)
    HWW_Simple_Mass   = self.HLL.HWW_met_simple(jet3.p4,jet4.p4,lepton.p4,self.corrMET.p4).M()                                                       if Hbb2Wj else op.c_float(0.)
    HWW_dR            = self.HLL.dR_Hww(jet3.p4,jet4.p4,lepton.p4,self.corrMET)                                                                      if Hbb2Wj else op.c_float(0.)
    cosThetaS_Wjj_simple    = op.extMethod("HHbbWWJPA::cosThetaS", returnType="float")(jet3.p4, jet4.p4)                                             if Hbb2Wj else op.c_float(0.)
    cosThetaS_WW_simple_met = op.extMethod("HHbbWWJPA::cosThetaS", returnType="float")(self.HLL.Wjj_simple(jet3.p4,jet4.p4), self.HLL.Wlep_met_simple(lepton.p4, self.corrMET.p4)) \
                              if Hbb2Wj else op.c_float(0.)
    cosThetaS_HH_simple_met = op.extMethod("HHbbWWJPA::cosThetaS", returnType="float")(fjet.subJet1.p4 + fjet.subJet2.p4, self.HLL.HWW_met_simple(jet3.p4,jet4.p4,lepton.p4,self.corrMET.p4)) \
                              if Hbb2Wj else op.c_float(0.)
    
    
    VBFj1Px        = op.switch(op.rng_len(VBFJetPairsJPA) > 0, 
                                 op.switch(VBFJetPairsJPA[0][0].p4.Px() > VBFJetPairsJPA[0][1].p4.Px(),
                                           VBFJetPairsJPA[0][0].p4.Px(), VBFJetPairsJPA[0][1].p4.Px()), op.c_float(0.)) 
    VBFj2Px        = op.switch(op.rng_len(VBFJetPairsJPA) > 0, 
                                 op.switch(VBFJetPairsJPA[0][0].p4.Px() > VBFJetPairsJPA[0][1].p4.Px(),
                                           VBFJetPairsJPA[0][1].p4.Px(), VBFJetPairsJPA[0][0].p4.Px()), op.c_float(0.)) 
    VBFj1Py        = op.switch(op.rng_len(VBFJetPairsJPA) > 0, 
                                 op.switch(VBFJetPairsJPA[0][0].p4.Py() > VBFJetPairsJPA[0][1].p4.Py(),
                                           VBFJetPairsJPA[0][0].p4.Py(), VBFJetPairsJPA[0][1].p4.Py()), op.c_float(0.)) 
    VBFj2Py        = op.switch(op.rng_len(VBFJetPairsJPA) > 0, 
                                 op.switch(VBFJetPairsJPA[0][0].p4.Py() > VBFJetPairsJPA[0][1].p4.Py(),
                                           VBFJetPairsJPA[0][1].p4.Py(), VBFJetPairsJPA[0][0].p4.Py()), op.c_float(0.)) 
    VBFj1Pz        = op.switch(op.rng_len(VBFJetPairsJPA) > 0, 
                                 op.switch(VBFJetPairsJPA[0][0].p4.Pz() > VBFJetPairsJPA[0][1].p4.Pz(),
                                           VBFJetPairsJPA[0][0].p4.Pz(), VBFJetPairsJPA[0][1].p4.Pz()), op.c_float(0.)) 
    VBFj2Pz        = op.switch(op.rng_len(VBFJetPairsJPA) > 0, 
                                 op.switch(VBFJetPairsJPA[0][0].p4.Pz() > VBFJetPairsJPA[0][1].p4.Pz(),
                                           VBFJetPairsJPA[0][1].p4.Pz(), VBFJetPairsJPA[0][0].p4.Pz()), op.c_float(0.)) 
    VBFj1E         = op.switch(op.rng_len(VBFJetPairsJPA) > 0, 
                                 op.switch(VBFJetPairsJPA[0][0].p4.E() > VBFJetPairsJPA[0][1].p4.E(),
                                           VBFJetPairsJPA[0][0].p4.E(), VBFJetPairsJPA[0][1].p4.E()), op.c_float(0.))   
    VBFj2E         = op.switch(op.rng_len(VBFJetPairsJPA) > 0, 
                                 op.switch(VBFJetPairsJPA[0][0].p4.E() > VBFJetPairsJPA[0][1].p4.E(), 
                                           VBFJetPairsJPA[0][1].p4.E(), VBFJetPairsJPA[0][0].p4.E()), op.c_float(0.))   
    VBFj1pt        = op.switch(op.rng_len(VBFJetPairsJPA) > 0, 
                                 op.switch(VBFJetPairsJPA[0][0].pt > VBFJetPairsJPA[0][1].pt,
                                           VBFJetPairsJPA[0][0].pt, VBFJetPairsJPA[0][1].pt), op.c_float(0.))    
    VBFj2pt        = op.switch(op.rng_len(VBFJetPairsJPA) > 0, 
                                 op.switch(VBFJetPairsJPA[0][0].pt > VBFJetPairsJPA[0][1].pt,
                                           VBFJetPairsJPA[0][1].pt, VBFJetPairsJPA[0][0].pt), op.c_float(0.))    
    VBFj1eta       = op.switch(op.rng_len(VBFJetPairsJPA) > 0, 
                                 op.switch(VBFJetPairsJPA[0][0].pt > VBFJetPairsJPA[0][1].pt,
                                           VBFJetPairsJPA[0][0].eta, VBFJetPairsJPA[0][1].eta), op.c_float(0.))  
    VBFj2eta       = op.switch(op.rng_len(VBFJetPairsJPA) > 0, 
                                 op.switch(VBFJetPairsJPA[0][0].pt > VBFJetPairsJPA[0][1].pt,
                                           VBFJetPairsJPA[0][1].eta, VBFJetPairsJPA[0][0].eta), op.c_float(0.))  
    
    VBFj1j2dEta    = op.switch(op.rng_len(VBFJetPairsJPA) > 0,op.abs(VBFJetPairsJPA[0][0].eta - VBFJetPairsJPA[0][1].eta), op.c_float(0.))
    VBFj1j2dPhi    = op.switch(op.rng_len(VBFJetPairsJPA) > 0,op.abs(op.deltaPhi(VBFJetPairsJPA[0][0].p4,VBFJetPairsJPA[0][1].p4)), op.c_float(0.))
    VBFj1j2invM    = op.switch(op.rng_len(VBFJetPairsJPA) > 0,op.invariant_mass(VBFJetPairsJPA[0][0].p4,  VBFJetPairsJPA[0][1].p4), op.c_float(0.))
    VBF_tag        = op.c_int(op.rng_len(VBFJetPairsJPA)>0)
    zeppenfeldVar  = op.switch(op.rng_len(VBFJetPairsJPA) > 0, op.abs(lepton.eta - (VBFJetPairsJPA[0][0].eta + VBFJetPairsJPA[0][1].eta)/2.0)/op.abs(VBFJetPairsJPA[0][0].eta - VBFJetPairsJPA[0][1].eta),op.c_float(0.))
    
    if Hbb2Wj:
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
        JPAcat      = op.c_int(1)                
        
    if Hbb1Wj:
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
        JPAcat      = op.c_int(2)                
        
    if Hbb0Wj:
        wj1wj2_pt   = op.c_float(0.)
        wj1wj2invM  = op.c_float(0.)   
        #jetMinDR    = op.deltaR(fjet.subJet1.p4, fjet.subJet2.p4)
        jetMinDR    = op.c_float(0.)
        jetMaxDR    = op.c_float(0.)
        jetLepMinDR = op.min(op.deltaR(fjet.subJet1.p4, lepton.p4), op.deltaR(fjet.subJet2.p4, lepton.p4)) 
        jetLepMaxDR = op.max(op.deltaR(fjet.subJet1.p4, lepton.p4), op.deltaR(fjet.subJet2.p4, lepton.p4)) 
        HT2         = self.HLL.HT2_2b0Wj(lepton,fjet.subJet1,fjet.subJet2,self.corrMET)
        HT2R        = self.HLL.HT2R_2b0Wj(lepton,fjet.subJet1,fjet.subJet2,self.corrMET)
        MT_W1W2     = op.c_float(0.)
        HT_taggedJets = fjet.pt
        M_taggedJetsAvg = fjet.msoftdrop
        E_taggedJetsAvg = fjet.p4.E()
        JPAcat      = op.c_int(3)                
        

    return{('L1_Hbb2Wj',               'L1_Hbb2Wj',               (25,0,1)):        L1out[0],
           ('L1_Hbb1Wj',               'L1_Hbb1Wj',               (25,0,1)):        L1out[1],
           ('L2_Hbb2Wj',               'L2_Hbb2Wj',               (25,0,1)):        L2out[0],
           ('L2_Hbb1Wj',               'L2_Hbb1Wj',               (25,0,1)):        L2out[1],
           ('L2_Hbb0Wj',               'L2_Hbb0Wj',               (25,0,1)):        L2out[2],
           ('leppdgId',                'leppdgId',                (45,-22.,22.)):   lepton.pdgId,
           ('lepcharge',               'lepcharge',               (3,-1,1.)):       lepton.charge,
           ('lepmet_DPhi',             'lepmet_DPhi',             (20,0,3.2)):      op.abs(self.HLL.SinglepMet_dPhi(lepton, self.corrMET)),
           ('lepmet_pt',               'lepmet_pt',               (50,0,300)):      self.HLL.SinglepMet_Pt(lepton, self.corrMET),
           ('lep_MT',                  'lep_MT',                  (40,0,200)):      self.HLL.MT(lepton, self.corrMET),
           ('lep_conept',              'lep_conept',              (40,0,200)):      self.HLL.lambdaConePt(lepton),
           ('minDR_lep_allJets',       'minDR_lep_allJets',       (5,0,5)):         self.HLL.MinDR_part1_partCont(lepton, self.ak4Jets),
           ('fatbj_Px',                'fatbj_Px',                (50,-250.,250.)): fatbj_Px,
           ('fatbj_Py',                'fatbj_Py',                (50,-250.,250.)): fatbj_Py,
           ('fatbj_Pz',                'fatbj_Py',                (50,-250.,250.)): fatbj_Pz,
           ('fatbj_E',                 'fatbj_E',                 (50,0.,500.)):    fatbj_E,
           ('fatbj_pt',                'fatbj_pt',                (50,0.,500.)):    fatbj_pt,
           ('fatbj_eta',               'fatbj_eta',               (50,0.,500.)):    fatbj_eta,
           ('fatbj_phi',               'fatbj_phi',               (50,0.,500.)):    fatbj_phi,
           ('fatbj_sub1pt',            'fatbj_sub1pt',            (50,0.,500.)):    fatbj_sub1pt,
           ('fatbj_sub1eta',           'fatbj_sub1eta',           (50,0.,500.)):    fatbj_sub1eta,
           ('fatbj_sub2pt',            'fatbj_sub2pt',            (50,0.,500.)):    fatbj_sub2pt,
           ('fatbj_sub2eta',           'fatbj_sub2eta',           (50,0.,500.)):    fatbj_sub2eta,
           ('fatbj_softdropMass',      'fatbj_softdropMass',      (50,0.,500.)):    fatbj_softdropMass,
           ('fatbj_tau1',              'fatbj_tau1',              (50,0.,500.)):    fatbj_tau1,
           ('fatbj_tau2',              'fatbj_tau2',              (50,0.,500.)):    fatbj_tau2,
           ('fatbj_tau3',              'fatbj_tau3',              (50,0.,500.)):    fatbj_tau3,
           ('fatbj_tau4',              'fatbj_tau4',              (50,0.,500.)):    fatbj_tau4,
           ('fatbj_tau43',             'fatbj_tau43',             (50,0.,500.)):    fatbj_tau43,
           ('fatbj_tau32',             'fatbj_tau32',             (50,0.,500.)):    fatbj_tau32,
           ('fatbj_tau21',             'fatbj_tau21',             (50,0.,500.)):    fatbj_tau21,
           ('fatbj_btagDDBvL',         'fatbj_btagDDBvL',         (50,0.,500.)):    fatbj_btagDDBvL,
           ('fatbj_btagDDBvL_noMD',    'fatbj_btagDDBvL_noMD',    (50,0.,500.)):    fatbj_btagDDBvL_noMD,
           ('fatbj_btagDDCvB',         'fatbj_btagDDCvB',         (50,0.,500.)):    fatbj_btagDDCvB,
           ('fatbj_btagDDCvB_noMD',    'fatbj_btagDDCvB_noMD',    (50,0.,500.)):    fatbj_btagDDCvB_noMD,
           ('fatbj_btagDDCvL',         'fatbj_btagDDCvL',         (50,0.,500.)):    fatbj_btagDDCvL,
           ('fatbj_btagDDCvL_noMD',    'fatbj_btagDDCvL_noMD',    (50,0.,500.)):    fatbj_btagDDCvL_noMD,
           ('fatbj_btagDeepB',         'fatbj_btagDeepB',         (50,0.,500.)):    fatbj_btagDeepB,
           ('cosThetaS_Hbb',           'cosThetaS_Hbb',           (50,0.,500.)):    cosThetaS_Hbb,
           ('mT_top_3particle',        'mT_top_3particle',        (50,0.,500.)):    mT_top_3particle,
           ('wj1_Px',                  'wj1_Px',                  (50,0.,500.)):    wj1_Px,           
           ('wj1_Py',                  'wj1_Py',                  (50,0.,500.)):    wj1_Py,           
           ('wj1_Pz',                  'wj1_Pz',                  (50,0.,500.)):    wj1_Pz,           
           ('wj1_E',                   'wj1_E',                   (50,0.,500.)):    wj1_E,           
           ('wj1_pt',                  'wj1_pt',                  (50,0.,500.)):    wj1_pt,           
           ('wj1_eta',                 'wj1_eta',                 (50,0.,500.)):    wj1_eta,           
           ('wj1_bTagDeepFlavB',       'wj1_bTagDeepFlavB',       (50,0.,500.)):    wj1_bTagDeepFlavB,
           ('wj2_Px',                  'wj2_Px',                  (50,0.,500.)):    wj2_Px,           
           ('wj2_Py',                  'wj2_Py',                  (50,0.,500.)):    wj2_Py,           
           ('wj2_Pz',                  'wj2_Pz',                  (50,0.,500.)):    wj2_Pz,           
           ('wj2_E',                   'wj2_E',                   (50,0.,500.)):    wj2_E,           
           ('wj2_pt',                  'wj2_pt',                  (50,0.,500.)):    wj2_pt,           
           ('wj2_eta',                 'wj2_eta',                 (50,0.,500.)):    wj2_eta,           
           ('wj2_bTagDeepFlavB',       'wj2_bTagDeepFlavB',       (50,0.,500.)):    wj2_bTagDeepFlavB,
           ('fatbj_lepDR',             'fatbj_lepDR',             (50,0.,500.)):    fatbj_lepDR,                     
           ('fatbjSub1_lepDR',         'fatbjSub1_lepDR',         (50,0.,500.)):    fatbjSub1_lepDR,                     
           ('fatbjSub2_lepDR',         'fatbjSub2_lepDR',         (50,0.,500.)):    fatbjSub2_lepDR,                     
           ('fatbj_lepDPhi',           'fatbj_lepDPhi',           (50,0.,500.)):    fatbj_lepDPhi,
           ('fatbjSub1_lepDPhi',       'fatbjSub1_lepDPhi',       (50,0.,500.)):    fatbjSub1_lepDPhi,
           ('minSubJetLepDR',          'minSubJetLepDR',          (50,0.,500.)):    minSubJetLepDR,
           ('fatbj_Wj1DR',             'fatbj_Wj1DR',             (50,0.,500.)):    fatbj_Wj1DR,
           ('fatbj_Wj1DPhi',           'fatbj_Wj1DPhi',           (50,0.,500.)):    fatbj_Wj1DPhi,
           ('fatbjSub1_Wj1DR',         'fatbjSub1_Wj1DR',         (50,0.,500.)):    fatbjSub1_Wj1DR,
           ('fatbjSub1_Wj1DPhi',       'fatbjSub1_Wj1DPhi',       (50,0.,500.)):    fatbjSub1_Wj1DPhi,           
           ('fatbjSub2_Wj1DR',         'fatbjSub2_Wj1DR',         (50,0.,500.)):    fatbjSub2_Wj1DR,
           ('fatbjSub2_Wj1DPhi',       'fatbjSub2_Wj1DPhi',       (50,0.,500.)):    fatbjSub2_Wj1DPhi,
           ('wj1_lepDR',               'wj1_lepDR',               (50,0.,500.)):    wj1_lepDR,
           ('wj1_lepDPhi',             'wj1_lepDPhi',             (50,0.,500.)):    wj1_lepDPhi,
           ('wj2_lepDR',               'wj2_lepDR',               (50,0.,500.)):    wj2_lepDR,
           ('wj2_lepDPhi',             'wj2_lepDPhi',             (50,0.,500.)):    wj2_lepDPhi,
           ('fatbj_wj2DR',             'fatbj_wj2DR',             (50,0.,500.)):    fatbj_wj2DR,
           ('fatbj_wj2DPhi',           'fatbj_wj2DPhi',           (50,0.,500.)):    fatbj_wj2DPhi,
           ('fatbjSub1_wj2DR',         'fatbjSub1_wj2DR',         (50,0.,500.)):    fatbjSub1_wj2DR,
           ('fatbjSub1_wj2DPhi',       'fatbjSub1_wj2DPhi',       (50,0.,500.)):    fatbjSub1_wj2DPhi,
           ('fatbjSub2_wj2DR',         'fatbjSub2_wj2DR',         (50,0.,500.)):    fatbjSub2_wj2DR,
           ('fatbjSub2_wj2DPhi',       'fatbjSub2_wj2DPhi',       (50,0.,500.)):    fatbjSub2_wj2DPhi,
           ('wj1wj2_pt',               'wj1wj2_pt',               (50,0.,500.)):    wj1wj2_pt,
           ('wj1wj2DR',                'wj1wj2DR',                (50,0.,500.)):    wj1wj2DR,
           ('wj1wj2DPhi',              'wj1wj2DPhi',              (50,0.,500.)):    wj1wj2DPhi,
           ('wj1wj2invM',              'wj1wj2invM',              (50,0.,500.)):    wj1wj2invM,
           ('WWplaneAngle',            'WWplaneAngle',            (50,0.,500.)):    WWplaneAngle,
           ('WWplaneAngle_withMET',    'WWplaneAngle_withMET',    (50,0.,500.)):    WWplaneAngle_withMET,
           ('HWplaneAngle',            'HWplaneAngle',            (50,0.,500.)):    HWplaneAngle,
           ('HWW_Mass',                'HWW_Mass',                (50,0.,500.)):    HWW_Mass,
           ('HWW_Simple_Mass',         'HWW_Simple_Mass',         (50,0.,500.)):    HWW_Simple_Mass,
           ('HWW_dR',                  'HWW_dR',                  (50,0.,500.)):    HWW_dR,
           ('cosThetaS_Wjj_simple',    'cosThetaS_Wjj_simple',    (50,0.,1.)):      cosThetaS_Wjj_simple,
           ('cosThetaS_WW_simple_met', 'cosThetaS_WW_simple_met', (50,0.,1.)):      cosThetaS_WW_simple_met,
           ('cosThetaS_HH_simple_met', 'cosThetaS_HH_simple_met', (50,0.,1.)):      cosThetaS_HH_simple_met,
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
           ('VBF_tag',                 'VBF_tag',                 (2, 0, 2)):       VBF_tag,
           ('zeppenfeldVar',           'zeppenfeldVar',           (50,0.,200)):     zeppenfeldVar,           
           ('jetMinDR',                'jetMinDR',                (50,0.,5)):       jetMinDR,           
           ('jetMaxDR',                'jetMaxDR',                (50,0.,5)):       jetMaxDR,           
           ('jetLepMinDR',             'jetLepMinDR',             (50,0.,5)):       jetLepMinDR,           
           ('jetLepMaxDR',             'jetLepMaxDR',             (50,0.,5)):       jetLepMaxDR,           
           ('HT2',                     'HT2',                     (50,0.,300)):     HT2,           
           ('HT2R',                    'HT2R',                    (50,0.,300)):     HT2R,
           ('MT_W1W2',                 'MT_W1W2',                 (50,0.,300)):     MT_W1W2,
           ('HT_taggedJets',           'HT_taggedJets',           (50,0.,300)):     HT_taggedJets,
           ('M_taggedJetsAvg',         'M_taggedJetsAvg',         (50,0.,300)):     M_taggedJetsAvg,
           ('E_taggedJetsAvg',         'E_taggedJetsAvg',         (50,0.,300)):     E_taggedJetsAvg,
           ('JPAcat',                  'JPAcat',                  (3,0.,3)):        JPAcat           
    }


def returnLBNInputs_Boosted(self,lep,jet1=None,jet2=None,jet3=None,jet4=None):
    wj1_Px = op.c_float(0.) if jet3 == None else jet3.p4.Px()
    wj1_Py = op.c_float(0.) if jet3 == None else jet3.p4.Py()
    wj1_Pz = op.c_float(0.) if jet3 == None else jet3.p4.Pz()
    wj1_E  = op.c_float(0.) if jet3 == None else jet3.p4.E()
    wj2_Px = op.c_float(0.) if jet4 == None else jet4.p4.Px()
    wj2_Py = op.c_float(0.) if jet4 == None else jet4.p4.Py()
    wj2_Pz = op.c_float(0.) if jet4 == None else jet4.p4.Pz()
    wj2_E  = op.c_float(0.) if jet4 == None else jet4.p4.E()
    return{('lep_E',   'lep_E',   (50,0.,500.)):    lep.p4.E(),
           ('lep_Px',  'lep_Px',  (50,-250.,250.)): lep.p4.Px(),
           ('lep_Py',  'lep_Py',  (50,-250.,250.)): lep.p4.Py(),
           ('lep_Pz',  'lep_Pz',  (50,-250.,250.)): lep.p4.Pz(),
           ('fatbj_E', 'fatbj_E', (50,0.,500.)):    jet1.p4.E(),
           ('fatbj_Px','bj1_Px',  (50,-250.,250.)): jet1.p4.Px(),
           ('fatbj_Py','bj1_Py',  (50,-250.,250.)): jet1.p4.Py(),
           ('fatbj_Pz','bj1_Pz',  (50,-250.,250.)): jet1.p4.Pz(),
           ('wj1_E',   'wj1_E',   (50,0.,500.)):    wj1_E,
           ('wj1_Px',  'wj1_Px',  (50,-250.,250.)): wj1_Px,
           ('wj1_Py',  'wj1_Py',  (50,-250.,250.)): wj1_Py,
           ('wj1_Pz',  'wj1_Pz',  (50,-250.,250.)): wj1_Pz,
           ('wj2_E',   'wj2_E',   (50,0.,500.)):    wj2_E,
           ('wj2_Px',  'wj2_Px',  (50,-250.,250.)): wj2_Px,
           ('wj2_Py',  'wj2_Py',  (50,-250.,250.)): wj2_Py,
           ('wj2_Pz',  'wj2_Pz',  (50,-250.,250.)): wj2_Pz
    }

### --------------------------------------------------------------------------------------------------- ###

def returnEventNr(self,t):
    return {('evennr','Event number',(100,0.,1e6)): t.event}
    
def inputStaticCast(inputDict,cast='float'):
    return [op.static_cast(cast,v) for v in inputDict.values()]
