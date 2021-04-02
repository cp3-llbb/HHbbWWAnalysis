import os
import sys
from copy import copy

from itertools import chain

import logging
logger = logging.getLogger(__name__) 

from bamboo.analysismodules import SkimmerModule
from bamboo import treefunctions as op
from bamboo.analysisutils import makePileupWeight

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)))) # Add scripts in this directory
from BaseHHtobbWW import BaseNanoHHtobbWW
from selectionDef import *
from highlevelLambdas import *
from JPA import *
from bamboo.root import gbl
import ROOT

#===============================================================================================#
#                                 SkimmerHHtobbWW                                               #
#===============================================================================================#
class SkimmerNanoHHtobbWWSL(BaseNanoHHtobbWW,SkimmerModule):
    """ Plotter module: HH->bbW(->e/µ nu)W(->e/µ nu) histograms from NanoAOD """
    def __init__(self, args):
        super(SkimmerNanoHHtobbWWSL, self).__init__(args)

    def initialize(self):
        super(SkimmerNanoHHtobbWWSL, self).initialize(True) # avoids doing the pseudo-data for skimmer

    def defineSkimSelection(self, t, noSel, sample=None, sampleCfg=None): 
        noSel = super(SkimmerNanoHHtobbWWSL,self).prepareObjects(t, noSel, sample, sampleCfg, "SL", forSkimmer=True)
        # For the Skimmer, SF must not use defineOnFirstUse -> segmentation fault

        era = sampleCfg['era'] 

        #self.datadrivenContributions = {} # Avoid all data-driven estimates

        # Initialize varsToKeep dict #
        varsToKeep = dict()  
        
        if not self.inclusive_sel:
            #----- Check arguments -----#
            jet_level    = ["Res2b2Wj","Res1b3Wj","Res2b1Wj","Res1b2Wj","Hbb2Wj","Hbb1Wj"]

            # Only one lepton_level must be in args and Only one jet_level must be in args
            if [boolean for (level,boolean) in self.args.__dict__.items() if level in jet_level].count(True) != 1:
                raise RuntimeError("Only one of the jet arguments must be used, check --help")

            if self.args.Channel not in ["El","Mu"]:
                raise RuntimeError("Channel must be either 'El' or 'Mu'")            
            
            #----- Lepton selection -----#
            # Args are passed within the self #
            ElSelObj, MuSelObj = makeSingleLeptonSelection(self,noSel,use_dd=False)
            if self.args.Channel == "El":
                selObj = ElSelObj
                lepton = self.electronsTightSel[0]
                
            if self.args.Channel == "Mu":
                selObj = MuSelObj
                lepton = self.muonsTightSel[0]
            
            # -------- Basic Reco Method --------- #
            def chooseWjj(bjet1, bjet2, lep, lightJetPairs, isResolved=False):
                if isResolved : 
                    mH1_mH2_diff = lambda dijet: op.abs((dijet[0].p4+dijet[1].p4+lep.p4+self.corrMET.p4).M() - (self.HLL.bJetCorrP4(bjet1) + self.HLL.bJetCorrP4(bjet2)).M())
                else:
                    mH1_mH2_diff = lambda dijet: op.abs((dijet[0].p4+dijet[1].p4+lep.p4+self.corrMET.p4).M() - self.HLL.bJetCorrP4(bjet1).M())
                findWjj = lambda jetPairs : op.sort(jetPairs, mH1_mH2_diff)
                return findWjj(lightJetPairs)

            # ----------- VBF ------------ #
            def findVBFJets_Resolved(jet1, jet2, jet3, jet4, nAk4JetsSelected):
                if nAk4JetsSelected > 3 : 
                    lambda_findVBFjets = lambda j : op.OR(op.deltaR(jet1.p4, j.p4) > 0.8,
                                                          op.deltaR(jet2.p4, j.p4) > 0.8,
                                                          op.deltaR(jet3.p4, j.p4) > 0.8,
                                                          op.deltaR(jet4.p4, j.p4) > 0.8)
                else :
                    lambda_findVBFjets = lambda j : op.OR(op.deltaR(jet1.p4, j.p4) > 0.8,
                                                          op.deltaR(jet2.p4, j.p4) > 0.8,
                                                          op.deltaR(jet3.p4, j.p4) > 0.8)
                    
                return op.sort(op.combine(op.select(self.VBFJets, lambda_findVBFjets), N=2, pred=self.lambda_VBFPair), 
                               lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))


            def findVBFJets_Boosted(fatJet, jet3, jet4, nAk4JetsSelected):
                if nAk4JetsSelected > 1 : 
                    lambda_findVBFjets = lambda j : op.AND(op.deltaR(fatJet.p4, j.p4) > 1.2, op.OR(op.deltaR(jet3.p4, j.p4) > 0.8,
                                                                                                   op.deltaR(jet4.p4, j.p4) > 0.8))
                else :
                    lambda_findVBFjets = lambda j : op.AND(op.deltaR(fatJet.p4, j.p4) > 1.2, op.deltaR(jet3.p4, j.p4) > 0.8)
                    
                return op.sort(op.combine(op.select(self.VBFJets, lambda_findVBFjets), N=2, pred=self.lambda_VBFPair), 
                               lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))


            #----- Jet selection -----#
            # Since the selections in one line, we can use the non copy option of the selection to modify the selection object internally
            if any([self.args.__dict__[item] for item in ["Res2b2Wj","Res1b3Wj"]]):
                makeCoarseResolvedSelection(self,selObj,4)
                if self.args.Res2b2Wj:
                    makeExclusiveResolvedJetComboSelection(self,selObj,2,4)
                    print('Resolved 2b2Wj')
                    jet1 = self.ak4BJets[0]
                    jet2 = self.ak4BJets[1]
                    # simple basic reco method
                    WjjPairs = chooseWjj(self.ak4BJets[0], self.ak4BJets[1], lepton, op.combine(self.ak4LightJetsByPt, N=2))
                    jet3 = WjjPairs[0][0]
                    jet4 = WjjPairs[0][1]
                    # VBF
                    VBFJetPairs = findVBFJets_Resolved(jet1, jet2, jet3, jet4, 4)
                    

                if self.args.Res1b3Wj:
                    makeExclusiveResolvedJetComboSelection(self,selObj,1,4)
                    print('Resolved 1b3Wj')
                    jet1 = self.ak4BJets[0]
                    jet2 = self.ak4LightJetsByBtagScore[0]
                    WjjPairs = chooseWjj(self.ak4BJets[0], self.ak4LightJetsByBtagScore[0], lepton, self.remainingJetPairs)
                    jet3 = WjjPairs[0][0]
                    jet4 = WjjPairs[0][1]
                    # VBF
                    VBFJetPairs = findVBFJets_Resolved(jet1, jet2, jet3, jet4, 4)


            if any([self.args.__dict__[item] for item in ["Res2b1Wj","Res1b2Wj"]]):
                makeCoarseResolvedSelection(self,selObj,3)
                if self.args.Res2b1Wj:
                    makeExclusiveResolvedJetComboSelection(self,selObj,2,3)
                    print('Resolved 2b1Wj')
                    jet1  = self.ak4BJets[0]
                    jet2  = self.ak4BJets[1]
                    jet3  = self.ak4LightJetsByPt[0]
                    #jet4 = None
                    # VBF
                    VBFJetPairs = findVBFJets_Resolved(jet1, jet2, jet3, None, 3)

                if self.args.Res1b2Wj:
                    makeExclusiveLooseResolvedJetComboSelection(self,selObj,1,3)
                    print('Resolved 1b2Wj')
                    jet1  = self.ak4BJets[0]
                    jet2  = self.ak4LightJetsByBtagScore[0]
                    jet3  = self.remainingJets[0]
                    #jet4 = None
                    # VBF
                    VBFJetPairs = findVBFJets_Resolved(jet1, jet2, jet3, None, 3)


            if any([self.args.__dict__[item] for item in ["Hbb2Wj","Hbb1Wj"]]):
                makeBoostedSelection(self,selObj)
                fatjet  = self.ak8BJets[0]
                if self.args.Hbb1Wj:
                    makeSemiBoostedHbbSelection(self,selObj,1)
                    print('Resolved Hbb1Wj')
                    jet1 = fatjet.subJet1
                    jet2 = fatjet.subJet2
                    jet3 = self.ak4JetsCleanedFromAk8b[0]
                    #jet4 = None
                    # VBF
                    VBFJetPairs = findVBFJets_Boosted(fatjet, jet3, None, 1)

                if self.args.Hbb2Wj:
                    makeSemiBoostedHbbSelection(self,selObj,2)
                    print('Resolved Hbb2Wj')
                    jet1 = fatjet.subJet1
                    jet2 = fatjet.subJet2
                    WjjPairs = chooseWjj(self.ak8BJets[0], None, lepton, op.combine(self.ak4JetsCleanedFromAk8b, N=2))
                    jet3 = WjjPairs[0][0]
                    jet4 = WjjPairs[0][1]
                    # VBF
                    VBFJetPairs = findVBFJets_Boosted(fatjet, jet3, jet4, 2)

        else:
            noSel = self.beforeJetselection(noSel)
            

        #---------------------------------------------------------------------------------------# 
        #                                    Selection tree                                     #
        #---------------------------------------------------------------------------------------#
        #----- EVT variables -----#
        varsToKeep["event"]     = None # Already in tree                                               
        varsToKeep["run"]       = None # Already in tree                                                                                                                                        
        varsToKeep["ls"]        = t.luminosityBlock

        #----- MET variables -----#
        MET = self.corrMET
        
        varsToKeep['METpt']       = MET.pt
        varsToKeep['METphi']      = MET.phi
        varsToKeep['METpx']       = MET.p4.Px()
        varsToKeep['METpy']       = MET.p4.Py()
        varsToKeep['METpz']       = MET.p4.Pz()
        varsToKeep['METenergy']   = MET.p4.E()
        
        varsToKeep['lep_Px']     = lepton.p4.Px()
        varsToKeep['lep_Py']     = lepton.p4.Py()
        varsToKeep['lep_Pz']     = lepton.p4.Pz()
        varsToKeep['lep_E']      = lepton.p4.E()
        varsToKeep['lep_pt']     = lepton.pt
        varsToKeep['lep_eta']    = lepton.eta
        varsToKeep['lep_phi']    = lepton.phi
        varsToKeep['lep_pdgId']  = lepton.pdgId
        varsToKeep['lep_charge'] = lepton.charge


        varsToKeep['lepmet_DPhi'] = op.abs(self.HLL.SinglepMet_dPhi(lepton,MET))
        varsToKeep['lepmet_pt']   = self.HLL.SinglepMet_Pt(lepton,MET)

        varsToKeep['lep_MT']      = self.HLL.MT(lepton,MET)
        ###varsToKeep['MET_LD']      = self.HLL.MET_LD(self.corrMET, self.ak4Jets, self.electronsFakeSel, self.muonsFakeSel)
        varsToKeep['MET_LD']      = self.HLL.MET_LD_DL(self.corrMET, self.ak4Jets, self.electronsFakeSel, self.muonsFakeSel) # highlevelLambdas.py : L-184
        varsToKeep['lep_conept']  = self.HLL.lambdaConePt(lepton)

        varsToKeep['nfakeableMuons']     = op.static_cast("UInt_t",op.rng_len(self.muonsFakeSel))
        varsToKeep['nfakeableElectrons'] = op.static_cast("UInt_t",op.rng_len(self.electronsFakeSel))
        varsToKeep['hT'] = op.static_cast("Float_t", self.HLL.HT_SL(self.ak4Jets)) # highlevelLambdas.py : L-53

        varsToKeep['nAk4Jets']  = op.static_cast("UInt_t",op.rng_len(self.ak4Jets))
        varsToKeep['nAk4BJets'] = op.static_cast("UInt_t",op.rng_len(self.ak4BJets))
        varsToKeep['nAk4BJetsLoose'] = op.static_cast("UInt_t",op.rng_len(self.ak4BJetsLoose))
        varsToKeep['nAk8Jets']  = op.static_cast("UInt_t",op.rng_len(self.ak8Jets))
        varsToKeep['nAk8BJets'] = op.static_cast("UInt_t",op.rng_len(self.ak8BJets))
        varsToKeep['nAk4JetsCleandWithAk8b'] = op.static_cast("UInt_t",op.rng_len(self.ak4JetsCleanedFromAk8b))
        
        #----- Jet variables -----#
        if any([self.args.__dict__[item] for item in ["Res2b2Wj","Res1b3Wj","Res2b1Wj","Res1b2Wj"]]):
            isFullReco = any([self.args.Res2b2Wj, self.args.Res1b3Wj])
            isMissReco = any([self.args.Res2b1Wj, self.args.Res1b2Wj])
            # Jet-1 variables :: Available for each category
            varsToKeep['bj1_Px']             = self.HLL.bJetCorrP4(jet1).Px()
            varsToKeep['bj1_Py']             = self.HLL.bJetCorrP4(jet1).Py()
            varsToKeep['bj1_Pz']             = self.HLL.bJetCorrP4(jet1).Pz()
            varsToKeep['bj1_E']              = self.HLL.bJetCorrP4(jet1).E() 
            varsToKeep['bj1_pt']             = self.HLL.bJetCorrP4(jet1).Pt()
            varsToKeep['bj1_eta']            = self.HLL.bJetCorrP4(jet1).Eta()
            varsToKeep['bj1_phi']            = self.HLL.bJetCorrP4(jet1).Phi()
            varsToKeep['bj1_bTagDeepFlavB']  = jet1.btagDeepFlavB             
            varsToKeep['bj1_bTagCSVV2']      = jet1.btagCSVV2                 
            varsToKeep['bj1LepDR']           = op.deltaR(jet1.p4, lepton.p4)  
            varsToKeep['bj1LepDPhi']         = op.abs(op.deltaPhi(jet1.p4, lepton.p4))
            varsToKeep['bj1MetDPhi']         = op.abs(self.HLL.SinglepMet_dPhi(jet1, MET))
            varsToKeep['minDR_lep_allJets']  = self.HLL.mindr_lep1_jet(lepton, self.ak4Jets)
            
            # Jet-2 variables :: Available for 2b2Wj, 2b1Wj. 2b0Wj
            varsToKeep['bj2_Px']            = self.HLL.bJetCorrP4(jet2).Px()
            varsToKeep['bj2_Py']            = self.HLL.bJetCorrP4(jet2).Py()
            varsToKeep['bj2_Pz']            = self.HLL.bJetCorrP4(jet2).Pz()
            varsToKeep['bj2_E']             = self.HLL.bJetCorrP4(jet2).E() 
            varsToKeep['bj2_pt']            = self.HLL.bJetCorrP4(jet2).Pt()
            varsToKeep['bj2_eta']           = self.HLL.bJetCorrP4(jet2).Eta()
            varsToKeep['bj2_phi']           = self.HLL.bJetCorrP4(jet2).Phi()
            varsToKeep['bj2_bTagDeepFlavB'] = jet2.btagDeepFlavB
            varsToKeep['bj2_bTagCSVV2']     = jet2.btagCSVV2
            varsToKeep['bj2LepDR']          = op.deltaR(jet2.p4, lepton.p4)
            varsToKeep['bj2LepDPhi']        = op.abs(op.deltaPhi(jet2.p4, lepton.p4))
            varsToKeep['bj2MetDPhi']        = op.abs(self.HLL.SinglepMet_dPhi(jet2, MET))

            # Di-Jet Variables
            varsToKeep['bj1bj2_pt']         = (self.HLL.bJetCorrP4(jet1)+self.HLL.bJetCorrP4(jet2)).Pt()
            varsToKeep['bj1bj2_M']          = op.invariant_mass(self.HLL.bJetCorrP4(jet1),self.HLL.bJetCorrP4(jet2))
            varsToKeep['cosThetaS_Hbb']     = op.static_cast("Float_t", self.HLL.comp_cosThetaS(self.HLL.bJetCorrP4(jet1), self.HLL.bJetCorrP4(jet2))) 
            varsToKeep['mT_top_3particle']  = op.min(self.HLL.mT2(self.HLL.bJetCorrP4(jet1),lepton.p4 ,MET.p4), 
                                                     self.HLL.mT2(self.HLL.bJetCorrP4(jet2), lepton.p4, MET.p4))
            
            # Jet-3 variables :: Available for 2b2Wj, 2b1Wj, 1b2Wj, 1b1Wj
            varsToKeep['wj1_Px']            = jet3.p4.Px()
            varsToKeep['wj1_Py']            = jet3.p4.Py()
            varsToKeep['wj1_Pz']            = jet3.p4.Pz()
            varsToKeep['wj1_E']             = jet3.p4.E()
            varsToKeep['wj1_pt']            = jet3.pt
            varsToKeep['wj1_eta']           = jet3.eta
            varsToKeep['wj1_phi']           = jet3.phi
            varsToKeep['wj1_bTagDeepFlavB'] = jet3.btagDeepFlavB
            varsToKeep['wj1_bTagCSVV2']     = jet3.btagCSVV2
            varsToKeep['wj1LepDR']          = op.deltaR(jet3.p4, lepton.p4)
            varsToKeep['wj1LepDPhi']        = op.abs(op.deltaPhi(jet3.p4, lepton.p4))
            varsToKeep['wj1MetDPhi']        = op.abs(self.HLL.SinglepMet_dPhi(jet3, MET))

            # Jet-4 variables :: Available for 2b2Wj, 1b2Wj
            varsToKeep['wj2_Px']  = jet4.p4.Px()                                            if isFullReco else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj2_Py']  = jet4.p4.Py()                                            if isFullReco else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj2_Pz']  = jet4.p4.Pz()                                            if isFullReco else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj2_E']   = jet4.p4.E()                                             if isFullReco else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj2_pt']  = jet4.pt                                                 if isFullReco else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj2_eta'] = jet4.eta                                                if isFullReco else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj2_phi'] = jet4.phi                                                if isFullReco else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj2_bTagDeepFlavB'] = jet4.btagDeepFlavB                            if isFullReco else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj2_bTagCSVV2']     = jet4.btagCSVV2                                if isFullReco else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj2LepDR']          = op.deltaR(jet4.p4, lepton.p4)                 if isFullReco else op.static_cast("Float_t",op.c_float(0.)) 
            varsToKeep['wj2LepDPhi']  = op.abs(op.deltaPhi(jet4.p4, lepton.p4))             if isFullReco else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj2MetDPhi']  = op.abs(self.HLL.SinglepMet_dPhi(jet4, MET))         if isFullReco else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj1wj2_pt']   = (jet3.p4 + jet4.p4).Pt()                            if isFullReco else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj1wj2_M']    = op.invariant_mass(jet3.p4, jet4.p4)                 if isFullReco else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['w1w2_MT']     = self.HLL.MT_W1W2_ljj(lepton,jet3,jet4,MET)          if isFullReco else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['HWW_Mass']    = self.HLL.HWW_simple(jet3.p4,jet4.p4,lepton.p4,MET).M() if isFullReco else op.static_cast("Float_t",op.c_float(0.)) # w.r.t. neutrino
            varsToKeep['HWW_Simple_Mass'] = self.HLL.HWW_met_simple(jet3.p4,jet4.p4,lepton.p4,MET.p4).M() if isFullReco else op.static_cast("Float_t",op.c_float(0.)) # not w.r.t. neutrino
            varsToKeep['HWW_dR'] = self.HLL.dR_Hww(jet3.p4,jet4.p4,lepton.p4,MET)           if isFullReco else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['mHH']    = (jet3.p4+jet4.p4+lepton.p4+MET.p4+self.HLL.bJetCorrP4(jet1)+self.HLL.bJetCorrP4(jet2)).M() if isFullReco else op.static_cast("Float_t",op.c_float(0.))


            varsToKeep['cosThetaS_Wjj_simple']    = op.static_cast("Float_t", self.HLL.comp_cosThetaS(jet3.p4, jet4.p4))  if isFullReco else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['cosThetaS_WW_simple_met'] = op.static_cast("Float_t", self.HLL.comp_cosThetaS(self.HLL.Wjj_simple(jet3.p4,jet4.p4), self.HLL.Wlep_met_simple(lepton.p4, MET.p4)))  \
                                                    if isFullReco else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['cosThetaS_HH_simple_met'] = op.static_cast("Float_t", self.HLL.comp_cosThetaS(self.HLL.bJetCorrP4(jet1)+self.HLL.bJetCorrP4(jet2), 
                                                                                                      self.HLL.HWW_met_simple(jet3.p4,jet4.p4,lepton.p4,MET.p4))) \
                if isFullReco else op.static_cast("Float_t",op.c_float(0.))

            # neu p4 only for 1b2Wj and 2b2Wj
            varsToKeep["neuPx"] = self.HLL.neuP4(jet3.p4+jet4.p4+lepton.p4, MET).Px()  if isFullReco else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep["neuPy"] = self.HLL.neuP4(jet3.p4+jet4.p4+lepton.p4, MET).Py()  if isFullReco else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep["neuPz"] = self.HLL.neuP4(jet3.p4+jet4.p4+lepton.p4, MET).Pz()  if isFullReco else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep["neuE"]  = self.HLL.neuP4(jet3.p4+jet4.p4+lepton.p4, MET).E()   if isFullReco else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep["neuPt"] = self.HLL.neuP4(jet3.p4+jet4.p4+lepton.p4, MET).Pt()  if isFullReco else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep["angleBetWWPlane"] = self.HLL.angleWWplane(lepton.p4, MET, jet3.p4, jet4.p4) if isFullReco else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep["angleBetHWPlane"] = self.HLL.angleBetPlanes(jet1.p4,jet2.p4,jet3.p4,jet4.p4) if isFullReco else op.static_cast("Float_t",op.c_float(0.))

            # diJet dR & dPhi
            varsToKeep['bj1bj2_DR']   = op.deltaR(jet1.p4,jet2.p4)
            varsToKeep['bj1bj2_DPhi'] = op.abs(op.deltaPhi(jet1.p4,jet2.p4))
            varsToKeep['bj2wj1_DR']   = op.deltaR(jet2.p4,jet3.p4)
            varsToKeep['bj2wj1_DPhi'] = op.abs(op.deltaPhi(jet2.p4,jet3.p4))
            varsToKeep['wj1wj2_DR']   = op.deltaR(jet3.p4,jet4.p4)             if isFullReco else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj1wj2_DPhi'] = op.abs(op.deltaPhi(jet3.p4,jet4.p4))   if isFullReco else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['bj1wj2_DR']   = op.deltaR(jet1.p4,jet4.p4)             if isFullReco else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['bj1wj2_DPhi'] = op.abs(op.deltaPhi(jet1.p4,jet4.p4))   if isFullReco else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['bj1wj1_DR']   = op.deltaR(jet1.p4,jet3.p4)         
            varsToKeep['bj1wj1_DPhi'] = op.abs(op.deltaPhi(jet1.p4,jet3.p4))

            # L-562, L-635 [VBF]
            varsToKeep['VBFj1pt']        = op.switch(op.rng_len(VBFJetPairs) > 0, 
                                                     op.switch(VBFJetPairs[0][0].pt > VBFJetPairs[0][1].pt,
                                                               VBFJetPairs[0][0].pt, VBFJetPairs[0][1].pt),
                                                     op.static_cast("Float_t",op.c_float(0.)))
            varsToKeep['VBFj2pt']        = op.switch(op.rng_len(VBFJetPairs) > 0, 
                                                     op.switch(VBFJetPairs[0][0].pt > VBFJetPairs[0][1].pt,
                                                               VBFJetPairs[0][1].pt, VBFJetPairs[0][0].pt),
                                                     op.static_cast("Float_t",op.c_float(0.)))
            varsToKeep['VBFj1eta']       = op.switch(op.rng_len(VBFJetPairs) > 0, 
                                                     op.switch(VBFJetPairs[0][0].pt > VBFJetPairs[0][1].pt,
                                                               VBFJetPairs[0][0].eta, VBFJetPairs[0][1].eta),
                                                     op.static_cast("Float_t",op.c_float(0.)))
            varsToKeep['VBFj2eta']       = op.switch(op.rng_len(VBFJetPairs) > 0, 
                                                     op.switch(VBFJetPairs[0][0].pt > VBFJetPairs[0][1].pt,
                                                               VBFJetPairs[0][1].eta, VBFJetPairs[0][0].eta),
                                                     op.static_cast("Float_t",op.c_float(0.)))
            varsToKeep['VBFj1j2dEta']    = op.switch(op.rng_len(VBFJetPairs) > 0,op.abs(VBFJetPairs[0][0].eta - VBFJetPairs[0][1].eta),
                                                     op.static_cast("Float_t",op.c_float(0.)))
            
            varsToKeep['VBFj1j2dPhi']    = op.switch(op.rng_len(VBFJetPairs) > 0,op.abs(op.deltaPhi(VBFJetPairs[0][0].p4,VBFJetPairs[0][1].p4)),
                                                     op.static_cast("Float_t",op.c_float(0.)))
            
            varsToKeep['VBFj1j2invM']    = op.switch(op.rng_len(VBFJetPairs) > 0,op.invariant_mass(VBFJetPairs[0][0].p4,  VBFJetPairs[0][1].p4),
                                                     op.static_cast("Float_t",op.c_float(0.)))
            
            varsToKeep['VBF_tag']        = op.c_int(op.rng_len(VBFJetPairs) > 0)
            varsToKeep['zeppenfeldVar']  = op.switch(op.rng_len(VBFJetPairs) > 0, 
                                                     op.abs(lepton.eta - (VBFJetPairs[0][0].eta + VBFJetPairs[0][1].eta)/2.0)/op.abs(VBFJetPairs[0][0].eta - VBFJetPairs[0][1].eta),
                                                     op.static_cast("Float_t",op.c_float(0.)))
                
            if self.args.Res2b2Wj or self.args.Res1b3Wj: 
                varsToKeep['minJetDR']       = self.HLL.MinDiJetDRTight(jet1,jet2,jet3,jet4)
                varsToKeep['minLepJetDR']    = self.HLL.MinDR_lep4j(lepton,jet1,jet2,jet3,jet4)
                varsToKeep['HT2_lepJetMet']  = self.HLL.HT2_2b2Wj(lepton,jet1,jet2,jet3,jet4,MET)
                varsToKeep['HT2R_lepJetMet'] = self.HLL.HT2R_2b2Wj(lepton,jet1,jet2,jet3,jet4,MET)

            if self.args.Res2b1Wj or self.args.Res1b2Wj: 
                varsToKeep['minJetDR']       = self.HLL.MinDiJetDRLoose(jet1,jet2,jet3)
                varsToKeep['minLepJetDR']    = self.HLL.MinDR_lep3j(lepton,jet1,jet2,jet3)
                varsToKeep['HT2_lepJetMet']  = self.HLL.HT2_2b1Wj(lepton,jet1,jet2,jet3,MET)
                varsToKeep['HT2R_lepJetMet'] = self.HLL.HT2R_2b1Wj(lepton,jet1,jet2,jet3,MET)


        #----- Fatjet variables -----#
        if any([self.args.__dict__[item] for item in ["Hbb2Wj","Hbb1Wj"]]):
                
            varsToKeep['fatbj_Px']             = fatjet.p4.Px()
            varsToKeep['fatbj_Py']             = fatjet.p4.Py()
            varsToKeep['fatbj_Pz']             = fatjet.p4.Pz()
            varsToKeep['fatbj_E']              = fatjet.p4.E()
            varsToKeep['fatbj_pt']             = fatjet.pt
            varsToKeep['fatbj_eta']            = fatjet.eta
            varsToKeep['fatbj_phi']            = fatjet.phi
            varsToKeep['fatbj_sub1pt']         = op.switch(jet1.btagDeepB > jet2.btagDeepB, jet1.pt, jet2.pt)
            varsToKeep['fatbj_sub1eta']        = op.switch(jet1.btagDeepB > jet2.btagDeepB, jet1.eta, jet2.eta)
            varsToKeep['fatbj_sub2pt']         = op.switch(jet1.btagDeepB > jet2.btagDeepB, jet2.pt, jet1.pt)
            varsToKeep['fatbj_sub2eta']        = op.switch(jet1.btagDeepB > jet2.btagDeepB, jet2.eta, jet1.eta)
            varsToKeep['fatbj_softdropMass']   = fatjet.msoftdrop
            # deepBtags
            varsToKeep['fatbj_btagDDBvL']      = fatjet.btagDDBvL
            varsToKeep['fatbj_btagDDBvL_noMD'] = fatjet.btagDDBvL_noMD
            varsToKeep['fatbj_btagDDCvB']      = fatjet.btagDDCvB
            varsToKeep['fatbj_btagDDCvB_noMD'] = fatjet.btagDDCvB_noMD
            varsToKeep['fatbj_btagDDCvL']      = fatjet.btagDDCvL
            varsToKeep['fatbj_btagDDCvL_noMD'] = fatjet.btagDDCvL_noMD
            varsToKeep['fatbj_btagDeepB']      = fatjet.btagDeepB

            
            varsToKeep['cosThetaS_Hbb']     = op.static_cast("Float_t", self.HLL.comp_cosThetaS(jet1.p4, jet2.p4))
            varsToKeep['mT_top_3particle']  = op.min(self.HLL.mT2(jet1.p4,lepton.p4,MET.p4), self.HLL.mT2(jet2.p4,lepton.p4,MET.p4))


            varsToKeep['wj1_Px']      = jet3.p4.Px()           
            varsToKeep['wj1_Py']      = jet3.p4.Py()           
            varsToKeep['wj1_Pz']      = jet3.p4.Pz()           
            varsToKeep['wj1_E']       = jet3.p4.E()            
            varsToKeep['wj1_pt']      = jet3.pt                
            varsToKeep['wj1_eta']     = jet3.eta               
            varsToKeep['wj1_bTagDeepFlavB'] = jet3.btagDeepFlavB
            varsToKeep['wj1_bTagCSVV2']     = jet3.btagCSVV2    
            varsToKeep['wj2_Px']      = jet4.p4.Px()               if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj2_Py']      = jet4.p4.Py()               if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj2_Pz']      = jet4.p4.Pz()               if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj2_E']       = jet4.p4.E()                if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj2_pt']      = jet4.pt                    if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj2_eta']     = jet4.eta                   if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj2_bTagDeepFlavB'] = jet4.btagDeepFlavB   if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj2_bTagCSVV2']     = jet4.btagCSVV2       if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))

            varsToKeep['fatbj_lepDR']       = op.deltaR(fatjet.p4, lepton.p4)
            varsToKeep['fatbjSub1_lepDR']   = op.switch(jet1.btagDeepB > jet2.btagDeepB, op.deltaR(jet1.p4, lepton.p4), op.deltaR(jet2.p4, lepton.p4))
            varsToKeep['fatbjSub2_lepDR']   = op.switch(jet1.btagDeepB > jet2.btagDeepB, op.deltaR(jet2.p4, lepton.p4), op.deltaR(jet1.p4, lepton.p4))
            varsToKeep['fatbj_lepDPhi']     = op.abs(op.deltaPhi(fatjet.p4, lepton.p4))
            varsToKeep['fatbjSub1_lepDPhi'] = op.switch(jet1.btagDeepB > jet2.btagDeepB, op.abs(op.deltaPhi(jet1.p4, lepton.p4)), op.abs(op.deltaPhi(jet2.p4, lepton.p4)))
            varsToKeep['fatbjSub2_lepDPhi'] = op.switch(jet1.btagDeepB > jet2.btagDeepB, op.abs(op.deltaPhi(jet2.p4, lepton.p4)), op.abs(op.deltaPhi(jet1.p4, lepton.p4)))
            varsToKeep['minSubJetLepDR']    = op.min(op.deltaR(jet1.p4, lepton.p4), op.deltaR(jet2.p4, lepton.p4))

            varsToKeep['fatbj_Wj1DR']       = op.deltaR(fatjet.p4, jet3.p4)           
            varsToKeep['fatbj_Wj1DPhi']     = op.abs(op.deltaPhi(fatjet.p4, jet3.p4)) 
            varsToKeep['fatbjSub1_Wj1DR']   = op.switch(jet1.btagDeepB > jet2.btagDeepB, 
                                                        op.deltaR(jet1.p4, jet3.p4), 
                                                        op.deltaR(jet2.p4, jet3.p4))  
            varsToKeep['fatbjSub1_Wj1DPhi'] = op.switch(jet1.btagDeepB > jet2.btagDeepB, 
                                                        op.abs(op.deltaPhi(jet1.p4, jet3.p4)), 
                                                        op.abs(op.deltaPhi(jet2.p4, jet3.p4)))
            varsToKeep['fatbjSub2_Wj1DR']   = op.switch(jet1.btagDeepB > jet2.btagDeepB,
                                                        op.deltaR(jet2.p4, jet3.p4),
                                                        op.deltaR(jet1.p4, jet3.p4))          
            varsToKeep['fatbjSub2_Wj1DPhi'] = op.switch(jet1.btagDeepB > jet2.btagDeepB,
                                                        op.abs(op.deltaPhi(jet2.p4, jet3.p4)),
                                                        op.abs(op.deltaPhi(jet1.p4, jet3.p4)))
            varsToKeep['wj1_lepDR']         = op.deltaR(jet3.p4, lepton.p4)                   
            varsToKeep['wj1_lepDPhi']       = op.abs(op.deltaPhi(jet3.p4, lepton.p4))         
            varsToKeep['wj2_lepDR']         = op.deltaR(jet4.p4, lepton.p4)             if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj2_lepDPhi']       = op.abs(op.deltaPhi(jet4.p4, lepton.p4))   if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))

            varsToKeep['fatbj_wj2DR']       = op.deltaR(fatjet.p4, jet4.p4)             if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.)) 
            varsToKeep['fatbj_wj2DPhi']     = op.abs(op.deltaPhi(fatjet.p4, jet4.p4))   if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['fatbjSub1_wj2DR']   = op.deltaR(jet1.p4, jet4.p4)               if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['fatbjSub1_wj2DPhi'] = op.abs(op.deltaPhi(jet1.p4, jet4.p4))     if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['fatbjSub2_wj2DR']   = op.deltaR(jet2.p4, jet4.p4)               if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['fatbjSub2_wj2DPhi'] = op.abs(op.deltaPhi(jet2.p4, jet4.p4))     if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))

            varsToKeep['wj1wj2_pt']    = (jet3.p4+jet4.p4).Pt()                    if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj1wj2DR']     = op.deltaR(jet3.p4, jet4.p4)               if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj1wj2DPhi']   = op.abs(op.deltaPhi(jet3.p4, jet4.p4))     if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj1wj2invM']   = op.invariant_mass(jet3.p4,jet4.p4)        if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))
            #varsToKeep["angleBetWWPlane"] = self.HLL.angleWWplane(lepton.p4, MET, jet3.p4, jet4.p4) if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep["HWplaneAngle"] = self.HLL.angleBetPlanes(jet1.p4,jet2.p4,jet3.p4,jet4.p4) if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))


            varsToKeep['HWW_Mass']        = self.HLL.HWW_simple(jet3.p4,jet4.p4,lepton.p4,MET).M()        if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['HWW_Simple_Mass'] = self.HLL.HWW_met_simple(jet3.p4,jet4.p4,lepton.p4,MET.p4).M() if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['HWW_dR']          = self.HLL.dR_Hww(jet3.p4,jet4.p4,lepton.p4,MET)                if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['cosThetaS_Wjj_simple']    = op.static_cast("Float_t", self.HLL.comp_cosThetaS(jet3.p4, jet4.p4))             if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['cosThetaS_WW_simple_met'] = op.static_cast("Float_t", self.HLL.comp_cosThetaS(self.HLL.Wjj_simple(jet3.p4,jet4.p4), 
                                                                                                      self.HLL.Wlep_met_simple(lepton.p4, MET.p4))) \
                if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['cosThetaS_HH_simple_met'] = op.static_cast("Float_t", self.HLL.comp_cosThetaS(jet1.p4 + jet2.p4, 
                                                                                                      self.HLL.HWW_met_simple(jet3.p4,jet4.p4,lepton.p4,MET.p4))) \
                if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))
                

            varsToKeep['VBFj1pt']        = op.switch(op.rng_len(VBFJetPairs) > 0, 
                                                     op.switch(VBFJetPairs[0][0].pt > VBFJetPairs[0][1].pt,
                                                               VBFJetPairs[0][0].pt, VBFJetPairs[0][1].pt),
                                                     op.static_cast("Float_t",op.c_float(0.)))
            varsToKeep['VBFj2pt']        = op.switch(op.rng_len(VBFJetPairs) > 0, 
                                                     op.switch(VBFJetPairs[0][0].pt > VBFJetPairs[0][1].pt,
                                                               VBFJetPairs[0][1].pt, VBFJetPairs[0][0].pt),
                                                     op.static_cast("Float_t",op.c_float(0.)))

            varsToKeep['VBFj1eta']       = op.switch(op.rng_len(VBFJetPairs) > 0, 
                                                     op.switch(VBFJetPairs[0][0].pt > VBFJetPairs[0][1].pt,
                                                               VBFJetPairs[0][0].eta, VBFJetPairs[0][1].eta),
                                                     op.static_cast("Float_t",op.c_float(0.)))
            varsToKeep['VBFj2eta']       = op.switch(op.rng_len(VBFJetPairs) > 0, 
                                                     op.switch(VBFJetPairs[0][0].pt > VBFJetPairs[0][1].pt,
                                                               VBFJetPairs[0][1].eta, VBFJetPairs[0][0].eta),
                                                     op.static_cast("Float_t",op.c_float(0.)))

            varsToKeep['VBFj1j2dEta']    = op.switch(op.rng_len(VBFJetPairs) > 0,op.abs(VBFJetPairs[0][0].eta - VBFJetPairs[0][1].eta),
                                                     op.static_cast("Float_t",op.c_float(0.)))

            varsToKeep['VBFj1j2dPhi']    = op.switch(op.rng_len(VBFJetPairs) > 0,op.abs(op.deltaPhi(VBFJetPairs[0][0].p4,VBFJetPairs[0][1].p4)),
                                                     op.static_cast("Float_t",op.c_float(0.)))

            varsToKeep['VBFj1j2invM']    = op.switch(op.rng_len(VBFJetPairs) > 0,op.invariant_mass(VBFJetPairs[0][0].p4,  VBFJetPairs[0][1].p4),
                                                     op.static_cast("Float_t",op.c_float(0.)))
            
            varsToKeep['VBF_tag']         = op.c_int(op.rng_len(VBFJetPairs)>0)
            varsToKeep['zeppenfeldVar']   = op.switch(op.rng_len(VBFJetPairs) > 0, 
                                                      op.abs(lepton.eta - (VBFJetPairs[0][0].eta + VBFJetPairs[0][1].eta)/2.0)/op.abs(VBFJetPairs[0][0].eta - VBFJetPairs[0][1].eta),
                                                      op.static_cast("Float_t",op.c_float(0.)))
            

            if self.args.Hbb2Wj:
                varsToKeep['jetMinDR']    = self.HLL.MinDiJetDRTight(jet1,jet2,jet3,jet4)
                varsToKeep['jetLepMinDR'] = self.HLL.MinDR_lep4j(lepton, jet1, jet2, jet3, jet4) 
                varsToKeep['HT2']         = self.HLL.HT2_2b2Wj(lepton, jet1, jet2,jet3,jet4,MET)
                varsToKeep['HT2R']        = self.HLL.HT2R_2b2Wj(lepton, jet1, jet2, jet3, jet4,MET)
                varsToKeep['MT_W1W2']     = self.HLL.MT_W1W2_ljj(lepton,jet3,jet4,MET)

            if self.args.Hbb1Wj:
                varsToKeep['jetMinDR']    = self.HLL.MinDiJetDRLoose(jet1,jet2,jet3)
                varsToKeep['jetLepMinDR'] = self.HLL.MinDR_lep3j(lepton,jet1,jet2,jet3) 
                varsToKeep['HT2']         = self.HLL.HT2_2b1Wj(lepton,jet1,jet2,jet3,MET)
                varsToKeep['HT2R']        = self.HLL.HT2R_2b1Wj(lepton,jet1,jet2,jet3,MET)
                varsToKeep['MT_W1W2']     = self.HLL.MT_W1W2_lj(lepton,jet3,MET)

                #----- Additional variables -----#
        varsToKeep["MC_weight"]    = t.genWeight
        varsToKeep['total_weight'] = selObj.sel.weight

        #return leptonSel.sel, varsToKeep
        return selObj.sel, varsToKeep
