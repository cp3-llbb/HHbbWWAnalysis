import os
import sys
from copy import copy
from operator import mul
from functools import reduce

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
from DDHelper import DataDrivenFake, DataDrivenDY
from variableMakerSL import *
#from bamboo.root import loadHeader
#loadHeader("jpa.h")

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

        # Initialize varsToKeep dict #
        varsToKeep = dict()  

        #---------------------------------------------------------------------------------------# 
        #                                     Selections                                        #
        #---------------------------------------------------------------------------------------#
        # keep the exact same order of nodes as mentioned in respective JPA model xml files
        ResolvedJPANodeList = ['2b2Wj','2b1Wj','1b2Wj','2b0Wj','1b1Wj','1b0Wj','0b']
        BoostedJPANodeList  = ['Hbb2Wj','Hbb1Wj','Hbb0Wj']

        # JPA Models
        basepath = os.path.join(os.path.abspath(os.path.dirname(__file__)),'MachineLearning','ml-models','JPA_Loose_ttH')
        resolvedModelDict = getResolvedJpaModelDict(basepath, ResolvedJPANodeList, era)
        boostedModelDict  = getBoostedJpaModelDict(basepath, BoostedJPANodeList, era)

        if not self.inclusive_sel:
            #----- Check arguments -----#
            jet_level    = ["Ak4","Ak8","Res2b2Wj","Res2b1Wj","Res2b0Wj","Res1b2Wj","Res1b1Wj","Res1b0Wj","Res0b","Hbb2Wj","Hbb1Wj","Hbb0Wj"]

            if [boolean for (level,boolean) in self.args.__dict__.items() if level in jet_level].count(True) != 1:
                raise RuntimeError("Only one of the jet arguments must be used, check --help")

            if self.args.Channel not in ["El","Mu"]:
                raise RuntimeError("Channel must be either 'El' or 'Mu'")            
            
            #----- Lepton selection -----#
            ElSelObj, MuSelObj = makeSingleLeptonSelection(self,noSel,use_dd=False,fake_selection=self.args.FakeCR)

            if self.args.Channel is None:
                raise RuntimeError("You need to specify --Channel")
            if self.args.Channel == "El":
                selObj = ElSelObj
                lep = self.electronsTightSel[0]
            if self.args.Channel == "Mu":
                selObj = MuSelObj
                lep = self.muonsTightSel[0]
            
            #----- Jet selection -----#
            if any([self.args.__dict__[item] for item in ["Ak4","Res2b2Wj","Res2b1Wj","Res2b0Wj","Res1b2Wj","Res1b1Wj","Res1b0Wj","Res0b"]]):
                makeResolvedSelection(self,selObj)
                if self.args.Channel == "El":
                    print('... Resolved :: El Channel')
                    L1out,L2out,selObjAndJetsPerJpaCatDict = findJPACategoryResolved (self, selObj, lep, 
                                                                                      self.muonsPreSel, self.electronsPreSel, 
                                                                                      self.ak4Jets, self.ak4BJetsLoose, self.ak4BJets, self.corrMET, 
                                                                                      resolvedModelDict, t.event, self.HLL, ResolvedJPANodeList, 
                                                                                      plot_yield=False)
                                        
                if self.args.Channel == "Mu":
                    print('... Resolved :: Mu Channel')
                    L1out,L2out,selObjAndJetsPerJpaCatDict = findJPACategoryResolved (self, selObj, lep, 
                                                                                      self.muonsPreSel, self.electronsPreSel, 
                                                                                      self.ak4Jets, self.ak4BJetsLoose, self.ak4BJets, self.corrMET, 
                                                                                      resolvedModelDict, t.event, self.HLL, ResolvedJPANodeList, 
                                                                                      plot_yield=False)
                    
            if any([self.args.__dict__[item] for item in ["Ak8","Hbb2Wj","Hbb1Wj","Hbb0Wj"]]):
                makeBoostedSelection(self,selObj)
                if self.args.Channel == "El":
                    print('... Boosted :: El Channel')
                    L1out,L2out,selObjAndJetsPerJpaCatDict = findJPACategoryBoosted (self, selObj, lep, 
                                                                                     self.muonsPreSel, self.electronsPreSel, 
                                                                                     self.ak8BJets, self.ak4JetsCleanedFromAk8b, self.ak4BJetsLoose, 
                                                                                     self.ak4BJets, self.corrMET, boostedModelDict, t.event, self.HLL, 
                                                                                     BoostedJPANodeList, plot_yield=False)

                if self.args.Channel == "Mu":
                    print('... Boosted :: Mu Channel')
                    L1out,L2out,selObjAndJetsPerJpaCatDict = findJPACategoryBoosted (self, selObj, lep, 
                                                                                     self.muonsPreSel, self.electronsPreSel, 
                                                                                     self.ak8BJets, self.ak4JetsCleanedFromAk8b, self.ak4BJetsLoose, 
                                                                                     self.ak4BJets, self.corrMET, boostedModelDict, t.event, self.HLL, 
                                                                                     BoostedJPANodeList, plot_yield=False)
                    
            if self.args.Res2b2Wj:
                print("...... 2b2Wj")
                selObj   = selObjAndJetsPerJpaCatDict.get('2b2Wj')[0]
                jpaJets  = selObjAndJetsPerJpaCatDict.get('2b2Wj')[1]
                jpaarg   = "Res2b2Wj"

            if self.args.Res2b1Wj:
                print("...... 2b1Wj")
                selObj  = selObjAndJetsPerJpaCatDict.get('2b1Wj')[0]
                jpaJets = selObjAndJetsPerJpaCatDict.get('2b1Wj')[1]
                jpaarg   = "Res2b1Wj"

            if self.args.Res1b2Wj:
                print("...... 1b2Wj")
                selObj  = selObjAndJetsPerJpaCatDict.get('1b2Wj')[0]
                jpaJets = selObjAndJetsPerJpaCatDict.get('1b2Wj')[1]
                jpaarg   = "Res1b2Wj"

            if self.args.Res2b0Wj:
                print("...... 2b0Wj")
                selObj  = selObjAndJetsPerJpaCatDict.get('2b0Wj')[0]
                jpaJets = selObjAndJetsPerJpaCatDict.get('2b0Wj')[1]
                jpaarg   = "Res2b0Wj"

            if self.args.Res1b1Wj:
                print("...... 1b1Wj")
                selObj  = selObjAndJetsPerJpaCatDict.get('1b1Wj')[0]
                jpaJets = selObjAndJetsPerJpaCatDict.get('1b1Wj')[1]
                jpaarg   = "Res1b1Wj"

            if self.args.Res1b0Wj:
                print("...... 1b0Wj")
                selObj  = selObjAndJetsPerJpaCatDict.get('1b0Wj')[0]
                jpaJets = selObjAndJetsPerJpaCatDict.get('1b0Wj')[1]
                jpaarg   = "Res1b0Wj"

            if self.args.Res0b:
                print("...... 0b")
                selObj  = selObjAndJetsPerJpaCatDict.get('0b')[0]
                jpaJets = None 
                jpaarg   = "Res0b"

            #######################################
            # Hbb2Wj : jet1 jet2   jet3   jet4
            # Hbb1Wj : jet1 jet2   jet3   jet4=0
            # Hbb0Wj : jet1 jet2   jet3=0 jet4=0
            #######################################
            if self.args.Hbb2Wj:
                print("...... Hbb2Wj")
                selObj  = selObjAndJetsPerJpaCatDict.get('Hbb2Wj')[0]
                jpaJets = selObjAndJetsPerJpaCatDict.get('Hbb2Wj')[1]
                jpaarg   = "Hbb2Wj"

            if self.args.Hbb1Wj:
                print("...... Hbb1Wj")
                selObj  = selObjAndJetsPerJpaCatDict.get('Hbb1Wj')[0]
                jpaJets = selObjAndJetsPerJpaCatDict.get('Hbb1Wj')[1]
                jpaarg   = "Hbb1Wj"

            if self.args.Hbb0Wj:
                print("...... Hbb0Wj")
                selObj  = selObjAndJetsPerJpaCatDict.get('Hbb0Wj')[0]
                jpaJets = None
                jpaarg   = "Hbb0Wj"

        else:
            noSel = self.beforeJetselection(noSel)
            
        #---------------------------------------------------------------------------------------# 
        #                                 Synchronization tree                                  #
        #---------------------------------------------------------------------------------------#
        if self.args.Synchronization:
            # Event variables #
            varsToKeep["event"]             = None # Already in tree
            varsToKeep["run"]               = None # Already in tree 
            varsToKeep["ls"]                = t.luminosityBlock
            varsToKeep["n_presel_mu"]       = op.static_cast("UInt_t",op.rng_len(self.muonsPreSel))
            varsToKeep["n_fakeablesel_mu"]  = op.static_cast("UInt_t",op.rng_len(self.muonsFakeSel))
            varsToKeep["n_mvasel_mu"]       = op.static_cast("UInt_t",op.rng_len(self.muonsTightSel))
            varsToKeep["n_presel_ele"]      = op.static_cast("UInt_t",op.rng_len(self.electronsPreSel))
            varsToKeep["n_fakeablesel_ele"] = op.static_cast("UInt_t",op.rng_len(self.electronsFakeSel))
            varsToKeep["n_mvasel_ele"]      = op.static_cast("UInt_t",op.rng_len(self.electronsTightSel))
            varsToKeep["n_presel_ak4Jet"]   = op.static_cast("UInt_t",op.rng_len(self.ak4Jets))
            varsToKeep["n_presel_ak8Jet"]   = op.static_cast("UInt_t",op.rng_len(self.ak8Jets))    
            varsToKeep["n_presel_ak8BJet"]  = op.static_cast("UInt_t",op.rng_len(self.ak8BJets))    
            varsToKeep["n_loose_ak4BJet"]   = op.static_cast("UInt_t",op.rng_len(self.ak4BJetsLoose))    
            varsToKeep["n_medium_ak4BJet"]  = op.static_cast("UInt_t",op.rng_len(self.ak4BJets))
            varsToKeep["n_ak4JetsCleanAk8b"]= op.static_cast("UInt_t",op.rng_len(self.ak4JetsCleanedFromAk8b))
            varsToKeep["n_presel_ak4JetVBF"]= op.static_cast("UInt_t",op.rng_len(self.VBFJetsPreSel))
 
            varsToKeep["is_SR"]             = op.static_cast("UInt_t",op.OR(op.rng_len(self.electronsTightSel)==1,
                                                                            op.rng_len(self.muonsTightSel)==1))

            varsToKeep["is_e"]          = op.c_float(True)  if self.args.Channel == 'El' else op.c_float(False)
            varsToKeep["is_m"]          = op.c_float(False) if self.args.Channel == 'El' else op.c_float(True)

            varsToKeep["is_resolved"]   = op.switch(op.AND(op.rng_len(self.ak4Jets)>=3,
                                                           op.rng_len(self.ak4BJets) >= 1,
                                                           op.rng_len(self.ak8BJets)==0), op.c_bool(True), op.c_bool(False))
            varsToKeep["is_boosted"]    = op.switch(op.AND(op.rng_len(self.ak8BJets)>=1,op.rng_len(self.ak4JetsCleanedFromAk8b)>=1), op.c_bool(True), op.c_bool(False))
            
            varsToKeep["n_tau"]         = op.static_cast("UInt_t", op.rng_len(self.tauCleanSel))

            varsToKeep['resolved_tag']  = op.static_cast("UInt_t",op.AND(op.rng_len(self.ak4Jets) >= 3,
                                                                         op.rng_len(self.ak4BJets) >= 1,
                                                                         op.rng_len(self.ak8BJets) == 0))
            varsToKeep['boosted_tag']   = op.static_cast("UInt_t",op.AND(op.rng_len(self.ak8BJets) >= 1,op.rng_len(self.ak4JetsCleanedFromAk8b) >= 1))

            # Triggers #
            '''
            varsToKeep["triggers"]                  = self.triggers
            varsToKeep["triggers_SingleElectron"]   = op.OR(*self.triggersPerPrimaryDataset['SingleElectron'])
            varsToKeep["triggers_SingleMuon"]       = op.OR(*self.triggersPerPrimaryDataset['SingleMuon'])
            varsToKeep["triggers_DoubleElectron"]   = op.OR(*self.triggersPerPrimaryDataset['DoubleEGamma'])
            varsToKeep["triggers_DoubleMuon"]       = op.OR(*self.triggersPerPrimaryDataset['DoubleMuon'])
            varsToKeep["triggers_MuonElectron"]     = op.OR(*self.triggersPerPrimaryDataset['MuonEG'])
            '''
            # Muons #
            
            for i in range(1,3): # 2 leading muons
                varsToKeep["mu{}_pt".format(i)]                    = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].pt, op.c_float(-9999., "float"))
                varsToKeep["mu{}_eta".format(i)]                   = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].eta, op.c_float(-9999.))
                varsToKeep["mu{}_phi".format(i)]                   = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].phi, op.c_float(-9999.))
                varsToKeep["mu{}_E".format(i)]                     = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].p4.E(), op.c_float(-9999., "float"))
                varsToKeep["mu{}_charge".format(i)]                = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].charge, op.c_int(-9999.))
                varsToKeep["mu{}_conept".format(i)]                = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muon_conept[self.muonsPreSel[i-1].idx], op.c_float(-9999.))
                varsToKeep["mu{}_miniRelIso".format(i)]            = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].miniPFRelIso_all, op.c_float(-9999.))
                varsToKeep["mu{}_PFRelIso04".format(i)]            = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].pfRelIso04_all, op.c_float(-9999.))
                varsToKeep["mu{}_jetNDauChargedMVASel".format(i)]  = op.c_float(-9999.)
                varsToKeep["mu{}_jetPtRel".format(i)]              = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].jetPtRelv2, op.c_float(-9999.))
                varsToKeep["mu{}_jetRelIso".format(i)]             = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].jetRelIso, op.c_float(-9999.))
                varsToKeep["mu{}_jetDeepJet".format(i)]            = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].jet.btagDeepFlavB, op.c_float(-9999.))
                varsToKeep["mu{}_sip3D".format(i)]                 = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].sip3d, op.c_float(-9999.))
                varsToKeep["mu{}_dxy".format(i)]                   = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].dxy, op.c_float(-9999.))
                varsToKeep["mu{}_dxyAbs".format(i)]                = op.switch(op.rng_len(self.muonsPreSel) >= i, op.abs(self.muonsPreSel[i-1].dxy), op.c_float(-9999.))
                varsToKeep["mu{}_dz".format(i)]                    = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].dz, op.c_float(-9999.))
                varsToKeep["mu{}_segmentCompatibility".format(i)]  = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].segmentComp, op.c_float(-9999.))
                varsToKeep["mu{}_leptonMVA".format(i)]             = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].mvaTTH, op.c_float(-9999.))
                varsToKeep["mu{}_mediumID".format(i)]              = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].mediumId, op.c_float(-9999.,"Bool_t"))
                varsToKeep["mu{}_dpt_div_pt".format(i)]            = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].tunepRelPt, op.c_float(-9999.))  # Not sure
                varsToKeep["mu{}_isfakeablesel".format(i)]         = op.switch(op.rng_len(self.muonsPreSel) >= i, op.switch(self.lambda_muonFakeSel(self.muonsPreSel[i-1]), op.c_int(1), op.c_int(0)), op.c_int(-9999))
                varsToKeep["mu{}_ismvasel".format(i)]              = op.switch(op.rng_len(self.muonsPreSel) >= i, op.switch(op.AND(self.lambda_muonTightSel(self.muonsPreSel[i-1]), self.lambda_muonFakeSel(self.muonsPreSel[i-1])), op.c_int(1), op.c_int(0)), op.c_int(-9999)) # mvasel encompasses fakeablesel
                varsToKeep["mu{}_isGenMatched".format(i)]          = op.switch(op.rng_len(self.muonsPreSel) >= i, op.switch(self.lambda_is_matched(self.muonsPreSel[i-1]), op.c_int(1), op.c_int(0)), op.c_int(-9999))
                varsToKeep["mu{}_genPartFlav".format(i)]           = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].genPartFlav, op.c_int(-9999))
                varsToKeep["mu{}_FR".format(i)]                    = op.switch(op.rng_len(self.muonsPreSel) >= i, self.lambda_FR_mu(self.muonsPreSel[i-1]), op.c_int(-9999))
                varsToKeep["mu{}_FRcorr".format(i)]                = op.switch(op.rng_len(self.muonsPreSel) >= i, self.lambda_FRcorr_mu(self.muonsPreSel[i-1]), op.c_int(-9999))
                varsToKeep["mu{}_FF".format(i)]                    = op.switch(op.rng_len(self.muonsPreSel) >= i, self.lambda_FF_mu(self.muonsPreSel[i-1]), op.c_int(-9999))
                varsToKeep["mu{}_looseSF".format(i)]               = op.switch(op.rng_len(self.muonsPreSel) >= i, reduce(mul,self.lambda_MuonLooseSF(self.muonsPreSel[i-1])), op.c_int(-9999))
                varsToKeep["mu{}_tightSF".format(i)]               = op.switch(op.rng_len(self.muonsPreSel) >= i, reduce(mul,self.lambda_MuonTightSF(self.muonsPreSel[i-1])), op.c_int(-9999))
            
            # Electrons #
            for i in range(1,3): # 2 leading electrons 
                varsToKeep["ele{}_pt".format(i)]                    = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].pt, op.c_float(-9999.))
                varsToKeep["ele{}_eta".format(i)]                   = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].eta, op.c_float(-9999.))
                varsToKeep["ele{}_phi".format(i)]                   = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].phi, op.c_float(-9999.))
                varsToKeep["ele{}_E".format(i)]                     = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].p4.E(), op.c_float(-9999.,))
                varsToKeep["ele{}_charge".format(i)]                = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].charge, op.c_int(-9999.))
                varsToKeep["ele{}_conept".format(i)]                = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electron_conept[self.electronsPreSel[i-1].idx], op.c_float(-9999.))
                varsToKeep["ele{}_miniRelIso".format(i)]            = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].miniPFRelIso_all, op.c_float(-9999.))
                varsToKeep["ele{}_PFRelIso03".format(i)]            = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].pfRelIso03_all, op.c_float(-9999.)) # Iso03, Iso04 not in NanoAOD
                varsToKeep["ele{}_jetNDauChargedMVASel".format(i)]  = op.c_float(-9999.)
                varsToKeep["ele{}_jetPtRel".format(i)]              = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].jetPtRelv2, op.c_float(-9999.))
                varsToKeep["ele{}_jetRelIso".format(i)]             = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].jetRelIso, op.c_float(-9999.))
                varsToKeep["ele{}_jetDeepJet".format(i)]            = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].jet.btagDeepFlavB, op.c_float(-9999.))
                varsToKeep["ele{}_sip3D".format(i)]                 = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].sip3d, op.c_float(-9999.))

                varsToKeep["ele{}_dxy".format(i)]                   = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].dxy, op.c_float(-9999.))
                varsToKeep["ele{}_dxyAbs".format(i)]                = op.switch(op.rng_len(self.electronsPreSel) >= i, op.abs(self.electronsPreSel[i-1].dxy), op.c_float(-9999.))
                varsToKeep["ele{}_dz".format(i)]                    = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].dz, op.c_float(-9999.))
                varsToKeep["ele{}_ntMVAeleID".format(i)]            = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].mvaFall17V2noIso, op.c_float(-9999.))
                varsToKeep["ele{}_leptonMVA".format(i)]             = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].mvaTTH, op.c_float(-9999.))
                varsToKeep["ele{}_passesConversionVeto".format(i)]  = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].convVeto, op.c_float(-9999.,"Bool_t"))
                varsToKeep["ele{}_nMissingHits".format(i)]          = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].lostHits, op.c_float(-9999.,"UChar_t"))
                varsToKeep["ele{}_sigmaEtaEta".format(i)]           = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].sieie, op.c_float(-9999.))
                varsToKeep["ele{}_HoE".format(i)]                   = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].hoe, op.c_float(-9999.))
                varsToKeep["ele{}_OoEminusOoP".format(i)]           = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].eInvMinusPInv, op.c_float(-9999.))
                varsToKeep["ele{}_isfakeablesel".format(i)]         = op.switch(op.rng_len(self.electronsPreSel) >= i, op.switch(self.lambda_electronFakeSel(self.electronsPreSel[i-1]), op.c_int(1), op.c_int(0)), op.c_int(-9999))
                varsToKeep["ele{}_ismvasel".format(i)]              = op.switch(op.rng_len(self.electronsPreSel) >= i, op.switch(op.AND(self.lambda_electronTightSel(self.electronsPreSel[i-1]), self.lambda_electronFakeSel(self.electronsPreSel[i-1])), op.c_int(1), op.c_int(0)), op.c_int(-9999)) # mvasel encompasses fakeablesel
                varsToKeep["ele{}_isGenMatched".format(i)]          = op.switch(op.rng_len(self.electronsPreSel) >= i, op.switch(self.lambda_is_matched(self.electronsPreSel[i-1]), op.c_int(1), op.c_int(0)), op.c_int(-9999))
                varsToKeep["ele{}_genPartFlav".format(i)]           = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].genPartFlav, op.c_int(-9999))
                varsToKeep["ele{}_deltaEtaSC".format(i)]            = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].deltaEtaSC, op.c_int(-9999))
                varsToKeep["ele{}_FR".format(i)]                    = op.switch(op.rng_len(self.electronsPreSel) >= i, self.lambda_FR_el(self.electronsPreSel[i-1]), op.c_int(-9999))
                varsToKeep["ele{}_FRcorr".format(i)]                = op.switch(op.rng_len(self.electronsPreSel) >= i, self.lambda_FRcorr_el(self.electronsPreSel[i-1]), op.c_int(-9999))
                varsToKeep["ele{}_FF".format(i)]                    = op.switch(op.rng_len(self.electronsPreSel) >= i, self.lambda_FF_el(self.electronsPreSel[i-1]), op.c_int(-9999))
                varsToKeep["ele{}_looseSF".format(i)]               = op.switch(op.rng_len(self.electronsPreSel) >= i, reduce(mul,self.lambda_ElectronLooseSF(self.electronsPreSel[i-1])), op.c_int(-9999))
                varsToKeep["ele{}_tightSF".format(i)]               = op.switch(op.rng_len(self.electronsPreSel) >= i, reduce(mul,self.lambda_ElectronTightSF(self.electronsPreSel[i-1])), op.c_int(-9999))

            # AK4 Jets #
            for i in range(1,5): # 4 leading jets 
                varsToKeep["ak4Jet{}_pt".format(i)]                 = op.switch(op.rng_len(self.ak4Jets) >= i, self.ak4Jets[i-1].pt, op.c_float(-9999.))
                varsToKeep["ak4Jet{}_eta".format(i)]                = op.switch(op.rng_len(self.ak4Jets) >= i, self.ak4Jets[i-1].eta, op.c_float(-9999.))
                varsToKeep["ak4Jet{}_phi".format(i)]                = op.switch(op.rng_len(self.ak4Jets) >= i, self.ak4Jets[i-1].phi, op.c_float(-9999.))
                varsToKeep["ak4Jet{}_E".format(i)]                  = op.switch(op.rng_len(self.ak4Jets) >= i, self.ak4Jets[i-1].p4.E(), op.c_float(-9999.))
                varsToKeep["ak4Jet{}_CSV".format(i)]                = op.switch(op.rng_len(self.ak4Jets) >= i, self.ak4Jets[i-1].btagDeepFlavB, op.c_float(-9999.))
                varsToKeep["ak4Jet{}_hadronFlavour".format(i)]      = op.switch(op.rng_len(self.ak4Jets) >= i, self.ak4Jets[i-1].hadronFlavour, op.c_float(-9999.))
                varsToKeep["ak4Jet{}_btagSF".format(i)]             = op.switch(op.rng_len(self.ak4Jets) >= i, self.DeepJetDiscReshapingSF(self.ak4Jets[i-1]), op.c_float(-9999.))
                varsToKeep["ak4Jet{}_puid_eff".format(i)]           = op.switch(op.rng_len(self.ak4Jets) >= i, self.jetpuid_mc_eff(self.ak4Jets[i-1]), op.c_float(-9999.))
                varsToKeep["ak4Jet{}_puid_sfeff".format(i)]         = op.switch(op.rng_len(self.ak4Jets) >= i, self.jetpuid_sf_eff(self.ak4Jets[i-1]), op.c_float(-9999.))
                varsToKeep["ak4Jet{}_puid_mis".format(i)]           = op.switch(op.rng_len(self.ak4Jets) >= i, self.jetpuid_mc_mis(self.ak4Jets[i-1]), op.c_float(-9999.))
                varsToKeep["ak4Jet{}_puid_sfmis".format(i)]         = op.switch(op.rng_len(self.ak4Jets) >= i, self.jetpuid_sf_mis(self.ak4Jets[i-1]), op.c_float(-9999.))

            # VBF Jets #
            for i in range(1,6): # 5 leading jets
                varsToKeep["ak4JetVBF{}_pt".format(i)]              = op.switch(op.rng_len(self.VBFJetsPreSel) >= i, self.VBFJetsPreSel[i-1].pt, op.c_float(-9999.))
                varsToKeep["ak4JetVBF{}_eta".format(i)]             = op.switch(op.rng_len(self.VBFJetsPreSel) >= i, self.VBFJetsPreSel[i-1].eta, op.c_float(-9999.))
                varsToKeep["ak4JetVBF{}_phi".format(i)]             = op.switch(op.rng_len(self.VBFJetsPreSel) >= i, self.VBFJetsPreSel[i-1].phi, op.c_float(-9999.))
                varsToKeep["ak4JetVBF{}_E".format(i)]               = op.switch(op.rng_len(self.VBFJetsPreSel) >= i, self.VBFJetsPreSel[i-1].p4.E(), op.c_float(-9999.))
                varsToKeep["ak4JetVBF{}_CSV".format(i)]             = op.switch(op.rng_len(self.VBFJetsPreSel) >= i, self.VBFJetsPreSel[i-1].btagDeepFlavB, op.c_float(-9999.))
                varsToKeep["ak4JetVBF{}_btagSF".format(i)]          = op.switch(op.rng_len(self.VBFJetsPreSel) >= i, self.DeepJetDiscReshapingSF(self.VBFJetsPreSel[i-1]), op.c_float(-9999.))

            if not self.inclusive_sel:
                if not any([self.args.Res0b, self.args.Ak4, self.args.Ak8]) : 
                    varsToKeep["ak4JetVBFPair1_pt"]              = op.switch(op.rng_len(VBFJetPairsJPA) >= 1, VBFJetPairsJPA[0][0].pt, op.c_float(-9999.))
                    varsToKeep["ak4JetVBFPair1_eta"]             = op.switch(op.rng_len(VBFJetPairsJPA) >= 1, VBFJetPairsJPA[0][0].eta, op.c_float(-9999.))
                    varsToKeep["ak4JetVBFPair1_phi"]             = op.switch(op.rng_len(VBFJetPairsJPA) >= 1, VBFJetPairsJPA[0][0].phi, op.c_float(-9999.))
                    varsToKeep["ak4JetVBFPair1_E"]               = op.switch(op.rng_len(VBFJetPairsJPA) >= 1, VBFJetPairsJPA[0][0].p4.E(), op.c_float(-9999.))
                    varsToKeep["ak4JetVBFPair1_CSV"]             = op.switch(op.rng_len(VBFJetPairsJPA) >= 1, VBFJetPairsJPA[0][0].btagDeepFlavB, op.c_float(-9999.))
                    varsToKeep["ak4JetVBFPair1_btagSF"]          = op.switch(op.rng_len(VBFJetPairsJPA) >= 1, self.DeepJetDiscReshapingSF(VBFJetPairsJPA[0][0]), op.c_float(-9999.))
                    varsToKeep["ak4JetVBFPair2_pt"]              = op.switch(op.rng_len(VBFJetPairsJPA) >= 1, VBFJetPairsJPA[0][1].pt, op.c_float(-9999.))
                    varsToKeep["ak4JetVBFPair2_eta"]             = op.switch(op.rng_len(VBFJetPairsJPA) >= 1, VBFJetPairsJPA[0][1].eta, op.c_float(-9999.))
                    varsToKeep["ak4JetVBFPair2_phi"]             = op.switch(op.rng_len(VBFJetPairsJPA) >= 1, VBFJetPairsJPA[0][1].phi, op.c_float(-9999.))
                    varsToKeep["ak4JetVBFPair2_E"]               = op.switch(op.rng_len(VBFJetPairsJPA) >= 1, VBFJetPairsJPA[0][1].p4.E(), op.c_float(-9999.))
                    varsToKeep["ak4JetVBFPair2_CSV"]             = op.switch(op.rng_len(VBFJetPairsJPA) >= 1, VBFJetPairsJPA[0][1].btagDeepFlavB, op.c_float(-9999.))
                    varsToKeep["ak4JetVBFPair2_btagSF"]          = op.switch(op.rng_len(VBFJetPairsJPA) >= 1, self.DeepJetDiscReshapingSF(VBFJetPairsJPA[0][1]), op.c_float(-9999.))
                else : 
                    varsToKeep["ak4JetVBFPair1_pt"]              = op.c_float(-9999.)
                    varsToKeep["ak4JetVBFPair1_eta"]             = op.c_float(-9999.)
                    varsToKeep["ak4JetVBFPair1_phi"]             = op.c_float(-9999.)
                    varsToKeep["ak4JetVBFPair1_E"]               = op.c_float(-9999.)
                    varsToKeep["ak4JetVBFPair1_CSV"]             = op.c_float(-9999.)
                    varsToKeep["ak4JetVBFPair1_btagSF"]          = op.c_float(-9999.)
                    varsToKeep["ak4JetVBFPair2_pt"]              = op.c_float(-9999.)
                    varsToKeep["ak4JetVBFPair2_eta"]             = op.c_float(-9999.)
                    varsToKeep["ak4JetVBFPair2_phi"]             = op.c_float(-9999.)
                    varsToKeep["ak4JetVBFPair2_E"]               = op.c_float(-9999.)
                    varsToKeep["ak4JetVBFPair2_CSV"]             = op.c_float(-9999.)
                    varsToKeep["ak4JetVBFPair2_btagSF"]          = op.c_float(-9999.)
                    
            # AK8 Jets #
            for i in range(1,3): # 2 leading fatjets 
                varsToKeep["ak8Jet{}_pt".format(i)]                 = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].pt, op.c_float(-9999.))
                varsToKeep["ak8Jet{}_eta".format(i)]                = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].eta, op.c_float(-9999.))
                varsToKeep["ak8Jet{}_phi".format(i)]                = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].phi, op.c_float(-9999.))
                varsToKeep["ak8Jet{}_E".format(i)]                  = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].p4.E(), op.c_float(-9999.))
                varsToKeep["ak8Jet{}_msoftdrop".format(i)]          = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].msoftdrop, op.c_float(-9999.))
                varsToKeep["ak8Jet{}_tau1".format(i)]               = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].tau1, op.c_float(-9999.))
                varsToKeep["ak8Jet{}_tau2".format(i)]               = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].tau2, op.c_float(-9999.))
                varsToKeep["ak8Jet{}_subjet0_pt".format(i)]         = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].subJet1.pt, op.c_float(-9999.))
                varsToKeep["ak8Jet{}_subjet0_eta".format(i)]        = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].subJet1.eta, op.c_float(-9999.))
                varsToKeep["ak8Jet{}_subjet0_phi".format(i)]        = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].subJet1.phi, op.c_float(-9999.))
                varsToKeep["ak8Jet{}_subjet0_CSV".format(i)]        = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].subJet1.btagDeepB, op.c_float(-9999.))
                varsToKeep["ak8Jet{}_subjet1_pt".format(i)]         = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].subJet2.pt, op.c_float(-9999.))
                varsToKeep["ak8Jet{}_subjet1_eta".format(i)]        = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].subJet2.eta, op.c_float(-9999.))
                varsToKeep["ak8Jet{}_subjet1_phi".format(i)]        = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].subJet2.phi, op.c_float(-9999.))
                varsToKeep["ak8Jet{}_subjet1_CSV".format(i)]        = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].subJet2.btagDeepB, op.c_float(-9999.))

            varsToKeep["PFMET"]    = self.corrMET.pt
            varsToKeep["PFMETphi"] = self.corrMET.phi
            varsToKeep["met1_E"]   = self.corrMET.p4.E()
            varsToKeep["met1_pt"]  = self.corrMET.pt
            varsToKeep["met1_eta"] = self.corrMET.eta
            varsToKeep["met1_phi"] = self.corrMET.phi

            # SF #
            electronMuon_cont        = op.combine((self.electronsFakeSel, self.muonsFakeSel))
            
            varsToKeep["trigger_SF"] = op.multiSwitch(
                (op.AND(op.rng_len(self.electronsTightSel)==1,op.rng_len(self.muonsTightSel)==0) , self.ttH_singleElectron_trigSF(self.electronsTightSel[0])),
                (op.AND(op.rng_len(self.electronsTightSel)==0,op.rng_len(self.muonsTightSel)==1) , self.ttH_singleMuon_trigSF(self.muonsTightSel[0])),
                (op.AND(op.rng_len(self.electronsTightSel)>=2,op.rng_len(self.muonsTightSel)==0) , self.lambda_ttH_doubleElectron_trigSF(self.electronsTightSel)),
                (op.AND(op.rng_len(self.electronsTightSel)==0,op.rng_len(self.muonsTightSel)>=2) , self.lambda_ttH_doubleMuon_trigSF(self.muonsTightSel)),
                (op.AND(op.rng_len(self.electronsTightSel)>=1,op.rng_len(self.muonsTightSel)>=1) , self.lambda_ttH_electronMuon_trigSF(electronMuon_cont[0])),
                op.c_float(1.))
            
            if not self.inclusive_sel:
                varsToKeep["weight_trigger_el_sf"] = op.switch(op.rng_len(self.electronsTightSel)>0, self.ttH_singleElectron_trigSF(lep),op.c_float(1.))
                varsToKeep["weight_trigger_mu_sf"] = op.switch(op.rng_len(self.muonsTightSel)>0, self.ttH_singleMuon_trigSF(lep),op.c_float(1.))
                
                varsToKeep["lepton_IDSF"] = op.rng_product(self.electronsFakeSel, lambda el : reduce(mul,self.lambda_ElectronLooseSF(el)+self.lambda_ElectronTightSF(el))) * \
                                            op.rng_product(self.muonsFakeSel, lambda mu : reduce(mul,self.lambda_MuonLooseSF(mu)+self.lambda_MuonTightSF(mu))) 
                
                varsToKeep["lepton_IDSF_recoToLoose"]  = op.rng_product(self.electronsFakeSel, lambda el : reduce(mul,self.lambda_ElectronLooseSF(el))) * \
                                                        op.rng_product(self.muonsFakeSel, lambda mu : reduce(mul,self.lambda_MuonLooseSF(mu)))
                varsToKeep["lepton_IDSF_looseToTight"] = op.rng_product(self.electronsFakeSel, lambda el : reduce(mul,self.lambda_ElectronTightSF(el))) * \
                                                         op.rng_product(self.muonsFakeSel, lambda mu : reduce(mul,self.lambda_MuonTightSF(mu)))
                if era == "2016" or era == "2017": 
                    if self.args.Channel == "El":
                        varsToKeep["weight_electron_reco_low"]    = op.switch(op.AND(self.lambda_is_matched(lep),lep.pt<=20.), self.elLooseRecoPtLt20(lep), op.c_float(1.))
                        varsToKeep["weight_electron_reco_high"]   = op.switch(op.AND(self.lambda_is_matched(lep),lep.pt>20.), self.elLooseRecoPtGt20(lep), op.c_float(1.))
                        varsToKeep["weight_muon_idiso_loose"]     = op.c_float(1.)
                        varsToKeep["weight_electron_id_loose_01"] = op.switch(self.lambda_is_matched(lep), self.elLooseEff(lep), op.c_float(1.))
                        varsToKeep["weight_electron_id_loose_02"] = op.switch(self.lambda_is_matched(lep), self.elLooseId(lep),  op.c_float(1.))
                        varsToKeep["weight_electron_tth_loose"]   = self.lambda_ElectronTightSF(lep)[0]
                        varsToKeep["weight_muon_tth_loose"]       = op.c_float(1.)

                    if self.args.Channel == "Mu":
                        varsToKeep["weight_muon_idiso_loose"]     = op.switch(self.lambda_is_matched(lep), self.muLooseId(lep), op.c_float(1.))
                        varsToKeep["weight_electron_reco_low"]    = op.c_float(1.)
                        varsToKeep["weight_electron_reco_high"]   = op.c_float(1.)
                        varsToKeep["weight_electron_id_loose_01"] = op.c_float(1.)
                        varsToKeep["weight_electron_id_loose_02"] = op.c_float(1.)
                        varsToKeep["weight_electron_tth_loose"]   = op.c_float(1.)
                        varsToKeep["weight_muon_tth_loose"]       = self.lambda_MuonTightSF(lep)[0]
                else:
                    raise NotImplementedError

            # L1 Prefire #
            if era in ["2016","2017"]:
                varsToKeep["L1prefire"] = self.L1Prefiring
                varsToKeep["weight_l1_ecal_prefiring"] = self.L1Prefiring
            else:
                varsToKeep["L1prefire"] = op.c_float(-9999.)
                varsToKeep["weight_l1_ecal_prefiring"] = op.c_float(-9999.)

            # Fake rate #
            if self.args.Channel == "El":
                varsToKeep["fakeRate"]                  = op.switch(self.lambda_electronTightSel(self.electronsFakeSel[0]), self.ElFakeFactor(self.electronsFakeSel[0]), op.c_float(1.))
                varsToKeep["weight_fake_electrons"]     = op.switch(self.lambda_electronTightSel(self.electronsFakeSel[0]), op.abs(self.ElFakeFactor(self.electronsFakeSel[0])), op.c_float(1.))
                varsToKeep["weight_fake_muons"]         = op.c_float(1.)
                varsToKeep["weight_fake_two_non_tight"] = op.c_float(999.0)
            if self.args.Channel == "Mu":
                varsToKeep["fakeRate"]                  = op.switch(self.lambda_muonTightSel(self.muonsFakeSel[0]), self.MuFakeFactor(self.muonsFakeSel[0]), op.c_float(1.))
                varsToKeep["weight_fake_electrons"]     = op.c_float(1.)
                varsToKeep["weight_fake_muons"]         = op.switch(self.lambda_muonTightSel(self.muonsFakeSel[0]), op.abs(self.MuFakeFactor(self.muonsFakeSel[0])), op.c_float(1.))
                varsToKeep["weight_fake_two_non_tight"] = op.c_float(999.0)

            if self.is_MC:
                varsToKeep["weight_fake_is_mc"] = op.c_float(-1.)
            else:
                varsToKeep["weight_fake_is_mc"] = op.c_float(1.)

            # PU ID SF #
            varsToKeep["PU_jetID_SF"]                = self.puid_reweighting
            varsToKeep["weight_jet_PUid_efficiency"] = self.puid_reweighting_efficiency
            varsToKeep["weight_jet_PUid_mistag"]     = self.puid_reweighting_mistag

            # Btagging SF #
            varsToKeep["btag_SF"]           = self.btagAk4SF
            varsToKeep["weight_btagWeight"] = self.btagAk4SF
            if "BtagRatioWeight" in self.__dict__.keys():
                varsToKeep["btag_ratio_SF"]   = self.BtagRatioWeight
                varsToKeep["weight_btagNorm"] = self.BtagRatioWeight

            # PS weights #
            varsToKeep["weight_PSWeight_ISR"] = self.psISRSyst
            varsToKeep["weight_PSWeight_FSR"] = self.psFSRSyst

            # ttbar PT reweighting #
            if "group" in sampleCfg and sampleCfg["group"] == 'ttbar':
                varsToKeep["topPt_wgt"] = self.ttbar_weight(self.genTop[0],self.genAntitop[0])

           # Event Weight #
            if self.is_MC:
                varsToKeep["MC_weight"] = t.genWeight
                puWeightsFile = os.path.join(os.path.dirname(__file__), "data", "pileup",sample+'_%s.json'%era)
                #puWeightsFile = os.path.join(os.path.dirname(__file__), "data" , "pileup", sampleCfg["pufile"])
                varsToKeep["PU_weight"] = makePileupWeight(puWeightsFile, t.Pileup_nTrueInt, nameHint=f"puweightFromFile{sample}".replace('-','_'))
                varsToKeep["eventWeight"] = noSel.weight if self.inclusive_sel else selObj.sel.weight

            if self.inclusive_sel:
                return noSel, varsToKeep
            else:
                return selObj.sel, varsToKeep

        #---------------------------------------------------------------------------------------# 
        #                                    Selection tree                                     #
        #---------------------------------------------------------------------------------------#
        #----- EVT variables -----#
        varsToKeep["event"]     = None # Already in tree                                               
        varsToKeep["run"]       = None # Already in tree             
        varsToKeep["ls"]        = t.luminosityBlock

        genVarsList = self.HLL.comp_cosThetaSbetBeamAndHiggs(t.GenPart)
        varsToKeep['mHH_gen']        = genVarsList[0]
        varsToKeep['consTheta1_gen'] = genVarsList[1]
        varsToKeep['consTheta2_gen'] = genVarsList[2]

        if any([self.args.__dict__[item] for item in ["Ak4","Res2b2Wj","Res2b1Wj","Res2b0Wj","Res1b2Wj","Res1b1Wj","Res1b0Wj","Res0b"]]):
            commonInputs_Resolved  = returnCommonInputs_Resolved  (self = self)
            classicInputs_Resolved = returnClassicInputs_Resolved (self = self, 
                                                                   lepton = lep, 
                                                                   jpaSelectedJets = jpaJets, 
                                                                   L1out = L1out, L2out = L2out, 
                                                                   jpaArg = jpaarg)
            inputs_resolved        = {**commonInputs_Resolved, **classicInputs_Resolved}

            for (varname,_,_),var in inputs_resolved.items():
                varsToKeep[varname] = var 

        #----- Fatjet variables -----#
        if any([self.args.__dict__[item] for item in ["Ak8","Hbb2Wj","Hbb1Wj","Hbb0Wj"]]):
            commonInputs_Boosted  = returnCommonInputs_Boosted  (self = self)
            classicInputs_Boosted = returnClassicInputs_Boosted (self = self, 
                                                                 lepton = lep, 
                                                                 jpaSelectedJets = jpaJets, 
                                                                 L1out = L1out, L2out = L2out,
                                                                 jpaArg = jpaarg)
            inputs_boosted        = {**commonInputs_Boosted, **classicInputs_Boosted}
            
            for (varname,_,_),var in inputs_boosted.items():
                varsToKeep[varname] = var 

        #----- Additional variables -----#
        varsToKeep["MC_weight"]    = t.genWeight
        varsToKeep['total_weight'] = selObj.sel.weight

        #return leptonSel.sel, varsToKeep
        return selObj.sel, varsToKeep
