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

if not hasattr(gbl, "bamboo_printEntry"):
    gbl.gInterpreter.Declare("""
    bool bamboo_printEntry(long entry, float val) {
      std::cout << "Processing entry #" << entry << ": val " << val << std::endl;
      return true;
    }""")
    gbl.gInterpreter.Declare("""
    bool bamboo_printComb3_(long entry, const rdfhelpers::Combination<3>& comb, std::size_t idx, double value, UInt_t nJets, UInt_t nComb, const ROOT::VecOps::RVec<rdfhelpers::Combination<3>>& combs, const ROOT::VecOps::RVec<std::size_t>& rng) {
      std::cout << "Combination #" << idx << " (value " << value << ") for " << entry << ": " << comb.get(0) << ", " << comb.get(1) << ", " << comb.get(2) << ", nJets=" << nJets << ", nSelJets=" << rng.size() << ", nCombinations=" << nComb << std::endl;
      if ( ( comb.get(0) >= nJets ) || ( comb.get(1) >= nJets ) || ( comb.get(2) >= nJets ) ) {
        std::cout << "All combinations: " << std::endl;
        for ( const auto& icomb : combs ) {
          std::cout << " - " << icomb.get(0) << ", " << icomb.get(1) << ", " << icomb.get(2) << std::endl;
        }
        std::cout << "From range with size " << rng.size() << ": ";
        for ( auto idx : rng) {
          std::cout << " " << idx;
        }
        std::cout << std::endl;
      }
      return true;
    }""")

def addPrint(nd, funName, *args):
    callPrint = op.extMethod(funName, returnType="bool")(*args)
    filterStr = nd(callPrint.op)
    logger.debug(f"Adding printout with {filterStr}")
    nd.df = nd.df.Filter(filterStr)

def addPrintout(selection, funName, *args):
    from bamboo import treefunctions as op
    from bamboo import treeproxies as _tp
    from bamboo import dataframebackend
    be = selection._fbe
    if selection.name not in be.selDFs:
        raise RuntimeError("This method will only work with a dynamically constructed RDataFrame")
    nd = be.selDFs[selection.name]
    callPrint = op.extMethod(funName, returnType=_tp.boolType)(*args)
    filterStr = nd(callPrint.op)
    logger.debug(f"Adding printout with {filterStr}")
    nd.df = nd.df.Filter(filterStr)

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

        #---------------------------------------------------------------------------------------# 
        #                                     Selections                                        #
        #---------------------------------------------------------------------------------------#
        # JPA Models
        basepath = os.path.join(os.path.abspath(os.path.dirname(__file__)),'MachineLearning','ml-models','JPA')
        resolvedModelDict = dict()
        resolvedModelDict['2b2Wj'] = [os.path.join(basepath, 'bb1l_jpa_4jet_resolved_even.xml'), 
                                      os.path.join(basepath, 'bb1l_jpa_4jet_resolved_odd.xml')]
        resolvedModelDict['2b1Wj'] = [os.path.join(basepath, 'bb1l_jpa_missingWJet_resolved_even.xml'), 
                                      os.path.join(basepath, 'bb1l_jpa_missingWJet_resolved_odd.xml')]
        resolvedModelDict['2b0Wj'] = [os.path.join(basepath, 'bb1l_jpa_missingAllWJet_even.xml'), 
                                      os.path.join(basepath, 'bb1l_jpa_missingAllWJet_odd.xml')]
        resolvedModelDict['1b2Wj'] = [os.path.join(basepath, 'bb1l_jpa_missingBJet_even.xml'), 
                                      os.path.join(basepath, 'bb1l_jpa_missingBJet_odd.xml')]
        resolvedModelDict['1b1Wj'] = [os.path.join(basepath, 'bb1l_jpa_missingBJet_missingWJet_even.xml'), 
                                      os.path.join(basepath, 'bb1l_jpa_missingBJet_missingWJet_odd.xml')]
        resolvedModelDict['1b0Wj'] = [os.path.join(basepath, 'bb1l_jpa_missingBJet_missingAllWJet_even.xml'), 
                                      os.path.join(basepath, 'bb1l_jpa_missingBJet_missingAllWJet_odd.xml')]
        resolvedModelDict['evCat'] = [os.path.join(basepath, 'bb1l_jpa_evtCat_resolved_even.xml'), 
                                      os.path.join(basepath, 'bb1l_jpa_evtCat_resolved_odd.xml')]
        
        boostedModelDict  = dict()
        boostedModelDict['Hbb2Wj'] = [os.path.join(basepath, 'bb1l_jpa_4jet_boosted_even.xml'), 
                                      os.path.join(basepath, 'bb1l_jpa_4jet_boosted_odd.xml')]
        boostedModelDict['Hbb1Wj'] = [os.path.join(basepath, 'bb1l_jpa_missingWJet_boosted_even.xml'), 
                                      os.path.join(basepath, 'bb1l_jpa_missingWJet_boosted_odd.xml')]
        boostedModelDict['evCat']  = [os.path.join(basepath, 'bb1l_jpa_evtCat_boosted_even.xml'), 
                                      os.path.join(basepath, 'bb1l_jpa_evtCat_boosted_odd.xml')]
        
        # ----------------------------------- NodeList ------------------------------------- #
        # keep the exact same order of nodes as mentioned in respective xml files
        ResolvedJPANodeList = ['2b2Wj','2b1Wj','1b2Wj','2b0Wj','1b1Wj','1b0Wj','0b']
        BoostedJPANodeList  = ['Hbb2Wj','Hbb1Wj','Hbb0Wj']
        
        
        if not self.inclusive_sel:
            #----- Check arguments -----#
            jet_level    = ["Ak4","Ak8","Res2b2Wj","Res2b1Wj","Res2b0Wj","Res1b2Wj","Res1b1Wj","Res1b0Wj","Res0b","Hbb2Wj","Hbb1Wj","Hbb0Wj"]

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
                lep = self.electronsTightSel[0]

            if self.args.Channel == "Mu":
                selObj = MuSelObj
                lep = self.muonsTightSel[0]

            
            #----- Jet selection -----#
            # Since the selections in one line, we can use the non copy option of the selection to modify the selection object internally
            if any([self.args.__dict__[item] for item in ["Ak4","Res2b2Wj","Res2b1Wj","Res2b0Wj","Res1b2Wj","Res1b1Wj","Res1b0Wj","Res0b"]]):
                selObj = makeResolvedSelection(self,selObj,copy_sel=True)
                if self.args.Channel == "El":
                    print('... Resolved :: El Channel')
                    L1out,L2out,selObjAndJetsPerJpaCatDict = findJPACategoryResolved (self, selObj, self.leadElectronTightSel[0], 
                                                                                      self.muonsPreSel, self.electronsPreSel, 
                                                                                      self.ak4Jets, self.ak4BJetsLoose, self.ak4BJets, self.corrMET, 
                                                                                      resolvedModelDict, t.event, self.HLL, ResolvedJPANodeList, 
                                                                                      plot_yield=False)
                    
                    Res2b2Wjets = selObjAndJetsPerJpaCatDict.get('2b2Wj')[1]
                    Res2b1Wjets = selObjAndJetsPerJpaCatDict.get('2b1Wj')[1]
                    Res1b2Wjets = selObjAndJetsPerJpaCatDict.get('1b2Wj')[1]
                    Res2b0Wjets = selObjAndJetsPerJpaCatDict.get('2b0Wj')[1]
                    Res1b1Wjets = selObjAndJetsPerJpaCatDict.get('1b1Wj')[1]
                    Res1b0Wjets = selObjAndJetsPerJpaCatDict.get('1b0Wj')[1]
                    
                if self.args.Channel == "Mu":
                    print('... Resolved :: Mu Channel')
                    L1out,L2out,selObjAndJetsPerJpaCatDict = findJPACategoryResolved (self, selObj, self.leadMuonTightSel[0], 
                                                                                      self.muonsPreSel, self.electronsPreSel, 
                                                                                      self.ak4Jets, self.ak4BJetsLoose, self.ak4BJets, self.corrMET, 
                                                                                      resolvedModelDict, t.event, self.HLL, ResolvedJPANodeList, 
                                                                                      plot_yield=False)
                    Res2b2Wjets = selObjAndJetsPerJpaCatDict.get('2b2Wj')[1]
                    Res2b1Wjets = selObjAndJetsPerJpaCatDict.get('2b1Wj')[1]
                    Res1b2Wjets = selObjAndJetsPerJpaCatDict.get('1b2Wj')[1]
                    Res2b0Wjets = selObjAndJetsPerJpaCatDict.get('2b0Wj')[1]
                    Res1b1Wjets = selObjAndJetsPerJpaCatDict.get('1b1Wj')[1]
                    Res1b0Wjets = selObjAndJetsPerJpaCatDict.get('1b0Wj')[1]
                    
            if any([self.args.__dict__[item] for item in ["Ak8","Hbb2Wj","Hbb1Wj","Hbb0Wj"]]):
                selObj = makeBoostedSelection(self,selObj,copy_sel=True)
                if self.args.Channel == "El":
                    print('... Boosted :: El Channel')
                    L1out,L2out,selObjAndJetsPerJpaCatDict = findJPACategoryBoosted (self, selObj, self.leadElectronTightSel[0], 
                                                                                     self.muonsPreSel, self.electronsPreSel, 
                                                                                     self.ak8BJets, self.ak4JetsCleanedFromAk8b, self.ak4BJetsLoose, 
                                                                                     self.ak4BJets, self.corrMET, boostedModelDict, t.event, self.HLL, 
                                                                                     BoostedJPANodeList, plot_yield=False)
                    Boo2b2Wjets = selObjAndJetsPerJpaCatDict.get('Hbb2Wj')[1]
                    Boo2b1Wjets = selObjAndJetsPerJpaCatDict.get('Hbb1Wj')[1]

                if self.args.Channel == "Mu":
                    print('... Boosted :: Mu Channel')
                    L1out,L2out,selObjAndJetsPerJpaCatDict = findJPACategoryBoosted (self, selObj, self.leadMuonTightSel[0], 
                                                                                     self.muonsPreSel, self.electronsPreSel, 
                                                                                     self.ak8BJets, self.ak4JetsCleanedFromAk8b, self.ak4BJetsLoose, 
                                                                                     self.ak4BJets, self.corrMET, boostedModelDict, t.event, self.HLL, 
                                                                                     BoostedJPANodeList, plot_yield=False)
                    Boo2b2Wjets = selObjAndJetsPerJpaCatDict.get('Hbb2Wj')[1]
                    Boo2b1Wjets = selObjAndJetsPerJpaCatDict.get('Hbb1Wj')[1]

                    
            if self.args.Res2b2Wj:
                print("...... 2b2Wj")
                selObj = selObjAndJetsPerJpaCatDict.get('2b2Wj')[0]
            if self.args.Res2b1Wj:
                print("...... 2b1Wj")
                selObj = selObjAndJetsPerJpaCatDict.get('2b1Wj')[0]
            if self.args.Res2b0Wj:
                print("...... 2b0Wj")
                selObj = selObjAndJetsPerJpaCatDict.get('2b0Wj')[0]
            if self.args.Res1b2Wj:
                print("...... 1b2Wj")
                selObj = selObjAndJetsPerJpaCatDict.get('1b2Wj')[0]
            if self.args.Res1b1Wj:
                print("...... 1b1Wj")
                selObj = selObjAndJetsPerJpaCatDict.get('1b1Wj')[0]
            if self.args.Res1b0Wj:
                print("...... 1b0Wj")
                selObj = selObjAndJetsPerJpaCatDict.get('1b0Wj')[0]
            if self.args.Res0b:
                print("...... 0b")
                selObj = selObjAndJetsPerJpaCatDict.get('0b')[0]
            if self.args.Hbb2Wj:
                print("...... Hbb2Wj")
                selObj = selObjAndJetsPerJpaCatDict.get('Hbb2Wj')[0]
            if self.args.Hbb1Wj:
                print("...... Hbb1Wj")
                selObj = selObjAndJetsPerJpaCatDict.get('Hbb1Wj')[0]
            if self.args.Hbb0Wj:
                print("...... Hbb0Wj")
                selObj = selObjAndJetsPerJpaCatDict.get('Hbb0Wj')[0]

            
        #---------------------------------------------------------------------------------------# 
        #                                 Synchronization tree                                  #
        #---------------------------------------------------------------------------------------#
        if self.args.Synchronization:


        #>>>> TODO : fix few variables 

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
            varsToKeep["has_lead_fakeable_mu"] = op.switch(op.rng_len(self.leadMuonFakeSel) == 1, op.c_bool(True), op.c_bool(False))
            varsToKeep["has_lead_fakeable_el"] = op.switch(op.rng_len(self.leadElectronFakeSel) == 1, op.c_bool(True), op.c_bool(False))
            varsToKeep["pt_leadfakeable_ele"]  = op.switch(op.rng_len(self.leadElectronFakeSel) == 1, 
                                                           op.static_cast("Float_t", self.electron_conept[0]), op.c_float(-1.0))
            varsToKeep["pt_leadfakeable_mu"]   = op.switch(op.rng_len(self.leadMuonFakeSel) == 1, 
                                                           op.static_cast("Float_t", self.muon_conept[0]), op.c_float(-1.0))
                                                          
            varsToKeep["n_presel_ak4Jet"]   = op.static_cast("UInt_t",op.rng_len(self.ak4Jets))
            varsToKeep["n_presel_ak8Jet"]   = op.static_cast("UInt_t",op.rng_len(self.ak8Jets))    
            varsToKeep["n_presel_ak8BJet"]  = op.static_cast("UInt_t",op.rng_len(self.ak8BJets))    
            varsToKeep["n_loose_ak4BJet"]   = op.static_cast("UInt_t",op.rng_len(self.ak4BJetsLoose))    
            varsToKeep["n_medium_ak4BJet"]  = op.static_cast("UInt_t",op.rng_len(self.ak4BJets))
            varsToKeep["n_ak4JetsCleanAk8b"]= op.static_cast("UInt_t",op.rng_len(self.ak4JetsCleanedFromAk8b)) 
            varsToKeep["is_SR"]             = op.static_cast("UInt_t",op.OR(op.rng_len(self.leadElectronTightSel)==1,
                                                                            op.rng_len(self.leadMuonTightSel)==1))
            varsToKeep["is_CR"]             = op.static_cast("UInt_t",op.OR(op.rng_len(self.leadElectronFakeExtrapolationSel)>=1,
                                                                            op.rng_len(self.leadMuonFakeExtrapolationSel)>=1))
            varsToKeep["is_e"]              = op.static_cast("UInt_t", op.rng_len(self.leadElectronTightSel)>=1)
            varsToKeep["is_m"]              = op.static_cast("UInt_t", op.rng_len(self.leadMuonTightSel)>=1)

            varsToKeep["has_leadingTightEl"]      = op.switch(op.rng_len(self.leadElectronTightSel)==1, op.c_bool(True), op.c_bool(False))
            varsToKeep["has_leadingTightMu"]      = op.switch(op.rng_len(self.leadMuonTightSel)==1, op.c_bool(True), op.c_bool(False))

            varsToKeep["is_resolved_UCL"]         = op.switch(op.AND(op.rng_len(self.ak4Jets)>=3,
                                                                     op.rng_len(self.ak4BJets) >= 1,
                                                                     op.rng_len(self.ak8Jets)==0), op.c_bool(True), op.c_bool(False))
            varsToKeep["is_resolved"]             = op.switch(op.AND(op.rng_len(self.ak4Jets)>=3,
                                                                     op.rng_len(self.ak4BJets) >= 1,
                                                                     op.rng_len(self.ak8BJets)==0), op.c_bool(True), op.c_bool(False))
            
            varsToKeep["is_semiboosted"]          = op.switch(op.AND(op.rng_len(self.ak8BJets)>=1,op.rng_len(self.ak4JetsCleanedFromAk8b)>=1), op.c_bool(True), op.c_bool(False))

            varsToKeep["n_tau"] = op.static_cast("UInt_t", op.rng_len(self.tauCleanSel))



            varsToKeep["lead_fakeable_mu_mvaTTH"] = op.switch(op.rng_len(self.leadMuonFakeSel) == 1, op.static_cast("Float_t", self.leadMuonFakeSel[0].mvaTTH), 
                                                              op.static_cast("Float_t", op.c_float(-999.9)))
            varsToKeep["lead_fakeable_el_mvaTTH"] = op.switch(op.rng_len(self.leadElectronFakeSel) == 1, op.static_cast("Float_t", self.leadElectronFakeSel[0].mvaTTH), 
                                                              op.static_cast("Float_t", op.c_float(-999.9)))





            varsToKeep["leadFakeableEl_isPrompt"]  = op.switch(op.rng_len(self.leadElectronFakeSel) >= 1, 
                                                               op.switch(self.lambda_is_matched(self.leadElectronFakeSel[0]), op.c_int(1), op.c_int(0)), op.c_int(-9999))
            varsToKeep["leadFakeableEl_mvaTTH"] = op.switch(op.rng_len(self.leadElectronFakeSel) >= 1, self.leadElectronFakeSel[0].mvaTTH, op.static_cast("Float_t", op.c_float(-9999.)))

            varsToKeep["leadFakeableMu_isPrompt"]  = op.switch(op.rng_len(self.leadMuonFakeSel) >= 1, 
                                                               op.switch(self.lambda_is_matched(self.leadMuonFakeSel[0]), op.c_int(1), op.c_int(0)), op.c_int(-9999))
            varsToKeep["leadFakeableMu_mvaTTH"] = op.switch(op.rng_len(self.leadMuonFakeSel) >= 1, self.leadMuonFakeSel[0].mvaTTH, op.static_cast("Float_t", op.c_float(-9999.)))








            # Triggers #
            #varsToKeep['n_leadfakeableSel_ele']     = op.static_cast("UInt_t",op.rng_len(self.leadElectronsFakeSel))
            #varsToKeep['n_leadfakeableSel_mu']      = op.static_cast("UInt_t",op.rng_len(self.leadMuonsFakeSel))
            varsToKeep["triggers"]                  = self.triggers
            varsToKeep["triggers_SingleElectron"]   = op.OR(*self.triggersPerPrimaryDataset['SingleElectron'])
            varsToKeep["triggers_SingleMuon"]       = op.OR(*self.triggersPerPrimaryDataset['SingleMuon'])
            varsToKeep["triggers_DoubleElectron"]   = op.OR(*self.triggersPerPrimaryDataset['DoubleEGamma'])
            varsToKeep["triggers_DoubleMuon"]       = op.OR(*self.triggersPerPrimaryDataset['DoubleMuon'])
            varsToKeep["triggers_MuonElectron"]     = op.OR(*self.triggersPerPrimaryDataset['MuonEG'])

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
                varsToKeep["mu{}_isfakeablesel".format(i)]         = op.switch(op.rng_len(self.muonsPreSel) >= i, op.switch(self.lambda_muonFakeSel(self.muonsPreSel[i-1]), 
                                                                                                                            op.c_int(1), op.c_int(0)), op.c_int(-9999))
                varsToKeep["mu{}_ismvasel".format(i)]              = op.switch(op.rng_len(self.muonsPreSel) >= i, op.switch(op.AND(self.lambda_muonTightSel(self.muonsPreSel[i-1]), 
                                                                                                                                   self.lambda_muonFakeSel(self.muonsPreSel[i-1])), 
                                                                                                                            op.c_int(1), op.c_int(0)), op.c_int(-9999)) # mvasel encompasses fakeablesel
                varsToKeep["mu{}_isGenMatched".format(i)]          = op.switch(op.rng_len(self.muonsPreSel) >= i, op.switch(self.lambda_is_matched(self.muonsPreSel[i-1]), 
                                                                                                                            op.c_int(1), op.c_int(0)), op.c_int(-9999))
                varsToKeep["mu{}_genPartFlav".format(i)]           = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonsPreSel[i-1].genPartFlav, op.c_int(-9999))
                #varsToKeep["mu{}_FR".format(i)]                    = op.switch(op.rng_len(self.muonsPreSel) >= i, self.muonFR(self.muonsPreSel[i-1]), op.c_int(-9999))
                #varsToKeep["mu{}_FRCorr".format(i)]                = op.switch(op.rng_len(self.muonsPreSel) >= i, self.lambda_FF_mu(self.muonsPreSel[i-1]), op.c_int(-9999))
                
                #varsToKeep["mu{}_looseSF".format(i)]               = op.switch(op.rng_len(self.muonsPreSel) >= i, reduce(mul,self.lambda_MuonLooseSF(self.muonsPreSel[i-1])), op.c_int(-9999))
                #varsToKeep["mu{}_tightSF".format(i)]               = op.switch(op.rng_len(self.muonsPreSel) >= i, reduce(mul,self.lambda_MuonTightSF(self.muonsPreSel[i-1])), op.c_int(-9999))
            
            for i in range(1,6):
                varsToKeep["tightEl{}_pt".format(i)]   = op.switch(op.rng_len(self.electronsTightSel) >= i, self.electronsTightSel[i-1].pt, op.static_cast("Float_t", op.c_float(-9999.)))
                varsToKeep["tightMu{}_pt".format(i)]   = op.switch(op.rng_len(self.muonsTightSel) >= i, self.muonsTightSel[i-1].pt, op.static_cast("Float_t", op.c_float(-9999.)))
                
                varsToKeep["dR_tau1_mu_{}".format(i)]   = op.switch(op.AND(op.rng_len(self.muonsFakeSel) >= i, op.rng_len(self.tauSel) >= 1), 
                                                                    op.deltaR(self.tauSel[0].p4, self.muonsFakeSel[i-1].p4),
                                                                    op.static_cast("Float_t", op.c_float(-9999.)))
                varsToKeep["dR_tau2_mu_{}".format(i)]   = op.switch(op.AND(op.rng_len(self.muonsFakeSel) >= i, op.rng_len(self.tauSel) >= 2), 
                                                                    op.deltaR(self.tauSel[1].p4, self.muonsFakeSel[i-1].p4),
                                                                    op.static_cast("Float_t", op.c_float(-9999.)))
                varsToKeep["dR_tau1_el_{}".format(i)]   = op.switch(op.AND(op.rng_len(self.electronsFakeSel) >= i, op.rng_len(self.tauSel) >= 1), 
                                                                    op.deltaR(self.tauSel[0].p4, self.electronsFakeSel[i-1].p4),
                                                                    op.static_cast("Float_t", op.c_float(-9999.)))
                varsToKeep["dR_tau2_el_{}".format(i)]   = op.switch(op.AND(op.rng_len(self.electronsFakeSel) >= i, op.rng_len(self.tauSel) >= 2), 
                                                                    op.deltaR(self.tauSel[1].p4, self.electronsFakeSel[i-1].p4),
                                                                    op.static_cast("Float_t", op.c_float(-9999.)))

                varsToKeep["fakeableEl_{}_pt".format(i)] = op.switch(op.rng_len(self.electronsFakeSel) >= i, self.electronsFakeSel[i-1].pt, op.static_cast("Float_t", op.c_float(-9999.)))
                varsToKeep["fakeableEl_{}_eta".format(i)] = op.switch(op.rng_len(self.electronsFakeSel) >= i, self.electronsFakeSel[i-1].eta, op.static_cast("Float_t", op.c_float(-9999.)))
                varsToKeep["fakeableEl_{}_phi".format(i)] = op.switch(op.rng_len(self.electronsFakeSel) >= i, self.electronsFakeSel[i-1].phi, op.static_cast("Float_t", op.c_float(-9999.)))
                varsToKeep["fakeableEl_{}_isPrompt".format(i)]  = op.switch(op.rng_len(self.electronsFakeSel) >= i, 
                                                                            op.switch(self.lambda_is_matched(self.electronsFakeSel[i-1]), op.c_int(1), op.c_int(0)), op.c_int(-9999))
                varsToKeep["fakeableEl_{}_mvaTTH".format(i)] = op.switch(op.rng_len(self.electronsFakeSel) >= i, self.electronsFakeSel[i-1].mvaTTH, op.static_cast("Float_t", op.c_float(-9999.)))


                varsToKeep["fakeableMu_{}_pt".format(i)] = op.switch(op.rng_len(self.muonsFakeSel) >= i, self.muonsFakeSel[i-1].pt, op.static_cast("Float_t", op.c_float(-9999.)))
                varsToKeep["fakeableMu_{}_eta".format(i)] = op.switch(op.rng_len(self.muonsFakeSel) >= i, self.muonsFakeSel[i-1].eta, op.static_cast("Float_t", op.c_float(-9999.)))
                varsToKeep["fakeableMu_{}_phi".format(i)] = op.switch(op.rng_len(self.muonsFakeSel) >= i, self.muonsFakeSel[i-1].phi, op.static_cast("Float_t", op.c_float(-9999.)))
                varsToKeep["fakeableMu_{}_isPrompt".format(i)]  = op.switch(op.rng_len(self.muonsFakeSel) >= i, 
                                                                            op.switch(self.lambda_is_matched(self.muonsFakeSel[i-1]), op.c_int(1), op.c_int(0)), op.c_int(-9999))
                varsToKeep["fakeableMu_{}_mvaTTH".format(i)] = op.switch(op.rng_len(self.muonsFakeSel) >= i, self.muonsFakeSel[i-1].mvaTTH, op.static_cast("Float_t", op.c_float(-9999.)))

                
                varsToKeep["tauVetoCand_{}_pt".format(i)] = op.switch(op.rng_len(self.tauCleanSel) >= i, self.tauCleanSel[i-1].pt, op.static_cast("Float_t", op.c_float(-9999.)))
                varsToKeep["tauVetoCand_{}_eta".format(i)] = op.switch(op.rng_len(self.tauCleanSel) >= i, self.tauCleanSel[i-1].eta, op.static_cast("Float_t", op.c_float(-9999.)))
                varsToKeep["tauVetoCand_{}_phi".format(i)] = op.switch(op.rng_len(self.tauCleanSel) >= i, self.tauCleanSel[i-1].phi, op.static_cast("Float_t", op.c_float(-9999.)))
                
                varsToKeep["dR_clnTau1_mu_{}".format(i)]   = op.switch(op.AND(op.rng_len(self.muonsFakeSel) >= i, op.rng_len(self.tauCleanSel) >= 1), 
                                                                       op.deltaR(self.tauCleanSel[0].p4, self.muonsFakeSel[i-1].p4),
                                                                       op.static_cast("Float_t", op.c_float(-9999.)))
                varsToKeep["dR_clnTau2_mu_{}".format(i)]   = op.switch(op.AND(op.rng_len(self.muonsFakeSel) >= i, op.rng_len(self.tauCleanSel) >= 2), 
                                                                       op.deltaR(self.tauCleanSel[1].p4, self.muonsFakeSel[i-1].p4),
                                                                       op.static_cast("Float_t", op.c_float(-9999.)))
                varsToKeep["dR_clnTau1_el_{}".format(i)]   = op.switch(op.AND(op.rng_len(self.electronsFakeSel) >= i, op.rng_len(self.tauCleanSel) >= 1), 
                                                                       op.deltaR(self.tauCleanSel[0].p4, self.electronsFakeSel[i-1].p4),
                                                                       op.static_cast("Float_t", op.c_float(-9999.)))
                varsToKeep["dR_clnTau2_el_{}".format(i)]   = op.switch(op.AND(op.rng_len(self.electronsFakeSel) >= i, op.rng_len(self.tauCleanSel) >= 2), 
                                                                       op.deltaR(self.tauCleanSel[1].p4, self.electronsFakeSel[i-1].p4),
                                                                       op.static_cast("Float_t", op.c_float(-9999.)))
                

            # Electrons #
            for i in range(1,3): # 2 leading electrons 
                varsToKeep["ele{}_pt".format(i)]                    = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].pt, op.c_float(-9999.))
                varsToKeep["ele{}_eta".format(i)]                   = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].eta, op.c_float(-9999.))
                varsToKeep["ele{}_phi".format(i)]                   = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].phi, op.c_float(-9999.))
                varsToKeep["ele{}_E".format(i)]                     = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].p4.E(), op.c_float(-9999., "float"))
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
                varsToKeep["ele{}_isfakeablesel".format(i)]         = op.switch(op.rng_len(self.electronsPreSel) >= i, op.switch(self.lambda_electronFakeSel(self.electronsPreSel[i-1]), 
                                                                                                                                 op.c_int(1), op.c_int(0)), op.c_int(-9999))
                varsToKeep["ele{}_ismvasel".format(i)]              = op.switch(op.rng_len(self.electronsPreSel) >= i, op.switch(op.AND(self.lambda_electronTightSel(self.electronsPreSel[i-1]), 
                                                                                                                                        self.lambda_electronFakeSel(self.electronsPreSel[i-1])), 
                                                                                                                                 op.c_int(1), op.c_int(0)), op.c_int(-9999)) # mvasel encompasses fakeablesel
                varsToKeep["ele{}_isGenMatched".format(i)]          = op.switch(op.rng_len(self.electronsPreSel) >= i, op.switch(self.lambda_is_matched(self.electronsPreSel[i-1]), op.c_int(1), op.c_int(0)), op.c_int(-9999))
                varsToKeep["ele{}_genPartFlav".format(i)]           = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].genPartFlav, op.c_int(-9999))
                varsToKeep["ele{}_deltaEtaSC".format(i)]            = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronsPreSel[i-1].deltaEtaSC, op.c_int(-9999))
                #varsToKeep["ele{}_FR".format(i)]                    = op.switch(op.rng_len(self.electronsPreSel) >= i, self.electronFR(self.electronsPreSel[i-1]), op.c_int(-9999))
                #varsToKeep["ele{}_FF".format(i)]                    = op.switch(op.rng_len(self.electronsPreSel) >= i, self.lambda_FF_el(self.electronsPreSel[i-1]), op.c_int(-9999))
              
                #varsToKeep["ele{}_looseSF".format(i)]               = op.switch(op.rng_len(self.electronsPreSel) >= i, reduce(mul,self.lambda_ElectronLooseSF(self.electronsPreSel[i-1])), op.c_int(-9999))
                #varsToKeep["ele{}_tightSF".format(i)]               = op.switch(op.rng_len(self.electronsPreSel) >= i, reduce(mul,self.lambda_ElectronTightSF(self.electronsPreSel[i-1])), op.c_int(-9999))

            # AK4 Jets #
            for i in range(1,5): # 4 leading jets 
                varsToKeep["ak4Jet{}_pt".format(i)]                 = op.switch(op.rng_len(self.ak4Jets) >= i, self.ak4Jets[i-1].pt, op.c_float(-9999.,"float"))
                varsToKeep["ak4Jet{}_eta".format(i)]                = op.switch(op.rng_len(self.ak4Jets) >= i, self.ak4Jets[i-1].eta, op.c_float(-9999.))
                varsToKeep["ak4Jet{}_phi".format(i)]                = op.switch(op.rng_len(self.ak4Jets) >= i, self.ak4Jets[i-1].phi, op.c_float(-9999.))
                varsToKeep["ak4Jet{}_E".format(i)]                  = op.switch(op.rng_len(self.ak4Jets) >= i, self.ak4Jets[i-1].p4.E(), op.c_float(-9999., "float"))
                varsToKeep["ak4Jet{}_CSV".format(i)]                = op.switch(op.rng_len(self.ak4Jets) >= i, self.ak4Jets[i-1].btagDeepFlavB, op.c_float(-9999.))
                varsToKeep["ak4Jet{}_hadronFlavour".format(i)]      = op.switch(op.rng_len(self.ak4Jets) >= i, self.ak4Jets[i-1].hadronFlavour, op.c_float(-9999.))
                varsToKeep["ak4Jet{}_btagSF".format(i)]             = op.switch(op.rng_len(self.ak4Jets) >= i, self.DeepJetDiscReshapingSF(self.ak4Jets[i-1]), op.c_float(-9999.))

            # AK8 Jets #
            for i in range(1,3): # 2 leading fatjets 
                varsToKeep["ak8Jet{}_pt".format(i)]                 = op.switch(op.rng_len(self.ak8Jets) >= i, self.ak8Jets[i-1].pt, op.c_float(-9999.))
                varsToKeep["ak8Jet{}_eta".format(i)]                = op.switch(op.rng_len(self.ak8Jets) >= i, self.ak8Jets[i-1].eta, op.c_float(-9999.))
                varsToKeep["ak8Jet{}_phi".format(i)]                = op.switch(op.rng_len(self.ak8Jets) >= i, self.ak8Jets[i-1].phi, op.c_float(-9999.))


                varsToKeep["ak8bJet{}_pt".format(i)]                 = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].pt, op.c_float(-9999.))
                varsToKeep["ak8bJet{}_eta".format(i)]                = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].eta, op.c_float(-9999.))
                varsToKeep["ak8bJet{}_phi".format(i)]                = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].phi, op.c_float(-9999.))
                varsToKeep["ak8bJet{}_E".format(i)]                  = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].p4.E(), op.c_float(-9999., "float"))
                varsToKeep["ak8bJet{}_msoftdrop".format(i)]          = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].msoftdrop, op.c_float(-9999.))
                varsToKeep["ak8bJet{}_tau1".format(i)]               = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].tau1, op.c_float(-9999.))
                varsToKeep["ak8bJet{}_tau2".format(i)]               = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].tau2, op.c_float(-9999.))
                varsToKeep["ak8bJet{}_subjet0_pt".format(i)]         = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].subJet1.pt, op.c_float(-9999.))
                varsToKeep["ak8bJet{}_subjet0_eta".format(i)]        = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].subJet1.eta, op.c_float(-9999.))
                varsToKeep["ak8bJet{}_subjet0_phi".format(i)]        = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].subJet1.phi, op.c_float(-9999.))
                varsToKeep["ak8bJet{}_subjet0_CSV".format(i)]        = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].subJet1.btagDeepB, op.c_float(-9999.))
                varsToKeep["ak8bJet{}_subjet1_pt".format(i)]         = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].subJet2.pt, op.c_float(-9999.))
                varsToKeep["ak8bJet{}_subjet1_eta".format(i)]        = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].subJet2.eta, op.c_float(-9999.))
                varsToKeep["ak8bJet{}_subjet1_phi".format(i)]        = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].subJet2.phi, op.c_float(-9999.))
                varsToKeep["ak8bJet{}_subjet1_CSV".format(i)]        = op.switch(op.rng_len(self.ak8BJets) >= i, self.ak8BJets[i-1].subJet2.btagDeepB, op.c_float(-9999.))

            # JPA #
            isResolved = any([self.args.__dict__[item] for item in ["Ak4","Res2b2Wj","Res2b1Wj","Res2b0Wj","Res1b2Wj","Res1b1Wj","Res1b0Wj","Res0b"]]) 
            isBoosted  = any([self.args.__dict__[item] for item in ["Ak8","Hbb2Wj","Hbb1Wj","Hbb0Wj"]]) 
            
            
            varsToKeep["jpa_1stLayer_2b2W_resolved"] = op.static_cast("Float_t",L1out[0]) if isResolved else op.static_cast("Float_t",op.c_float(-1.)) 
            varsToKeep["jpa_1stLayer_2b1W_resolved"] = op.static_cast("Float_t",L1out[1]) if isResolved else op.static_cast("Float_t",op.c_float(-1.)) 
            varsToKeep["jpa_1stLayer_1b2W_resolved"] = op.static_cast("Float_t",L1out[2]) if isResolved else op.static_cast("Float_t",op.c_float(-1.)) 
            varsToKeep["jpa_1stLayer_2b0W_resolved"] = op.static_cast("Float_t",L1out[3]) if isResolved else op.static_cast("Float_t",op.c_float(-1.)) 
            varsToKeep["jpa_1stLayer_1b1W_resolved"] = op.static_cast("Float_t",L1out[4]) if isResolved else op.static_cast("Float_t", op.c_float(-1.)) 
            varsToKeep["jpa_1stLayer_1b0W_resolved"] = op.static_cast("Float_t",L1out[5]) if isResolved else op.static_cast("Float_t",op.c_float(-1.)) 
            varsToKeep["jpa_1stLayer_0b_resolved"]   = op.static_cast("Float_t",op.c_float(-1.0)) if isResolved else op.static_cast("Float_t",op.c_float(-1.))
            
            varsToKeep["jpa_2ndLayer_2b2W_resolved"] = op.static_cast("Float_t",L2out[0]) if isResolved else op.static_cast("Float_t",op.c_float(-1.)) 
            varsToKeep["jpa_2ndLayer_2b1W_resolved"] = op.static_cast("Float_t",L2out[1]) if isResolved else op.static_cast("Float_t",op.c_float(-1.)) 
            varsToKeep["jpa_2ndLayer_1b2W_resolved"] = op.static_cast("Float_t",L2out[2]) if isResolved else op.static_cast("Float_t",op.c_float(-1.)) 
            varsToKeep["jpa_2ndLayer_2b0W_resolved"] = op.static_cast("Float_t",L2out[3]) if isResolved else op.static_cast("Float_t",op.c_float(-1.)) 
            varsToKeep["jpa_2ndLayer_1b1W_resolved"] = op.static_cast("Float_t",L2out[4]) if isResolved else op.static_cast("Float_t",op.c_float(-1.)) 
            varsToKeep["jpa_2ndLayer_1b0W_resolved"] = op.static_cast("Float_t",L2out[5]) if isResolved else op.static_cast("Float_t",op.c_float(-1.)) 
            varsToKeep["jpa_2ndLayer_0b_resolved"]   = op.static_cast("Float_t",L2out[6]) if isResolved else op.static_cast("Float_t",op.c_float(-1.))
            
            varsToKeep["jpa_1stLayer_2b2W_boosted"] = op.static_cast("Float_t",L1out[0]) if isBoosted else op.static_cast("Float_t",op.c_float(-1.)) 
            varsToKeep["jpa_1stLayer_2b1W_boosted"] = op.static_cast("Float_t",L1out[1]) if isBoosted else op.static_cast("Float_t",op.c_float(-1.)) 
            varsToKeep["jpa_1stLayer_2b0W_boosted"] = op.static_cast("Float_t",op.c_float(-1.0)) if isBoosted else op.static_cast("Float_t",op.c_float(-1.))
            
            varsToKeep["jpa_2ndLayer_2b2W_boosted"] = op.static_cast("Float_t",L2out[0]) if isBoosted else op.static_cast("Float_t",op.c_float(-1.)) 
            varsToKeep["jpa_2ndLayer_2b1W_boosted"] = op.static_cast("Float_t",L2out[1]) if isBoosted else op.static_cast("Float_t",op.c_float(-1.)) 
            varsToKeep["jpa_2ndLayer_2b0W_boosted"] = op.static_cast("Float_t",L2out[2]) if isBoosted else op.static_cast("Float_t",op.c_float(-1.)) 
            
            
            # JPA-jet variables
            
            # 2b1Wj
            
            varsToKeep["jpa_1stLayer_2b1W_resolved_bjet1_btagCSV"] = op.static_cast("Float_t", Res2b1Wjets[0].btagCSVV2) if isResolved else op.static_cast("Float_t",op.c_float(-9999.))
            varsToKeep["jpa_1stLayer_2b1W_resolved_bjet2_ptReg"]   = op.static_cast("Float_t", bJetCorrPT(Res2b1Wjets[1])) if isResolved else op.static_cast("Float_t",op.c_float(-9999.))
            varsToKeep["jpa_1stLayer_2b1W_resolved_bjet2_btagCSV"] = op.static_cast("Float_t", Res2b1Wjets[1].btagCSVV2) if isResolved else op.static_cast("Float_t",op.c_float(-9999.))
            varsToKeep["jpa_1stLayer_2b1W_resolved_bjet2_qgDiscr"] = op.static_cast("Float_t", Res2b1Wjets[1].qgl) if isResolved else op.static_cast("Float_t",op.c_float(-9999.))
            varsToKeep["jpa_1stLayer_2b1W_resolved_dEta_bjet2_lep"]= op.static_cast("Float_t", op.abs(lep.eta - Res2b1Wjets[1].eta)) if isResolved else op.static_cast("Float_t",op.c_float(-9999.))
            varsToKeep["jpa_1stLayer_2b1W_resolved_wjet1_ptReg"]   = op.static_cast("Float_t", bJetCorrPT(Res2b1Wjets[2])) if isResolved else op.static_cast("Float_t",op.c_float(-9999.))
            varsToKeep["jpa_1stLayer_2b1W_resolved_wjet1_btagCSV"] = op.static_cast("Float_t", Res2b1Wjets[2].btagCSVV2) if isResolved else op.static_cast("Float_t",op.c_float(-9999.))
            varsToKeep["jpa_1stLayer_2b1W_resolved_wjet1_qgDiscr"] = op.static_cast("Float_t", Res2b1Wjets[2].qgl) if isResolved else op.static_cast("Float_t",op.c_float(-9999.))
            varsToKeep["jpa_1stLayer_2b1W_resolved_dEta_wjet1_lep"]= op.static_cast("Float_t", op.abs(lep.eta - Res2b1Wjets[2].eta)) if isResolved else op.static_cast("Float_t",op.c_float(-9999.))
            varsToKeep["jpa_1stLayer_2b1W_resolved_dR_bjet1bjet2"] = op.static_cast("Float_t", op.deltaR(Res2b1Wjets[0].p4, Res2b1Wjets[1].p4)) if isResolved else op.static_cast("Float_t",op.c_float(-9999.))
            varsToKeep["jpa_1stLayer_2b1W_resolved_nJet"]          = op.static_cast("Float_t", op.rng_len(self.ak4Jets)) if isResolved else op.static_cast("Float_t",op.c_float(-9999.))
            varsToKeep["jpa_1stLayer_2b1W_resolved_nBJetLoose"]    = op.static_cast("Float_t", op.rng_len(self.ak4BJetsLoose)) if isResolved else op.static_cast("Float_t",op.c_float(-9999.))
            varsToKeep["jpa_1stLayer_2b1W_resolved_nBJetMedium"]   = op.static_cast("Float_t", op.rng_len(self.ak4BJets)) if isResolved else op.static_cast("Float_t",op.c_float(-9999.))
            
            # 1b2Wj
            


            # MET #
             
            varsToKeep["PFMET"]    = self.corrMET.pt
            varsToKeep["PFMETphi"] = self.corrMET.phi

            # HME #

            # SF #
            electronMuon_cont = op.combine((self.electronsFakeSel, self.muonsFakeSel))
            varsToKeep["trigger_SF"] = op.multiSwitch(
                    (op.AND(op.rng_len(self.electronsTightSel)==1,op.rng_len(self.muonsTightSel)==0) , self.ttH_singleElectron_trigSF(self.electronsTightSel[0])),
                    (op.AND(op.rng_len(self.electronsTightSel)==0,op.rng_len(self.muonsTightSel)==1) , self.ttH_singleMuon_trigSF(self.muonsTightSel[0])),
                    (op.AND(op.rng_len(self.electronsTightSel)>=2,op.rng_len(self.muonsTightSel)==0) , self.lambda_ttH_doubleElectron_trigSF(self.electronsTightSel)),
                    (op.AND(op.rng_len(self.electronsTightSel)==0,op.rng_len(self.muonsTightSel)>=2) , self.lambda_ttH_doubleMuon_trigSF(self.muonsTightSel)),
                    (op.AND(op.rng_len(self.electronsTightSel)>=1,op.rng_len(self.muonsTightSel)>=1) , self.lambda_ttH_electronMuon_trigSF(electronMuon_cont[0])),
                     op.c_float(1.))

            #varsToKeep["lepton_IDSF"] = op.rng_product(self.electronsTightSel, lambda el : reduce(mul,self.lambda_ElectronLooseSF(el)+self.lambda_ElectronTightSF(el))) * \
            #                            op.rng_product(self.muonsTightSel, lambda mu : reduce(mul,self.lambda_MuonLooseSF(mu)+self.lambda_MuonTightSF(mu))) 

            #varsToKeep["lepton_IDSF_recoToLoose"] = op.rng_product(self.electronsTightSel, lambda el : reduce(mul,self.lambda_ElectronLooseSF(el))) * \
            #                                        op.rng_product(self.muonsTightSel, lambda mu : reduce(mul,self.lambda_MuonLooseSF(mu)))
            #varsToKeep["lepton_IDSF_looseToTight"] = op.rng_product(self.electronsTightSel, lambda el : reduce(mul,self.lambda_ElectronTightSF(el))) * \
            #                                         op.rng_product(self.muonsTightSel, lambda mu : reduce(mul,self.lambda_MuonTightSF(mu)))

            # L1 Prefire #
            if era in ["2016","2017"]:
                varsToKeep["L1prefire"] = self.L1Prefiring
            else:
                varsToKeep["L1prefire"] = op.c_float(-9999.)

            # Fake rate #
            if self.args.FakeExtrapolation:
                varsToKeep["fakeRate"] = op.multiSwitch((op.rng_len(self.ElElDileptonFakeExtrapolationSel)>=1,self.ElElFakeFactor(self.ElElDileptonFakeExtrapolationSel[0])),
                                                        (op.rng_len(self.MuMuDileptonFakeExtrapolationSel)>=1,self.MuMuFakeFactor(self.MuMuDileptonFakeExtrapolationSel[0])),
                                                        (op.rng_len(self.ElMuDileptonFakeExtrapolationSel)>=1,self.ElMuFakeFactor(self.ElMuDileptonFakeExtrapolationSel[0])),
                                                        op.c_float(0.))
            else:
                varsToKeep["fakeRate"] = op.c_float(-9999.)

            # Btagging SF #
            varsToKeep["btag_SF"] = self.btagAk4SF
            if "BtagRatioWeight" in self.__dict__.keys():
                varsToKeep["btag_reweighting"] = self.BtagRatioWeight
                varsToKeep["btag_reweighting_SF"] = self.btagAk4SF * self.BtagRatioWeight

            # ttbar PT reweighting #
            if "group" in sampleCfg and sampleCfg["group"] == 'ttbar':
                varsToKeep["topPt_wgt"] = self.ttbar_weight(self.genTop[0],self.genAntitop[0])

           # Event Weight #
            if self.is_MC:
                #varsToKeep["MC_weight"] = op.sign(t.genWeight)
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

        #----- MET variables -----#
        MET = self.corrMET
        
        varsToKeep['METpt']       = MET.pt
        varsToKeep['METphi']      = MET.phi
        varsToKeep['METpx']       = MET.p4.Px()
        varsToKeep['METpy']       = MET.p4.Py()
        varsToKeep['METpz']       = MET.p4.Pz()
        varsToKeep['METenergy']   = MET.p4.E()
        

        #----- Lepton variables -----#
        if self.args.Channel is None:
            raise RuntimeError("You need to specify --Channel")
        if self.args.Channel   == "El": 
            lepton = self.electronsTightSel[0] 
        elif self.args.Channel == "Mu": 
            lepton = self.muonsTightSel[0] 

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
        #varsToKeep['MET_LD']      = self.HLL.MET_LD(self.corrMET, self.ak4Jets, self.electronsFakeSel) if self.args.Channel == "El" else self.HLL.MET_LD(self.corrMET, self.ak4Jets, self.muonsFakeSel)
        varsToKeep['MET_LD']      = self.HLL.MET_LD_DL(self.corrMET, self.ak4Jets, self.electronsFakeSel, self.muonsFakeSel)
        varsToKeep['lep_conept']  = self.HLL.lambdaConePt(lepton)

        #----- Jet variables -----#
        if any([self.args.__dict__[item] for item in ["Ak4","Res2b2Wj","Res2b1Wj","Res2b0Wj","Res1b2Wj","Res1b1Wj","Res1b0Wj","Res0b"]]):
            ##########################################################################################################
            # Res2b2Wj : jet1   jet2     jet3     jet4
            # Res2b1Wj : jet1   jet2     jet3     jet4=0
            # Res2b0Wj : jet1   jet2     jet3=0   jet4=0
            # Res1b2Wj : jet1   jet2=0   jet3     jet4
            # Res1b1Wj : jet1   jet2=0   jet3     jet4=0
            # Res1b0Wj : jet1   jet2=0   jet3=0   jet4=0
            # Res0b    : all jet = 0 because, this is not a JPA category. Only lepton and MET variables will be saved.
            ##########################################################################################################
            
            if self.args.Res2b2Wj:
                print('||| Resolved 2b2Wj category')
                jet1 = Res2b2Wjets[0]
                jet2 = Res2b2Wjets[1]
                jet3 = Res2b2Wjets[2]
                jet4 = Res2b2Wjets[3]

                lambda_cleanVBFresolved = lambda j : op.NOT(op.OR(*(op.deltaR(Res2b2Wjets[i].p4, j.p4) <= 0.8 for i in range(3))))
                VBFJetsJPA              = op.select(self.VBFJets, lambda_cleanVBFresolved)
                VBFJetPairsJPA          = op.sort(op.combine(VBFJetsJPA, N=2, pred=self.lambda_VBFPair), lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))

            if self.args.Res2b1Wj:
                print('||| Resolved 2b1Wj category')
                jet1 = Res2b1Wjets[0]
                jet2 = Res2b1Wjets[1]
                jet3 = Res2b1Wjets[2]
                # jet4 missing
                
                lambda_cleanVBFresolved = lambda j : op.NOT(op.OR(*(op.deltaR(Res2b1Wjets[i].p4, j.p4) <= 0.8 for i in range(2))))
                VBFJetsJPA              = op.select(self.VBFJets, lambda_cleanVBFresolved)
                VBFJetPairsJPA          = op.sort(op.combine(VBFJetsJPA, N=2, pred=self.lambda_VBFPair), lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))

            if self.args.Res2b0Wj:
                print('||| Resolved 2b0Wj category')
                jet1 = Res2b0Wjets[0]
                jet2 = Res2b1Wjets[1]
                # jet3 missing 
                # jet4 missing 

                lambda_cleanVBFresolved = lambda j : op.NOT(op.OR(*(op.deltaR(Res2b0Wjets[i].p4, j.p4) <= 0.8 for i in range(1))))
                VBFJetsJPA              = op.select(self.VBFJets, lambda_cleanVBFresolved)
                VBFJetPairsJPA          = op.sort(op.combine(VBFJetsJPA, N=2, pred=self.lambda_VBFPair), lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))


            if self.args.Res1b2Wj:
                print('||| Resolved 1b2Wj category')
                jet1 = Res1b2Wjets[0]
                # jet2 missing 
                jet3 = Res1b2Wjets[1]
                jet4 = Res1b2Wjets[2]

                lambda_cleanVBFresolved = lambda j : op.NOT(op.OR(*(op.deltaR(Res1b2Wjets[i].p4, j.p4) <= 0.8 for i in range(2))))
                VBFJetsJPA              = op.select(self.VBFJets, lambda_cleanVBFresolved)
                VBFJetPairsJPA          = op.sort(op.combine(VBFJetsJPA, N=2, pred=self.lambda_VBFPair), lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))


            if self.args.Res1b1Wj:
                print('||| Resolved 1b1Wj category')
                jet1 = Res1b1Wjets[0]
                # jet2 missing 
                jet3 = Res1b1Wjets[1]
                # jet4 missing 

                lambda_cleanVBFresolved = lambda j : op.NOT(op.OR(*(op.deltaR(Res1b1Wjets[i].p4, j.p4) <= 0.8 for i in range(1))))
                VBFJetsJPA              = op.select(self.VBFJets, lambda_cleanVBFresolved)
                VBFJetPairsJPA          = op.sort(op.combine(VBFJetsJPA, N=2, pred=self.lambda_VBFPair), lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))

                
            if self.args.Res1b0Wj:
                print('||| Resolved 1b0Wj category')
                jet1 = Res1b0Wjets[0]
                # jet2 missing 
                # jet3 missing 
                # jet4 missing 

                lambda_cleanVBFresolved = lambda j : op.NOT(op.OR(*(op.deltaR(Res2b1Wjets[i].p4, j.p4) <= 0.8 for i in range(0))))
                VBFJetsJPA              = op.select(self.VBFJets, lambda_cleanVBFresolved)
                VBFJetPairsJPA          = op.sort(op.combine(VBFJetsJPA, N=2, pred=self.lambda_VBFPair), lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))

            if self.args.Res0b:
                print('||| Resolved 0b category')
                # no jpa jets 


            varsToKeep['nAk4Jets']  = op.static_cast("UInt_t",op.rng_len(self.ak4Jets))
            varsToKeep['nAk4BJets'] = op.static_cast("UInt_t",op.rng_len(self.ak4BJets))

            varsToKeep['L1_2b2Wj'] = L1out[0]
            varsToKeep['L1_2b1Wj'] = L1out[1]
            varsToKeep['L1_1b2Wj'] = L1out[2]
            varsToKeep['L1_2b0Wj'] = L1out[3]
            varsToKeep['L1_1b1Wj'] = L1out[4]
            varsToKeep['L1_1b0Wj'] = L1out[5]

            varsToKeep['L2_2b2Wj'] = L2out[0]
            varsToKeep['L2_2b1Wj'] = L2out[1]
            varsToKeep['L2_1b2Wj'] = L2out[2]
            varsToKeep['L2_2b0Wj'] = L2out[3]
            varsToKeep['L2_1b1Wj'] = L2out[4]
            varsToKeep['L2_1b0Wj'] = L2out[5]
            varsToKeep['L2_0b'] = L2out[6]

            # Jet-1 variables :: Available for each category
            varsToKeep['bj1_Px']             = self.HLL.bJetCorrP4(jet1).Px()                if not self.args.Res0b else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['bj1_Py']             = self.HLL.bJetCorrP4(jet1).Py()                if not self.args.Res0b else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['bj1_Pz']             = self.HLL.bJetCorrP4(jet1).Pz()                if not self.args.Res0b else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['bj1_E']              = self.HLL.bJetCorrP4(jet1).E()                 if not self.args.Res0b else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['bj1_pt']             = self.HLL.bJetCorrP4(jet1).Pt()                if not self.args.Res0b else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['bj1_eta']            = self.HLL.bJetCorrP4(jet1).Eta()               if not self.args.Res0b else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['bj1_phi']            = self.HLL.bJetCorrP4(jet1).Phi()               if not self.args.Res0b else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['bj1_bTagDeepFlavB']  = jet1.btagDeepFlavB                            if not self.args.Res0b else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['bj1_bTagCSVV2']      = jet1.btagCSVV2                                if not self.args.Res0b else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['bj1LepDR']           = op.deltaR(jet1.p4, lepton.p4)                 if not self.args.Res0b else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['bj1LepDPhi']         = op.abs(op.deltaPhi(jet1.p4, lepton.p4))       if not self.args.Res0b else op.static_cast("Float_t",op.c_float(0.))  
            varsToKeep['bj1MetDPhi']         = op.abs(self.HLL.SinglepMet_dPhi(jet1, MET))   if not self.args.Res0b else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['minDR_lep_allJets']  = self.HLL.mindr_lep1_jet(lepton, self.ak4Jets) if not self.args.Res0b else op.static_cast("Float_t",op.c_float(0.))
            
            # Jet-2 variables :: Available for 2b2Wj, 2b1Wj. 2b0Wj
            has2b = any([self.args.Res2b2Wj, self.args.Res2b1Wj, self.args.Res2b0Wj])
            varsToKeep['bj2_Px']            = self.HLL.bJetCorrP4(jet2).Px()                 if has2b else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['bj2_Py']            = self.HLL.bJetCorrP4(jet2).Py()                 if has2b else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['bj2_Pz']            = self.HLL.bJetCorrP4(jet2).Pz()                 if has2b else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['bj2_E']             = self.HLL.bJetCorrP4(jet2).E()                  if has2b else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['bj2_pt']            = self.HLL.bJetCorrP4(jet2).Pt()                 if has2b else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['bj2_eta']           = self.HLL.bJetCorrP4(jet2).Eta()                if has2b else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['bj2_phi']           = self.HLL.bJetCorrP4(jet2).Phi()                if has2b else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['bj2_bTagDeepFlavB'] = jet2.btagDeepFlavB                             if has2b else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['bj2_bTagCSVV2']     = jet2.btagCSVV2                                 if has2b else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['bj2LepDR']          = op.deltaR(jet2.p4, lepton.p4)                  if has2b else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['bj2LepDPhi']        = op.abs(op.deltaPhi(jet2.p4, lepton.p4))        if has2b else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['bj2MetDPhi']        = op.abs(self.HLL.SinglepMet_dPhi(jet2, MET))    if has2b else op.static_cast("Float_t",op.c_float(0.))
            # Di-Jet Variables
            varsToKeep['bj1bj2_pt']         = (self.HLL.bJetCorrP4(jet1)+self.HLL.bJetCorrP4(jet2)).Pt()                        if has2b else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['bj1bj2_M']          = op.invariant_mass(self.HLL.bJetCorrP4(jet1),self.HLL.bJetCorrP4(jet2))            if has2b else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['cosThetaS_Hbb']     = self.HLL.comp_cosThetaS(self.HLL.bJetCorrP4(jet1), self.HLL.bJetCorrP4(jet2))     if has2b else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['mT_top_3particle']  = op.min(self.HLL.mT2(self.HLL.bJetCorrP4(jet1),lepton.p4 ,MET.p4), 
                                                     self.HLL.mT2(self.HLL.bJetCorrP4(jet2), lepton.p4, MET.p4))                if has2b else op.static_cast("Float_t",op.c_float(0.))
            
            # Jet-3 variables :: Available for 2b2Wj, 2b1Wj, 1b2Wj, 1b1Wj
            has1stWjet = any([self.args.Res2b2Wj, self.args.Res2b1Wj, self.args.Res1b2Wj, self.args.Res1b1Wj])
            varsToKeep['wj1_Px']            = jet3.p4.Px()                                if has1stWjet else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj1_Py']            = jet3.p4.Py()                                if has1stWjet else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj1_Pz']            = jet3.p4.Pz()                                if has1stWjet else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj1_E']             = jet3.p4.E()                                 if has1stWjet else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj1_pt']            = jet3.pt                                     if has1stWjet else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj1_eta']           = jet3.eta                                    if has1stWjet else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj1_phi']           = jet3.phi                                    if has1stWjet else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj1_bTagDeepFlavB'] = jet3.btagDeepFlavB                          if has1stWjet else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj1_bTagCSVV2']     = jet3.btagCSVV2                              if has1stWjet else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj1LepDR']          = op.deltaR(jet3.p4, lepton.p4)               if has1stWjet else op.static_cast("Float_t",op.c_float(0.)) 
            varsToKeep['wj1LepDPhi']        = op.abs(op.deltaPhi(jet3.p4, lepton.p4))     if has1stWjet else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj1MetDPhi']        = op.abs(self.HLL.SinglepMet_dPhi(jet3, MET)) if has1stWjet else op.static_cast("Float_t",op.c_float(0.))

            # Jet-4 variables :: Available for 2b2Wj, 1b2Wj
            has2ndWjet = any([self.args.Res2b2Wj, self.args.Res1b2Wj])
            varsToKeep['wj2_Px']  = jet4.p4.Px()                                            if has2ndWjet else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj2_Py']  = jet4.p4.Py()                                            if has2ndWjet else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj2_Pz']  = jet4.p4.Pz()                                            if has2ndWjet else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj2_E']   = jet4.p4.E()                                             if has2ndWjet else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj2_pt']  = jet4.pt                                                 if has2ndWjet else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj2_eta'] = jet4.eta                                                if has2ndWjet else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj2_phi'] = jet4.phi                                                if has2ndWjet else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj2_bTagDeepFlavB'] = jet4.btagDeepFlavB                            if has2ndWjet else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj2_bTagCSVV2']     = jet4.btagCSVV2                                if has2ndWjet else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj2LepDR']          = op.deltaR(jet4.p4, lepton.p4)                 if has2ndWjet else op.static_cast("Float_t",op.c_float(0.)) 
            varsToKeep['wj2LepDPhi']  = op.abs(op.deltaPhi(jet4.p4, lepton.p4))             if has2ndWjet else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj2MetDPhi']  = op.abs(self.HLL.SinglepMet_dPhi(jet4, MET))         if has2ndWjet else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj1wj2_pt']   = (jet3.p4 + jet4.p4).Pt()                            if has2ndWjet else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj1wj2_M']    = op.invariant_mass(jet3.p4, jet4.p4)                 if has2ndWjet else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['w1w2_MT']     = self.HLL.MT_W1W2_ljj(lepton,jet3,jet4,MET)          if has2ndWjet else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['HWW_Mass']    = self.HLL.HWW_simple(jet3.p4,jet4.p4,lepton.p4,MET).M() if has2ndWjet else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['HWW_Simple_Mass'] = self.HLL.HWW_met_simple(jet3.p4,jet4.p4,lepton.p4,MET.p4).M() if has2ndWjet else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['HWW_dR'] = self.HLL.dR_Hww(jet3.p4,jet4.p4,lepton.p4,MET)           if has2ndWjet else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['cosThetaS_Wjj_simple'] = self.HLL.comp_cosThetaS(jet3.p4, jet4.p4)  if has2ndWjet else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['cosThetaS_WW_simple_met'] = self.HLL.comp_cosThetaS(self.HLL.Wjj_simple(jet3.p4,jet4.p4), 
                                                                            self.HLL.Wlep_met_simple(lepton.p4, MET.p4))               if has2ndWjet else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['cosThetaS_HH_simple_met'] = self.HLL.comp_cosThetaS(self.HLL.bJetCorrP4(jet1)+self.HLL.bJetCorrP4(jet2), 
                                                                            self.HLL.HWW_met_simple(jet3.p4,jet4.p4,lepton.p4,MET.p4)) if self.args.Res2b2Wj else op.static_cast("Float_t",op.c_float(0.))

            # neu p4 only for 1b2Wj and 2b2Wj
            varsToKeep["neuPx"] = self.HLL.neuP4(jet3.p4+jet4.p4+lepton.p4, MET).Px()  if any([self.args.Res2b2Wj, self.args.Res1b2Wj]) else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep["neuPy"] = self.HLL.neuP4(jet3.p4+jet4.p4+lepton.p4, MET).Py()  if any([self.args.Res2b2Wj, self.args.Res1b2Wj]) else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep["neuPz"] = self.HLL.neuP4(jet3.p4+jet4.p4+lepton.p4, MET).Pz()  if any([self.args.Res2b2Wj, self.args.Res1b2Wj]) else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep["neuE"] = self.HLL.neuP4(jet3.p4+jet4.p4+lepton.p4, MET).E()    if any([self.args.Res2b2Wj, self.args.Res1b2Wj]) else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep["neuPt"] = self.HLL.neuP4(jet3.p4+jet4.p4+lepton.p4, MET).Pt()  if any([self.args.Res2b2Wj, self.args.Res1b2Wj]) else op.static_cast("Float_t",op.c_float(0.))


            # diJet dR & dPhi
            hasbj1wj1 = any([self.args.Res2b2Wj, self.args.Res2b1Wj, self.args.Res1b2Wj, self.args.Res1b1Wj])
            varsToKeep['bj1bj2_DR']   = op.deltaR(jet1.p4,jet2.p4)             if has2b else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['bj1bj2_DPhi'] = op.abs(op.deltaPhi(jet1.p4,jet2.p4))   if has2b else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['bj2wj1_DR']   = op.deltaR(jet2.p4,jet3.p4)             if any([self.args.Res2b2Wj, self.args.Res2b1Wj]) else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['bj2wj1_DPhi'] = op.abs(op.deltaPhi(jet2.p4,jet3.p4))   if any([self.args.Res2b2Wj, self.args.Res2b1Wj]) else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj1wj2_DR']   = op.deltaR(jet3.p4,jet4.p4)             if any([self.args.Res2b2Wj, self.args.Res1b2Wj]) else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj1wj2_DPhi'] = op.abs(op.deltaPhi(jet3.p4,jet4.p4))   if any([self.args.Res2b2Wj, self.args.Res1b2Wj]) else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['bj1wj2_DR']   = op.deltaR(jet1.p4,jet4.p4)             if any([self.args.Res2b2Wj, self.args.Res1b2Wj]) else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['bj1wj2_DPhi'] = op.abs(op.deltaPhi(jet1.p4,jet4.p4))   if any([self.args.Res2b2Wj, self.args.Res1b2Wj]) else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['bj1wj1_DR']   = op.deltaR(jet1.p4,jet3.p4)             if hasbj1wj1 else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['bj1wj1_DPhi'] = op.abs(op.deltaPhi(jet1.p4,jet3.p4))   if hasbj1wj1 else op.static_cast("Float_t",op.c_float(0.))

            varsToKeep['VBFj1pt']        = op.switch(op.rng_len(VBFJetPairsJPA) > 0, 
                                                     op.switch(VBFJetPairsJPA[0][0].pt > VBFJetPairsJPA[0][1].pt,
                                                               VBFJetPairsJPA[0][0].pt, VBFJetPairsJPA[0][1].pt),
                                                     op.static_cast("Float_t",op.c_float(0.))) if not self.args.Res0b else op.c_float(0.)
            varsToKeep['VBFj2pt']        = op.switch(op.rng_len(VBFJetPairsJPA) > 0, 
                                                     op.switch(VBFJetPairsJPA[0][0].pt > VBFJetPairsJPA[0][1].pt,
                                                               VBFJetPairsJPA[0][1].pt, VBFJetPairsJPA[0][0].pt),
                                                     op.static_cast("Float_t",op.c_float(0.))) if not self.args.Res0b else op.c_float(0.)

            varsToKeep['VBFj1eta']       = op.switch(op.rng_len(VBFJetPairsJPA) > 0, 
                                                     op.switch(VBFJetPairsJPA[0][0].pt > VBFJetPairsJPA[0][1].pt,
                                                               VBFJetPairsJPA[0][0].eta, VBFJetPairsJPA[0][1].eta),
                                                     op.static_cast("Float_t",op.c_float(0.))) if not self.args.Res0b else op.c_float(0.)
            varsToKeep['VBFj2eta']       = op.switch(op.rng_len(VBFJetPairsJPA) > 0, 
                                                     op.switch(VBFJetPairsJPA[0][0].pt > VBFJetPairsJPA[0][1].pt,
                                                               VBFJetPairsJPA[0][1].eta, VBFJetPairsJPA[0][0].eta),
                                                     op.static_cast("Float_t",op.c_float(0.))) if not self.args.Res0b else op.c_float(0.)

            varsToKeep['VBFj1j2dEta']    = op.switch(op.rng_len(VBFJetPairsJPA) > 0,op.abs(VBFJetPairsJPA[0][0].eta - VBFJetPairsJPA[0][1].eta),
                                                     op.static_cast("Float_t",op.c_float(0.))) if not self.args.Res0b else op.c_float(0.)

            varsToKeep['VBFj1j2dPhi']    = op.switch(op.rng_len(VBFJetPairsJPA) > 0,op.abs(op.deltaPhi(VBFJetPairsJPA[0][0].p4,VBFJetPairsJPA[0][1].p4)),
                                                     op.static_cast("Float_t",op.c_float(0.))) if not self.args.Res0b else op.c_float(0.)

            varsToKeep['VBFj1j2invM']    = op.switch(op.rng_len(VBFJetPairsJPA) > 0,op.invariant_mass(VBFJetPairsJPA[0][0].p4,  VBFJetPairsJPA[0][1].p4),
                                                     op.static_cast("Float_t",op.c_float(0.))) if not self.args.Res0b else op.c_float(0.)
            
            varsToKeep['VBF_tag']         = op.c_int(op.rng_len(VBFJetPairsJPA) > 0) if not self.args.Res0b else op.c_float(0.) 

            if self.args.Res2b2Wj : 
                varsToKeep['minJetDR']       = self.HLL.MinDiJetDRTight(jet1,jet2,jet3,jet4)
                varsToKeep['minLepJetDR']    = self.HLL.MinDR_lep4j(lepton,jet1,jet2,jet3,jet4)
                varsToKeep['HT2_lepJetMet']  = self.HLL.HT2_l4jmet(lepton,jet1,jet2,jet3,jet4,MET)
                varsToKeep['HT2R_lepJetMet'] = self.HLL.HT2R_l4jmet(lepton,jet1,jet2,jet3,jet4,MET)
                varsToKeep['JPAcat']         = op.c_int(11)

            if self.args.Res2b1Wj : 
                varsToKeep['minJetDR']       = self.HLL.MinDiJetDRLoose(jet1,jet2,jet3)
                varsToKeep['minLepJetDR']    = self.HLL.MinDR_lep3j(lepton,jet1,jet2,jet3)
                varsToKeep['HT2_lepJetMet']  = self.HLL.HT2_l3jmet(lepton,jet1,jet2,jet3,MET)
                varsToKeep['HT2R_lepJetMet'] = self.HLL.HT2R_l3jmet(lepton,jet1,jet2,jet3,MET)
                varsToKeep['JPAcat']         = op.c_int(12)

            if self.args.Res2b0Wj : 
                varsToKeep['minJetDR']       = op.deltaR(jet1.p4,jet2.p4)
                varsToKeep['minLepJetDR']    = self.HLL.MinDR_lep2j(lepton,jet1,jet2)
                varsToKeep['HT2_lepJetMet']  = self.HLL.HT2_l2jmet(lepton,jet1,jet2,MET)
                varsToKeep['HT2R_lepJetMet'] = self.HLL.HT2R_l2jmet(lepton,jet1,jet2,MET)
                varsToKeep['JPAcat']         = op.c_int(13)

            if self.args.Res1b2Wj : 
                varsToKeep['minJetDR']       = self.HLL.MinDiJetDRLoose(jet1,jet3,jet4)
                varsToKeep['minLepJetDR']    = self.HLL.MinDR_lep3j(lepton,jet1,jet3,jet4)
                varsToKeep['HT2_lepJetMet']  = self.HLL.HT2_l3jmet(lepton,jet1,jet3,jet4,MET)
                varsToKeep['HT2R_lepJetMet'] = self.HLL.HT2R_l3jmet(lepton,jet1,jet3,jet4,MET)
                varsToKeep['JPAcat']         = op.c_int(14)

            if self.args.Res1b1Wj : 
                varsToKeep['minJetDR']       = op.deltaR(jet1.p4,jet3.p4)
                varsToKeep['minLepJetDR']    = self.HLL.MinDR_lep2j(lepton,jet1,jet3)
                varsToKeep['HT2_lepJetMet']  = self.HLL.HT2_l2jmet(lepton,jet1,jet3,MET)
                varsToKeep['HT2R_lepJetMet'] = self.HLL.HT2R_l2jmet(lepton,jet1,jet3,MET)
                varsToKeep['JPAcat']         = op.c_int(15)

            if self.args.Res1b0Wj : 
                varsToKeep['minJetDR']       = op.static_cast("Float_t",op.c_float(0.))
                varsToKeep['minLepJetDR']    = op.deltaR(lepton.p4, jet1.p4)
                varsToKeep['HT2_lepJetMet']  = self.HLL.HT2_l1jmet(lepton,jet1,MET)
                varsToKeep['HT2R_lepJetMet'] = self.HLL.HT2R_l1jmet(lepton,jet1,MET)
                varsToKeep['JPAcat']         = op.c_int(16)

            if self.args.Res0b : 
                varsToKeep['minJetDR']       = op.static_cast("Float_t",op.c_float(0.))
                varsToKeep['minLepJetDR']    = op.static_cast("Float_t",op.c_float(0.))
                varsToKeep['HT2_lepJetMet']  = op.static_cast("Float_t",op.c_float(0.))
                varsToKeep['HT2R_lepJetMet'] = op.static_cast("Float_t",op.c_float(0.))
                varsToKeep['JPAcat']         = op.c_int(17)

            

        #----- Fatjet variables -----#
        if any([self.args.__dict__[item] for item in ["Ak8","Hbb2Wj","Hbb1Wj","Hbb0Wj"]]):
            #######################################
            # Hbb2Wj : jet1 jet2   jet3   jet4
            # Hbb1Wj : jet1 jet2   jet3   jet4=0
            # Hbb0Wj : jet1 jet2   jet3=0 jet4=0
            #######################################
            fatjet  = self.ak8BJets[0]
            if self.args.Hbb2Wj:
                jet1    = fatjet.subJet1
                jet2    = fatjet.subJet2
                jet3    = Boo2b2Wjets[0]
                jet4    = Boo2b2Wjets[1]

                lambda_cleanVBFresolved = lambda j : op.NOT(op.OR(op.OR(*(op.deltaR(Res2b0Wjets[i].p4, j.p4) <= 0.8 for i in range(1)))),
                                                            op.deltaR(fatjet.p4, j.p4) <= 1.2)
                VBFJetsJPA              = op.select(self.VBFJets, lambda_cleanVBFresolved)
                VBFJetPairsJPA          = op.sort(op.combine(VBFJetsJPA, N=2, pred=self.lambda_VBFPair), lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))



            if self.args.Hbb1Wj:
                jet1    = fatjet.subJet1
                jet2    = fatjet.subJet2
                jet3    = Boo2b1Wjets[0]

                lambda_cleanVBFresolved = lambda j : op.NOT(op.OR(op.OR(*(op.deltaR(Res2b0Wjets[i].p4, j.p4) <= 0.8 for i in range(0)))),
                                                            op.deltaR(fatjet.p4, j.p4) <= 1.2)
                VBFJetsJPA              = op.select(self.VBFJets, lambda_cleanVBFresolved)
                VBFJetPairsJPA          = op.sort(op.combine(VBFJetsJPA, N=2, pred=self.lambda_VBFPair), lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))


            if self.args.Hbb0Wj:
                jet1    = fatjet.subJet1
                jet2    = fatjet.subJet2

                lambda_cleanVBFresolved = lambda j : op.deltaR(fatjet.p4, j.p4) > 1.2
                VBFJetsJPA              = op.select(self.VBFJets, lambda_cleanVBFresolved)
                VBFJetPairsJPA          = op.sort(op.combine(VBFJetsJPA, N=2, pred=self.lambda_VBFPair), lambda dijet : -op.invariant_mass(dijet[0].p4,dijet[1].p4))

                

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

            
            varsToKeep['cosThetaS_Hbb']     = self.HLL.comp_cosThetaS(jet1.p4, jet2.p4)
            varsToKeep['mT_top_3particle']  = op.min(self.HLL.mT2(jet1.p4,lepton.p4,MET.p4), self.HLL.mT2(jet2.p4,lepton.p4,MET.p4))


            varsToKeep['wj1_Px']      = jet3.p4.Px()                 if not self.args.Hbb0Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj1_Py']      = jet3.p4.Py()                 if not self.args.Hbb0Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj1_Pz']      = jet3.p4.Pz()                 if not self.args.Hbb0Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj1_E']       = jet3.p4.E()                  if not self.args.Hbb0Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj1_pt']      = jet3.pt                      if not self.args.Hbb0Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj1_eta']     = jet3.eta                     if not self.args.Hbb0Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj1_bTagDeepFlavB'] = jet3.btagDeepFlavB     if not self.args.Hbb0Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj1_bTagCSVV2']     = jet3.btagCSVV2         if not self.args.Hbb0Wj else op.static_cast("Float_t",op.c_float(0.))
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

            varsToKeep['fatbj_Wj1DR']       = op.deltaR(fatjet.p4, jet3.p4)                    if not self.args.Hbb0Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['fatbj_Wj1DPhi']     = op.abs(op.deltaPhi(fatjet.p4, jet3.p4))          if not self.args.Hbb0Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['fatbjSub1_Wj1DR']   = op.switch(jet1.btagDeepB > jet2.btagDeepB, 
                                                        op.deltaR(jet1.p4, jet3.p4), 
                                                        op.deltaR(jet2.p4, jet3.p4))           if not self.args.Hbb0Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['fatbjSub1_Wj1DPhi'] = op.switch(jet1.btagDeepB > jet2.btagDeepB, 
                                                        op.abs(op.deltaPhi(jet1.p4, jet3.p4)), 
                                                        op.abs(op.deltaPhi(jet2.p4, jet3.p4))) if not self.args.Hbb0Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['fatbjSub2_Wj1DR']   = op.switch(jet1.btagDeepB > jet2.btagDeepB,
                                                        op.deltaR(jet2.p4, jet3.p4),
                                                        op.deltaR(jet1.p4, jet3.p4))           if not self.args.Hbb0Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['fatbjSub2_Wj1DPhi'] = op.switch(jet1.btagDeepB > jet2.btagDeepB,
                                                        op.abs(op.deltaPhi(jet2.p4, jet3.p4)),
                                                        op.abs(op.deltaPhi(jet1.p4, jet3.p4)))  if not self.args.Hbb0Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj1_lepDR']         = op.deltaR(jet3.p4, lepton.p4)                     if not self.args.Hbb0Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj1_lepDPhi']       = op.abs(op.deltaPhi(jet3.p4, lepton.p4))           if not self.args.Hbb0Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj2_lepDR']         = op.deltaR(jet4.p4, lepton.p4)                     if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj2_lepDPhi']       = op.abs(op.deltaPhi(jet4.p4, lepton.p4))           if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))


            varsToKeep['fatbj_wj2DR']       = op.deltaR(fatjet.p4, jet4.p4)             if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.)) 
            varsToKeep['fatbj_wj2DPhi']     = op.abs(op.deltaPhi(fatjet.p4, jet4.p4))   if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['fatbjSub1_wj2DR']   = op.deltaR(jet1.p4, jet4.p4)               if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['fatbjSub1_wj2DPhi'] = op.abs(op.deltaPhi(jet1.p4, jet4.p4))     if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['fatbjSub2_wj2DR']   = op.deltaR(jet2.p4, jet4.p4)               if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['fatbjSub2_wj2DPhi'] = op.abs(op.deltaPhi(jet2.p4, jet4.p4))     if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))

            varsToKeep['wj1wj2_pt']   = (jet3.p4+jet4.p4).Pt()                    if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj1wj2DR']    = op.deltaR(jet3.p4, jet4.p4)               if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj1wj2DPhi']  = op.abs(op.deltaPhi(jet3.p4, jet4.p4))     if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['wj1wj2invM']  = op.invariant_mass(jet3.p4,jet4.p4)        if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))

            varsToKeep['HWW_Mass']        = self.HLL.HWW_simple(jet3.p4,jet4.p4,lepton.p4,MET).M()        if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['HWW_Simple_Mass'] = self.HLL.HWW_met_simple(jet3.p4,jet4.p4,lepton.p4,MET.p4).M() if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['HWW_dR']          = self.HLL.dR_Hww(jet3.p4,jet4.p4,lepton.p4,MET)                if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['cosThetaS_Wjj_simple']    = self.HLL.comp_cosThetaS(jet3.p4, jet4.p4)             if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['cosThetaS_WW_simple_met'] = self.HLL.comp_cosThetaS(self.HLL.Wjj_simple(jet3.p4,jet4.p4), 
                                                                            self.HLL.Wlep_met_simple(lepton.p4, MET.p4))               if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))
            varsToKeep['cosThetaS_HH_simple_met'] = self.HLL.comp_cosThetaS(jet1.p4 + jet2.p4, 
                                                                            self.HLL.HWW_met_simple(jet3.p4,jet4.p4,lepton.p4,MET.p4)) if self.args.Hbb2Wj else op.static_cast("Float_t",op.c_float(0.))
                

            varsToKeep['VBFj1pt']        = op.switch(op.rng_len(VBFJetPairsJPA) > 0, 
                                                     op.switch(VBFJetPairsJPA[0][0].pt > VBFJetPairsJPA[0][1].pt,
                                                               VBFJetPairsJPA[0][0].pt, VBFJetPairsJPA[0][1].pt),
                                                     op.static_cast("Float_t",op.c_float(0.)))
            varsToKeep['VBFj2pt']        = op.switch(op.rng_len(VBFJetPairsJPA) > 0, 
                                                     op.switch(VBFJetPairsJPA[0][0].pt > VBFJetPairsJPA[0][1].pt,
                                                               VBFJetPairsJPA[0][1].pt, VBFJetPairsJPA[0][0].pt),
                                                     op.static_cast("Float_t",op.c_float(0.)))

            varsToKeep['VBFj1eta']       = op.switch(op.rng_len(VBFJetPairsJPA) > 0, 
                                                     op.switch(VBFJetPairsJPA[0][0].pt > VBFJetPairsJPA[0][1].pt,
                                                               VBFJetPairsJPA[0][0].eta, VBFJetPairsJPA[0][1].eta),
                                                     op.static_cast("Float_t",op.c_float(0.)))
            varsToKeep['VBFj2eta']       = op.switch(op.rng_len(VBFJetPairsJPA) > 0, 
                                                     op.switch(VBFJetPairsJPA[0][0].pt > VBFJetPairsJPA[0][1].pt,
                                                               VBFJetPairsJPA[0][1].eta, VBFJetPairsJPA[0][0].eta),
                                                     op.static_cast("Float_t",op.c_float(0.)))

            varsToKeep['VBFj1j2dEta']    = op.switch(op.rng_len(VBFJetPairsJPA) > 0,op.abs(VBFJetPairsJPA[0][0].eta - VBFJetPairsJPA[0][1].eta),
                                                     op.static_cast("Float_t",op.c_float(0.)))

            varsToKeep['VBFj1j2dPhi']    = op.switch(op.rng_len(VBFJetPairsJPA) > 0,op.abs(op.deltaPhi(VBFJetPairsJPA[0][0].p4,VBFJetPairsJPA[0][1].p4)),
                                                     op.static_cast("Float_t",op.c_float(0.)))

            varsToKeep['VBFj1j2invM']    = op.switch(op.rng_len(VBFJetPairsJPA) > 0,op.invariant_mass(VBFJetPairsJPA[0][0].p4,  VBFJetPairsJPA[0][1].p4),
                                                     op.static_cast("Float_t",op.c_float(0.)))
            
            varsToKeep['VBF_tag']         = op.c_int(op.rng_len(self.VBFJetPairs)>0)


            if self.args.Hbb2Wj:
                varsToKeep['jetMinDR']    = self.HLL.MinDiJetDRTight(jet1,jet2,jet3,jet4)
                varsToKeep['jetLepMinDR'] = self.HLL.MinDR_lep4j(lepton, jet1, jet2, jet3, jet4) 
                varsToKeep['HT2']         = self.HLL.HT2_l4jmet(lepton, jet1, jet2,jet3,jet4,MET)
                varsToKeep['HT2R']        = self.HLL.HT2R_l4jmet(lepton, jet1, jet2, jet3, jet4,MET)
                varsToKeep['MT_W1W2']     = self.HLL.MT_W1W2_ljj(lepton,jet3,jet4,MET)
                varsToKeep['JPAcat']      = op.c_int(21)                

            if self.args.Hbb1Wj:
                varsToKeep['jetMinDR']    = self.HLL.MinDiJetDRLoose(jet1,jet2,jet3)
                varsToKeep['jetLepMinDR'] = self.HLL.MinDR_lep3j(lepton,jet1,jet2,jet3) 
                varsToKeep['HT2']         = self.HLL.HT2_l3jmet(lepton,jet1,jet2,jet3,MET)
                varsToKeep['HT2R']        = self.HLL.HT2R_l3jmet(lepton,jet1,jet2,jet3,MET)
                varsToKeep['MT_W1W2']     = self.HLL.MT_W1W2_lj(lepton,jet3,MET)
                varsToKeep['JPAcat']      = op.c_int(22)                

            if self.args.Hbb0Wj:
                varsToKeep['jetMinDR']    = op.static_cast("Float_t",op.c_float(0.))
                varsToKeep['jetLepMinDR'] = op.static_cast("Float_t",op.c_float(0.)) 
                varsToKeep['HT2']         = op.static_cast("Float_t",op.c_float(0.))
                varsToKeep['HT2R']        = op.static_cast("Float_t",op.c_float(0.))
                varsToKeep['MT_W1W2']     = op.static_cast("Float_t",op.c_float(0.))
                varsToKeep['JPAcat']      = op.c_int(23)                


                #----- Additional variables -----#
        varsToKeep["MC_weight"]    = t.genWeight
        varsToKeep['total_weight'] = selObj.sel.weight

        #return leptonSel.sel, varsToKeep
        return selObj.sel, varsToKeep
