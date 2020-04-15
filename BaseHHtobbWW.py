import os
import sys
from itertools import chain

from bamboo import treefunctions as op
from bamboo.analysismodules import NanoAODModule, NanoAODHistoModule, NanoAODSkimmerModule
from bamboo.analysisutils import makeMultiPrimaryDatasetTriggerSelection
from bamboo.scalefactors import binningVariables_nano

from METScripts import METFilter, METcorrection
from scalefactorsbbWW import ScaleFactorsbbWW

#===============================================================================================#
#                                  BaseHHtobbWW                                                 #
#===============================================================================================#
class BaseNanoHHtobbWW(NanoAODModule):
    """ Base module: HH->bbW(->e/µ nu)W(->e/µ nu) histograms from NanoAOD """
    def __init__(self, args):
            super(BaseNanoHHtobbWW, self).__init__(args)
            # Set plots options #
            self.plotDefaults = {"show-ratio": True,
                                 "y-axis-show-zero" : True,
                                 #"normalized": True,
                                 "y-axis": "Events",
                                 "log-y"  : "both",
                                 "ratio-y-axis-range" : [0.8,1.2],
                                 "ratio-y-axis" : 'Ratio Data/MC',
                                 "sort-by-yields" : False}

    #-------------------------------------------------------------------------------------------#
    #                                       addArgs                                             #
    #-------------------------------------------------------------------------------------------#
    def addArgs(self,parser):
        super(BaseNanoHHtobbWW, self).addArgs(parser)

        #----- Technical arguments -----#
        parser.add_argument("--backend", 
                            type=str, 
                            default="dataframe", 
                            help="Backend to use, 'dataframe' (default) or 'lazy'")

        parser.title("""

Arguments for the HH->bbWW analysis on bamboo framework

----- Argument groups -----
    * Lepton arguments *
        --Preselected
        --Fakeable
        --Tight
        --FakeExtrapolation
        --NoZVeto

    * Jet arguments *
        --Ak4
        --Ak8
        --Resolved0Btag
        --Resolved1Btag
        --Resolved2Btag
        --Boosted

    * Skimmer arguments *
        --Synchronization
        --Channel

    * Plotter arguments *
        --OnlyYield

----- Plotter mode -----
Every combination of lepton and jet arguments can be used, if none specified they will all be plotted

----- Skimmer mode -----
One lepton and and one jet argument must be specified in addition to the required channel
(this is because only one selection at a time can produce a ntuple)

                    """)

        #----- Lepton selection arguments -----#
        parser.add_argument("--Preselected", 
                            action      = "store_true",
                            default     = False,
                            help        = "Produce the plots/skim for a preselected pair of leptons")
        parser.add_argument("--Fakeable", 
                            action      = "store_true",
                            default     = False,
                            help        = "Produce the plots/skim for a fakeable pair of leptons")
        parser.add_argument("--Tight", 
                            action      = "store_true",
                            default     = False,
                            help        = "Produce the plots/skim for a tight pair of leptons (MC matched)")
        parser.add_argument("--FakeExtrapolation", 
                            action      = "store_true",
                            default     = False,
                            help        = "Produce the plots/skim for a pair of leptons passing the fakeable selection but not the tight (MC matched)")
        parser.add_argument("--NoZVeto", 
                            action      = "store_true",
                            default     = False,
                            help        = "Remove the cut of preselected leptons |M_ll-M_Z|<10 GeV")

        #----- Jet selection arguments -----#
        parser.add_argument("--Ak4", 
                                action      = "store_true",
                                default     = False,
                                help        = "Produce the plots/skim for two Ak4 jets passing the selection criteria")
        parser.add_argument("--Ak8", 
                                action      = "store_true",
                                default     = False,
                                help        = "Produce the plots/skim for one Ak8 jet passing the selection criteria")
        parser.add_argument("--Resolved0Btag", 
                                action      = "store_true",
                                default     = False,
                                help        = "Produce the plots/skim for the exclusive resolved category with no btagged jet")
        parser.add_argument("--Resolved1Btag", 
                                action      = "store_true",
                                default     = False,
                                help        = "Produce the plots/skim for the exclusive resolved category with only one btagged jet")
        parser.add_argument("--Resolved2Btag", 
                                action      = "store_true",
                                default     = False,
                                help        = "Produce the plots/skim for the exclusive resolved category with two btagged jets")
        parser.add_argument("--Boosted", 
                                action      = "store_true",
                                default     = False,
                                help        = "Produce the plots/skim for the inclusive boosted category")

        #----- Skimmer Arguments -----#
        parser.add_argument("--Synchronization", 
                            action      = "store_true",
                            default     = False,
                            help        = "Produce the skims for the synchronization (without triggers, corrections of flags)")
        parser.add_argument("--Channel", 
                            action      = "store",
                            type        = str,
                            help        = "Specify the channel for the tuple : ElEl, ElMu, MuEl")

        #----- Plotter arguments -----#
        parser.add_argument("--OnlyYield", 
                            action      = "store_true",
                            default     = False,
                            help        = "Only produce the yield plots")

    def prepareTree(self, tree, sample=None, sampleCfg=None):
        from bamboo.treedecorators import NanoAODDescription, nanoRochesterCalc, nanoJetMETCalc
        # JEC's Recommendation for Full RunII: https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC
        # JER : -----------------------------: https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
        # Get base aguments #
        era = sampleCfg['era']
        self.is_MC = self.isMC(sample) # no confusion between boolean is_MC and method isMC()
        metName = "METFixEE2017" if era == "2017" else "MET"
        tree,noSel,be,lumiArgs = super(BaseNanoHHtobbWW,self).prepareTree(tree          = tree, 
                                                                          sample        = sample, 
                                                                          sampleCfg     = sampleCfg, 
                                                                          description   = NanoAODDescription.get("v5", year=(era if era else "2016"),
                                                                          isMC          = self.is_MC,
                                                                          systVariations= [ nanoRochesterCalc, (nanoJetMETCalc_METFixEE2017 if era == "2017" else nanoJetMETCalc) ]), 
                                                                                # will do Jet and MET variations, and the Rochester correction
                                                                          lazyBackend   = (self.args.backend == "lazy"))
    
        self.triggersPerPrimaryDataset = {}
        from bamboo.analysisutils import configureJets ,configureRochesterCorrection, configureType1MET 

        ## Check distributed option #
        isNotWorker = (self.args.distributed != "worker") 

        # Check era #
        if era != "2016" and era != "2017" and era != "2018":
            raise RuntimeError("Unknown era {0}".format(era))

        # Rochester and JEC corrections (depends on era) #     
        cachJEC_dir = '/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/cacheJEC'
        ############################################################################################
        # ERA 2016 #
        ############################################################################################
        if era == "2016" and not self.args.Synchronization:
            # Rochester corrections #
            configureRochesterCorrection(variProxy  = tree._Muon,
                                         paramsFile = os.path.join(os.path.dirname(__file__), "data", "RoccoR2016.txt"),
                                         isMC       = self.is_MC,
                                         backend    = be, 
                                         uName      = sample)

            # Trigger efficiencies #
            # tree.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_(DZ) 
            # tree.HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_(DZ) 
            # -> can be present with, without and without both the DZ
            if self.is_MC or "2016F" in sample or "2016G" in sample:# or "2016H" in sample:
                # Found in 2016F : both
                # Found in 2016G : both
                self.triggersPerPrimaryDataset = {
                    "SingleMuon" :  [ tree.HLT.IsoMu24],
                    "SingleElectron":  [ tree.HLT.Ele27_WPTight_Gsf],
                    "DoubleMuon" :  [ tree.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL,
                                      tree.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ],
                    "DoubleEGamma": [ tree.HLT.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ],
                    "MuonEG":       [ tree.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,
                                      tree.HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ,
                                      tree.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL,
                                      tree.HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL]}
            elif "2016H" in sample:
                # Found in 2016H : has DZ but not without
                self.triggersPerPrimaryDataset = {
                    "SingleMuon" :  [ tree.HLT.IsoMu24],
                    "SingleElectron":  [ tree.HLT.Ele27_WPTight_Gsf],
                    "DoubleMuon" :  [ tree.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL,
                                      tree.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ],
                    "DoubleEGamma": [ tree.HLT.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ],
                    "MuonEG":       [ tree.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,
                                      tree.HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ]}

            elif "2016B" in sample or "2016C" in sample or "2016D" in sample or "2016E" in sample : 
                # Found in 2016B : only without DZ
                # Found in 2016C : only without DZ
                # Found in 2016D : only without DZ
                # Found in 2016E : only without DZ

                self.triggersPerPrimaryDataset = {
                    "SingleMuon" :  [ tree.HLT.IsoMu24],
                    "SingleElectron":  [ tree.HLT.Ele27_WPTight_Gsf],
                    "DoubleMuon" :  [ tree.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL,
                                      tree.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ],
                    "DoubleEGamma": [ tree.HLT.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ],
                    "MuonEG":       [ tree.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL, 
                                      tree.HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL]} 

            # JetMET treatment #
            if self.is_MC:   # if MC -> needs smearing
                configureJets(variProxy             = tree._Jet, 
                              jetType               = "AK4PFchs",
                              jec                   = "Summer16_07Aug2017_V20_MC",
                              smear                 = "Summer16_25nsV1_MC",
                              jesUncertaintySources = ["Total"],
                              mayWriteCache         = isNotWorker,
                              isMC                  = self.is_MC,
                              backend               = be, 
                              uName                 = sample,
                              cachedir              = cachJEC_dir)
                configureType1MET(variProxy             = tree._MET, 
                                  jec                   = "Summer16_07Aug2017_V20_MC",
                                  smear                 = "Summer16_25nsV1_MC",
                                  jesUncertaintySources = ["Total"],
                                  mayWriteCache         = isNotWorker,
                                  isMC                  = self.is_MC,
                                  backend               = be,
                                  uName                 = sample,
                                  cachedir              = cachJEC_dir)

            else:                   # If data -> extract info from config 
                jecTag = None
                if "2016B" in sample or "2016C" in sample or "2016D" in sample:
                    jecTag = "Summer16_07Aug2017BCD_V11_DATA"
                elif "2016E" in sample or "2016F" in sample:
                    jecTag = "Summer16_07Aug2017EF_V11_DATA"
                elif "2016G" in sample or "2016H" in sample:
                    jecTag = "Summer16_07Aug2017GH_V11_DATA"
                configureJets(variProxy             = tree._Jet, 
                              jetType               = "AK4PFchs",
                              jec                   = jecTag,
                              mayWriteCache         = isNotWorker,
                              isMC                  = self.is_MC,
                              backend               = be, 
                              uName                 = sample,
                              cachedir              = cachJEC_dir)
                configureType1MET(variProxy         = tree._MET, 
                                  jec               = jecTag,
                                  mayWriteCache     = isNotWorker,
                                  isMC              = self.is_MC,
                                  backend           = be, 
                                  uName             = sample,
                                  cachedir          = cachJEC_dir)

        ############################################################################################
        # ERA 2017 #
        ############################################################################################
        #elif era == "2017" and not self.args.Synchronization:
        #    configureRochesterCorrection(tree._Muon.calc,os.path.join(os.path.dirname(__file__), "data", "RoccoR2017.txt"))
        #    self.triggersPerPrimaryDataset = {
        #        "DoubleMuon" : [ tree.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL,
        #                         tree.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ,
        #                         tree.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8,
        #                         tree.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 ],
        #        "DoubleEG"   : [ tree.HLT.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL ],
        #        # it's recommended to not use the DZ  version  for 2017 and 2018, it would be a needless efficiency loss
        #        #---> https://twiki.cern.ch/twiki/bin/view/CMS/EgHLTRunIISummary
        #        "MuonEG"     : [ tree.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL,
        #                         tree.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,
        #                         tree.HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL,
        #                         tree.HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ ]
        #        }
        #    if self.is_MC:
        #        configureJets(tree, "Jet", "AK4PFchs",
        #            jec="Fall17_17Nov2017_V32_MC",
        #            smear="Fall17_V3_MC",
        #            jesUncertaintySources=["Total"], mayWriteCache=isNotWorker)
        #    else:
        #        if "2017B" in sample:
        #            configureJets(tree, "Jet", "AK4PFchs", jec="Fall17_17Nov2017B_V32_DATA", mayWriteCache=isNotWorker)
        #        elif "2017C" in sample:
        #            configureJets(tree, "Jet", "AK4PFchs", jec="Fall17_17Nov2017C_V32_DATA", mayWriteCache=isNotWorker)
        #        elif "2017D" in sample or "2017E" in sample:
        #            configureJets(tree, "Jet", "AK4PFchs", jec="Fall17_17Nov2017DE_V32_DATA", mayWriteCache=isNotWorker)
        #        elif "2017F" in sample:
        #            configureJets(tree, "Jet", "AK4PFchs", jec="Fall17_17Nov2017F_V32_DATA", mayWriteCache=isNotWorker)

        ############################################################################################
        # ERA 2018 #
        ############################################################################################
        #elif era == "2018" and not self.args.Synchronization:
        #    configureRochesterCorrection(tree._Muon.calc,os.path.join(os.path.dirname(__file__), "data", "RoccoR2018.txt"))
        #    self.triggersPerPrimaryDataset = {
        #        "DoubleMuon" : [ tree.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL,
        #                         tree.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ,
        #                         tree.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8,
        #                         tree.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 ],
        #        "EGamma"     : [ tree.HLT.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL ], 
        #        "MuonEG"     : [ tree.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL,
        #                         tree.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,
        #                         tree.HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL,
        #                         tree.HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ ]
        #        }
        #    if self.is_MC:
        #        configureJets(tree, "Jet", "AK4PFchs",
        #            jec="Autumn18_V8_MC",
        #            smear="Autumn18_V1_MC",
        #            jesUncertaintySources=["Total"], mayWriteCache=isNotWorker)
        #    else:
        #        if "2018A" in sample:
        #            configureJets(tree, "Jet", "AK4PFchs",
        #                jec="Autumn18_RunA_V8_DATA", mayWriteCache=isNotWorker)
        #        elif "2018B" in sample:
        #            configureJets(tree, "Jet", "AK4PFchs",
        #                jec="Autumn18_RunB_V8_DATA", mayWriteCache=isNotWorker)
        #        elif "2018C" in sample:
        #            configureJets(tree, "Jet", "AK4PFchs",
        #                jec="Autumn18_RunC_V8_DATA", mayWriteCache=isNotWorker)
        #        elif "2018D" in sample:
        #            configureJets(tree, "Jet", "AK4PFchs",
        #                jec="Autumn18_RunD_V8_DATA", mayWriteCache=isNotWorker)

        # Triggers and Gen Weight #
        if not self.args.Synchronization:
            if self.is_MC:
                noSel = noSel.refine("genWeight", weight=tree.genWeight, cut=op.OR(*chain.from_iterable(self.triggersPerPrimaryDataset.values())))
            else:
                noSel = noSel.refine("withTrig", cut=makeMultiPrimaryDatasetTriggerSelection(sample, self.triggersPerPrimaryDataset))

        return tree,noSel,be,lumiArgs

    def prepareObjects(self, t, noSel, sample, sampleCfg):
        # Some imports #
        from bamboo.analysisutils import forceDefine

        era = sampleCfg['era']

        # Forcedefine : calculate once per event (for every event) #
        if not self.args.Synchronization:
            forceDefine(t._Muon.calcProd, noSel) # Muons for Rochester corrections
            forceDefine(t._Jet.calcProd, noSel)  # Jets for configureJets
            forceDefine(getattr(t, "_{0}".format("MET" if era != "2017" else "METFixEE2017")).calcProd,noSel) # MET for configureMET

        ###########################################################################
        #                           TTbar reweighting                             #
        ###########################################################################
        #if self.is_MC and sample.startswith("TT"):
        #    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting#Use_case_3_ttbar_MC_is_used_to_m
        #    # Get tops #
        #    genTop_all = op.select(t.GenPart,lambda g : g.pdgId==6)
        #    genTop = op.select(genTop_all,lambda g : g.statusFlags & ( 0x1 << 13))
        #    genAntitop_all = op.select(t.GenPart,lambda g : g.pdgId==-6)
        #    genAntitop = op.select(genAntitop_all,lambda g : g.statusFlags & ( 0x1 << 13))
        #        # statusFlags==13 : isLastCopy
        #        # Pdgid == 6 : top
        #    #hasttbar = noSel.refine("hasttbar",cut=[op.rng_len(genTop)>=1,op.rng_len(genAntitop)>=1])
        #    # Lambda to compute weight if there is a ttbar #
        #    ttbar_SF = lambda t : op.exp(0.0615-0.0005*t.pt)
        #    ttbar_weight = lambda t,tbar : op.sqrt(ttbar_SF(t)*ttbar_SF(tbar))
        #   # Apply correction to TT #
        #    noSel = noSel.refine("ttbarWeight",weight=ttbar_weight(genTop[0],genAntitop[0]))

        #############################################################################
        #                             Pile-up                                       #
        #############################################################################
        puWeightsFile = None
        # Select sfTaf #
        if era == "2016":
            self.sfTag="94X"
        elif era == "2017":
            self.sfTag="94X"     
        elif era == "2018":
            self.sfTag="102X"

        # Get MC PU weight file #
        if self.is_MC and not self.args.Synchronization:
            if "pufile" not in sampleCfg:
                raise KeyError("Could not find 'pufile' entry for sample %s in the YAML file"%sampleCfg["sample"])
            puWeightsFile = os.path.join(os.path.dirname(__file__), "data", "pileup",sampleCfg["pufile"])
            if not os.path.exists(puWeightsFile):
                raise RuntimeError("Could not find pileup file %s"%puWeightsFile)
            from bamboo.analysisutils import makePileupWeight
            PUWeight = makePileupWeight(puWeightsFile, t.Pileup_nTrueInt, systName="pileup",nameHint=f"puweightFromFile{sample}".replace('-','_'))
            noSel = noSel.refine("puWeight", weight=PUWeight)

        #############################################################################
        #                                 MET                                       #
        #############################################################################
        # MET filter #
        if not self.args.Synchronization:
            noSel = noSel.refine("passMETFlags", cut=METFilter(t.Flag, era, self.is_MC) )

        # MET corrections #
        MET = t.MET if era != "2017" else t.METFixEE2017
        #if not self.args.Synchronization:
        #    self.corrMET = METcorrection(MET,t.PV,sample,era,self.is_MC) # Flatness correction might not be needed
        #else:
        #    self.corrMET = MET
        self.corrMET = MET

        #############################################################################
        #                      Lepton Lambdas Variables                             #
        #############################################################################
        # Associated jet Btagging #
        #self.lambda_hasAssociatedJet = lambda lep : op.AND(lep.jet.idx != -1 , op.deltaR(lep.p4,lep.jet.p4) <= 0.4)
        self.lambda_hasAssociatedJet = lambda lep : lep.jet.idx != -1
        if era == "2016": 
            #self.lambda_lepton_associatedJetNoBtag = lambda lep : lep.jet.btagDeepFlavB < 0.3093
            self.lambda_lepton_associatedJetNoBtag = lambda lep : op.OR(op.NOT(self.lambda_hasAssociatedJet(lep)),
                                                                        lep.jet.btagDeepFlavB <= 0.3093)
        elif era =="2017":
            self.lambda_lepton_associatedJetNoBtag = lambda lep : op.OR(op.NOT(self.lambda_hasAssociatedJet(lep)),
                                                                        lep.jet.btagDeepFlavB < 0.3033)
        elif era == "2018":
            self.lambda_lepton_associatedJetNoBtag = lambda lep : op.OR(op.NOT(self.lambda_hasAssociatedJet(lep)),
                                                                        lep.jet.btagDeepFlavB < 0.2770)
                     # If no associated jet, isolated lepton : cool !  

        # Cone pt #
                    # Def conept : https://github.com/CERN-PH-CMG/cmgtools-lite/blob/f8a34c64a4489d94ff9ac4c0d8b0b06dad46e521/TTHAnalysis/python/tools/conept.py#L74
        self.lambda_conept_electron = lambda lep : op.multiSwitch((op.AND(op.abs(lep.pdgId)!=11 , op.abs(lep.pdgId)!=13) , op.static_cast("Float_t",lep.pt)),
                                                                  # if (abs(lep.pdgId)!=11 and abs(lep.pdgId)!=13): return lep.pt : anything that is not muon or electron
                                                                  (op.AND(op.abs(lep.pdgId)==11 , lep.mvaTTH > 0.80) , op.static_cast("Float_t",lep.pt)),
                                                                  # if (abs(lep.pdgId)!=13 or lep.mediumMuonId>0) and lep.mvaTTH > 0.90: return lep.pt 
                                                                    # if electron, check above MVA 
                                                                    # If muon, check that passes medium and above MVA
                                                                   op.static_cast("Float_t",0.9*lep.pt*(1.+lep.jetRelIso)))
                                                                  # else: return 0.90 * lep.pt / jetPtRatio where jetPtRatio = 1./(Electron_jetRelIso + 1.)

                    # Def conept : https://github.com/CERN-PH-CMG/cmgtools-lite/blob/f8a34c64a4489d94ff9ac4c0d8b0b06dad46e521/TTHAnalysis/python/tools/conept.py#L74
        self.lambda_conept_muon = lambda lep : op.multiSwitch((op.AND(op.abs(lep.pdgId)!=11 , op.abs(lep.pdgId)!=13) , op.static_cast("Float_t",lep.pt)),
                                                               # if (abs(lep.pdgId)!=11 and abs(lep.pdgId)!=13): return lep.pt : anything that is not muon or electron
                                                              (op.AND(op.abs(lep.pdgId)==13 , lep.mediumId ,lep.mvaTTH > 0.85) , op.static_cast("Float_t",lep.pt)),
                                                               # if (abs(lep.pdgId)!=13 or lep.mediumMuonId>0) and lep.mvaTTH > 0.90: return lep.pt 
                                                                    # if electron, check above MVA 
                                                                    # If muon, check that passes medium and above MVA
                                                               op.static_cast("Float_t",0.9*lep.pt*(1.+lep.jetRelIso)))
                                                               # else: return 0.90 * lep.pt / lep.jetPtRatiov2

        # Btag interpolation #
                    # https://indico.cern.ch/event/812025/contributions/3475878/attachments/1867083/3070589/gp-fr-run2b.pdf (slide 7)
        self.lambda_muon_x = lambda mu : op.min(op.max(0.,(0.9*mu.pt*(1+mu.jetRelIso))-20.)/(45.-20.), 1.)
                    # x = min(max(0, jet_pt-PT_min)/(PT_max-PT_min), 1) where jet_pt = 0.9*PT_muon*(1+MuonJetRelIso), PT_min=25, PT_max=40
        if era == "2016": 
            self.lambda_muon_btagInterpolation = lambda mu : self.lambda_muon_x(mu)*0.0614 + (1-self.lambda_muon_x(mu))*0.3093
        elif era =="2017":
            self.lambda_muon_btagInterpolation = lambda mu : self.lambda_muon_x(mu)*0.0521 + (1-self.lambda_muon_x(mu))*0.3033
        elif era == "2018":
            self.lambda_muon_btagInterpolation = lambda mu : self.lambda_muon_x(mu)*0.0494 + (1-self.lambda_muon_x(mu))*0.2770
            # return x*WP_loose+(1-x)*WP_medium
        #self.lambda_muon_deepJetInterpIfMvaFailed = lambda mu : mu.jet.btagDeepFlavB < self.lambda_muon_btagInterpolation(mu)
        self.lambda_muon_deepJetInterpIfMvaFailed = lambda mu : op.OR(op.NOT(self.lambda_hasAssociatedJet(mu)),     # If no associated jet, isolated lepton : cool !
                                                                      mu.jet.btagDeepFlavB < self.lambda_muon_btagInterpolation(mu))

        # Dilepton lambdas #
        self.lambda_leptonOS  = lambda l1,l2 : l1.charge != l2.charge
        if self.is_MC:
            self.lambda_is_matched = lambda lep : op.OR(lep.genPartFlav==1,  # Prompt muon or electron
                                                        lep.genPartFlav==15, # From tau decay
                                                        lep.genPartFlav==22) # From photon conversion (only available for electrons)
        else:
            self.lambda_is_matched = lambda lep : op.c_bool(True)
        lambda_dilepton_OS_matched  = lambda l1,l2 : op.AND(self.lambda_leptonOS(l1,l2),
                                                            self.lambda_is_matched(l1),
                                                            self.lambda_is_matched(l2))

        #############################################################################
        #                                 Muons                                     #
        #############################################################################

        # Preselection #
        muonsByPt = op.sort(t.Muon, lambda mu : -self.lambda_conept_muon(mu)) # Ordering done by conept
        self.lambda_muonPreSel = lambda mu : op.AND(
                                                    mu.p4.Pt() >= 5.,
                                                    op.abs(mu.p4.Eta()) <= 2.4,
                                                    op.abs(mu.dxy) <= 0.05,
                                                    op.abs(mu.dz) <= 0.1,
                                                    mu.miniPFRelIso_all <= 0.4, # mini PF relative isolation, total (with scaled rho*EA PU corrections)
                                                    mu.sip3d <= 8,
                                                    mu.looseId,
                                              )
        self.muonsPreSel = op.select(muonsByPt, self.lambda_muonPreSel)
        # Fakeable selection #
        self.lambda_muonFakeSel = lambda mu : op.AND(
                                                    self.lambda_conept_muon(mu) >= op.c_float(10.),
                                                    self.lambda_lepton_associatedJetNoBtag(mu),
                                                    op.OR(mu.mvaTTH >= 0.85, op.AND(mu.jetRelIso<0.5 , self.lambda_muon_deepJetInterpIfMvaFailed(mu))), 
                                                        # If mvaTTH < 0.85 : jetRelIso <0.5 and < deepJet medium with interpolation
                                                )
        self.muonsFakeSel = op.select(self.muonsPreSel, self.lambda_muonFakeSel)
        # Tight selection #
        self.lambda_muonTightSel = lambda mu : op.AND(mu.mvaTTH >= 0.85, # Lepton MVA id from ttH
                                                      mu.mediumId)
        self.muonsTightSel = op.select(self.muonsFakeSel, self.lambda_muonTightSel)

        #############################################################################
        #                              Electrons                                    #
        #############################################################################
        # Preselection #
        electronsByPt = op.sort(t.Electron, lambda ele : -self.lambda_conept_electron(ele)) # Ordering done by conept
        self.lambda_electronPreSel = lambda ele : op.AND(
                                                        ele.p4.Pt() >= 7.,
                                                        op.abs(ele.p4.Eta()) <= 2.5,
                                                        op.abs(ele.dxy) <= 0.05,
                                                        op.abs(ele.dz) <= 0.1,
                                                        ele.miniPFRelIso_all <= 0.4, # mini PF relative isolation, total (with scaled rho*EA PU corrections)

                                                        ele.sip3d <= 8,
                                                        ele.mvaFall17V2noIso_WPL, 
                                                        ele.lostHits <=1    # number of missing inner hits
                                                        )
        self.electronsPreSelInclu = op.select(electronsByPt, self.lambda_electronPreSel) # can include a muon in cone
        self.lambda_cleanElectron = lambda ele : op.NOT(op.rng_any(self.muonsPreSel, lambda mu : op.deltaR(mu.p4, ele.p4) <= 0.3 ))
            # No overlap between electron and muon in cone of DR<=0.3
        self.electronsPreSel = op.select(self.electronsPreSelInclu, self.lambda_cleanElectron)
        # Fakeable selection #
        self.lambda_electronFakeSel = lambda ele : op.AND(
                                                        self.lambda_conept_electron(ele) >= 10,
                                                        op.OR(
                                                                op.AND(op.abs(ele.eta+ele.deltaEtaSC)<=1.479, ele.sieie<=0.011), 
                                                                #op.AND(op.AND(op.abs(ele.eta+ele.deltaEtaSC)>=1.479,op.abs(ele.eta+ele.deltaEtaSC)<=2.5), ele.sieie<=0.030)),
                                                                op.AND(op.abs(ele.eta+ele.deltaEtaSC)>1.479, ele.sieie<=0.030)),
                                                        self.lambda_lepton_associatedJetNoBtag(ele),
                                                        ele.hoe <= 0.10,
                                                        ele.eInvMinusPInv >= -0.04,
                                                        op.OR(ele.mvaTTH >= 0.80, op.AND(ele.jetRelIso<0.7, ele.mvaFall17V2noIso_WP80)), # Lepton MVA id from ttH
                                                            # If mvaTTH < 0.80 : jetRelIso <0.5 and Fall17V2noIso_WP80
                                                        ele.lostHits == 0,    # number of missing inner hits
                                                        ele.convVeto        # Passes conversion veto
                                                        )
        self.electronsFakeSel = op.select(self.electronsPreSel, self.lambda_electronFakeSel)
        # Tight selection #
        self.lambda_electronTightSel = lambda ele : ele.mvaTTH >= 0.80
        self.electronsTightSel = op.select(self.electronsFakeSel, self.lambda_electronTightSel)

        #############################################################################
        #                               Dileptons                                   #
        #############################################################################
        # Presecleted dilepton #
        self.OSElElDileptonPreSel = op.combine(self.electronsPreSel, N=2, pred=self.lambda_leptonOS)
        self.OSMuMuDileptonPreSel = op.combine(self.muonsPreSel, N=2, pred=self.lambda_leptonOS)
        self.OSElMuDileptonPreSel = op.combine((self.electronsPreSel, self.muonsPreSel), pred=self.lambda_leptonOS)

        # Fakeable dilepton #
        self.OSElElDileptonFakeSel = op.combine(self.electronsFakeSel, N=2, pred=self.lambda_leptonOS)
        self.OSMuMuDileptonFakeSel = op.combine(self.muonsFakeSel, N=2, pred=self.lambda_leptonOS)
        self.OSElMuDileptonFakeSel = op.combine((self.electronsFakeSel, self.muonsFakeSel), pred=self.lambda_leptonOS)

        # Tight Dilepton : must also be Gen matched if MC #
        self.OSElElDileptonTightSel = op.combine(self.electronsTightSel, N=2, 
                                                 pred=lambda_dilepton_OS_matched)
        self.OSMuMuDileptonTightSel = op.combine(self.muonsTightSel, N=2, 
                                                 pred=lambda_dilepton_OS_matched)
        self.OSElMuDileptonTightSel = op.combine((self.electronsTightSel, self.muonsTightSel), 
                                                 pred=lambda_dilepton_OS_matched)
        # Fake extrapolated Dilepton : must also be Gen matched if MC #
        self.OSElElDileptonFakeExtrapolationSel = op.combine(self.electronsFakeSel, N=2, 
                                                 pred=lambda_dilepton_OS_matched)
        self.OSMuMuDileptonFakeExtrapolationSel = op.combine(self.muonsFakeSel, N=2, 
                                                 pred=lambda_dilepton_OS_matched)
        self.OSElMuDileptonFakeExtrapolationSel = op.combine((self.electronsFakeSel, self.muonsFakeSel), 
                                                 pred=lambda_dilepton_OS_matched)

        #############################################################################
        #                                AK4 Jets                                   #
        #############################################################################
        self.ak4JetsByPt = op.sort(t.Jet, lambda jet : -jet.p4.Pt())
        # Preselection #
        if era == "2016":
            self.lambda_ak4JetsPreSel = lambda j : op.AND(
                                                        j.jetId & 1, # Jet ID flags bit1 is loose, bit2 is tight, bit3 is tightLepVeto
                                                        j.p4.Pt() >= 25.,
                                                        op.abs(j.p4.Eta())<= 2.4
                                                    )
        elif  era == "2017" or era == "2018":
            self.lambda_ak4JetsPreSel = lambda j : op.AND(
                                                        j.jetId & 2, # Jet ID flags bit1 is loose, bit2 is tight, bit3 is tightLepVeto
                                                        j.p4.Pt() >= 25.,
                                                        op.abs(j.p4.Eta())<= 2.4
                                                    )
        self.ak4JetsPreSel = op.select(self.ak4JetsByPt, self.lambda_ak4JetsPreSel)
        # Cleaning #
        self.lambda_cleanAk4Jets = lambda j : op.AND(
                                                op.NOT(op.rng_any(self.electronsFakeSel, lambda ele : op.deltaR(j.p4, ele.p4) <= 0.4 )), 
                                                op.NOT(op.rng_any(self.muonsFakeSel, lambda mu : op.deltaR(j.p4, mu.p4) <= 0.4 )))
            # remove jets within cone of DR<0.4 of preselected electrons and muons
        self.ak4Jets = op.select(self.ak4JetsPreSel,self.lambda_cleanAk4Jets)
    
        ############     Btagging     #############
        # The pfDeepFlavour (DeepJet) algorithm is used
        if era == "2016": 
            self.lambda_ak4Btag =   lambda jet    : jet.btagDeepFlavB > 0.3093
            self.lambda_ak4NoBtag = lambda jet    : jet.btagDeepFlavB <= 0.3093
        elif era =="2017":
            self.lambda_ak4Btag =   lambda jet    : jet.btagDeepFlavB > 0.3033
            lambda_ak4NoBtag = lambda jet    : jet.btagDeepFlavB <= 0.3033
        elif era == "2018":
            self.lambda_ak4Btag =   lambda jet    : jet.btagDeepFlavB > 0.2770
            self.lambda_ak4NoBtag = lambda jet    : jet.btagDeepFlavB <= 0.2770

        self.ak4BJets     = op.select(self.ak4Jets, self.lambda_ak4Btag)
        self.ak4LightJets = op.select(self.ak4Jets, self.lambda_ak4NoBtag)

        #############################################################################
        #                                AK8 Jets                                   #
        #############################################################################
        self.ak8JetsByPt = op.sort(t.FatJet, lambda jet : -jet.p4.Pt())
        # Preselection #
        if era == "2016":
            self.lambda_ak8JetsPreSel = lambda j : op.AND(
                                                        j.jetId & 1, # Jet ID flags bit1 is loose, bit2 is tight, bit3 is tightLepVeto
                                                        j.p4.Pt() >= 200.,
                                                        op.abs(j.p4.Eta())<= 2.4,
                                                        op.AND(j.subJet1._idx.result != -1, j.subJet1.pt >= 20. , op.abs(j.subJet1.eta)<=2.4,
                                                               j.subJet2._idx.result != -1, j.subJet2.pt >= 20. , op.abs(j.subJet2.eta)<=2.4),
                                                               # Fatjet subjets must exist before checking Pt and eta 
                                                        op.AND(j.msoftdrop >= 30, j.msoftdrop <= 210),
                                                        j.tau2/j.tau1 <= 0.75
                                                    )
        elif  era == "2017" or era == "2018":
            self.lambda_ak8JetsPreSel = lambda j : op.AND(
                                                        j.jetId & 2, # Jet ID flags bit1 is loose, bit2 is tight, bit3 is tightLepVeto
                                                        j.p4.Pt() >= 200.,
                                                        op.abs(j.p4.Eta())<= 2.4,
                                                        op.AND(j.subJet1._idx.result != -1, j.subJet1.pt >= 20. , op.abs(j.subJet1.eta)<=2.4,
                                                               j.subJet2._idx.result != -1, j.subJet2.pt >= 20. , op.abs(j.subJet2.eta)<=2.4),
                                                               # Fatjet subjets must exist before checking Pt and eta 
                                                        op.AND(j.msoftdrop >= 30, j.msoftdrop <= 210),
                                                        j.tau2/j.tau1 <= 0.75
                                                    )
        self.ak8JetsPreSel = op.select(self.ak8JetsByPt, self.lambda_ak8JetsPreSel)
        # Cleaning #
        self.lambda_cleanAk8Jets = lambda j : op.AND(
                                                op.NOT(op.rng_any(self.electronsFakeSel, lambda ele : op.deltaR(j.p4, ele.p4) <= 0.8 )), 
                                                op.NOT(op.rng_any(self.muonsFakeSel, lambda mu : op.deltaR(j.p4, mu.p4) <= 0.8 ))
                                            )
            # remove jets within cone of DR<0.8 of preselected electrons and muons
        self.ak8Jets = op.select(self.ak8JetsPreSel,self.lambda_cleanAk8Jets)

        ############     Btagging     #############
        # The DeepCSV b-tagging algorithm is used on subjets 
        if era == "2016": 
            self.lambda_ak8Btag = lambda fatjet : op.OR(op.AND(fatjet.subJet1.pt >= 30, fatjet.subJet1.btagDeepB > 0.6321),
                                                        op.AND(fatjet.subJet2.pt >= 30, fatjet.subJet2.btagDeepB > 0.6321))
        elif era =="2017":
            self.lambda_ak8Btag = lambda fatjet : op.OR(op.AND(fatjet.subJet1.pt >= 30, fatjet.subJet1.btagDeepB > 0.4941),
                                                        op.AND(fatjet.subJet2.pt >= 30, fatjet.subJet2.btagDeepB > 0.4941))
        elif era == "2018":
            self.lambda_ak8Btag = lambda fatjet : op.OR(op.AND(fatjet.subJet1.pt >= 30, fatjet.subJet1.btagDeepB > 0.4184),
                                                        op.AND(fatjet.subJet2.pt >= 30, fatjet.subJet2.btagDeepB > 0.4184))

        self.ak8BJets = op.select(self.ak8Jets, self.lambda_ak8Btag)

        #############################################################################
        #                             Scalefactors                                  #
        #############################################################################
        if self.is_MC:
            # Initialize scalefactors class #
            SF = ScaleFactorsbbWW()    
            ####  Muons ####
            self.muLooseId = SF.get_scalefactor("lepton", 'muon_loose_{}'.format(era), combine="weight", systName="mu_loose")
            self.muTightMVA = SF.get_scalefactor("lepton", 'muon_tightMVA_{}'.format(era), combine="weight", systName="mu_tightmva")
   
            ####  Electrons ####
            self.elLooseRecoPtLt20 = SF.get_scalefactor("lepton", ('electron_loosereco_{}'.format(era) , 'electron_loosereco_ptgt20'), combine="weight", systName="el_looserecoptlt20")
            self.elLooseRecoPtGt20 = SF.get_scalefactor("lepton", ('electron_loosereco_{}'.format(era) , 'electron_loosereco_ptlt20'), combine="weight", systName="el_looserecoptgt20")
            self.elLooseId = SF.get_scalefactor("lepton", 'electron_looseid_{}'.format(era) , combine="weight", systName="el_looseid")
            self.elLooseEff = SF.get_scalefactor("lepton", 'electron_looseeff_{}'.format(era) , combine="weight", systName="el_looseeff")
            self.elTightMVA = SF.get_scalefactor("lepton", 'electron_tightMVA_{}'.format(era) , combine="weight", systName="el_tightmva")
 
            #### Triggers (from TTH) ####
            self.ttH_doubleMuon_trigSF = op.systematic(op.c_float(1.010), name="ttH_doubleMuon_trigSF", up=op.c_float(1.020), down=op.c_float(1.000))
            self.ttH_doubleElectron_trigSF = op.systematic(op.c_float(1.020), name="ttH_doubleElectron_trigSF", up=op.c_float(1.040), down=op.c_float(1.000))
            self.ttH_electronMuon_trigSF = op.systematic(op.c_float(1.020), name="ttH_electronMuon_trigSF", up=op.c_float(1.030), down=op.c_float(1.010))
            #### Ak4 Btagging ####
            DeepJetTag_discriVar = {"BTagDiscri": lambda j : j.btagDeepFlavB}
            self.DeepJetMediumSF = SF.get_scalefactor("jet", ("btag_"+era+"_"+self.sfTag, "DeepJet_medium"), additionalVariables=DeepJetTag_discriVar, systName="deepjet") # For RESOLVED
            #### Ak8 Btagging ####
 #           DeepCSVTag_discriVar = {"BTagDiscri": lambda j : j.btagDeepB}
 #           self.DeepCSVMediumSF = SF.get_scalefactor("jet", ("subjet_btag_"+era+"_"+self.sfTag, "DeepCSV_medium"), additionalVariables=DeepCSVTag_discriVar, systName="deepcsv") # For BOOSTED (btag o    n subjet)
                # TODO : add Boosted SF (need NanoAODv7)
        



        # Return #
        return noSel
