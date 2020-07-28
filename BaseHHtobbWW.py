import os
import sys
from itertools import chain

from bamboo import treefunctions as op
from bamboo.analysismodules import NanoAODModule, NanoAODHistoModule, NanoAODSkimmerModule
from bamboo.analysisutils import makeMultiPrimaryDatasetTriggerSelection
from bamboo.scalefactors import binningVariables_nano, BtagSF

from METScripts import METFilter, METcorrection
from scalefactorsbbWW import ScaleFactorsbbWW
from btagHelper import makeBtagRatioReweighting
from triggers import returnTriggerRanges

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
                                 "ratio-y-axis-range" : [0.5,1.5],
                                 "ratio-y-axis" : 'Ratio Data/MC',
                                 "sort-by-yields" : False}

    #-------------------------------------------------------------------------------------------#
    #                                       addArgs                                             #
    #-------------------------------------------------------------------------------------------#
    def addArgs(self,parser):
        super(BaseNanoHHtobbWW, self).addArgs(parser)

        parser.title = """

Arguments for the HH->bbWW analysis on bamboo framework

----- Argument groups -----
    * Lepton arguments *
        --Preselected
        --Fakeable
        --Tight
        --FakeExtrapolation
        --NoZVeto
        --ZPeak

    * Jet arguments *
        --Ak4
        --Ak8
        --Resolved0Btag
        --Resolved1Btag
        --Resolved2Btag
        --Boosted
        --TTBarCR
        --BtagReweightingOn
        --BtagReweightingOff

    * Skimmer arguments *
        --Synchronization
        --Channel

    * Plotter arguments *
        --OnlyYield

    * DY correction *
        --DYDataEstimation1Btag
        --DYDataEstimation2Btag

    * Technical *
        --backend
        --NoSystematics

----- Plotter mode -----
Every combination of lepton and jet arguments can be used, if none specified they will all be plotted

----- Skimmer mode -----
One lepton and and one jet argument must be specified in addition to the required channel
(this is because only one selection at a time can produce a ntuple)

----- Detailed explanation -----

                    """
        #----- Technical arguments -----#
        parser.add_argument("--backend", 
                            type=str, 
                            default="dataframe", 
                            help="Backend to use, 'dataframe' (default) or 'lazy'")
        parser.add_argument("--NoSystematics", 
                            action      = "store_true",
                            default     = False,
                            help="Disable all systematic variations (default=False)")


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
                            help        = "Remove the cut of preselected leptons |M_ll-M_Z|>10 GeV")
        parser.add_argument("--ZPeak", 
                            action      = "store_true",
                            default     = False,
                            help        = "Select the Z peak at tight level |M_ll-M_Z|<10 GeV (must be used with --NoZVeto, only effective with --Tight)")

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
        parser.add_argument("--TTBarCR", 
                                action      = "store_true",
                                default     = False,
                                help        = "Apply cut on Mbb for ttbar CR (only effective with --Resolved2Btag)")
        parser.add_argument("--BtagReweightingOn", 
                                action      = "store_true",
                                default     = False,
                                help        = "Btag ratio study : Btag SF applied (without the ratio), will only do the plots for reweighting (jets and leptons args are ignored)")
        parser.add_argument("--BtagReweightingOff", 
                                action      = "store_true",
                                default     = False,
                                help        = "Btag ratio study : Btag Sf not applied (without the ratio), will only do the plots for reweighting (jets and leptons args are ignored)")
        parser.add_argument("--DYStitching", 
                                action      = "store_true",
                                default     = False,
                                help        = "DY stitching studies : only produce LHE jet multiplicities (inclusive analysis, only on DY events, rest of plots ignored)")


        #----- Skimmer arguments -----#
        parser.add_argument("--Synchronization", 
                            action      = "store_true",
                            default     = False,
                            help        = "Produce the skims for the synchronization (without triggers, corrections of flags) if alone. If sync for specific selection, lepton, jet and channel arguments need to be used")
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
        from bamboo.treedecorators import NanoAODDescription, nanoRochesterCalc, nanoJetMETCalc, nanoJetMETCalc_METFixEE2017
        # JEC's Recommendation for Full RunII: https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC
        # JER : -----------------------------: https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
        # Get base aguments #
        era = sampleCfg['era']
        self.is_MC = self.isMC(sample) # no confusion between boolean is_MC and method isMC()
        metName = "METFixEE2017" if era == "2017" else "MET"
        tree,noSel,be,lumiArgs = super(BaseNanoHHtobbWW,self).prepareTree(tree          = tree, 
                                                                          sample        = sample, 
                                                                          sampleCfg     = sampleCfg, 
                                                                          description   = NanoAODDescription.get(
                                                                                            tag             = "v5", 
                                                                                            year            = (era if era else "2016"),
                                                                                            isMC            = self.is_MC,
                                                                                            #systVariations  = [ nanoRochesterCalc, (nanoJetMETCalc_METFixEE2017 if era == "2017" else nanoJetMETCalc) ]), 
                                                                                            systVariations  = [ (nanoJetMETCalc_METFixEE2017 if era == "2017" else nanoJetMETCalc)]),
                                                                                            #systVariations = None ),
                                                                                            # will do Jet and MET variations, and the Rochester correction
                                                                          lazyBackend   = (self.args.backend == "lazy" or self.args.onlypost))
    
        self.triggersPerPrimaryDataset = {}
        from bamboo.analysisutils import configureJets ,configureRochesterCorrection, configureType1MET 

        # Check distributed option #
        isNotWorker = (self.args.distributed != "worker") 

        # Turn off systs #
        if self.args.NoSystematics:
            noSel = noSel.refine('SystOff',autoSyst=False)

        # Check era #
        if era != "2016" and era != "2017" and era != "2018":
            raise RuntimeError("Unknown era {0}".format(era))

        # Rochester and JEC corrections (depends on era) #     
        cachJEC_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),'cacheJEC')

        # Check if basic synchronization is required (no corrections and triggers) #
        self.inclusive_sel = ((self.args.Synchronization 
                                  and not any([self.args.__dict__[key] for key in['Preselected', 'Fakeable', 'Tight', 'FakeExtrapolation']]) \
                                  and not any([self.args.__dict__[key] for key in['Ak4', 'Ak8', 'Resolved0Btag', 'Resolved1Btag', 'Resolved2Btag', 'Boosted']]) \
                                  and (self.args.Channel is None or self.args.Channel=='None')) \
                                  # No channel selection
                                  # None is local mode, "None" in distributed mode
                              or self.args.BtagReweightingOn 
                              or self.args.BtagReweightingOff
                              or self.args.DYStitching)
                                # Inclusive plots
            # If no lepton, jet and channel selection : basic object selection (no trigger nor corrections)
        if self.inclusive_sel:
            print ("Inclusive analysis, no selections applied")

        #----- Theory uncertainties -----#
#        if self.is_MC:
#            qcdScaleVariations = { f"qcdScalevar{i}": tree.LHEScaleWeight[i] for i in [0, 1, 3, 5, 7, 8] }
#            qcdScaleSyst = op.systematic(op.c_float(1.), name="qcdScale", **qcdScaleVariations)
#            psISRSyst = op.systematic(op.c_float(1.), name="psISR", up=tree.PSWeight[2], down=tree.PSWeight[0])
#            psFSRSyst = op.systematic(op.c_float(1.), name="psFSR", up=tree.PSWeight[3], down=tree.PSWeight[1])
#            pdfsWeight = op.systematic(op.c_float(1.), name="pdfsWgt", up=tree.LHEPdfWeight, down=tree.LHEPdfWeight)
#            noSel = noSel.refine("theorySystematics", weight = [qcdScaleSyst,psISRSyst,psFSRSyst,pdfsWeight])

        #----- Triggers and Corrections -----#
        self.triggersPerPrimaryDataset = {}
        def addHLTPath(key,HLT):
            if key not in self.triggersPerPrimaryDataset.keys():
                self.triggersPerPrimaryDataset[key] = []
            try:
                self.triggersPerPrimaryDataset[key].append(getattr(tree.HLT,HLT))
            except AttributeError:
                print ("Could not find branch tree.HLT.%s, will omit it"%HLT)

        ############################################################################################
        # ERA 2016 #
        ############################################################################################
        if era == "2016":
#            # Rochester corrections #
#            configureRochesterCorrection(variProxy  = tree._Muon,
#                                         paramsFile = os.path.join(os.path.dirname(__file__), "data", "RoccoR2016.txt"),
#                                         isMC       = self.is_MC,
#                                         backend    = be, 
#                                         uName      = sample)

##                self.triggersPerPrimaryDataset = {
##                    "SingleMuon" :  [ tree.HLT.IsoMu24],
##                    "SingleElectron":  [ tree.HLT.Ele27_WPTight_Gsf],
##                    "DoubleMuon" :  [ tree.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL,
##                                      tree.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ],
##                    "DoubleEGamma": [ tree.HLT.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ],
##                    "MuonEG":       [ tree.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,
##                                      tree.HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ,
##                                      tree.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL,
##                                      tree.HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL]}
 
            # SingleMuon #
            addHLTPath("SingleMuon","IsoMu22")
            addHLTPath("SingleMuon","IsoTkMu22")
            addHLTPath("SingleMuon","IsoMu22_eta2p1")
            addHLTPath("SingleMuon","IsoTkMu22_eta2p1")
            addHLTPath("SingleMuon","IsoMu24")
            addHLTPath("SingleMuon","IsoTkMu24")
            # SingleElectron #
            addHLTPath("SingleElectron","Ele27_WPTight_Gsf")
            addHLTPath("SingleElectron","Ele25_eta2p1_WPTight_Gsf")
            addHLTPath("SingleElectron","Ele27_eta2p1_WPLoose_Gsf")
            # DoubleMuon #
            addHLTPath("DoubleMuon","Mu17_TrkIsoVVL_Mu8_TrkIsoVVL")
            addHLTPath("DoubleMuon","Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ")
            addHLTPath("DoubleMuon","Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL")
            addHLTPath("DoubleMuon","Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ")
            # DoubleEGamma #
            addHLTPath("DoubleEGamma","Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ")
            # MuonEG #
            addHLTPath("MuonEG","Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL")
            addHLTPath("MuonEG","Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ")
            addHLTPath("MuonEG","Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL")
            addHLTPath("MuonEG","Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ")

            # Links : 
            # JEC : https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECDataMC
            # JER (smear) : https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution

            # JetMET treatment #
            if not self.args.Synchronization:
                if self.is_MC:   # if MC -> needs smearing
                    configureJets(variProxy             = tree._Jet, 
                                  jetType               = "AK4PFchs",
                                  jec                   = "Summer16_07Aug2017_V11_MC",
                                  smear                 = "Summer16_25nsV1_MC",
                                  jesUncertaintySources = ["Total"],
                                  mayWriteCache         = isNotWorker,
                                  isMC                  = self.is_MC,
                                  backend               = be, 
                                  uName                 = sample,
                                  cachedir              = cachJEC_dir)
                    configureType1MET(variProxy             = getattr(tree, f"_{metName}"),
                                      jec                   = "Summer16_07Aug2017_V11_MC",
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
                    else:
                        raise RuntimeError("Could not find appropriate JEC tag for data")
                    configureJets(variProxy             = tree._Jet, 
                                  jetType               = "AK4PFchs",
                                  jec                   = jecTag,
                                  mayWriteCache         = isNotWorker,
                                  isMC                  = self.is_MC,
                                  backend               = be, 
                                  uName                 = sample,
                                  cachedir              = cachJEC_dir)
                    configureType1MET(variProxy         = getattr(tree, f"_{metName}"),
                                      jec               = jecTag,
                                      mayWriteCache     = isNotWorker,
                                      isMC              = self.is_MC,
                                      backend           = be, 
                                      uName             = sample,
                                      cachedir          = cachJEC_dir)

        ############################################################################################
        # ERA 2017 #
        ############################################################################################
        elif era == "2017":
            #configureRochesterCorrection(variProxy  = tree._Muon,
            #                             paramsFile = os.path.join(os.path.dirname(__file__), "data", "RoccoR2017.txt"),
            #                             isMC       = self.is_MC,
            #                             backend    = be, 
            #                             uName      = sample)

            # SingleMuon #
            addHLTPath("SingleMuon","IsoMu24")
            addHLTPath("SingleMuon","IsoMu27")
            # SingleElectron #
            addHLTPath("SingleElectron","Ele35_WPTight_Gsf")
            addHLTPath("SingleElectron","Ele32_WPTight_Gsf")
            # DoubleMuon #
            addHLTPath("DoubleMuon","Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8")
            addHLTPath("DoubleMuon","Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8")
            # DoubleEGamma #
            addHLTPath("DoubleEGamma","Ele23_Ele12_CaloIdL_TrackIdL_IsoVL")
            # MuonEG #
            addHLTPath("MuonEG","Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ")
            addHLTPath("MuonEG","Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ")
            addHLTPath("MuonEG","Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL")
            addHLTPath("MuonEG","Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ")
            
            # JetMET treatment #
            if not self.args.Synchronization:
                if self.is_MC:   # if MC -> needs smearing
                    configureJets(variProxy             = tree._Jet, 
                                  jetType               = "AK4PFchs",
                                  jec                   = "Fall17_17Nov2017_V32_MC",
                                  smear                 = "Fall17_V3b_MC",
                                  jesUncertaintySources = ["Total"],
                                  mayWriteCache         = isNotWorker,
                                  isMC                  = self.is_MC,
                                  backend               = be, 
                                  uName                 = sample,
                                  cachedir              = cachJEC_dir)
                    configureType1MET(variProxy             = getattr(tree, f"_{metName}"),
                                      jec                   = "Fall17_17Nov2017_V32_MC",
                                      smear                 = "Fall17_V3b_MC",
                                      jesUncertaintySources = ["Total"],
                                      mayWriteCache         = isNotWorker,
                                      isMC                  = self.is_MC,
                                      backend               = be,
                                      uName                 = sample,
                                      cachedir              = cachJEC_dir)

                else:                   # If data -> extract info from config 
                    jecTag = None
                    if "2017B" in sample:
                        jecTag = "Fall17_17Nov2017B_V32_DATA"
                    elif "2017C" in sample:
                        jecTag = "Fall17_17Nov2017C_V32_DATA"
                    elif "2017D" in sample or "2017E" in sample:
                        jecTag = "Fall17_17Nov2017DE_V32_DATA"
                    elif "2017F" in sample:
                        jecTag = "Fall17_17Nov2017F_V32_DATA"
                    else:
                        raise RuntimeError("Could not find appropriate JEC tag for data")
                    configureJets(variProxy             = tree._Jet, 
                                  jetType               = "AK4PFchs",
                                  jec                   = jecTag,
                                  mayWriteCache         = isNotWorker,
                                  isMC                  = self.is_MC,
                                  backend               = be, 
                                  uName                 = sample,
                                  cachedir              = cachJEC_dir)
                    configureType1MET(variProxy         = getattr(tree, f"_{metName}"),
                                      jec               = jecTag,
                                      mayWriteCache     = isNotWorker,
                                      isMC              = self.is_MC,
                                      backend           = be, 
                                      uName             = sample,
                                      cachedir          = cachJEC_dir)

        ############################################################################################
        # ERA 2018 #
        ############################################################################################
        elif era == "2018":
            #configureRochesterCorrection(variProxy  = tree._Muon,
            #                             paramsFile = os.path.join(os.path.dirname(__file__), "data", "RoccoR2018.txt"),
            #                             isMC       = self.is_MC,
            #                             backend    = be, 
            #                             uName      = sample)

            # SingleMuon #
            addHLTPath("SingleMuon","IsoMu24")
            addHLTPath("SingleMuon","IsoMu27")
            # SingleElectron #
            addHLTPath("SingleElectron","Ele32_WPTight_Gsf")
            # DoubleMuon #
            addHLTPath("DoubleMuon","Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8")
            # DoubleEGamma #
            addHLTPath("DoubleEGamma","Ele23_Ele12_CaloIdL_TrackIdL_IsoVL")
            # MuonEG #
            addHLTPath("MuonEG","Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ")
            addHLTPath("MuonEG","Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ")
            addHLTPath("MuonEG","Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL")

            # JetMET treatment #
            if not self.args.Synchronization:
                if self.is_MC:   # if MC -> needs smearing
                    configureJets(variProxy             = tree._Jet, 
                                  jetType               = "AK4PFchs",
                                  jec                   = "Autumn18_V19_MC",
                                  smear                 = "Autumn18_V7b_MC",
                                  jesUncertaintySources = ["Total"],
                                  mayWriteCache         = isNotWorker,
                                  isMC                  = self.is_MC,
                                  backend               = be, 
                                  uName                 = sample,
                                  cachedir              = cachJEC_dir)
                    configureType1MET(variProxy             = getattr(tree, f"_{metName}"),
                                      jec                   = "Autumn18_V19_MC",
                                      smear                 = "Autumn18_V7b_MC",
                                      jesUncertaintySources = ["Total"],
                                      mayWriteCache         = isNotWorker,
                                      isMC                  = self.is_MC,
                                      backend               = be,
                                      uName                 = sample,
                                      cachedir              = cachJEC_dir)

                else:                   # If data -> extract info from config 
                    jecTag = None
                    if "2018A" in sample:
                        jecTag = "Autumn18_RunA_V19_DATA"
                    elif "2018B" in sample:
                        jecTag = "Autumn18_RunB_V19_DATA"
                    elif "2018C" in sample:
                        jecTag = "Autumn18_RunC_V19_DATA"
                    elif "2018D" in sample:
                        jecTag = "Autumn18_RunD_V19_DATA"
                    else:
                        raise RuntimeError("Could not find appropriate JEC tag for data")
                    configureJets(variProxy             = tree._Jet, 
                                  jetType               = "AK4PFchs",
                                  jec                   = jecTag,
                                  mayWriteCache         = isNotWorker,
                                  isMC                  = self.is_MC,
                                  backend               = be, 
                                  uName                 = sample,
                                  cachedir              = cachJEC_dir)
                    configureType1MET(variProxy         = getattr(tree, f"_{metName}"),
                                      jec               = jecTag,
                                      mayWriteCache     = isNotWorker,
                                      isMC              = self.is_MC,
                                      backend           = be, 
                                      uName             = sample,
                                      cachedir          = cachJEC_dir)

        return tree,noSel,be,lumiArgs

    def prepareObjects(self, t, noSel, sample, sampleCfg, channel, forSkimmer=False):
        # Some imports #
        from bamboo.analysisutils import forceDefine

        if channel not in ["DL","SL"]:
            raise RuntimeError('Channel %s not understood'%channel)

        era = sampleCfg['era']
        self.era = era

        # Forcedefine : calculate once per event (for every event) #
        if not self.args.Synchronization:
            #forceDefine(t._Muon.calcProd, noSel) # Muons for Rochester corrections
            forceDefine(t._Jet.calcProd, noSel)  # Jets for configureJets
            forceDefine(getattr(t, "_{0}".format("MET" if era != "2017" else "METFixEE2017")).calcProd,noSel) # MET for configureMET
        else:
            print ("No jet corrections applied")

        ###########################################################################
        #                           TTbar reweighting                             #
        ###########################################################################
        if "group" in sampleCfg and sampleCfg["group"] == 'ttbar': 
            print ('Applied TT top reweighting')
            # https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting#Use_case_3_ttbar_MC_is_used_to_m -> do not use when ttbar is background
            # Correct : https://indico.cern.ch/event/904971/contributions/3857701/attachments/2036949/3410728/TopPt_20.05.12.pdf
            #       -> Weight formula : slide 2
            #       -> top/antitop SF : slide 12 bottom left 
            # Get tops #
            self.genTop = op.select(t.GenPart,lambda g : op.AND(g.pdgId==6, g.statusFlags & ( 0x1 << 13)))
            self.genAntitop = op.select(t.GenPart,lambda g : op.AND(g.pdgId==-6, g.statusFlags & ( 0x1 << 13)))
                # statusFlags==13 : isLastCopy
                # Pdgid == 6 : top
            # Lambda to compute weight if there is a ttbar #
            self.t_SF = lambda t : op.exp(-2.02274e-01 + 1.09734e-04*t.pt + -1.30088e-07*t.pt**2 + (5.83494e+01/(t.pt+1.96252e+02)))
            self.ttbar_weight = lambda t,tbar : op.sqrt(self.t_SF(t)*self.t_SF(tbar))
            self.ttbar_sys = op.systematic(self.ttbar_weight(self.genTop[0],self.genAntitop[0]),
                                           name = "ttbarweightsyst",
                                           up   = op.pow(self.ttbar_weight(self.genTop[0],self.genAntitop[0]),2), # Up = apply twice
                                           down = op.c_float(1))                                                  # Down = not apply
            # Apply correction to TT #
            noSel = noSel.refine("ttbarWeight",weight=self.ttbar_sys)

        ###########################################################################
        #                            DY stitching                                 #
        ###########################################################################
#       Still experimental
#        if "group" in sampleCfg and sampleCfg["group"] == 'DY':
#            if sample == "DYJetsToLL_M-50":
#                if era == "2016":
#                    stitching_weights = [0.739,0.434,0.329,0.338]
#                if era == "2017":
#                    stitching_weights = [0.815,0.685,0.604,0.582]
#                if era == "2018":
#                    stitching_weights = [0.801,0.655,0.542,0.533]
#            elif sample == "DYToLL_0J":
#                if era == "2016":
#                    stitching_weights = [0.261,0.153,0.116,0.119]
#                if era == "2017":
#                    stitching_weights = [0.146,0.123,0.109,0.105]
#                if era == "2018":
#                    stitching_weights = [0.165,0.135,0.111,0.110]
#            elif sample == "DYToLL_1J":
#                if era == "2016":
#                    stitching_weights = [1.083,0.635,0.481,0.495]
#                if era == "2017":
#                    stitching_weights = [0.652,0.548,0.483,0.465] 
#                if era == "2018":
#                    stitching_weights = [0.697,0.570,0.471,0.464]
#            elif sample == "DYToLL_2J":
#                if era == "2016":
#                    stitching_weights = [1.480,0.868,0.658,0.676]
#                if era == "2017":
#                    stitching_weights = [0.547,0.460,0.406,0.391]
#                if era == "2018":
#                    stitching_weights = [0.682,0.557,0.461,0.454]
#            else:
#                stitching_weights = [1.,1.,1.,1.]
#            print ('Stitching weights :',stitching_weights)
#            noSel = noSel.refine("DYStitching",weight = op.multiSwitch((t.LHE.Njets==0,op.c_float(stitching_weights[0])),
#                                                                       (t.LHE.Njets==1,op.c_float(stitching_weights[1])),
#                                                                       (t.LHE.Njets==2,op.c_float(stitching_weights[2])),
#                                                                       (t.LHE.Njets==3,op.c_float(stitching_weights[3])),
#                                                                       op.c_float(1.)))
                
        ###########################################################################
        #                               tH samples                                #
        ###########################################################################
        if sample.startswith('TH'):
            print ('Applied tH LHE weights')
            noSel = noSel.refine("tHWeight",weight=t.LHEReweightingWeight[12])

        #############################################################################
        #                            Pre-firing rates                               #
        #############################################################################
        if era in ["2016","2017"]:
            self.L1Prefiring = op.systematic(t.L1PreFiringWeight_Nom,
                                             name = "L1PreFiring",
                                             up   = t.L1PreFiringWeight_Up,
                                             down = t.L1PreFiringWeight_Dn)
            
            noSel = noSel.refine("L1PreFiringRate", weight = self.L1Prefiring)

        #############################################################################
        #                             Pile-up                                       #
        #############################################################################
        # Get MC PU weight file #
        puWeightsFile = None
        if self.is_MC:
            puWeightsFile = os.path.join(os.path.dirname(__file__), "data", "pileup",sample+'_%s.json'%era)
            if not os.path.exists(puWeightsFile):
                raise RuntimeError("Could not find pileup file %s"%puWeightsFile)
            from bamboo.analysisutils import makePileupWeight
            PUWeight = makePileupWeight(puWeightsFile, t.Pileup_nTrueInt, systName="pileup",nameHint=f"puweightFromFile{sample}".replace('-','_'))
            noSel = noSel.refine("puWeight", weight = PUWeight)

        #############################################################################
        #                                 MET                                       #
        #############################################################################
        # MET filter #
        if not self.inclusive_sel:
            noSel = noSel.refine("passMETFlags", cut=METFilter(t.Flag, era, self.is_MC) )

        # MET corrections #
        MET = t.MET if era != "2017" else t.METFixEE2017
        #if not self.inclusive_sel:
        #    self.corrMET = METcorrection(MET,t.PV,sample,era,self.is_MC) # Flatness correction might not be needed
        #else:
        #    self.corrMET = MET
        self.corrMET = MET

        #############################################################################
        #                      Lepton Lambdas Variables                             #
        #############################################################################
        # Associated jet Btagging #
        self.lambda_hasAssociatedJet = lambda lep : lep.jet.idx != -1
        if era == "2016": 
            self.lambda_lepton_associatedJetNoBtag = lambda lep : op.OR(op.NOT(self.lambda_hasAssociatedJet(lep)),
                                                                        lep.jet.btagDeepFlavB <= 0.3093)
        elif era =="2017":
            self.lambda_lepton_associatedJetNoBtag = lambda lep : op.OR(op.NOT(self.lambda_hasAssociatedJet(lep)),
                                                                        lep.jet.btagDeepFlavB <= 0.3033)
        elif era == "2018":
            self.lambda_lepton_associatedJetNoBtag = lambda lep : op.OR(op.NOT(self.lambda_hasAssociatedJet(lep)),
                                                                        lep.jet.btagDeepFlavB <= 0.2770)
                     # If no associated jet, isolated lepton : cool !  

        # Cone pt #
                    # Def conept : https://github.com/CERN-PH-CMG/cmgtools-lite/blob/f8a34c64a4489d94ff9ac4c0d8b0b06dad46e521/TTHAnalysis/python/tools/conept.py#L74
        self.lambda_conept_electron = lambda lep : op.multiSwitch((op.AND(op.abs(lep.pdgId)!=11 , op.abs(lep.pdgId)!=13) , op.static_cast("Float_t",lep.pt)),
                                                                  # if (abs(lep.pdgId)!=11 and abs(lep.pdgId)!=13): return lep.pt : anything that is not muon or electron
                                                                  (op.AND(op.abs(lep.pdgId)==11 , lep.mvaTTH > 0.80) , op.static_cast("Float_t",lep.pt)),
                                                                  # if electron, check above MVA 
                                                                   op.static_cast("Float_t",0.9*lep.pt*(1.+lep.jetRelIso)))
                                                                  # else: return 0.90 * lep.pt / jetPtRatio where jetPtRatio = 1./(Electron_jetRelIso + 1.)
        self.electron_conept = op.map(t.Electron, self.lambda_conept_electron)

                    # Def conept : https://github.com/CERN-PH-CMG/cmgtools-lite/blob/f8a34c64a4489d94ff9ac4c0d8b0b06dad46e521/TTHAnalysis/python/tools/conept.py#L74
        self.lambda_conept_muon = lambda lep : op.multiSwitch((op.AND(op.abs(lep.pdgId)!=11 , op.abs(lep.pdgId)!=13) , op.static_cast("Float_t",lep.pt)),
                                                               # if (abs(lep.pdgId)!=11 and abs(lep.pdgId)!=13): return lep.pt : anything that is not muon or electron
                                                              (op.AND(op.abs(lep.pdgId)==13 , lep.mediumId ,lep.mvaTTH > 0.85) , op.static_cast("Float_t",lep.pt)),
                                                               # if muon, check that passes medium and above MVA
                                                               op.static_cast("Float_t",0.9*lep.pt*(1.+lep.jetRelIso)))
                                                               # else: return 0.90 * lep.pt / lep.jetPtRatiov2
        self.muon_conept = op.map(t.Muon, self.lambda_conept_muon)

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
                                                        lep.genPartFlav==15) # From tau decay
                                                        #lep.genPartFlav==22) # From photon conversion (only available for electrons)
        else:
            self.lambda_is_matched = lambda lep : op.c_bool(True)
        
        #############################################################################
        #                                 Muons                                     #
        #############################################################################
        # Preselection #
        muonsByPt = op.sort(t.Muon, lambda mu : -self.muon_conept[mu.idx]) # Ordering done by conept
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
                                                    self.muon_conept[mu.idx] >= op.c_float(10.),
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
        electronsByPt = op.sort(t.Electron, lambda ele : -self.electron_conept[ele.idx]) # Ordering done by conept
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
                                                        self.electron_conept[ele.idx] >= 10,
                                                        op.OR(
                                                                op.AND(op.abs(ele.eta+ele.deltaEtaSC)<=1.479, ele.sieie<=0.011), 
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

        #----- Needed for both SL and DL -----#
        # Preselected dilepton -> Mll>12 cut #
        self.ElElDileptonPreSel = op.combine(self.electronsPreSelInclu, N=2)
        self.MuMuDileptonPreSel = op.combine(self.muonsPreSel, N=2)
        self.ElMuDileptonPreSel = op.combine((self.electronsPreSelInclu, self.muonsPreSel))
            # Need pairs without sign for one of the selection
            # We use electrons before muon cleaning : electron-muon pair that has low mass are typically close to each other
        # OS Preselected dilepton -> Z Veto #
        self.OSElElDileptonPreSel = op.combine(self.electronsPreSel, N=2, pred=self.lambda_leptonOS)
        self.OSMuMuDileptonPreSel = op.combine(self.muonsPreSel, N=2, pred=self.lambda_leptonOS)
        self.OSElMuDileptonPreSel = op.combine((self.electronsPreSel, self.muonsPreSel), pred=self.lambda_leptonOS) 

        #----- Needed only for DL -----#
        if channel == "DL":
            # We need to only select the leading (in conept) fakeable leptons and only work with those
            # First issue : electrons and muons and in separate containers that cannot be appended
            # Second issue : before we compare conept of first or second lepton in a container ... we need to check that it is there !
            # Solution : slice the containers of fakeable leptons but need to separate each case
            self.leadElectronsFakeSel = self.electronsFakeSel[:op.multiSwitch(
                                (op.rng_len(self.muonsFakeSel)==0 , op.min(op.static_cast("std::size_t", op.c_int(2)), op.rng_len(self.electronsFakeSel))),
                                    # No muon -> take at most the two first electrons
                                (op.rng_len(self.electronsFakeSel)==0 , op.static_cast("std::size_t", op.c_int(0))),
                                    # No electron -> slice at 0 
                                (op.AND(op.rng_len(self.electronsFakeSel)==1 , op.rng_len(self.muonsFakeSel)==1) , op.static_cast("std::size_t", op.c_int(1))),
                                    # One electron and one muon : slice electrons at 1
                                (op.AND(op.rng_len(self.electronsFakeSel)>=2 , op.rng_len(self.muonsFakeSel)==1), 
                                    # At least 2 electrons and 1 muon : If both electrons have higher conept, leading leptons will be ele+ele (slice at 2)
                                    #                                   If muon has higher conept than one of them, leading leptons will be ele+mu (slice at 1)
                                        op.switch(op.AND(self.electron_conept[self.electronsFakeSel[0].idx]>self.muon_conept[self.muonsFakeSel[0].idx],
                                                         self.electron_conept[self.electronsFakeSel[1].idx]>self.muon_conept[self.muonsFakeSel[0].idx]),
                                                  op.static_cast("std::size_t", op.c_int(2)),
                                                  op.static_cast("std::size_t", op.c_int(1)))),
                                (op.AND(op.rng_len(self.electronsFakeSel)==1 , op.rng_len(self.muonsFakeSel)>=2),
                                    # 1 electron and at least 2 muons: If electron had higher conept than one of the leading muons : leading will be ele+mu (slice at 1)
                                    #                                  If electron has lower conept than both muons : leading leptons will be mu+mu (slice at 0) 
                                        op.switch(op.OR(self.electron_conept[self.electronsFakeSel[0].idx]>self.muon_conept[self.muonsFakeSel[0].idx],
                                                        self.electron_conept[self.electronsFakeSel[0].idx]>self.muon_conept[self.muonsFakeSel[1].idx]),
                                                  op.static_cast("std::size_t", op.c_int(1)),
                                                  op.static_cast("std::size_t", op.c_int(0)))),
                                op.switch(self.electron_conept[self.electronsFakeSel[0].idx]>self.muon_conept[self.muonsFakeSel[0].idx],
                                    # At least 2 electrons and 2 muons : If ele1 > mu1 : if ele2>mu1 : leading is ele+ele (slice at 2)
                                    #                                                    if ele2<mu1 : leading is ele+mu (slice at 1)
                                    #                                    If ele1 < mu1 : if ele1>mu2 : leading is ele+mu (slice at 1)
                                    #                                                    if ele1<mu2 : leading is mu+mu (slice at 0)
                                          op.switch(self.electron_conept[self.electronsFakeSel[1].idx]>self.muon_conept[self.muonsFakeSel[0].idx],
                                                    op.static_cast("std::size_t", op.c_int(2)),
                                                    op.static_cast("std::size_t", op.c_int(1))),
                                          op.switch(self.electron_conept[self.electronsFakeSel[0].idx]>self.muon_conept[self.muonsFakeSel[1].idx],
                                                    op.static_cast("std::size_t", op.c_int(1)),
                                                    op.static_cast("std::size_t", op.c_int(0))))
                                                                                )]

            self.leadMuonsFakeSel = self.muonsFakeSel[:op.multiSwitch(  
                                # Same logic as for electrons 
                                (op.rng_len(self.electronsFakeSel)==0 , op.min(op.static_cast("std::size_t", op.c_int(2)), op.rng_len(self.muonsFakeSel))),
                                (op.rng_len(self.muonsFakeSel)==0 , op.static_cast("std::size_t", op.c_int(0))),
                                (op.AND(op.rng_len(self.muonsFakeSel)==1 , op.rng_len(self.electronsFakeSel)==1) , op.static_cast("std::size_t", op.c_int(1))),
                                (op.AND(op.rng_len(self.muonsFakeSel)>=2 , op.rng_len(self.electronsFakeSel)==1), 
                                        op.switch(op.AND(self.muon_conept[self.muonsFakeSel[0].idx]>self.electron_conept[self.electronsFakeSel[0].idx],
                                                         self.muon_conept[self.muonsFakeSel[1].idx]>self.electron_conept[self.electronsFakeSel[0].idx]),
                                                  op.static_cast("std::size_t", op.c_int(2)),
                                                  op.static_cast("std::size_t", op.c_int(1)))),
                                (op.AND(op.rng_len(self.muonsFakeSel)==1 , op.rng_len(self.electronsFakeSel)>=2),
                                        op.switch(op.OR(self.muon_conept[self.muonsFakeSel[0].idx]>self.electron_conept[self.electronsFakeSel[0].idx],
                                                        self.muon_conept[self.muonsFakeSel[0].idx]>self.electron_conept[self.electronsFakeSel[1].idx]),
                                                  op.static_cast("std::size_t", op.c_int(1)),
                                                  op.static_cast("std::size_t", op.c_int(0)))),
                                op.switch(self.muon_conept[self.muonsFakeSel[0].idx]>self.electron_conept[self.electronsFakeSel[0].idx],
                                          op.switch(self.muon_conept[self.muonsFakeSel[1].idx]>self.electron_conept[self.electronsFakeSel[0].idx],
                                                    op.static_cast("std::size_t", op.c_int(2)),
                                                    op.static_cast("std::size_t", op.c_int(1))),
                                          op.switch(self.muon_conept[self.muonsFakeSel[0].idx]>self.electron_conept[self.electronsFakeSel[1].idx],
                                                    op.static_cast("std::size_t", op.c_int(1)),
                                                    op.static_cast("std::size_t", op.c_int(0))))
                                                                        )]


            # Fakeable dilepton #
            self.ElElDileptonFakeSel = op.combine(self.leadElectronsFakeSel, N=2)
            self.MuMuDileptonFakeSel = op.combine(self.leadMuonsFakeSel, N=2)
            self.ElMuDileptonFakeSel = op.combine((self.leadElectronsFakeSel,self.leadMuonsFakeSel))
                # Only use the containers with the leading fakeable leptons 

            # Tight Dilepton : must also be Gen matched if MC #
            lambda_dilepton_matched = lambda dilep : op.AND(self.lambda_is_matched(dilep[0]),self.lambda_is_matched(dilep[1]))
            lambda_tightpair_ElEl = lambda dilep : op.AND(lambda_dilepton_matched(dilep),
                                                          self.lambda_electronTightSel(dilep[0]),
                                                          self.lambda_electronTightSel(dilep[1]))
            lambda_tightpair_MuMu = lambda dilep : op.AND(lambda_dilepton_matched(dilep),
                                                          self.lambda_muonTightSel(dilep[0]),
                                                          self.lambda_muonTightSel(dilep[1]))
            lambda_tightpair_ElMu = lambda dilep : op.AND(lambda_dilepton_matched(dilep),
                                                          self.lambda_electronTightSel(dilep[0]),
                                                          self.lambda_muonTightSel(dilep[1]))

            self.ElElDileptonTightSel = op.select(self.ElElDileptonFakeSel, pred=lambda_tightpair_ElEl)
            self.MuMuDileptonTightSel = op.select(self.MuMuDileptonFakeSel, pred=lambda_tightpair_MuMu)
            self.ElMuDileptonTightSel = op.select(self.ElMuDileptonFakeSel, pred=lambda_tightpair_ElMu)
                # Because self.*DileptonFakeSel has len <= 1 , self.*DileptonTightSel will have length of 0 or 1

            # Fake Extrapolation dilepton #
            lambda_fakepair_ElEl = lambda dilep : op.AND(lambda_dilepton_matched(dilep),
                                                         op.NOT(op.AND(self.lambda_electronTightSel(dilep[0]),
                                                                       self.lambda_electronTightSel(dilep[1]))))
            lambda_fakepair_MuMu = lambda dilep : op.AND(lambda_dilepton_matched(dilep),
                                                         op.NOT(op.AND(self.lambda_muonTightSel(dilep[0]),
                                                                       self.lambda_muonTightSel(dilep[1]))))
            lambda_fakepair_ElMu = lambda dilep : op.AND(lambda_dilepton_matched(dilep),
                                                         op.NOT(op.AND(self.lambda_electronTightSel(dilep[0]),
                                                                       self.lambda_muonTightSel(dilep[1]))))
        
            self.ElElDileptonFakeExtrapolationSel = op.select(self.ElElDileptonFakeSel, pred=lambda_fakepair_ElEl)
            self.MuMuDileptonFakeExtrapolationSel = op.select(self.MuMuDileptonFakeSel, pred=lambda_fakepair_MuMu)
            self.ElMuDileptonFakeExtrapolationSel = op.select(self.ElMuDileptonFakeSel, pred=lambda_fakepair_ElMu)

        #----- Needed only for SL -----#
        if channel == "SL":
            # Same issues as DL 
            # We want self.leadElectronsFakeSel to have rng_len of [0,1] and self.leadMuonssFakeSel have rng_len of [0,1]
            # Such that the two containers have a combined rng_len of at most 1 (0 when no fakeable lepton)
            self.leadElectronFakeSel = self.electronsFakeSel[:op.multiSwitch(
                            (op.rng_len(self.muonsFakeSel)==0 , op.min(op.static_cast("std::size_t", op.c_int(1)), op.rng_len(self.electronsFakeSel))),
                                # No muon -> take the first electron, or zero if empty
                            (op.rng_len(self.electronsFakeSel)==0 , op.static_cast("std::size_t", op.c_int(0))),
                                # No electron -> slice at 0
                            # else : at least one electron and one muon, check the leading in conept
                            op.switch(self.electron_conept[self.electronsFakeSel[0].idx]>self.muon_conept[self.muonsFakeSel[0].idx],
                                      op.static_cast("std::size_t", op.c_int(1)), # electron is leading, take the first
                                      op.static_cast("std::size_t", op.c_int(0))) # muon is leading, take no electron
                                                                        )]
            self.leadMuonFakeSel = self.muonsFakeSel[:op.multiSwitch(
                            (op.rng_len(self.electronsFakeSel)==0 , op.min(op.static_cast("std::size_t", op.c_int(1)), op.rng_len(self.muonsFakeSel))),
                                # No electron -> take the first muon, or zero if empty
                            (op.rng_len(self.muonsFakeSel)==0 , op.static_cast("std::size_t", op.c_int(0))),
                                # No muon -> slice at 0
                            # else : at least one electron and one muon, check the leading in conept
                            op.switch(self.muon_conept[self.muonsFakeSel[0].idx]>self.electron_conept[self.electronsFakeSel[0].idx],
                                      op.static_cast("std::size_t", op.c_int(1)), # muon is leading, take the first
                                      op.static_cast("std::size_t", op.c_int(0))) # electron is leading, take no muon
                                                                        )]

            # Tight selection #
            # Only one the fakeable leading lepton
            lambda_tight_ele = lambda ele : op.AND(self.lambda_is_matched(ele) , self.lambda_electronTightSel(ele))
            lambda_tight_mu  = lambda mu  : op.AND(self.lambda_is_matched(mu)  , self.lambda_muonTightSel(mu))

            self.leadElectronTightSel = op.select(self.leadElectronFakeSel, pred=lambda_tight_ele)
            self.leadMuonTightSel     = op.select(self.leadMuonFakeSel, pred=lambda_tight_mu)

            # Fake Extrapolation selection #

            lambda_fake_ele = lambda ele : op.AND(self.lambda_is_matched(ele) , op.NOT(self.lambda_electronTightSel(ele)))
            lambda_fake_mu  = lambda mu  : op.AND(self.lambda_is_matched(mu)  , op.NOT(self.lambda_muonTightSel(mu)))

            self.leadElectronFakeExtrapolationSel = op.select(self.leadElectronFakeSel, pred=lambda_tight_ele)
            self.leadMuonFakeExtrapolationSel     = op.select(self.leadMuonFakeSel, pred=lambda_tight_mu)
                # NOTE 
                # Given that self.lead*FakeSel has len of 0 or 1, by definition we have :
                #       len(self.lead*TightSel) + len(self.lead*FakeExtrapolationSel) <= 1
                #           (can be 0 for MC where the lepton is not genMatched)
                # SR : len(self.lead*TightSel) == 1
                # Fake CR : len(self.lead*FakeExtrapolationSel) == 1 

        #############################################################################
        #                                AK4 Jets                                   #
        #############################################################################
        self.ak4JetsByPt = op.sort(t.Jet, lambda jet : -jet.pt)
        #self.ak4JetsByPt = op.sort(t.Jet, lambda jet : -jet.btagDeepFlavB)
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
            self.lambda_ak4NoBtag = lambda jet    : jet.btagDeepFlavB <= 0.3033
        elif era == "2018":
            self.lambda_ak4Btag =   lambda jet    : jet.btagDeepFlavB > 0.2770
            self.lambda_ak4NoBtag = lambda jet    : jet.btagDeepFlavB <= 0.2770

        self.ak4BJets     = op.select(self.ak4Jets, self.lambda_ak4Btag)
        self.ak4LightJetsByPt = op.select(self.ak4Jets, self.lambda_ak4NoBtag)
        self.ak4LightJets = op.sort(self.ak4LightJetsByPt, lambda jet : -jet.btagDeepFlavB)
            # Sorted by btag score because for 0 and 1 btag categories, 


        #############################################################################
        #                                AK8 Jets                                   #
        #############################################################################
        self.ak8JetsByPt = op.sort(t.FatJet, lambda jet : -jet.pt)
        #self.ak8JetsByPt = op.sort(t.FatJet, lambda jet : -jet.btagDeepB)
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
        #                               Triggers                                    #
        #############################################################################

        #----- Genweight -----#
        if self.is_MC:
            noSel = noSel.refine("genWeight", weight=t.genWeight)
            #noSel = noSel.refine("genWeight", weight=sign(t.genWeight))

        triggerRanges = returnTriggerRanges(era)
#        #----- Select triggers -----#
#       Still experimental 
#        def returnTriggers(keys):
#            # MC : just OR the different paths 
#            if self.is_MC:
#                return op.OR(*[trig for k in keys for trig in self.triggersPerPrimaryDataset[k]])
#            # Data : due to bug in NanoAOD production, check that the event run number is inside the ranges computed by bricalc
#            else:
#                conversion = {'SingleElectron':'1e','SingleMuon':'1mu','MuonEG':'1e1mu','DoubleEGamma':'2e','DoubleMuon':'2mu'}
#                list_cond = []
#                for trigKey in keys:
#                    rangekey = conversion[trigKey]
#                    rangeDict = {triggerRanges[rangekey][i]['name']:triggerRanges[rangekey][i]['runs'] for i in range(len(triggerRanges[rangekey]))}
#                    listTrig = self.triggersPerPrimaryDataset[trigKey]
#                    trigNames = [trig._parent.name for trig in listTrig]
#                    for trig,trigName in zip(listTrig,trigNames):
#                        trigRanges = rangeDict[trigName]
#                        list_cond.append(op.AND(trig,op.OR(*[op.in_range(r[0],t.run,r[1]) for r in trigRanges])))
#                return op.OR(*list_cond)

        def returnTriggers(keys):
            return op.OR(*[trig for k in keys for trig in self.triggersPerPrimaryDataset[k]])

        # NOTE 
        # For DL the naming convention is self.leadElectronsFakeSel and self.leadMuonsFakeSel (with an "s")
        # For SL it is self.leadElectronFakeSel and self.leadMuonFakeSel (without "s")

        if channel == "DL":
            # Select correct triggers #
            self.triggers = op.multiSwitch((op.rng_len(self.leadElectronsFakeSel)>=2 , returnTriggers(["SingleElectron","DoubleEGamma"])),
                                                  # 2 fakeable electrons as leading
                                            (op.AND(op.rng_len(self.leadElectronsFakeSel)>=1,op.rng_len(self.leadMuonsFakeSel)>=1) , returnTriggers(["SingleMuon","SingleElectron","MuonEG"])),
                                                  # 1 fakeable electron + 1 fakeable muon as leading
                                            (op.rng_len(self.leadMuonsFakeSel)>=2 , returnTriggers(["SingleMuon","DoubleMuon"])),
                                                  # 2 fakeable muons as leading
                                            op.c_bool(False))
                                                  # Number of fakeable leptons < 2 : not selected
        if channel == "SL":
            # Select correct triggers #
            self.triggers = op.multiSwitch((op.rng_len(self.leadElectronFakeSel)==1 , returnTriggers(["SingleElectron"])),
                                           (op.rng_len(self.leadMuonFakeSel)==1 , returnTriggers(["SingleMuon"])),
                                           op.c_bool(False)) # In principle does not happen
                             # Keep in mind either op.rng_len(self.leadElectronFakeSel)==1 or op.rng_len(self.leadMuonFakeSel)==1 but not both
                             # If both ... I made a mistake

        #----- Apply triggers -----#
        if not self.inclusive_sel:
            # Fakeable cut #
            if channel == "DL":
                # Make sure to have at least two fakeable leptons #
                noSel = noSel.refine("hasAtLeast2FakeableLeptons", cut = [op.rng_len(self.electronsFakeSel)+op.rng_len(self.muonsFakeSel) >= 2])
            if channel == "SL":
                # Make sure to have at least one fakeable lepton #
                noSel = noSel.refine("hasAtLeast1FakeableLepton", cut = [op.rng_len(self.electronsFakeSel)+op.rng_len(self.muonsFakeSel) >= 1])
            # Triggers #
            if self.is_MC:
                noSel = noSel.refine("withTrig", cut=self.triggers)
            else:
                if era == "2018": # For 2018 the electron samples are merged, need to build according dict for primary dataset
                    trigger_per_pd = { pd : trig for pd, trig in self.triggersPerPrimaryDataset.items() if pd not in ("SingleElectron", "DoubleEGamma") } 
                    trigger_per_pd["EGamma"] = self.triggersPerPrimaryDataset["SingleElectron"] + self.triggersPerPrimaryDataset["DoubleEGamma"]
                    noSel = noSel.refine("withTrig", cut=[makeMultiPrimaryDatasetTriggerSelection(sample, trigger_per_pd), self.triggers])
                else:
                    noSel = noSel.refine("withTrig", cut=[makeMultiPrimaryDatasetTriggerSelection(sample, self.triggersPerPrimaryDataset), self.triggers])
                # makeMultiPrimaryDatasetTriggerSelection must be done first, to make sure an event is not passed several times from several datasets
                # Then the trigger selection is based on fakeable leptons

        #############################################################################
        #                             Scalefactors                                  #
        #############################################################################
        SF = ScaleFactorsbbWW()    
        if self.is_MC:
            #---- Object SF -----# (Will take as argument the era)
            ####  Muons ####
            self.muLooseId = SF.get_scalefactor("lepton", 'muon_loose_{}'.format(era), combine="weight", systName="mu_loose", defineOnFirstUse=(not forSkimmer))
            self.lambda_MuonLooseSF = lambda mu : [self.muLooseId(mu)]
                # Need to be defined here to be accessible to Skimmer (check on self.is_MC must be applied there)
            self.muTightMVA = SF.get_scalefactor("lepton", 'muon_tightMVA_{}'.format(era), combine="weight", systName="mu_tightmva", defineOnFirstUse=(not forSkimmer))
            self.lambda_MuonTightSF = lambda mu : [self.muTightMVA(mu)] 

            ####  Electrons ####
            if era == "2016" or era == "2017": # Electron reco eff depend on Pt for 2016 and 2017
                self.elLooseRecoPtLt20 = SF.get_scalefactor("lepton", ('electron_loosereco_{}'.format(era) , 'electron_loosereco_ptlt20'), combine="weight", systName="el_looserecoptlt20", defineOnFirstUse=(not forSkimmer))
                self.elLooseRecoPtGt20 = SF.get_scalefactor("lepton", ('electron_loosereco_{}'.format(era) , 'electron_loosereco_ptgt20'), combine="weight", systName="el_looserecoptgt20", defineOnFirstUse=(not forSkimmer))
                # /!\ In analysis YAML file 2016 and 2017 for systematics : must use el_looserecoptlt20 and el_looserecoptgt20

            elif era == "2018": # Does not depend on pt for 2018
                self.elLooseReco = SF.get_scalefactor("lepton", 'electron_loosereco_{}'.format(era), combine="weight", systName="el_loosereco", defineOnFirstUse=(not forSkimmer))
                # /!\ In analysis YAML file 2018 for systematics : must use el_loosereco

            self.elLooseId = SF.get_scalefactor("lepton", 'electron_looseid_{}'.format(era) , combine="weight", systName="el_looseid", defineOnFirstUse=(not forSkimmer))
            self.elLooseEff = SF.get_scalefactor("lepton", 'electron_looseeff_{}'.format(era) , combine="weight", systName="el_looseeff", defineOnFirstUse=(not forSkimmer))

            if era == "2016" or era == "2017":
                # Need to be defined here to be accessible to Skimmer (check on self.is_MC must be applied there)
                self.lambda_ElectronLooseSF = lambda el : [self.elLooseId(el) , self.elLooseEff(el), op.switch(el.pt>20 , self.elLooseRecoPtGt20(el) , self.elLooseRecoPtLt20(el))]
            elif era == "2018": 
                # Need to be defined here to be accessible to Skimmer (check on self.is_MC must be applied there)
                self.lambda_ElectronLooseSF = lambda el : [self.elLooseId(el) , self.elLooseEff(el), self.elLooseReco(el)] 

            self.elTightMVA = SF.get_scalefactor("lepton", 'electron_tightMVA_{}'.format(era) , combine="weight", systName="el_tightmva", defineOnFirstUse=(not forSkimmer))
            self.lambda_ElectronTightSF = lambda el : [self.elTightMVA(el)] 
 
            #### Ak4 Btagging ####
            #DeepJetTag_discriVar = {"BTagDiscri": lambda j : j.btagDeepFlavB}
            #self.DeepJetMediumSF = SF.get_scalefactor("jet", ("btag_"+era, "DeepJet_medium"), additionalVariables=DeepJetTag_discriVar, systName="deepjet", defineOnFirstUse=(not forSkimmer)) # For RESOLVED

            ## Split by flavour #
            #self.DeepJetMediumSF_lightjets = SF.get_scalefactor("lepton", ("btag_lightjets_"+era, "DeepJet_medium"), additionalVariables=DeepJetTag_discriVar, systName="deepjet_lightjets", defineOnFirstUse=(not forSkimmer))
            #self.DeepJetMediumSF_cjets = SF.get_scalefactor("lepton", ("btag_cjets_"+era, "DeepJet_medium"), additionalVariables=DeepJetTag_discriVar, systName="deepjet_cjets", defineOnFirstUse=(not forSkimmer))
            #self.DeepJetMediumSF_bjets = SF.get_scalefactor("lepton", ("btag_bjets_"+era, "DeepJet_medium"), additionalVariables=DeepJetTag_discriVar, systName="deepjet_bjets", defineOnFirstUse=(not forSkimmer))

            #### Ak8 Btagging ####
            #DeepCSVTag_discriVar = {"BTagDiscri": lambda j : j.btagDeepB}
            #self.DeepCSVMediumSF = SF.get_scalefactor("jet", ("subjet_btag_"+era, "DeepCSV_medium"), additionalVariables=DeepCSVTag_discriVar, systName="deepcsv", defineOnFirstUse=(not forSkimmer)) # For BOOSTED (btag on subjet)
                # TODO : add Boosted SF (need NanoAODv7)

            #----- Triggers -----# 
            
            #### Single lepton triggers ####
            self.ttH_singleElectron_trigSF = SF.get_scalefactor("lepton", 'singleTrigger_electron_{}'.format(era) , combine="weight", systName="ttH_singleElectron_trigSF", defineOnFirstUse=(not forSkimmer))
            self.ttH_singleMuon_trigSF = SF.get_scalefactor("lepton", 'singleTrigger_muon_{}'.format(era) , combine="weight", systName="ttH_singleMuon_trigSF", defineOnFirstUse=(not forSkimmer))
            
            #### Double lepton triggers #### (Need to split according to era) 
                # https://gitlab.cern.ch/ttH_leptons/doc/-/blob/master/Legacy/data_to_mc_corrections.md#trigger-efficiency-scale-factors -> deprecated
                # New ref (more up to date) : https://cernbox.cern.ch/index.php/s/lW2BiTli5tJR0MN 
                # Same but clearer  : http://cms.cern.ch/iCMS/jsp/iCMS.jsp?mode=single&part=publications (AN-2019/111) at Table 30, page 43
                # -> Based on subleading conept lepton (except mumu 2018), use lambda 
                # Lambdas return list to be easily concatenated with other SF in lists as well (the weight argument of op.refine requires a list)
            if era == "2016":
                # Double Electron #
                self.lambda_ttH_doubleElectron_trigSF = lambda dilep : op.multiSwitch(
                    (self.electron_conept[dilep[1].idx]<25 , op.systematic(op.c_float(0.98), name="ttH_doubleElectron_trigSF", up=op.c_float(0.98*1.02), down=op.c_float(0.98*0.98))), 
                    op.systematic(op.c_float(1.), name="ttH_doubleElectron_trigSF", up=op.c_float(1.02), down=op.c_float(0.98)))
                # Double Muon #
                self.lambda_ttH_doubleMuon_trigSF = lambda dilep : op.systematic(op.c_float(0.99), name="ttH_doubleMuon_trigSF", up=op.c_float(0.99*1.01), down=op.c_float(0.99*0.99))
                # Electron Muon #
                self.lambda_ttH_electronMuon_trigSF = lambda dilep : op.systematic(op.c_float(1.00), name="ttH_electronMuon_trigSF", up=op.c_float(1.01), down=op.c_float(0.99))
            elif era == "2017":
                # Double Electron #
                self.lambda_ttH_doubleElectron_trigSF = lambda dilep : op.multiSwitch(
                    (self.electron_conept[dilep[1].idx]<40 , op.systematic(op.c_float(0.98), name="ttH_doubleElectron_trigSF", up=op.c_float(0.98*1.01), down=op.c_float(0.98*0.99))), 
                    op.systematic(op.c_float(1.), name="ttH_doubleElectron_trigSF", up=op.c_float(1.01), down=op.c_float(0.99)))
                # Double Muon #
                self.lambda_ttH_doubleMuon_trigSF = lambda dilep : op.multiSwitch(
                    (self.muon_conept[dilep[1].idx]<40 , op.systematic(op.c_float(0.97), name="ttH_doubleMuon_trigSF", up=op.c_float(0.97*1.02), down=op.c_float(0.97*0.98))),
                    (self.muon_conept[dilep[1].idx]<55 , op.systematic(op.c_float(0.995), name="ttH_doubleMuon_trigSF", up=op.c_float(0.995*1.02), down=op.c_float(0.995*0.98))),
                    (self.muon_conept[dilep[1].idx]<70 , op.systematic(op.c_float(0.96), name="ttH_doubleMuon_trigSF", up=op.c_float(0.96*1.02), down=op.c_float(0.96*0.98))),
                    op.systematic(op.c_float(0.94), name="ttH_doubleMuon_trigSF", up=op.c_float(0.94*1.02), down=op.c_float(0.94*0.98)))
                # Electron Muon : /!\ While ElEl and MuMu is conept ordered, ElMu is not #
                self.lambda_ttH_electronMuon_trigSF = lambda dilep : op.multiSwitch(
                    (op.min(self.electron_conept[dilep[0].idx],self.muon_conept[dilep[1].idx])<40, op.systematic(op.c_float(0.98), name="ttH_electronMuon_trigSF", up=op.c_float(0.98*1.01), down=op.c_float(0.98*0.99))),
                    op.systematic(op.c_float(0.99), name="ttH_electronMuon_trigSF", up=op.c_float(0.99*1.01), down=op.c_float(0.99*0.99)))
            elif era == "2018":
                # Double Electron #
                self.lambda_ttH_doubleElectron_trigSF = lambda dilep : op.multiSwitch(
                    (self.electron_conept[dilep[1].idx]<25 , op.systematic(op.c_float(0.98), name="ttH_doubleElectron_trigSF", up=op.c_float(0.98*1.01), down=op.c_float(0.98*0.99))), 
                    op.systematic(op.c_float(1.), name="ttH_doubleElectron_trigSF", up=op.c_float(1.01), down=op.c_float(0.99)))
                # Double Muon /!\ only this one using leading conept lepton #
                self.lambda_ttH_doubleMuon_trigSF = lambda dilep : op.multiSwitch(
                    (self.muon_conept[dilep[0].idx]<40 , op.systematic(op.c_float(1.01), name="ttH_doubleMuon_trigSF", up=op.c_float(1.01*1.01), down=op.c_float(1.01*0.99))),
                    (self.muon_conept[dilep[0].idx]<70 , op.systematic(op.c_float(0.995), name="ttH_doubleMuon_trigSF", up=op.c_float(0.995*1.01), down=op.c_float(0.995*0.99))),
                    op.systematic(op.c_float(0.98), name="ttH_doubleMuon_trigSF", up=op.c_float(0.98*1.01), down=op.c_float(0.98*0.99)))
                # Electron Muon : /!\ While ElEl and MuMu is conept ordered, ElMu is not #
                self.lambda_ttH_electronMuon_trigSF = lambda dilep : op.multiSwitch(
                    (op.min(self.electron_conept[dilep[0].idx],self.muon_conept[dilep[1].idx])<25, op.systematic(op.c_float(0.98), name="ttH_electronMuon_trigSF", up=op.c_float(0.98*1.01), down=op.c_float(0.98*0.99))),
                    op.systematic(op.c_float(1.00), name="ttH_electronMuon_trigSF", up=op.c_float(1.01), down=op.c_float(0.99)))
        #----- DY reweighting -----#
        if era == "2016":
            self.DYReweighting1bElEl = SF.get_scalefactor("lepton", ('DY_{}'.format(era),'ElEl_data_1b'), combine="weight", systName="dy_reweighting_1b_elel", defineOnFirstUse=(not forSkimmer))
            self.DYReweighting1bMuMu = SF.get_scalefactor("lepton", ('DY_{}'.format(era),'MuMu_data_1b'), combine="weight", systName="dy_reweighting_1b_mumu", defineOnFirstUse=(not forSkimmer))
            self.DYReweighting2bElEl = SF.get_scalefactor("lepton", ('DY_{}'.format(era),'ElEl_data_2b'), combine="weight", systName="dy_reweighting_2b_elel", defineOnFirstUse=(not forSkimmer))
            self.DYReweighting2bMuMu = SF.get_scalefactor("lepton", ('DY_{}'.format(era),'MuMu_data_2b'), combine="weight", systName="dy_reweighting_2b_mumu", defineOnFirstUse=(not forSkimmer))
        else:
            self.DYReweighting1bElEl = lambda dilep : None
            self.DYReweighting1bMuMu = lambda dilep : None
            self.DYReweighting2bElEl = lambda dilep : None
            self.DYReweighting2bMuMu = lambda dilep : None

        #----- Fake rates -----#
        FRSysts = ['pt_syst','barrel_syst','norm_syst']
        self.electronFRList = [SF.get_scalefactor("lepton", ('electron_fakerates_'+era, syst), combine="weight", systName="el_FR_"+syst, defineOnFirstUse=(not forSkimmer),
                                             additionalVariables={'Pt' : lambda obj : self.electron_conept[obj.idx]}) for syst in FRSysts]
        self.muonFRList = [SF.get_scalefactor("lepton", ('muon_fakerates_'+era, syst), combine="weight", systName="mu_FR_"+syst, defineOnFirstUse=(not forSkimmer),
                                         additionalVariables={'Pt' : lambda obj : self.muon_conept[obj.idx]}) for syst in FRSysts ] 

        def returnFFSF(obj,list_SF,systName):
            """ Helper when several systematics are present  """
            args = [ a(obj) for a in list_SF[0]._args ] ## get the things the SF depends on
            systArgs = {'nominal':list_SF[0].sfOp.get(*(args+[op.extVar("int", "Nominal")])),'name':systName}
            systArgs.update({SF._systName+'up':SF.sfOp.get(*(args+[op.extVar("int", "Up")])) for SF in list_SF})
            systArgs.update({SF._systName+'down':SF.sfOp.get(*(args+[op.extVar("int", "Down")])) for SF in list_SF})
            return op.systematic(**systArgs)
        
        #############################################################################
        #                             Fake Factors                                  #
        #############################################################################
        # Non closure between ttbar and QCD correction #
        # https://gitlab.cern.ch/cms-hh-bbww/cms-hh-to-bbww/-/blob/master/Legacy/backgroundEstimation.md
        if era == '2016':
            self.electronCorrFR = op.systematic(op.c_float(1.376), name="electronCorrFR",up=op.c_float(1.376*1.376),down=op.c_float(1.))
            self.muonCorrFR     = op.systematic(op.c_float(1.050), name="muonCorrFR",up=op.c_float(1.050*1.050),down=op.c_float(1.))
        elif era == '2017':
            self.electronCorrFR = op.systematic(op.c_float(1.252), name="electronCorrFR",up=op.c_float(1.252*1.252),down=op.c_float(1.))
            self.muonCorrFR     = op.systematic(op.c_float(1.157), name="muonCorrFR",up=op.c_float(1.157*1.157),down=op.c_float(1.))
        elif era == '2018':
            self.electronCorrFR = op.systematic(op.c_float(1.325), name="electronCorrFR",up=op.c_float(1.325*1.325),down=op.c_float(1.))
            self.muonCorrFR     = op.systematic(op.c_float(1.067), name="muonCorrFR",up=op.c_float(1.067*1.067),down=op.c_float(1.))

        self.lambda_FF_el = lambda el : self.electronCorrFR*returnFFSF(el,self.electronFRList,"el_FR")/(1-self.electronCorrFR*returnFFSF(el,self.electronFRList,"el_FR"))
        self.lambda_FF_mu = lambda mu : self.muonCorrFR*returnFFSF(mu,self.muonFRList,"mu_FR")/(1-self.muonCorrFR*returnFFSF(mu,self.muonFRList,"mu_FR"))

        if channel == "DL":
            self.ElElFakeFactor = lambda dilep : op.multiSwitch((op.AND(op.NOT(self.lambda_electronTightSel(dilep[0])),op.NOT(self.lambda_electronTightSel(dilep[1]))),
                                                                 # Both electrons fail tight -> -F1*F2
                                                                 -self.lambda_FF_el(dilep[0])*self.lambda_FF_el(dilep[1])),
                                                                (op.AND(op.NOT(self.lambda_electronTightSel(dilep[0])),self.lambda_electronTightSel(dilep[1])),
                                                                 # Only leading electron fails tight -> F1
                                                                 self.lambda_FF_el(dilep[0])),
                                                                (op.AND(self.lambda_electronTightSel(dilep[0]),op.NOT(self.lambda_electronTightSel(dilep[1]))),
                                                                 # Only subleading electron fails tight -> F2
                                                                 self.lambda_FF_el(dilep[1])),
                                                                 op.c_float(1.)) # Should not happen
            self.MuMuFakeFactor = lambda dilep : op.multiSwitch((op.AND(op.NOT(self.lambda_muonTightSel(dilep[0])),op.NOT(self.lambda_muonTightSel(dilep[1]))),
                                                                 # Both muons fail tight -> -F1*F2
                                                                 -self.lambda_FF_mu(dilep[0])*self.lambda_FF_mu(dilep[1])),
                                                                (op.AND(op.NOT(self.lambda_muonTightSel(dilep[0])),self.lambda_muonTightSel(dilep[1])),
                                                                 # Only leading muon fails tight -> F1
                                                                 self.lambda_FF_mu(dilep[0])),
                                                                (op.AND(self.lambda_muonTightSel(dilep[0]),op.NOT(self.lambda_muonTightSel(dilep[1]))),
                                                                 # Only subleading muon fails tight -> F2
                                                                 self.lambda_FF_mu(dilep[1])),
                                                                 op.c_float(1.)) # Should not happen
            self.ElMuFakeFactor = lambda dilep : op.multiSwitch((op.AND(op.NOT(self.lambda_electronTightSel(dilep[0])),op.NOT(self.lambda_muonTightSel(dilep[1]))),
                                                                 # Both electron and muon fail tight -> -F1*F2
                                                                 -self.lambda_FF_el(dilep[0])*self.lambda_FF_mu(dilep[1])),
                                                                (op.AND(op.NOT(self.lambda_electronTightSel(dilep[0])),self.lambda_muonTightSel(dilep[1])),
                                                                 # Only electron fails tight -> F1
                                                                 self.lambda_FF_el(dilep[0])),
                                                                (op.AND(self.lambda_electronTightSel(dilep[0]),op.NOT(self.lambda_muonTightSel(dilep[1]))),
                                                                 # Only subleading electron fails tight -> F2
                                                                 self.lambda_FF_mu(dilep[1])),
                                                                 op.c_float(1.)) # Should not happen


        ###########################################################################
        #                    b-tagging efficiency scale factors                   #
        ###########################################################################
        if self.is_MC and not self.args.DYStitching:
            #----- Method 1.d -----#
            # See https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagShapeCalibration
            # W_btag = Π_i(all jets) SD(jet_i)  which must be multiplied by r = Σ w(before)/Σ w(after) (before/after using the btag weight, no btag selection for both)
                
            if era == '2016':
                if ('db' in sampleCfg.keys() and 'TuneCP5' in sampleCfg['db']) or ('files' in sampleCfg.keys() and all(['TuneCP5' in f for f in sampleCfg['files']])):
                    csvFileName = os.path.join(os.path.abspath(os.path.dirname(__file__)), "data", "ScaleFactors_POG" , "DeepJet_2016LegacySF_V1_TuneCP5.csv")
                else: # With CUETP8M1
                    csvFileName = os.path.join(os.path.abspath(os.path.dirname(__file__)), "data", "ScaleFactors_POG" , "DeepJet_2016LegacySF_V1.csv")
            if era == '2017':
                csvFileName = os.path.join(os.path.abspath(os.path.dirname(__file__)), "data", "ScaleFactors_POG" , "DeepFlavour_94XSF_V4_B_F.csv")
            if era == '2018':
                csvFileName = os.path.join(os.path.abspath(os.path.dirname(__file__)), "data", "ScaleFactors_POG" , "DeepJet_102XSF_V2.csv")
                

            print ('Btag CSV file',csvFileName)
    
            systTypes = ["jes", "lf", "hf", "hfstats1", "hfstats2","lfstats1", "lfstats2", "cferr1", "cferr2"]
                # From https://twiki.cern.ch/twiki/bin/view/CMS/BTagShapeCalibration#Systematic_uncertainties
            self.DeepJetDiscReshapingSF = BtagSF(taggerName       = "deepjet", 
                                                 csvFileName      = csvFileName,
                                                 wp               = "reshaping",  # "loose", "medium", "tight" or "reshaping"
                                                 sysType          = "central", 
                                                 otherSysTypes    = ['up_'+s for s in systTypes]+['down_'+s for s in systTypes], 
                                                 measurementType  = "iterativefit", 
                                                 sel              = noSel, 
                                                 getters          = {'Discri':lambda j : j.btagDeepFlavB},
                                                 uName            = sample)

            self.btagSF = op.rng_product(self.ak4Jets , lambda j : self.DeepJetDiscReshapingSF(j))
            if self.args.BtagReweightingOn and self.args.BtagReweightingOff:
                raise RuntimeError("Reweighting cannot be both on and off") 
            if self.args.BtagReweightingOn: # Do not apply the ratio
                noSel = noSel.refine("BtagSF" , weight = self.btagSF)
            elif self.args.BtagReweightingOff:
                pass # Do not apply any SF
            else:
                #ReweightingFileName = os.path.join(os.path.dirname(os.path.abspath(__file__)),'data','ScaleFactors_Btag','BtagReweightingRatio_jetN_{}_{}.json'.format(sample,era))
                #if not os.path.exists(ReweightingFileName):
                #    raise RuntimeError("Could not find reweighting file %s"%ReweightingFileName)
                #print ('Reweighting file',ReweightingFileName)

                #self.BtagRatioWeight = makeBtagRatioReweighting(jsonFile = ReweightingFileName,
                #                                                numJets  = op.rng_len(self.ak4Jets),
                #                                                systName = 'btag_ratio',
                #                                                nameHint = f"bamboo_nJetsWeight{sample}".replace('-','_'))
                #noSel = noSel.refine("BtagSF" , weight = [self.btagSF,self.BtagRatioWeight])
                noSel = noSel.refine("BtagSF" , weight = self.btagSF)

        # Return #
        return noSel
