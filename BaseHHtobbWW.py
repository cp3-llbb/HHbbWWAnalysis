import os
import sys
from itertools import chain

from bamboo import treefunctions as op
from bamboo.analysismodules import NanoAODHistoModule
from bamboo.analysisutils import makeMultiPrimaryDatasetTriggerSelection
from bamboo.scalefactors import binningVariables_nano

from scalefactorsbbWW import ScaleFactorsbbWW
from METScripts import METFilter, METcorrection

#===============================================================================================#
#                                  BaseHHtobbWW                                                 #
#===============================================================================================#
class BaseNanoHHtobbWW(NanoAODHistoModule):
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

    def prepareTree(self, tree, sample=None, sampleCfg=None, enableSystematics=None):
        # JEC's Recommendation for Full RunII: https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC
        # JER : -----------------------------: https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
        # Get base aguments #
        era = sampleCfg['era']
        isMC = self.isMC(sample)
        metName = "METFixEE2017" if era == "2017" else "MET"
        tree,noSel,be,lumiArgs = super(BaseNanoHHtobbWW,self).prepareTree(tree, sample=sample, sampleCfg=sampleCfg, calcToAdd=["nJet", metName, "nMuon"])
        triggersPerPrimaryDataset = {}
        from bamboo.analysisutils import configureJets ,configureRochesterCorrection

        ## Check distributed option #
        isNotWorker = (self.args.distributed != "worker") 

        # Rochester and JEC corrections (depends on era) #     

        ############################################################################################
        # ERA 2016 #
        ############################################################################################
        if era == "2016":
            # Rochester corrections #
            configureRochesterCorrection(variProxy  = tree._Muon,
                                         paramsFile = os.path.join(os.path.dirname(__file__), "data", "RoccoR2016.txt"),
                                         isMC       = self.isMC(sample),
                                         backend    = be, 
                                         uName      = sample)

            # Trigger efficiencies #
            # tree.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_(DZ) 
            # tree.HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_(DZ) 
            # -> can be present with, without and without both the DZ
            if self.isMC(sample) or "2016F" in sample or "2016G" in sample:# or "2016H" in sample:
                # Found in 2016F : both
                # Found in 2016G : both
                triggersPerPrimaryDataset = {
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
                triggersPerPrimaryDataset = {
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

                triggersPerPrimaryDataset = {
                    "SingleMuon" :  [ tree.HLT.IsoMu24],
                    "SingleElectron":  [ tree.HLT.Ele27_WPTight_Gsf],
                    "DoubleMuon" :  [ tree.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL,
                                      tree.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ],
                    "DoubleEGamma": [ tree.HLT.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ],
                    "MuonEG":       [ tree.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL, 
                                      tree.HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL]} 

            # Jet treatment #
            cachJEC_dir = '/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/cacheJEC'
            if self.isMC(sample):   # if MC -> needs smearing
                configureJets(variProxy             = tree._Jet, 
                              jetType               = "AK4PFchs",
                              jec                   = "Summer16_07Aug2017_V20_MC",
                              smear                 = "Summer16_25nsV1_MC",
                              jesUncertaintySources = ["Total"],
                              mayWriteCache         = isNotWorker,
                              isMC                  = self.isMC(sample),
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
                              isMC                  = self.isMC(sample),
                              backend               = be, 
                              uName                 = sample,
                              cachedir              = cachJEC_dir)

        ############################################################################################
        # ERA 2017 #
        ############################################################################################
        #elif era == "2017":
        #    configureRochesterCorrection(tree._Muon.calc,os.path.join(os.path.dirname(__file__), "data", "RoccoR2017.txt"))
        #    triggersPerPrimaryDataset = {
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
        #    if self.isMC(sample):
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
        #elif era == "2018":
        #    configureRochesterCorrection(tree._Muon.calc,os.path.join(os.path.dirname(__file__), "data", "RoccoR2018.txt"))
        #    triggersPerPrimaryDataset = {
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
        #    if self.isMC(sample):
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
        else:
            raise RuntimeError("Unknown era {0}".format(era))

        # Get weights #
        if self.isMC(sample):
            noSel = noSel.refine("genWeight", weight=tree.genWeight, cut=op.OR(*chain.from_iterable(triggersPerPrimaryDataset.values())))
            #noSel = noSel.refine("genWeight", weight=op.abs(tree.genWeight), cut=op.OR(*chain.from_iterable(triggersPerPrimaryDataset.values())))
            #noSel = noSel.refine("negWeight", cut=[tree.genWeight<0])
        else:
            noSel = noSel.refine("withTrig", cut=makeMultiPrimaryDatasetTriggerSelection(sample, triggersPerPrimaryDataset))

        return tree,noSel,be,lumiArgs

    def preparePlots(self, t, noSel, sample, sampleCfg):
        # Some imports #
        from bamboo.analysisutils import forceDefine

        era = sampleCfg['era']

        isMC = self.isMC(sample)

        # Forcedefne #
        forceDefine(t._Muon.calcProd, noSel)
        forceDefine(t._Jet.calcProd, noSel) # calculate once per event (for every event)

        # Initialize scalefactors class #
        SF = ScaleFactorsbbWW()

        ###########################################################################
        #                           TTbar reweighting                             #
        ###########################################################################
        if self.isMC(sample) and sample.startswith("TT"):
            # https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting#Use_case_3_ttbar_MC_is_used_to_m
            # Get tops #
            genTop_all = op.select(t.GenPart,lambda g : g.pdgId==6)
            genTop = op.select(genTop_all,lambda g : g.statusFlags & ( 0x1 << 13))
            genAntitop_all = op.select(t.GenPart,lambda g : g.pdgId==-6)
            genAntitop = op.select(genAntitop_all,lambda g : g.statusFlags & ( 0x1 << 13))
                # statusFlags==13 : isLastCopy
                # Pdgid == 6 : top
            hasttbar = noSel.refine("hasttbar",cut=[op.rng_len(genTop)>=1,op.rng_len(genAntitop)>=1])

            # Check plots #
            #plots.append(Plot.make1D("N_tops",
            #                        op.rng_len(genTop),
            #                        noSel,
            #                        EquidistantBinning(5,0.,5.),
            #                        title='N tops',
            #                        xTitle='N tops'))
            #plots.append(Plot.make1D("N_antitops",
            #                        op.rng_len(genAntitop),
            #                        noSel,
            #                        EquidistantBinning(5,0.,5.),
            #                        title='N antitops',
            #                        xTitle='N antitops'))
            ## Top pt #
            #plots.append(Plot.make1D("top_pt",
            #                        genTop[0].pt,
            #                        hasttbar,
            #                        EquidistantBinning(50,0,500),
            #                        xTitle='P_{T} top'))
            ## Antitop Pt #
            #plots.append(Plot.make1D("antitop_pt",
            #                        genAntitop[0].pt,
            #                        hasttbar,
            #                        EquidistantBinning(50,0,500),
            #                        xTitle='P_{T} lead antitop'))
            ## Top - Antitop plots #
            #plots.append(Plot.make2D("top_pt_vs_antitop_pt",
            #                        [genTop[0].pt,genAntitop[0].pt],
            #                        hasttbar,
            #                        [EquidistantBinning(50,0,500),EquidistantBinning(50,0,500)],
            #                        xTitle='P_{T} lead top',
            #                        yTitle='P_{T} lead antitop'))
 

            # Compute weight if there is a ttbar #
            ttbar_SF = lambda t : op.exp(0.0615-0.0005*t.pt)
            ttbar_weight = lambda t,tbar : op.sqrt(ttbar_SF(t)*ttbar_SF(tbar))
            #plots.append(Plot.make1D("ttbar_weight",
            #                        ttbar_weight(genTop[0],genAntitop[0]),
            #                        hasttbar,
            #                        EquidistantBinning(100,0.,2.),
            #                        title='ttbar weight',
            #                        xTitle='ttbar weight'))
            #plots.append(Plot.make3D("ttbar_weight_vs_pt",
            #                        [genTop[0].pt,genAntitop[0].pt,ttbar_weight(genTop[0],genAntitop[0])],
            #                        hasttbar,
            #                        [EquidistantBinning(50,0,500),EquidistantBinning(50,0,500),EquidistantBinning(100,0.,2.)],
            #                        title='ttbar weight',
            #                        xTitle='top P_{T}',
            #                        yTitle='antitop P_{T}',
            #                        zTitle='ttbar weight'))
            # Apply correction to TT #
            noSel = noSel.refine("ttbarWeight",weight=ttbar_weight(genTop[0],genAntitop[0]))

        #############################################################################
        #                             Pile-up                                       #
        #############################################################################
        puWeightsFile = None
        if era == "2016":
            sfTag="94X"
            puWeightsFile = os.path.join(os.path.dirname(__file__), "data", "puweights2016.json")
        elif era == "2017":
            sfTag="94X"     
            puWeightsFile = os.path.join(os.path.dirname(__file__), "data", "puweights2017.json")
        elif era == "2018":
            sfTag="102X"
            puWeightsFile = os.path.join(os.path.dirname(__file__), "data", "puweights2018.json")
        if self.isMC(sample) and puWeightsFile is not None:
            from bamboo.analysisutils import makePileupWeight
            noSel = noSel.refine("puWeight", weight=makePileupWeight(puWeightsFile, t.Pileup_nTrueInt, systName="pileup"))

        #############################################################################
        #                                 MET                                       #
        #############################################################################
        # MET filter #
        noSel = noSel.refine("passMETFlags", cut=METFilter(t.Flag, era, isMC) )

        # MET corrections #
        MET = t.MET if era != "2017" else t.METFixEE2017
        self.corrMET = METcorrection(MET,t.PV,sample,era,self.isMC(sample))

        #############################################################################
        #                                 Muons                                     #
        #############################################################################
        # Wp // 2016- 2017 -2018   // https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Muon_Isolation
        muonsByPt = op.sort(t.Muon, lambda mu : -mu.p4.Pt())
        #muons = op.select(muonsByPt, lambda mu : op.AND(mu.p4.Pt() > 15., op.abs(mu.p4.Eta()) < 2.4, mu.mediumId, mu.pfRelIso04_all<0.15)) # MEDIUM
        self.muons = op.select(muonsByPt, lambda mu : op.AND(mu.p4.Pt() > 15., op.abs(mu.p4.Eta()) < 2.4, mu.tightId, mu.pfRelIso04_all<0.15)) # TIGHT
            # Subleading lepton pt cut is at 15 GeV so better start at that point
            # isolation : tight (Muon::PFIsoTight) cut value = 0.15 (ε~0.95)
      
        # Scalefactors #
        if self.isMC(sample):
            muTightIDSF = SF.get_scalefactor("lepton", ("muon_{0}_{1}".format(era, sfTag), "id_tight"), combine="weight", systName="muid")
            muTightISOSF = SF.get_scalefactor("lepton", ("muon_{0}_{1}".format(era, sfTag), "iso_tight_id_medium"), combine="weight", systName="muiso")
#            if era=="2016":
#                muTightIDSF = SF.get_scalefactor("lepton", ("muon_{0}_{1}".format(era, sfTag), "id_tight"), combine="weight", systName="muid")
#                muTightISOSF = SF.get_scalefactor("lepton", ("muon_{0}_{1}".format(era, sfTag), "iso_tight_id_medium"), combine="weight", systName="muiso")
#                TrkIDSF = SF.get_scalefactor("lepton", ("muon_{0}_{1}".format(era, sfTag), "idtrk_highpt"), combine="weight")         # Need to ask what it is
#                TrkISOSF = SF.get_scalefactor("lepton", ("muon_{0}_{1}".format(era, sfTag), "isotrk_loose_idtrk_highptidandipcut"), combine="weight") # Need to ask what it is
#
#            else:
#                muTightIDSF = SF.get_scalefactor("lepton", ("muon_{0}_{1}".format(era, sfTag), "id_tight"), systName="muid")
#                muTightISOSF = SF.get_scalefactor("lepton", ("muon_{0}_{1}".format(era, sfTag), "iso_tight_id_medium"), systName="muiso") 

        #############################################################################
        #                              Electrons                                    #
        #############################################################################
        #Wp  // 2016: Electron_cutBased_Sum16==3  -> medium     // 2017 -2018  : Electron_cutBased ==3   --> medium ( Fall17_V2)
        # asking for electrons to be in the Barrel region with dz<1mm & dxy< 0.5mm   //   Endcap region dz<2mm & dxy< 0.5mm 
        electronsByPt = op.sort(t.Electron, lambda ele : -ele.p4.Pt())
        #electrons = op.select(electronsByPt, lambda ele : op.AND(ele.p4.Pt() > 15., op.abs(ele.p4.Eta()) < 2.5 , ele.cutBased>=3 )) # //cut-based ID Fall17 V2 the recommended one from POG for the FullRunII MEDIUM
        self.electrons = op.select(electronsByPt, lambda ele : op.AND(ele.p4.Pt() > 15., op.abs(ele.p4.Eta()) < 2.5 , ele.cutBased>=4 )) # //cut-based ID Fall17 V2 the recommended one from POG for the FullRunII TIGHT 
            # electron cut based ID = 0:fail, 1: veto, 2:loose, 3:medium, 4:tight (From root file -> Events tree)
            # Subleading lepton pt cut is at 15 GeV so better start at that point

        # Scalefactors #
        if self.isMC(sample):
            elTightIDSF = SF.get_scalefactor("lepton", ("electron_{0}_{1}".format(era,sfTag), "id_tight"), systName="elid")

        #############################################################################
        #                               Dilepton                                    #
        #############################################################################
        # Combine dileptons #
        # Get OS leptons #
        lambdaOS  = lambda l1,l2 : l1.charge != l2.charge
        preOsElEl = op.combine(self.electrons, N=2, pred=lambdaOS)
        preOsMuMu = op.combine(self.muons, N=2, pred=lambdaOS)
        preOsElMu = op.combine((self.electrons, self.muons), pred=lambdaOS)
    
        # PT cut on leading > 25 GeV #
        lambdaPTCut = lambda dilep : op.OR( dilep[0].p4.Pt() > 25 , dilep[1].p4.Pt() > 25)
        self.OsElEl = op.select(preOsElEl, pred=lambdaPTCut)
        self.OsMuMu = op.select(preOsMuMu, pred=lambdaPTCut)
        self.OsElMu = op.select(preOsElMu, pred=lambdaPTCut)
 
        # Scalefactors #
        if self.isMC(sample):
            doubleEleTrigSF = SF.get_scalefactor("dilepton", ("doubleEleLeg_HHMoriond17_2016"), systName="eleltrig")     
            doubleMuTrigSF  = SF.get_scalefactor("dilepton", ("doubleMuLeg_HHMoriond17_2016"), systName="mumutrig")    
            elemuTrigSF     = SF.get_scalefactor("dilepton", ("elemuLeg_HHMoriond17_2016"), systName="elmutrig")
            mueleTrigSF     = SF.get_scalefactor("dilepton", ("mueleLeg_HHMoriond17_2016"), systName="mueltrig")

            # From https://gitlab.cern.ch/ttH_leptons/doc/blob/master/Legacy/data_to_mc_corrections.md#trigger-efficiency-scale-factors
        ttH_doubleMuon_trigSF = op.systematic(op.c_float(1.010), name="ttH_doubleMuon_trigSF", up=op.c_float(1.020), down=op.c_float(1.000))
        ttH_doubleElectron_trigSF = op.systematic(op.c_float(1.020), name="ttH_doubleElectron_trigSF", up=op.c_float(1.040), down=op.c_float(1.000))
        ttH_electronMuon_trigSF = op.systematic(op.c_float(1.020), name="ttH_electronMuon_trigSF", up=op.c_float(1.030), down=op.c_float(1.010))

        llSF =  {
            "ElEl" : (lambda ll : [ elTightIDSF(ll[0][0]),                                                                     # First lepton SF
                                    elTightIDSF(ll[0][1]),                                                                     # Second lepton SF
                                    #doubleEleTrigSF(ll[0]),                                                                  # Dilepton SF
                                    ttH_doubleElectron_trigSF,
                                  ]),
            "MuMu" : (lambda ll : [ muTightIDSF(ll[0][0]), muTightISOSF(ll[0][0]),                                             # First lepton SF
                                    muTightIDSF(ll[0][1]), muTightISOSF(ll[0][1]),                                             # Second lepton SF
                                    #doubleMuTrigSF(ll[0]),                                                                   # Dilepton SF
                                    ttH_doubleMuon_trigSF,
                                  ]),
            "ElMu" : (lambda ll : [ elTightIDSF(ll[0][0]),                                                                     # First lepton SF
                                    muTightIDSF(ll[0][1]), muTightISOSF(ll[0][1]),                                             # Second lepton SF
                                    #elemuTrigSF(ll[0]),                                                                      # Dilepton SF
                                    ttH_electronMuon_trigSF,
                                  ]),
                # ll is a proxy list of dileptons 
                # ll[0] is the first dilepton 
                # ll[0][0] is the first lepton and ll[0][1] the second in the dilepton
                }


        self.llSFApplied = {
            "ElEl": llSF["ElEl"](self.OsElEl) if isMC else None,
            "MuMu": llSF["MuMu"](self.OsMuMu) if isMC else None,
            "ElMu": llSF["ElMu"](self.OsElMu) if isMC else None,
                      }

        #############################################################################
        #                                 Jets                                      #
        #############################################################################
        # select jets   // 2016 (medium) - 2017 - 2018 (tight)  ( j.jetId &2) ->      tight jet ID
        jetsByPt = op.sort(t.Jet, lambda jet : -jet.p4.Pt())
        fatjetsByPt = op.sort(t.FatJet, lambda fatjet : -fatjet.p4.Pt())

        if era == "2016":
            self.jetsSel = op.select(jetsByPt, lambda j : op.AND(j.p4.Pt() > 25., op.abs(j.p4.Eta())< 2.4, (j.jetId &1)))        # Jets = AK4 jets
            self.fatjetsSel = op.select(fatjetsByPt, lambda j : op.AND(j.p4.Pt() > 200., op.abs(j.p4.Eta())< 2.4, (j.jetId &1)))        # FatJets = AK8 jets
        elif era == "2017" or era == "2018":
            self.jetsSel = op.select(jetsByPt, lambda j : op.AND(j.p4.Pt() > 25., op.abs(j.p4.Eta())< 2.4, (j.jetId &2)))        # Jets = AK4 jets
            self.fatjetsSel = op.select(fatjetsByPt, lambda j : op.AND(j.p4.Pt() > 200., op.abs(j.p4.Eta())< 2.4, (j.jetId &2)))        # FatJets = AK8 jets

        # exclude from the jetsSel any jet that happens to include within its reconstruction cone a muon or an electron.
        self.jets = op.select(self.jetsSel, lambda j : op.AND(op.NOT(op.rng_any(self.electrons, lambda ele : op.deltaR(j.p4, ele.p4) < 0.3 )), op.NOT(op.rng_any(self.muons, lambda mu : op.deltaR(j.p4, mu.p4) < 0.3 ))))
        self.fatjets = op.select(self.fatjetsSel, lambda j : op.AND(op.NOT(op.rng_any(self.electrons, lambda ele : op.deltaR(j.p4, ele.p4) < 0.3 )), op.NOT(op.rng_any(self.muons, lambda mu : op.deltaR(j.p4, mu.p4) < 0.3 ))))

       
        #############################################################################
        #                              BTagging                                     #
        #############################################################################
        # Boosted and resolved jets categories #
        # Boosted -> at least one AK8 jet (fatjet) with at least one subjet passing medium working point of DeepCSV (btagDeepB branch)
        # Resolved -> at least two Ak4 jets (jet) with at least one passing the medium working point of DeepJet (btagDeepFlavB branch)
        if era == "2016": # Must check that subJet exists before looking at the btag
            lambda_boosted  = lambda fatjet : op.OR(op.AND(fatjet.subJet1._idx.result != -1,fatjet.subJet1.btagDeepB > 0.6321), op.AND(fatjet.subJet2._idx.result != -1,fatjet.subJet2.btagDeepB > 0.6321))
            lambda_resolved = lambda jet    : jet.btagDeepFlavB > 0.3093
            lambda_notResolved = lambda jet : jet.btagDeepFlavB <= 0.3093

        elif era =="2017":
            lambda_boosted  = lambda fatjet : op.OR(op.AND(fatjet.subJet1._idx.result != -1,fatjet.subJet1.btagDeepB > 0.4941), op.AND(fatjet.subJet2._idx.result != -1,fatjet.subJet2.btagDeepB > 0.4941))
            lambda_resolved = lambda jet    : jet.btagDeepFlavB > 0.3033
            lambda_notResolved = lambda jet : jet.btagDeepFlavB <= 0.3033
        elif era == "2018":
            lambda_boosted  = lambda fatjet : op.OR(op.AND(fatjet.subJet1._idx.result != -1,fatjet.subJet1.btagDeepB > 0.4184), op.AND(fatjet.subJet2._idx.result != -1,fatjet.subJet2.btagDeepB > 0.4184))
            lambda_resolved = lambda jet    : jet.btagDeepFlavB > 0.2770
            lambda_notResolved = lambda jet : jet.btagDeepFlavB <= 0.2770

        # Select the bjets we want #
        self.bjetsBoosted        = op.select(self.fatjets, lambda_boosted)
        self.bjetsResolved       = op.select(self.jets, lambda_resolved)
        self.lightjetsResolved   = op.select(self.jets, lambda_notResolved) # To plot the case when 1 bjet + 1 lightjet

        # Scalefactors : 2017 and 2018 not yet present in the dict #
        self.DeepCSVMediumSFApplied = None
        self.DeepJetMediumSFApplied = None
        if self.isMC(sample):
            DeepJetTag_discriVar = {"BTagDiscri": lambda j : j.btagDeepFlavB}
            self.DeepJetMediumSF = SF.get_scalefactor("jet", ("btag_"+era+"_"+sfTag, "DeepJet_medium"), additionalVariables=DeepJetTag_discriVar, systName="deepjet") # For RESOLVED
            DeepJetMediumSFApplied = [self.DeepJetMediumSF(self.bjetsResolved[0])] # TODO : check if more than one bjet and apply to all
#            DeepCSVTag_discriVar = {"BTagDiscri": lambda j : j.btagDeepB}
#            DeepCSVMediumSF = SF.get_scalefactor("jet", ("subjet_btag_"+era+"_"+sfTag, "DeepCSV_medium"), additionalVariables=DeepCSVTag_discriVar, systName="deepcsv") # For BOOSTED (btag on subjet)
#            DeepCSVMediumSFApplied = [DeepCSVMediumSF(bjetsBoosted[0].subJet1)] # Must be applied on subjets : need to check each time which one has been btagged
 
        return noSel
