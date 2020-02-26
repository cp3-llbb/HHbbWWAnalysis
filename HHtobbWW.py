import sys
import os
from itertools import chain

from bamboo.analysismodules import NanoAODHistoModule
from bamboo.analysisutils import makeMultiPrimaryDatasetTriggerSelection
from bamboo.scalefactors import binningVariables_nano

from bamboo import treefunctions as op
from bamboo.plots import Plot, EquidistantBinning, SummedPlot

sys.path.append('/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/') # Add scripts in this directory -- TODO : make cleaner
from plotDef import makeDileptonPlots, makeJetsPlots, makeFatJetPlots, makeMETPlots, makeDeltaRPlots, makeYieldPlot, makeHighLevelQuantities
from scalefactorsbbWW import ScaleFactorsbbWW
from METScripts import METFilter, METcorrection

class NanoHHTobbWW(NanoAODHistoModule):
    """ Example module: HH->bbW(->e/µ nu)W(->e/µ nu) histograms from NanoAOD """
    def __init__(self, args):
            super(NanoHHTobbWW, self).__init__(args)
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
        tree,noSel,be,lumiArgs = super(NanoHHTobbWW,self).prepareTree(tree, sample=sample, sampleCfg=sampleCfg, calcToAdd=["nJet", metName, "nMuon"])
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
        
    def definePlots(self, t, noSel, sample=None, sampleCfg=None): 
        # Some imports #
        from bamboo.analysisutils import forceDefine

        era = sampleCfg['era']

        isMC = self.isMC(sample)
        plots = []

        # Forcedefne #
        forceDefine(t._Muon.calcProd, noSel)
        forceDefine(t._Jet.calcProd, noSel) # calculate once per event (for every event)

        # Initialize scalefactors class #
        SF = ScaleFactorsbbWW()

        ############################################################################
        ########################### TTbar reweighting #############################
        ############################################################################
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
        ##########################    Pile-up    ####################################
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
        ################################# MET #######################################
        #############################################################################
        # MET filter #
        noSel = noSel.refine("passMETFlags", cut=METFilter(t.Flag, era, isMC) )

        # MET corrections #
        MET = t.MET if era != "2017" else t.METFixEE2017
        corrMET = METcorrection(MET,t.PV,sample,era,self.isMC(sample))

        #############################################################################
        ################################  Muons #####################################
        #############################################################################
        # Wp // 2016- 2017 -2018   // https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Muon_Isolation
        muonsByPt = op.sort(t.Muon, lambda mu : -mu.p4.Pt())
        #muons = op.select(muonsByPt, lambda mu : op.AND(mu.p4.Pt() > 15., op.abs(mu.p4.Eta()) < 2.4, mu.mediumId, mu.pfRelIso04_all<0.15)) # MEDIUM
        muons = op.select(muonsByPt, lambda mu : op.AND(mu.p4.Pt() > 15., op.abs(mu.p4.Eta()) < 2.4, mu.tightId, mu.pfRelIso04_all<0.15)) # TIGHT
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
        #############################  Electrons  ###################################
        #############################################################################
        #Wp  // 2016: Electron_cutBased_Sum16==3  -> medium     // 2017 -2018  : Electron_cutBased ==3   --> medium ( Fall17_V2)
        # asking for electrons to be in the Barrel region with dz<1mm & dxy< 0.5mm   //   Endcap region dz<2mm & dxy< 0.5mm 
        electronsByPt = op.sort(t.Electron, lambda ele : -ele.p4.Pt())
        #electrons = op.select(electronsByPt, lambda ele : op.AND(ele.p4.Pt() > 15., op.abs(ele.p4.Eta()) < 2.5 , ele.cutBased>=3 )) # //cut-based ID Fall17 V2 the recommended one from POG for the FullRunII MEDIUM
        electrons = op.select(electronsByPt, lambda ele : op.AND(ele.p4.Pt() > 15., op.abs(ele.p4.Eta()) < 2.5 , ele.cutBased>=4 )) # //cut-based ID Fall17 V2 the recommended one from POG for the FullRunII TIGHT 
            # electron cut based ID = 0:fail, 1: veto, 2:loose, 3:medium, 4:tight (From root file -> Events tree)
            # Subleading lepton pt cut is at 15 GeV so better start at that point

        # Scalefactors #
        if self.isMC(sample):
            elTightIDSF = SF.get_scalefactor("lepton", ("electron_{0}_{1}".format(era,sfTag), "id_tight"), systName="elid")

        #############################################################################
        #############################   Dilepton  ###################################
        #############################################################################
        # Combine dileptons #
        # Get OS leptons #
        lambdaOS  = lambda l1,l2 : l1.charge != l2.charge
        preOsElEl = op.combine(electrons, N=2, pred=lambdaOS)
        preOsMuMu = op.combine(muons, N=2, pred=lambdaOS)
        preOsElMu = op.combine((electrons, muons), pred=lambdaOS)
    
        # PT cut on leading > 25 GeV #
        lambdaPTCut = lambda dilep : op.OR( dilep[0].p4.Pt() > 25 , dilep[1].p4.Pt() > 25)
        OsElEl = op.select(preOsElEl, pred=lambdaPTCut)
        OsMuMu = op.select(preOsMuMu, pred=lambdaPTCut)
        OsElMu = op.select(preOsElMu, pred=lambdaPTCut)


        # Plots Numbers of dilepton in each channel #
        plots.append(Plot.make1D("ElEl_channel",op.rng_len(OsElEl),noSel,EquidistantBinning(5,0,5.),title='Number of dilepton events in ElEl channel',xTitle='N_{dilepton} (ElEl channel)'))
        plots.append(Plot.make1D("MuMu_channel",op.rng_len(OsMuMu),noSel,EquidistantBinning(5,0,5.),title='Number of dilepton events in MuMu channel',xTitle='N_{dilepton} (MuMu channel)'))
        plots.append(Plot.make1D("ElMu_channel",op.rng_len(OsElMu),noSel,EquidistantBinning(5,0,5.),title='Number of dilepton events in ElMu channel',xTitle='N_{dilepton} (ElMu channel)'))

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


        llSFApplied = {
            "ElEl": llSF["ElEl"](OsElEl) if isMC else None,
            "MuMu": llSF["MuMu"](OsMuMu) if isMC else None,
            "ElMu": llSF["ElMu"](OsElMu) if isMC else None,
                      }

        # Selection #
        hasOsElEl = noSel.refine("hasOsElEl",
                                 cut = [op.rng_len(OsElEl) >= 1,                        # Require at least one dilepton ElEl
                                        (op.rng_len(electrons)+op.rng_len(muons))<=2],  # Not more than two tight leptons
                                 weight = llSFApplied["ElEl"])
        hasOsMuMu = noSel.refine("hasOsMuMu",
                                 cut = [op.rng_len(OsMuMu) >= 1,                        # Require at least one dilepton MuMu
                                        (op.rng_len(electrons)+op.rng_len(muons))<=2],  # Not more than two tight leptons
                                 weight = llSFApplied["MuMu"])
        hasOsElMu = noSel.refine("hasOsElMu",
                                 cut = [op.rng_len(OsElMu) >= 1,                        # Require at least one dilepton ElMu
                                        (op.rng_len(electrons)+op.rng_len(muons))<=2],  # Not more than two tight leptons
                                 weight = llSFApplied["ElMu"])

        # Yield plots #
        plots.append(makeYieldPlot(self,hasOsElEl,"OSElEl","OS leptons (channel $e^+e^-$)",0))
        plots.append(makeYieldPlot(self,hasOsMuMu,"OSMuMu","OS leptons (channel : $\mu^+\mu^-$)",1))
        plots.append(makeYieldPlot(self,hasOsElMu,"OSMuEl","OS leptons (channel $e^{\pm}\mu^{\mp}$)",2))

        # Dilepton channel plot #

        hasOsCutChannelList = [
                             {'channel':'ElEl','sel':hasOsElEl,'dilepton':OsElEl[0],'suffix':'hasOsElEl'},
                             {'channel':'MuMu','sel':hasOsMuMu,'dilepton':OsMuMu[0],'suffix':'hasOsMuMu'},
                             {'channel':'ElMu','sel':hasOsElMu,'dilepton':OsElMu[0],'suffix':'hasOsElMu'},
                                       ]
        for channelDict in hasOsCutChannelList:
            # Dilepton plots #
            plots.extend(makeDileptonPlots(self, **channelDict))
            # MET plots #
            plots.extend(makeMETPlots(self      = self,
                                      sel       = channelDict['sel'],
                                      met       = corrMET,
                                      suffix    = channelDict['suffix'],
                                      channel   = channelDict['channel']))

        # Dilepton Z peak exclusion (charge already done in previous selection) # 
        lambda_lowMllCut    = lambda dilep: op.invariant_mass(dilep[0].p4, dilep[1].p4)>12.
        lambda_outZ         = lambda dilep: op.NOT(op.in_range(80.,op.invariant_mass(dilep[0].p4, dilep[1].p4),100.))

        hasOsElElLowMllCutOutZ = hasOsElEl.refine("hasOsElElLowMllCutOutZ",cut=[lambda_lowMllCut(OsElEl[0]),lambda_outZ(OsElEl[0])])
        hasOsMuMuLowMllCutOutZ = hasOsMuMu.refine("hasOsMuMuLowMllCutOutZ",cut=[lambda_lowMllCut(OsMuMu[0]),lambda_outZ(OsMuMu[0])])
        hasOsElMuLowMllCut     = hasOsElMu.refine("hasOsElMuLowMllCut",cut=[lambda_lowMllCut(OsElMu[0])]) # Z peak cut not needed because Opposite Flavour

        # Yield plots #
        plots.append(makeYieldPlot(self,hasOsElElLowMllCutOutZ,"OSElElMllCutOutZ","OS leptons + $M_{ll}$ (channel : $e^+e^-$)",3))
        plots.append(makeYieldPlot(self,hasOsMuMuLowMllCutOutZ,"OSMuMuMllCutOutZ","OS leptons + $M_{ll}$ (channel : $\mu^+\mu^-$)",4))
        plots.append(makeYieldPlot(self,hasOsElMuLowMllCut,"OSMuMuMllCut","OS leptons + $M_{ll}$ (channel : $e^{\pm}\mu^{\mp}$)",5))

        # Dilepton plots #

        hasOsMllCutChannelList = [
                             {'channel':'ElEl','sel':hasOsElElLowMllCutOutZ,'dilepton':OsElEl[0],'suffix':'hasOsElElLowMllCutOutZ'},
                             {'channel':'MuMu','sel':hasOsMuMuLowMllCutOutZ,'dilepton':OsMuMu[0],'suffix':'hasOsMuMuLowMllCutOutZ'},
                             {'channel':'ElMu','sel':hasOsElMuLowMllCut,    'dilepton':OsElMu[0],'suffix':'hasOsElMuLowMllCut'},
                           ]

        for channelDict in hasOsMllCutChannelList:
            # Dilepton plots #
            plots.extend(makeDileptonPlots(self, **channelDict))
            # MET plots #
            plots.extend(makeMETPlots(self      = self,
                                      sel       = channelDict['sel'],
                                      met       = corrMET,
                                      suffix    = channelDict['suffix'],
                                      channel   = channelDict['channel']))

        #############################################################################
        ################################  Jets  #####################################
        #############################################################################
        # select jets   // 2016 (medium) - 2017 - 2018 (tight)  ( j.jetId &2) ->      tight jet ID
        jetsByPt = op.sort(t.Jet, lambda jet : -jet.p4.Pt())
        fatjetsByPt = op.sort(t.FatJet, lambda fatjet : -fatjet.p4.Pt())
        if era == "2016":
            jetsSel = op.select(jetsByPt, lambda j : op.AND(j.p4.Pt() > 25., op.abs(j.p4.Eta())< 2.4, (j.jetId &1)))        # Jets = AK4 jets
            fatjetsSel = op.select(fatjetsByPt, lambda j : op.AND(j.p4.Pt() > 200., op.abs(j.p4.Eta())< 2.4, (j.jetId &1)))        # FatJets = AK8 jets
        elif era == "2017" or era == "2018":
            jetsSel = op.select(jetsByPt, lambda j : op.AND(j.p4.Pt() > 25., op.abs(j.p4.Eta())< 2.4, (j.jetId &2)))        # Jets = AK4 jets
            fatjetsSel = op.select(fatjetsByPt, lambda j : op.AND(j.p4.Pt() > 200., op.abs(j.p4.Eta())< 2.4, (j.jetId &2)))        # FatJets = AK8 jets

        # Plot lepton/jet angle separartion
        DeltaRChannelList = [
                     {'isMC':isMC,'channel':'ElEl','sel':hasOsElElLowMllCutOutZ,'cont1':electrons,'cont2':jetsSel,'suffix':'hasOsElElLowMllCutOutZ_ElectronJet'},
                     {'isMC':isMC,'channel':'MuMu','sel':hasOsMuMuLowMllCutOutZ,'cont1':electrons,'cont2':jetsSel,'suffix':'hasOsMuMuLowMllCutOutZ_ElectronJet'},
                     {'isMC':isMC,'channel':'ElMu','sel':hasOsElMuLowMllCut,    'cont1':electrons,'cont2':jetsSel,'suffix':'hasOsElMuLowMllCut_ElectronJet'},
                     {'isMC':isMC,'channel':'ElEl','sel':hasOsElElLowMllCutOutZ,'cont1':electrons,'cont2':fatjetsSel,'suffix':'hasOsElElLowMllCutOutZ_ElectronFatjet'},
                     {'isMC':isMC,'channel':'MuMu','sel':hasOsMuMuLowMllCutOutZ,'cont1':electrons,'cont2':fatjetsSel,'suffix':'hasOsMuMuLowMllCutOutZ_ElectronFatjet'},
                     {'isMC':isMC,'channel':'ElMu','sel':hasOsElMuLowMllCut,    'cont1':electrons,'cont2':fatjetsSel,'suffix':'hasOsElMuLowMllCut_ElectronFatjet'},
                     {'isMC':isMC,'channel':'ElEl','sel':hasOsElElLowMllCutOutZ,'cont1':muons,'cont2':jetsSel,'suffix':'hasOsElElLowMllCutOutZ_MuonJet'},
                     {'isMC':isMC,'channel':'MuMu','sel':hasOsMuMuLowMllCutOutZ,'cont1':muons,'cont2':jetsSel,'suffix':'hasOsMuMuLowMllCutOutZ_MuonJet'},
                     {'isMC':isMC,'channel':'ElMu','sel':hasOsElMuLowMllCut,    'cont1':muons,'cont2':jetsSel,'suffix':'hasOsElMuLowMllCut_MuonJet'},
                     {'isMC':isMC,'channel':'ElEl','sel':hasOsElElLowMllCutOutZ,'cont1':muons,'cont2':fatjetsSel,'suffix':'hasOsElElLowMllCutOutZ_MuonFatjet'},
                     {'isMC':isMC,'channel':'MuMu','sel':hasOsMuMuLowMllCutOutZ,'cont1':muons,'cont2':fatjetsSel,'suffix':'hasOsMuMuLowMllCutOutZ_MuonFatjet'},
                     {'isMC':isMC,'channel':'ElMu','sel':hasOsElMuLowMllCut,    'cont1':muons,'cont2':fatjetsSel,'suffix':'hasOsElMuLowMllCut_MuonFatjet'},
                   ]

        for channelDict in DeltaRChannelList:
            plots.extend(makeDeltaRPlots(self, **channelDict))

        # exclude from the jetsSel any jet that happens to include within its reconstruction cone a muon or an electron.
        jets = op.select(jetsSel, lambda j : op.AND(op.NOT(op.rng_any(electrons, lambda ele : op.deltaR(j.p4, ele.p4) < 0.3 )), op.NOT(op.rng_any(muons, lambda mu : op.deltaR(j.p4, mu.p4) < 0.3 ))))
        fatjets = op.select(fatjetsSel, lambda j : op.AND(op.NOT(op.rng_any(electrons, lambda ele : op.deltaR(j.p4, ele.p4) < 0.3 )), op.NOT(op.rng_any(muons, lambda mu : op.deltaR(j.p4, mu.p4) < 0.3 ))))
        # Yield plots #
        hasOsElElLowMllCutOutZAtLeast2Jets = hasOsElElLowMllCutOutZ.refine("hasOsElElLowMllCutOutZAtLeast2Jets",cut=[op.rng_len(jets)>=2])
        hasOsMuMuLowMllCutOutZAtLeast2Jets = hasOsMuMuLowMllCutOutZ.refine("hasOsMuMuLowMllCutOutZAtLeast2Jets",cut=[op.rng_len(jets)>=2])
        hasOsElMuLowMllCutAtLeast2Jets     = hasOsElMuLowMllCut.refine("hasOsElMuLowMllCutOutZAtLeast2Jets",cut=[op.rng_len(jets)>=2])

        hasOsElElLowMllCutOutZAtLeast1Fatjet = hasOsElElLowMllCutOutZ.refine("hasOsElElLowMllCutOutZAtLeast1Fatjet",cut=[op.rng_len(fatjets)>=1])
        hasOsMuMuLowMllCutOutZAtLeast1Fatjet = hasOsMuMuLowMllCutOutZ.refine("hasOsMuMuLowMllCutOutZAtLeast1Fatjet",cut=[op.rng_len(fatjets)>=1])
        hasOsElMuLowMllCutAtLeast1Fatjet     = hasOsElMuLowMllCut.refine("hasOsElMuLowMllCutOutZAtLeast1Fatjet",cut=[op.rng_len(fatjets)>=1])
        
        plots.append(makeYieldPlot(self,hasOsElElLowMllCutOutZAtLeast2Jets,"OSElElMllCutOutZAtLeast2Jets","OS leptons + $M_{ll}$ + Jets $\geq 2$ (channel : $e^+e^-$)",6))
        plots.append(makeYieldPlot(self,hasOsMuMuLowMllCutOutZAtLeast2Jets,"OSMuMuMllCutOutZAtLeast2Jets","OS leptons + $M_{ll}$ + Jets $\geq 2$ (channel : $\mu^+\mu^-$)",7))
        plots.append(makeYieldPlot(self,hasOsElMuLowMllCutAtLeast2Jets,"OSElMuMllCutAtLeast2Jets","OS leptons + $M_{ll}$ + Jets $\geq 2$ (channel : $e^{\pm}\mu^{\mp}$",8))

        plots.append(makeYieldPlot(self,hasOsElElLowMllCutOutZAtLeast1Fatjet,"OSElElMllCutOutZAtLeast1Fatjet","OS leptons + $M_{ll}$ + Fatjets $\geq 1$ (channel : $e^+e^-$)",9))
        plots.append(makeYieldPlot(self,hasOsMuMuLowMllCutOutZAtLeast1Fatjet,"OSMuMuMllCutOutZAtLeast1Fatjet","OS leptons + $M_{ll}$ + Fatjets $\geq 1$ (channel : $\mu^+\mu^-$)",10))
        plots.append(makeYieldPlot(self,hasOsElMuLowMllCutAtLeast1Fatjet,"OSElMuMllCutAtLeast1Fatjet","OS leptons + $M_{ll}$ + Fatjets $\geq 1$ (channel : $e^{\pm}\mu^{\mp}$)",11))

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
        bjetsBoosted        = op.select(fatjets, lambda_boosted)
        bjetsResolved       = op.select(jets, lambda_resolved)
        lightjetsResolved   = op.select(jets, lambda_notResolved) # To plot the case when 1 bjet + 1 lightjet

        # Scalefactors : 2017 and 2018 not yet present in the dict #
        DeepCSVMediumSFApplied = None
        DeepJetMediumSFApplied = None
        if self.isMC(sample):
            DeepJetTag_discriVar = {"BTagDiscri": lambda j : j.btagDeepFlavB}
            DeepJetMediumSF = SF.get_scalefactor("jet", ("btag_"+era+"_"+sfTag, "DeepJet_medium"), additionalVariables=DeepJetTag_discriVar, systName="deepjet") # For RESOLVED
            DeepJetMediumSFApplied = [DeepJetMediumSF(bjetsResolved[0])] # TODO : check if more than one bjet and apply to all
#            DeepCSVTag_discriVar = {"BTagDiscri": lambda j : j.btagDeepB}
#            DeepCSVMediumSF = SF.get_scalefactor("jet", ("subjet_btag_"+era+"_"+sfTag, "DeepCSV_medium"), additionalVariables=DeepCSVTag_discriVar, systName="deepcsv") # For BOOSTED (btag on subjet)
#            DeepCSVMediumSFApplied = [DeepCSVMediumSF(bjetsBoosted[0].subJet1)] # Must be applied on subjets : need to check each time which one has been btagged
    
        # Define the boosted and Resolved (+exclusive) selections #
        hasBoostedJets = noSel.refine("hasBoostedJets",
                               cut=[op.rng_len(bjetsBoosted)>=1],
                               weight = DeepCSVMediumSFApplied)
        hasResolvedJets = noSel.refine("hasResolvedJets",
                               cut=[op.rng_len(jets)>=2,op.rng_len(bjetsResolved)>=1],
                               weight = DeepJetMediumSFApplied)
        hasExclusiveResolvedJets = noSel.refine("hasExclusiveResolved",
                               cut=[op.rng_len(jets)>=2, op.rng_len(bjetsResolved)>=1,op.rng_len(bjetsBoosted)==0],
                               weight = DeepJetMediumSFApplied)
        hasExclusiveBoostedJets = noSel.refine("hasExclusiveBoostedJets",
                               cut=[op.rng_len(bjetsBoosted)>=1,op.OR(op.rng_len(jets)<=1,op.rng_len(bjetsResolved)==0)],
                               weight = DeepCSVMediumSFApplied)

        hasNotBoostedJets = noSel.refine("hasNotBoostedJets",
                               cut=[op.rng_len(bjetsBoosted)==0])
        hasNotResolvedJets = noSel.refine("hasNotResolvedJets",
                               cut=[op.OR(op.rng_len(jets)<=1,op.rng_len(bjetsResolved)==0)])
        hasBoostedAndResolvedJets = noSel.refine("hasBoostedAndResolvedJets", 
                               cut=[op.rng_len(bjetsBoosted)>=1,op.rng_len(jets)>=2,op.rng_len(bjetsResolved)>=1])
        hasNotBoostedAndResolvedJets = noSel.refine("hasNotBoostedAndResolvedJets", 
                               cut=[op.OR(op.rng_len(bjetsBoosted)==0,op.rng_len(jets)<=1,op.rng_len(bjetsResolved)==0)])
        hasNotExclusiveResolvedJets = noSel.refine("hasNotExclusiveResolved",
                               cut=[op.OR(op.OR(op.rng_len(jets)<=1,op.rng_len(bjetsResolved)==0),op.AND(op.rng_len(bjetsBoosted)>=1,op.rng_len(jets)>=2,op.rng_len(bjetsResolved)>=1))])
        hasNotExclusiveBoostedJets = noSel.refine("hasNotExclusiveBoostedJets",
                               cut=[op.OR(op.rng_len(bjetsBoosted)==0,op.AND(op.rng_len(jets)>=2,op.rng_len(bjetsResolved)>=1))])
        # Note that these selection should be done the same (and SF) when combining the OS lepton selection #
        # Because so far there is no way to "concatenate" two refine's #


        # Scalefactors plot #
#        plots.append(Plot.make1D("DeepJetMediumSF_ResolvedJets",
#                                DeepJetMediumSF(bjetsResolved[0]),
#                                hasResolvedJets,
#                                EquidistantBinning(100,0.,2.),
#                                title='DeepJetMediumSF',
#                                xTitle='DeepJetMediumSF'))
#        plots.append(Plot.make1D("DeepJetMediumSF_ExclusiveResolvedJets",
#                                DeepJetMediumSF(bjetsResolved[0]),
#                                hasExclusiveResolvedJets,
#                                EquidistantBinning(100,0.,2.),
#                                title='DeepJetMediumSF',
#                                xTitle='DeepJetMediumSF'))

        # Counting events from different selections for debugging #
        # Passing Boosted selection #
        PassedBoosted = Plot.make1D("PassedBoosted",
                                    op.c_int(1),
                                    hasBoostedJets,
                                    EquidistantBinning(2,0.,2.),
                                    title='Passed Boosted',
                                    xTitle='Passed Boosted')
        FailedBoosted = Plot.make1D("FailedBoosted",
                                    op.c_int(0),
                                    hasNotBoostedJets,
                                    EquidistantBinning(2,0.,2.),
                                    title='Failed Boosted',
                                    xTitle='Failed Boosted')
        plots.append(SummedPlot("BoostedCase",
                                [FailedBoosted,PassedBoosted],
                                xTitle="Boosted selection"))

        # Passing Resolved selection #
        PassedResolved = Plot.make1D("PassedResolved",
                                     op.c_int(1),
                                     hasResolvedJets,
                                     EquidistantBinning(2,0.,2.),
                                     title='Passed Resolved',
                                     xTitle='Passed Resolved')
        FailedResolved = Plot.make1D("FailedResolved",
                                     op.c_int(0),
                                     hasNotResolvedJets,
                                     EquidistantBinning(2,0.,2.),
                                     title='Failed Resolved',
                                     xTitle='Failed Resolved')
        plots.append(SummedPlot("ResolvedCase",
                                [FailedResolved,PassedResolved],
                                xTitle="Resolved selection"))

        # Passing Exclusive Resolved (Resolved AND NOT Boosted) #
        PassedExclusiveResolved = Plot.make1D("PassedExclusiveResolved",
                                              op.c_int(1),
                                              hasExclusiveResolvedJets,
                                              EquidistantBinning(2,0.,2.),
                                              title='Passed Exclusive Resolved',
                                              xTitle='Passed Exclusive Resolved')
        FailedExclusiveResolved = Plot.make1D("FailedExclusiveResolved",
                                              op.c_int(0),
                                              hasNotExclusiveResolvedJets,
                                              EquidistantBinning(2,0.,2.),
                                              title='Failed Exclusive Resolved',
                                              xTitle='Failed Exclusive Resolved')
        plots.append(SummedPlot("ExclusiveResolvedCase",
                                [FailedExclusiveResolved,PassedExclusiveResolved],
                                xTitle="Exclusive Resolved selection"))

        # Passing Exclusive Boosted (Boosted AND NOT Resolved) #
        PassedExclusiveBoosted = Plot.make1D("PassedExclusiveBoosted",
                                             op.c_int(1),
                                             hasExclusiveBoostedJets,
                                             EquidistantBinning(2,0.,2.),
                                             title='Passed Exclusive Boosted',
                                             xTitle='Passed Exclusive Boosted')
        FailedExclusiveBoosted = Plot.make1D("FailedExclusiveBoosted",
                                             op.c_int(0),
                                             hasNotExclusiveBoostedJets,
                                             EquidistantBinning(2,0.,2.),
                                             title='Failed Exclusive Boosted',
                                             xTitle='Failed Exclusive Boosted')
        plots.append(SummedPlot("ExclusiveBoostedCase",
                                [FailedExclusiveBoosted,PassedExclusiveBoosted],
                                xTitle="Exclusive Boosted selection"))

        # Passing Boosted AND Resolved #
        PassedBoth = Plot.make1D("PassedBoth",
                                 op.c_int(1),
                                 hasBoostedAndResolvedJets,
                                 EquidistantBinning(2,0.,2.),
                                 title='Passed Both Boosted and Resolved',
                                 xTitle='Passed Boosted and Resolved')
        FailedBoth = Plot.make1D("FailedBoth", # Means failed the (Boosted AND Resolved) = either one or the other 
                                 op.c_int(0),
                                 hasNotBoostedAndResolvedJets,
                                 EquidistantBinning(2,0.,2.),
                                 title='Failed combination Boosted and Resolved',
                                 xTitle='Failed combination')
        plots.append(SummedPlot("BoostedAndResolvedCase",[FailedBoth,PassedBoth],xTitle="Boosted and Resolved selection"))

        # Plot number of subjets in the boosted fatjets #
        lambda_noSubjet  = lambda fatjet : op.AND(fatjet.subJet1._idx.result == -1, op.AND(fatjet.subJet2._idx.result == -1 )) 
        lambda_oneSubjet = lambda fatjet : op.AND(fatjet.subJet1._idx.result != -1, op.AND(fatjet.subJet2._idx.result == -1 ))
        lambda_twoSubjet = lambda fatjet : op.AND(fatjet.subJet1._idx.result != -1, op.AND(fatjet.subJet2._idx.result != -1 ))

        hasNoSubjet   = hasBoostedJets.refine("hasNoSubjet",
                       cut=[lambda_noSubjet(bjetsBoosted[0])])
        hasOneSubjet  = hasBoostedJets.refine("hasOneSubjet",
                       cut=[lambda_oneSubjet(bjetsBoosted[0])])
        hasTwoSubjet  = hasBoostedJets.refine("hasTwoSubjet",
                       cut=[lambda_twoSubjet(bjetsBoosted[0])])

        plot_hasNoSubjet = Plot.make1D("plot_hasNoSubjet", # Fill bin 0
                                       op.c_int(0),
                                       hasNoSubjet,
                                       EquidistantBinning(3,0.,3.),
                                       title='Boosted jet without subjet')
        plot_hasOneSubjet = Plot.make1D("plot_hasOneSubjet", # Fill bin 1
                                        op.c_int(1),
                                        hasOneSubjet,
                                        EquidistantBinning(3,0.,3.),
                                        title='Boosted jet with one subjet')
        plot_hasTwoSubjet = Plot.make1D("plot_hasTwoSubjet", # Fill bin 2
                                        op.c_int(2),
                                        hasTwoSubjet,
                                        EquidistantBinning(3,0.,3.),title='Boosted jet with two subjets')
        plots.append(SummedPlot("NumberOfSubjets",
                                [plot_hasNoSubjet,plot_hasOneSubjet,plot_hasTwoSubjet],
                                xTitle="Number of subjets in boosted jet"))


        #############################################################################
        ##################### Jets + Dilepton combination ###########################
        #############################################################################

        ##### BOOSTED #####
        # Combine dilepton and Boosted selections #
        hasOsElElLowMllCutOutZBoostedJets = hasOsElElLowMllCutOutZ.refine("hasOsElElLowMllCutOutZBoostedJets", 
                                      cut=[op.rng_len(bjetsBoosted)>=1],
                                      weight = DeepCSVMediumSFApplied)
        hasOsMuMuLowMllCutOutZBoostedJets = hasOsMuMuLowMllCutOutZ.refine("hasOsMuMuLowMllCutOutZBoostedJets", 
                                      cut=[op.rng_len(bjetsBoosted)>=1],
                                      weight = DeepCSVMediumSFApplied)
        hasOsElMuLowMllCutBoostedJets = hasOsElMuLowMllCut.refine("hasOsElMuLowMllCutBoostedJets", 
                                      cut=[op.rng_len(bjetsBoosted)>=1],
                                      weight = DeepCSVMediumSFApplied)

        hasOsElElLowMllCutOutZExclusiveBoostedJets = hasOsElElLowMllCutOutZ.refine("hasOsElElLowMllCutOutZExclusiveBoostedJets", 
                                      cut=[op.rng_len(bjetsBoosted)>=1,op.OR(op.rng_len(jets)<=1,op.rng_len(bjetsResolved)==0)],
                                      weight = DeepCSVMediumSFApplied)
        hasOsMuMuLowMllCutOutZExclusiveBoostedJets = hasOsMuMuLowMllCutOutZ.refine("hasOsMuMuLowMllCutOutZExclusiveBoostedJets", 
                                      cut=[op.rng_len(bjetsBoosted)>=1,op.OR(op.rng_len(jets)<=1,op.rng_len(bjetsResolved)==0)],
                                      weight = DeepCSVMediumSFApplied)
        hasOsElMuLowMllCutExclusiveBoostedJets = hasOsElMuLowMllCut.refine("hasOsElMuLowMllCutExclusiveBoostedJets", 
                                      cut=[op.rng_len(bjetsBoosted)>=1,op.OR(op.rng_len(jets)<=1,op.rng_len(bjetsResolved)==0)],
                                      weight = DeepCSVMediumSFApplied)

        # Yield plots #
        plots.append(makeYieldPlot(self,hasOsElElLowMllCutOutZBoostedJets,"OSElElMllCutOutZBoosted","OS leptons + $M_{ll}$ + Boosted (channel : $e^+e^-$)",12))
        plots.append(makeYieldPlot(self,hasOsMuMuLowMllCutOutZBoostedJets,"OSMuMuMllCutOutZBoosted","OS leptons + $M_{ll}$ + Boosted (channel : $\mu^+\mu^-$)",13))
        plots.append(makeYieldPlot(self,hasOsElMuLowMllCutBoostedJets,"OSMuMuMllCutBoosted","OS leptons + $M_{ll}$ + Boosted (channel : $e^{\pm}\mu^{\mp}$)",14))

        plots.append(makeYieldPlot(self,hasOsElElLowMllCutOutZExclusiveBoostedJets,"OSElElMllCutOutZExclusiveBoosted","OS leptons + $M_{ll}$ + Exclusive Boosted (channel : $e^+e^-$)",15))
        plots.append(makeYieldPlot(self,hasOsMuMuLowMllCutOutZExclusiveBoostedJets,"OSMuMuMllCutOutZExclusiveBoosted","OS leptons + $M_{ll}$ + Exclusive Boosted (channel : $\mu^+\mu^-$)",16))
        plots.append(makeYieldPlot(self,hasOsElMuLowMllCutExclusiveBoostedJets,"OSMuMuMllCutExclusiveBoosted","OS leptons + $M_{ll}$ + Exclusive Boosted (channel : $e^{\pm}\mu^{\mp}$)",17))

        # Boosted + OS dilepton plots #
        hasOsMllCutBoostedChannelList = [
                             {'channel':'ElEl','sel':hasOsElElLowMllCutOutZBoostedJets,'dilepton':OsElEl[0],'suffix':'hasOsElElLowMllCutOutZBoostedJets'},
                             {'channel':'MuMu','sel':hasOsMuMuLowMllCutOutZBoostedJets,'dilepton':OsMuMu[0],'suffix':'hasOsMuMuLowMllCutOutZBoostedJets'},
                             {'channel':'ElMu','sel':hasOsElMuLowMllCutBoostedJets,    'dilepton':OsElMu[0],'suffix':'hasOsElMuLowMllCutBoostedJets'},

        #                     {'channel':'ElEl','sel':hasOsElElLowMllCutOutZExclusiveBoostedJets,'dilepton':OsElEl[0],'suffix':'hasOsElElLowMllCutOutZExclusiveBoostedJets'},
        #                     {'channel':'MuMu','sel':hasOsMuMuLowMllCutOutZExclusiveBoostedJets,'dilepton':OsMuMu[0],'suffix':'hasOsMuMuLowMllCutOutZExclusiveBoostedJets'},
        #                     {'channel':'ElMu','sel':hasOsElMuLowMllCutExclusiveBoostedJets,    'dilepton':OsElMu[0],'suffix':'hasOsElMuLowMllCutExclusiveBoostedJets'},
                           ]

        for channelDict in hasOsMllCutBoostedChannelList:
            # Dilepton plots #
            plots.extend(makeDileptonPlots(self, **channelDict))
            # MET plots #
            plots.extend(makeMETPlots(self      = self,
                                      sel       = channelDict['sel'],
                                      met       = corrMET,
                                      suffix    = channelDict['suffix'],
                                      channel   = channelDict['channel']))
            # Fatjet plots #
            plots.extend(makeFatJetPlots(self,
                                         sel       = channelDict['sel'],
                                         fatjets   = bjetsBoosted,
                                         suffix    = channelDict['suffix'],
                                         channel   = channelDict['channel']))
  
        ##### RESOLVED #####
        # Combine dilepton and Exclusive Resolved (Exclusive = NOT Boosted) selections #
        hasOsElElLowMllCutOutZResolvedJets = hasOsElElLowMllCutOutZ.refine("hasOsElElLowMllCutOutZResolvedJets", 
                                      cut=[op.rng_len(jets)>=2, op.rng_len(bjetsResolved)>=1],
                                      weight = DeepJetMediumSFApplied)
        hasOsMuMuLowMllCutOutZResolvedJets = hasOsMuMuLowMllCutOutZ.refine("hasOsMuMuLowMllCutOutZResolvedJets", 
                                      cut=[op.rng_len(jets)>=2, op.rng_len(bjetsResolved)>=1], 
                                      weight = DeepJetMediumSFApplied)
        hasOsElMuLowMllCutResolvedJets = hasOsElMuLowMllCut.refine("hasOsElMuLowMllCutResolvedJets", 
                                      cut=[op.rng_len(jets)>=2, op.rng_len(bjetsResolved)>=1],
                                      weight = DeepJetMediumSFApplied)

        hasOsElElLowMllCutOutZExclusiveResolvedJets = hasOsElElLowMllCutOutZ.refine("hasOsElElLowMllCutOutZExclusiveResolvedJets", 
                                      cut=[op.rng_len(jets)>=2, op.rng_len(bjetsResolved)>=1,op.rng_len(bjetsBoosted)==0],
                                      weight = DeepJetMediumSFApplied)
        hasOsMuMuLowMllCutOutZExclusiveResolvedJets = hasOsMuMuLowMllCutOutZ.refine("hasOsMuMuLowMllCutOutZExclusiveResolvedJets", 
                                      cut=[op.rng_len(jets)>=2, op.rng_len(bjetsResolved)>=1,op.rng_len(bjetsBoosted)==0], 
                                      weight = DeepJetMediumSFApplied)
        hasOsElMuLowMllCutExclusiveResolvedJets = hasOsElMuLowMllCut.refine("hasOsElMuLowMllCutExclusiveResolvedJets", 
                                      cut=[op.rng_len(jets)>=2, op.rng_len(bjetsResolved)>=1,op.rng_len(bjetsBoosted)==0],
                                      weight = DeepJetMediumSFApplied)

        # Yield plots #
        plots.append(makeYieldPlot(self,hasOsElElLowMllCutOutZResolvedJets,"OSElElMllCutOutZResolved","OS leptons + $M_{ll}$ + Resolved (channel : $e^+e^-$)",18))
        plots.append(makeYieldPlot(self,hasOsMuMuLowMllCutOutZResolvedJets,"OSMuMuMllCutOutZResolved","OS leptons + $M_{ll}$ + Resolved (channel : $\mu^+\mu^-$)",19))
        plots.append(makeYieldPlot(self,hasOsElMuLowMllCutResolvedJets,"OSMuMuMllCutResolved","OS leptons + $M_{ll}$ + Resolved (channel : $e^{\pm}\mu^{\mp}$)",20))

        plots.append(makeYieldPlot(self,hasOsElElLowMllCutOutZExclusiveResolvedJets,"OSElElMllCutOutZExclusiveResolved","OS leptons + $M_{ll}$ + Exclusive Resolved (channel : $e^+e^-$)",21))
        plots.append(makeYieldPlot(self,hasOsMuMuLowMllCutOutZExclusiveResolvedJets,"OSMuMuMllCutOutZExclusiveResolved","OS leptons + $M_{ll}$ + Exclusive Resolved (channel : $\mu^+\mu^-$)",22))
        plots.append(makeYieldPlot(self,hasOsElMuLowMllCutExclusiveResolvedJets,"OSMuMuMllCutExclusiveResolved","OS leptons + $M_{ll}$ + Exclusive Resolved (channel : $e^{\pm}\mu^{\mp}$)",23))

        # ExclusiveResolved + OS dilepton plots #
        hasOsMllCutExclusiveResolvedChannelList = [
        #                     {'channel':'ElEl','sel':hasOsElElLowMllCutOutZResolvedJets,'dilepton':OsElEl[0],'suffix':'hasOsElElLowMllCutOutZResolvedJets'},
        #                     {'channel':'MuMu','sel':hasOsMuMuLowMllCutOutZResolvedJets,'dilepton':OsMuMu[0],'suffix':'hasOsMuMuLowMllCutOutZResolvedJets'},
        #                     {'channel':'ElMu','sel':hasOsElMuLowMllCutResolvedJets,    'dilepton':OsElMu[0],'suffix':'hasOsElMuLowMllCutResolvedJets'},

                             {'channel':'ElEl','sel':hasOsElElLowMllCutOutZExclusiveResolvedJets,'dilepton':OsElEl[0],'suffix':'hasOsElElLowMllCutOutZExclusiveResolvedJets'},
                             {'channel':'MuMu','sel':hasOsMuMuLowMllCutOutZExclusiveResolvedJets,'dilepton':OsMuMu[0],'suffix':'hasOsMuMuLowMllCutOutZExclusiveResolvedJets'},
                             {'channel':'ElMu','sel':hasOsElMuLowMllCutExclusiveResolvedJets,    'dilepton':OsElMu[0],'suffix':'hasOsElMuLowMllCutExclusiveResolvedJets'},
                           ]

        for channelDict in hasOsMllCutExclusiveResolvedChannelList:
            # Dilepton plots #
            plots.extend(makeDileptonPlots(self, **channelDict))
            # MET plots #
            plots.extend(makeMETPlots(self      = self,
                                      sel       = channelDict['sel'],
                                      met       = corrMET,
                                      suffix    = channelDict['suffix'],
                                      channel   = channelDict['channel']))
            # Jets plots #
            plots.extend(makeJetsPlots(self,
                                       sel         = channelDict['sel'],
                                       bjets       = bjetsResolved,
                                       lightjets   = lightjetsResolved,
                                       alljets     = jets,
                                       suffix      = channelDict['suffix'],
                                       channel     = channelDict['channel']))

        # High level quantities #
        highLevelBaseQuantities = { # No need to repeat for each channelDict
             'met'          : corrMET,
             'jets'         : jets, 
             'resolvedjets' : bjetsResolved,
             'lightjets'    : lightjetsResolved,
             'boostedjets'  : bjetsBoosted,
                                  }
 
        hasOsMllCutHighLevelVariablesChannelList = [
            # Boosted #
            {'channel'      :'ElEl',
             'dilepton'     : OsElEl[0],
             'sel'          : hasOsElElLowMllCutOutZBoostedJets,
             'suffix'       : 'hasOsElElLowMllCutOutZBoostedJets'},
            {'channel'      :'MuMu',
             'dilepton'     : OsMuMu[0],
             'sel'          : hasOsMuMuLowMllCutOutZBoostedJets,
             'suffix'       : 'hasOsMuMuLowMllCutOutZBoostedJets'},
            {'channel'      :'ElMu',
             'dilepton'     : OsElMu[0],
             'sel'          : hasOsElMuLowMllCutBoostedJets,
             'suffix'       : 'hasOsElMuLowMllCutBoostedJets'},
            # Resolved #
            {'channel'      :'ElEl',
             'dilepton'     : OsElEl[0],
             'sel'          : hasOsElElLowMllCutOutZExclusiveResolvedJets,
             'suffix'       : 'hasOsElElLowMllCutOutZExclusiveResolvedJets'},
            {'channel'      :'MuMu',
             'dilepton'     : OsMuMu[0],
             'sel'          : hasOsMuMuLowMllCutOutZExclusiveResolvedJets,
             'suffix'       : 'hasOsMuMuLowMllCutOutZExclusiveResolvedJets'},
            {'channel'      :'ElMu',
             'dilepton'     : OsElMu[0],
             'sel'          : hasOsElMuLowMllCutExclusiveResolvedJets,
             'suffix'       : 'hasOsElMuLowMllCutExclusiveResolvedJets'},
                                                   ]

        for channelDict in hasOsMllCutHighLevelVariablesChannelList:
            plots.extend(makeHighLevelQuantities(self,**highLevelBaseQuantities,**channelDict))


        ## helper selection (OR) to make sure jet calculations are only done once
        #hasOSLL = noSel.refine("hasOSLL", cut=op.OR(*( hasOSLL_cmbRng(rng) for rng in osLLRng.values())))
        #forceDefine(t._Jet.calcProd, hasOSLL)
        #for varNm in t._Jet.available:
        #    forceDefine(t._Jet[varNm], hasOSLL)
        return plots

