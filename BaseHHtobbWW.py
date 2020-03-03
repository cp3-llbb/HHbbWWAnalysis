import os
import sys
from itertools import chain

from bamboo import treefunctions as op
from bamboo.analysismodules import NanoAODModule, NanoAODHistoModule, NanoAODSkimmerModule
from bamboo.analysisutils import makeMultiPrimaryDatasetTriggerSelection
from bamboo.scalefactors import binningVariables_nano

from scalefactorsbbWW import ScaleFactorsbbWW
from METScripts import METFilter, METcorrection

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

    def prepareTree(self, tree, sample=None, sampleCfg=None):
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
        cachJEC_dir = '/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/cacheJEC'
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

    def prepareObjects(self, t, noSel, sample, sampleCfg):
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
#        if self.isMC(sample) and sample.startswith("TT"):
#            # https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting#Use_case_3_ttbar_MC_is_used_to_m
#            # Get tops #
#            genTop_all = op.select(t.GenPart,lambda g : g.pdgId==6)
#            genTop = op.select(genTop_all,lambda g : g.statusFlags & ( 0x1 << 13))
#            genAntitop_all = op.select(t.GenPart,lambda g : g.pdgId==-6)
#            genAntitop = op.select(genAntitop_all,lambda g : g.statusFlags & ( 0x1 << 13))
#                # statusFlags==13 : isLastCopy
#                # Pdgid == 6 : top
#            hasttbar = noSel.refine("hasttbar",cut=[op.rng_len(genTop)>=1,op.rng_len(genAntitop)>=1])
#
#            # Check plots #
#            #plots.append(Plot.make1D("N_tops",
#            #                        op.rng_len(genTop),
#            #                        noSel,
#            #                        EquidistantBinning(5,0.,5.),
#            #                        title='N tops',
#            #                        xTitle='N tops'))
#            #plots.append(Plot.make1D("N_antitops",
#            #                        op.rng_len(genAntitop),
#            #                        noSel,
#            #                        EquidistantBinning(5,0.,5.),
#            #                        title='N antitops',
#            #                        xTitle='N antitops'))
#            ## Top pt #
#            #plots.append(Plot.make1D("top_pt",
#            #                        genTop[0].pt,
#            #                        hasttbar,
#            #                        EquidistantBinning(50,0,500),
#            #                        xTitle='P_{T} top'))
#            ## Antitop Pt #
#            #plots.append(Plot.make1D("antitop_pt",
#            #                        genAntitop[0].pt,
#            #                        hasttbar,
#            #                        EquidistantBinning(50,0,500),
#            #                        xTitle='P_{T} lead antitop'))
#            ## Top - Antitop plots #
#            #plots.append(Plot.make2D("top_pt_vs_antitop_pt",
#            #                        [genTop[0].pt,genAntitop[0].pt],
#            #                        hasttbar,
#            #                        [EquidistantBinning(50,0,500),EquidistantBinning(50,0,500)],
#            #                        xTitle='P_{T} lead top',
#            #                        yTitle='P_{T} lead antitop'))
# 
#
#            # Compute weight if there is a ttbar #
#            ttbar_SF = lambda t : op.exp(0.0615-0.0005*t.pt)
#            ttbar_weight = lambda t,tbar : op.sqrt(ttbar_SF(t)*ttbar_SF(tbar))
#            #plots.append(Plot.make1D("ttbar_weight",
#            #                        ttbar_weight(genTop[0],genAntitop[0]),
#            #                        hasttbar,
#            #                        EquidistantBinning(100,0.,2.),
#            #                        title='ttbar weight',
#            #                        xTitle='ttbar weight'))
#            #plots.append(Plot.make3D("ttbar_weight_vs_pt",
#            #                        [genTop[0].pt,genAntitop[0].pt,ttbar_weight(genTop[0],genAntitop[0])],
#            #                        hasttbar,
#            #                        [EquidistantBinning(50,0,500),EquidistantBinning(50,0,500),EquidistantBinning(100,0.,2.)],
#            #                        title='ttbar weight',
#            #                        xTitle='top P_{T}',
#            #                        yTitle='antitop P_{T}',
#            #                        zTitle='ttbar weight'))
#            # Apply correction to TT #
#            noSel = noSel.refine("ttbarWeight",weight=ttbar_weight(genTop[0],genAntitop[0]))

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
            self.PUWeight = makePileupWeight(puWeightsFile, t.Pileup_nTrueInt, systName="pileup")
            #noSel = noSel.refine("puWeight", weight=makePileupWeight(puWeightsFile, t.Pileup_nTrueInt, systName="pileup"))
            noSel = noSel.refine("puWeight", weight=self.PUWeight)

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
        # lepton lambdas #
        self.lambda_hasAssociatedJet = lambda lep : op.AND(lep.jet.idx != -1 , op.deltaR(lep.p4,lep.jet.p4) <= 0.4)
        if era == "2016": 
            self.lambda_lepton_associatedJetNoBtag = lambda lep : op.OR(op.NOT(self.lambda_hasAssociatedJet(lep)),
                                                                        lep.jet.btagDeepFlavB < 0.3093)
        elif era =="2017":
            self.lambda_lepton_associatedJetNoBtag = lambda lep : op.OR(op.NOT(self.lambda_hasAssociatedJet(lep)),
                                                                        lep.jet.btagDeepFlavB < 0.3033)
        elif era == "2018":
            self.lambda_lepton_associatedJetNoBtag = lambda lep : op.OR(op.NOT(self.lambda_hasAssociatedJet(lep)),
                                                                        lep.jet.btagDeepFlavB < 0.2770)

        # Preselection #
        muonsByPt = op.sort(t.Muon, lambda mu : -mu.p4.Pt())
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
        #self.lambda_muon_conept = lambda mu : op.switch(mu.mvaTTH >= 0.85, mu.pt, op.static_cast("float",0.9*mu.jet.pt))
        self.lambda_muon_conept = lambda mu : op.multiSwitch((mu.mvaTTH >= 0.85, mu.pt>=10),
                                                             (self.lambda_hasAssociatedJet(mu), 0.9*mu.jet.pt>=10),
                                                             op.c_bool(False))
        self.lambda_muon_x = lambda mu : op.min(op.max(0.,(0.9*mu.pt*(1+mu.jetRelIso))-20.)/(45.-20.), 1.)
                    # x = min(max(0, jet_pt-PT_min)/(PT_max-PT_min), 1) where jet_pt = 0.9*PT_muon*(1+MuonJetRelIso), PT_min=25, PT_max=40
        if era == "2016": 
            self.lambda_muon_btagInterpolation = lambda mu : self.lambda_muon_x(mu)*0.0614 + (1-self.lambda_muon_x(mu))*0.3093
        elif era =="2017":
            self.lambda_muon_btagInterpolation = lambda mu : self.lambda_muon_x(mu)*0.0521 + (1-self.lambda_muon_x(mu))*0.3033
        elif era == "2018":
            self.lambda_muon_btagInterpolation = lambda mu : self.lambda_muon_x(mu)*0.0494 + (1-self.lambda_muon_x(mu))*0.2770
            # return x*WP_loose+(1-x)*WP_medium
        self.lambda_muon_deepJetInterpIfMvaFailed = lambda mu : op.AND(self.lambda_hasAssociatedJet(mu),
                                                                      mu.jet.btagDeepFlavB < self.lambda_muon_btagInterpolation(mu))
        self.lambda_muonFakeSel = lambda mu : op.AND(
                                                    self.lambda_muon_conept(mu),
                                                    self.lambda_lepton_associatedJetNoBtag(mu),
                                                    op.OR(mu.mvaTTH >= 0.85, op.AND(mu.jetRelIso<0.5 , self.lambda_muon_deepJetInterpIfMvaFailed(mu))), 
                                                        # If mvaTTH < 0.85 : jetRelIso <0.5 and < deepJet medium with interpolation
                                                        # https://indico.cern.ch/event/812025/contributions/3475878/attachments/1867083/3070589/gp-fr-run2b.pdf (slide 7)
                                                )
        self.muonsFakeSel = op.select(self.muonsPreSel, self.lambda_muonFakeSel)
        # Tight selection #
        self.lambda_muonTightSel = lambda mu : op.AND(
                                                    mu.mvaTTH >= 0.85, # Lepton MVA id from ttH
                                                    mu.mediumId
                                                )
        self.muonsTightSel = op.select(self.muonsFakeSel, self.lambda_muonTightSel)

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
        # Preselection #
        electronsByPt = op.sort(t.Electron, lambda ele : -ele.p4.Pt())
        self.lambda_electronPreSel = lambda ele : op.AND(
                                                        ele.p4.Pt() >= 7.,
                                                        op.abs(ele.p4.Eta()) <= 2.5,
                                                        op.abs(ele.dxy) <= 0.05,
                                                        op.abs(ele.dz) <= 0.1,
                                                        ele.miniPFRelIso_all <= 0.4, # mini PF relative isolation, total (with scaled rho*EA PU corrections)
                                                        ele.sip3d <= 8,
                                                        ele.cutBased>=2,     # cut-based ID Fall17 V2 (0:fail, 1:veto, 2:loose, 3:medium, 4:tight)
                                                        ele.lostHits <=1    # number of missing inner hits
                                                        )
        self.electronsPreSelInclu = op.select(electronsByPt, self.lambda_electronPreSel) # can include a muon in cone
        self.lambda_cleanElectron = lambda ele : op.NOT(op.rng_any(self.muonsPreSel, lambda mu : op.deltaR(mu.p4, ele.p4) < 0.3 ))
            # No overlap between electron and muon in cone of DR<=0.3
        self.electronsPreSel = op.select(self.electronsPreSelInclu, self.lambda_cleanElectron)
        # Fakeable selection #
        #self.lambda_electron_conept = lambda ele : op.switch(ele.mvaTTH >= 0.80, ele.pt, op.static_cast("Float_t",0.9*ele.jet.pt))
        self.lambda_electron_conept = lambda ele : op.multiSwitch((ele.mvaTTH >= 0.80, ele.pt>=10),
                                                                   (self.lambda_hasAssociatedJet(ele), 0.9*ele.jet.pt>=10),
                                                                   op.c_bool(False))
        self.lambda_electronFakeSel = lambda ele : op.AND(
                                                        self.lambda_electron_conept(ele),
                                                        op.OR(
                                                                op.AND(op.abs(ele.eta)<1.479, ele.sieie<=0.011), 
                                                                op.AND(op.AND(op.abs(ele.eta)>=1.479,op.abs(ele.eta)<=2.5), ele.sieie<=0.030)),
                                                        self.lambda_lepton_associatedJetNoBtag(ele),
                                                        ele.hoe <= 0.10,
                                                        ele.eInvMinusPInv >= -0.04,
                                                        op.OR(ele.mvaTTH >= 0.80, op.AND(ele.jetRelIso<0.7, ele.mvaFall17V2noIso_WP80)), # Lepton MVA id from ttH
                                                            # If mvaTTH < 0.80 : jetRelIso <0.5 (and deepJet medium ?) and Fall17V2noIso_WP80
                                                        ele.lostHits == 0,    # number of missing inner hits
                                                        ele.convVeto        # Passes conversion veto
                                                        )
        self.electronsFakeSel = op.select(self.electronsPreSel, self.lambda_electronFakeSel)
        # Tight selection #
        self.lambda_electronTightSel = lambda ele : op.AND(
                                                    ele.mvaTTH >= 0.80
                                                )
        self.electronsTightSel = op.select(self.electronsFakeSel, self.lambda_electronTightSel)

        # Scalefactors #
        if self.isMC(sample):
            elTightIDSF = SF.get_scalefactor("lepton", ("electron_{0}_{1}".format(era,sfTag), "id_tight"), systName="elid")

        #############################################################################
        #                                AK4 Jets                                   #
        #############################################################################
        self.ak4JetsByPt = op.sort(t.Jet, lambda jet : -jet.p4.Pt())
        # Preselection #
        if era == "2016":
            self.lambda_ak4JetsPreSel = lambda j : op.AND(
                                                        j.jetId & 1, # Jet ID flags bit1 is loose, bit2 is tight, bit3 is tightLepVeto
                                                        j.p4.Pt() > 25.,
                                                        op.abs(j.p4.Eta())< 2.4
                                                    )
        elif  era == "2017" or era == "2018":
            self.lambda_ak4JetsPreSel = lambda j : op.AND(
                                                        j.jetId & 2, # Jet ID flags bit1 is loose, bit2 is tight, bit3 is tightLepVeto
                                                        j.p4.Pt() > 25.,
                                                        op.abs(j.p4.Eta())< 2.4
                                                    )
        self.ak4JetsPreSel = op.select(self.ak4JetsByPt, self.lambda_ak4JetsPreSel)
        # Cleaning #
        self.lambda_cleanAk4Jets = lambda j : op.AND(
                                                op.NOT(op.rng_any(self.electronsPreSel, lambda ele : op.deltaR(j.p4, ele.p4) < 0.4 )), 
                                                op.NOT(op.rng_any(self.muonsPreSel, lambda mu : op.deltaR(j.p4, mu.p4) < 0.4 )))
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

        ############     ScaleFactors    #############
        self.DeepJetMediumSFApplied = None
        if self.isMC(sample):
            DeepJetTag_discriVar = {"BTagDiscri": lambda j : j.btagDeepFlavB}
            self.DeepJetMediumSF = SF.get_scalefactor("jet", ("btag_"+era+"_"+sfTag, "DeepJet_medium"), additionalVariables=DeepJetTag_discriVar, systName="deepjet") # For RESOLVED
            DeepJetMediumSFApplied = [self.DeepJetMediumSF(self.ak4BJets[0])] # TODO : check if more than one bjet and apply to all

        #############################################################################
        #                                AK8 Jets                                   #
        #############################################################################
        self.ak8JetsByPt = op.sort(t.FatJet, lambda jet : -jet.p4.Pt())
        # Preselection #
        if era == "2016":
            self.lambda_ak8JetsPreSel = lambda j : op.AND(
                                                        j.jetId & 1, # Jet ID flags bit1 is loose, bit2 is tight, bit3 is tightLepVeto
                                                        j.p4.Pt() > 200.,
                                                        op.abs(j.p4.Eta())< 2.4,
                                                        op.AND(j.msoftdrop > 30, j.msoftdrop < 210),
                                                        j.tau2/j.tau1 < 0.75
                                                    )
        elif  era == "2017" or era == "2018":
            self.lambda_ak8JetsPreSel = lambda j : op.AND(
                                                        j.jetId & 2, # Jet ID flags bit1 is loose, bit2 is tight, bit3 is tightLepVeto
                                                        j.p4.Pt() > 200.,
                                                        op.abs(j.p4.Eta())< 2.4,
                                                        op.AND(j.msoftdrop > 30, j.msoftdrop < 210),
                                                        j.tau2/j.tau1 < 0.75
                                                    )
        self.ak8JetsPreSel = op.select(self.ak8JetsByPt, self.lambda_ak8JetsPreSel)
        # Cleaning #
        self.lambda_cleanAk8Jets = lambda j : op.AND(
                                                op.NOT(op.rng_any(self.electronsPreSel, lambda ele : op.deltaR(j.p4, ele.p4) < 0.8 )), 
                                                op.NOT(op.rng_any(self.muonsPreSel, lambda mu : op.deltaR(j.p4, mu.p4) < 0.8 ))
                                            )
            # remove jets within cone of DR<0.8 of preselected electrons and muons
        self.ak8Jets = op.select(self.ak8JetsPreSel,self.lambda_cleanAk8Jets)

        ############     Btagging     #############
        # The DeepCSV b-tagging algorithm is used on subjets 
        # If must check that subJet exists before looking at the btag : fatjet.subJet1._idx.result != -1
                                                        #op.AND(fatjet.subJet1._idx.result != -1,fatjet.subJet1.btagDeepB > 0.6321), 
                                                        #op.AND(fatjet.subJet2._idx.result != -1,fatjet.subJet2.btagDeepB > 0.6321)))
        if era == "2016": 
            self.lambda_ak8Btag = lambda fatjet : op.OR(op.AND(fatjet.subJet1.pt >30, fatjet.subJet1.btagDeepB > 0.6321),
                                                        op.AND(fatjet.subJet2.pt >30, fatjet.subJet2.btagDeepB > 0.6321))
        elif era =="2017":
            self.lambda_ak8Btag = lambda fatjet : op.OR(op.AND(fatjet.subJet1.pt >30, fatjet.subJet1.btagDeepB > 0.4941),
                                                        op.AND(fatjet.subJet2.pt >30, fatjet.subJet2.btagDeepB > 0.4941))
        elif era == "2018":
            self.lambda_ak8Btag = lambda fatjet : op.OR(op.AND(fatjet.subJet1.pt >30, fatjet.subJet1.btagDeepB > 0.4184),
                                                        op.AND(fatjet.subJet2.pt >30, fatjet.subJet2.btagDeepB > 0.4184))

        self.ak8BJets = op.select(self.ak8Jets, self.lambda_ak8Btag)

        ############     ScaleFactors    #############
        if self.isMC(sample): # Needs new nanoAOD version 
            pass
#            DeepCSVTag_discriVar = {"BTagDiscri": lambda j : j.btagDeepB}
#            DeepCSVMediumSF = SF.get_scalefactor("jet", ("subjet_btag_"+era+"_"+sfTag, "DeepCSV_medium"), additionalVariables=DeepCSVTag_discriVar, systName="deepcsv") # For BOOSTED (btag on subjet)
#            DeepCSVMediumSFApplied = [DeepCSVMediumSF(self.ak8BJets[0].subJet1)] # Must be applied on subjets : need to check each time which one has been btagged

        # Return #
        return noSel
