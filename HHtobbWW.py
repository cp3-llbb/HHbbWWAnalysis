import sys
import os
from itertools import chain

from bamboo.analysismodules import NanoAODHistoModule
from bamboo.analysisutils import makeMultiPrimaryDatasetTriggerSelection
from bamboo.scalefactors import binningVariables_nano

from bamboo import treefunctions as op
from bamboo.plots import Plot, EquidistantBinning, SummedPlot
from bamboo import treefunctions as op

sys.path.append('/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/') # Add scripts in this directory -- TODO : make cleaner
from plotDef import makeDileptonPlots, makeJetsPlots, makeFatJetPlots

class NanoHHTobbWW(NanoAODHistoModule):
    """ Example module: HH->bbW(->e/µ nu)W(->e/µ nu) histograms from NanoAOD """
    def __init__(self, args):
            super(NanoHHTobbWW, self).__init__(args)
            self.calcToAdd += ["nJet", "nMuon"]

    def prepareTree(self, tree, sample=None, sampleCfg=None, enableSystematics=None):
        # JEC's Recommendation for Full RunII: https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC
        # JER : -----------------------------: https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
        # Get base aguments #
        tree,noSel,be,lumiArgs = super(NanoHHTobbWW,self).prepareTree(tree, sample=sample, sampleCfg=sampleCfg)
        era = sampleCfg['era']
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
            configureRochesterCorrection(tree._Muon.calc,os.path.join(os.path.dirname(__file__), "data", "RoccoR2016.txt"))

            # Trigger efficiencies #
            if self.isMC(sample):
                triggersPerPrimaryDataset = {
                    "SingleMuon" :  [ tree.HLT.IsoMu24],
                    "SingleEG"   :  [ tree.HLT.Ele27_WPTight_Gsf],
                    "DoubleMuon" :  [tree.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL,
                                     tree.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ],
                    "DoubleEGamma": [tree.HLT.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ],
                    "MuonEG":       [tree.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,
                                     tree.HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ]}
            else: # Some branches have different names in data
                triggersPerPrimaryDataset = {
                    "SingleMuon" :  [ tree.HLT.IsoMu24],
                    "SingleEG"   :  [ tree.HLT.Ele27_WPTight_Gsf],
                    "DoubleMuon" :  [tree.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL,
                                     tree.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ],
                    "DoubleEGamma": [tree.HLT.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ],
                    "MuonEG":       [tree.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL, 
                                     tree.HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL]}


            # Jet treatment #
            cachJEC_dir = '/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/cacheJEC'
            if self.isMC(sample):   # if MC -> needs smearing
                configureJets(tree, "Jet", "AK4PFchs",
                    jec="Summer16_07Aug2017_V20_MC",
                    #smear="Summer16_25nsV1_MC",
                    jesUncertaintySources=["Total"],
                    mayWriteCache=isNotWorker,
                    cachedir=cachJEC_dir)
            else:                   # If data -> extract info from config 
                if "2016B" in sample or "2016C" in sample or "2016D" in sample:
                    configureJets(tree,"Jet","AK4PFchs",
                        jec="Summer16_07Aug2017BCD_V11_DATA", mayWriteCache=isNotWorker,cachedir=cachJEC_dir)
                elif "2016E" in sample or "2016F" in sample:
                    configureJets(tree,"Jet","AK4PFchs",
                        jec="Summer16_07Aug2017EF_V11_DATA", mayWriteCache=isNotWorker,cachedir=cachJEC_dir)
                elif "2016G" in sample or "2016H" in sample:
                    configureJets(tree,"Jet","AK4PFchs",
                        jec="Summer16_07Aug2017GH_V11_DATA", mayWriteCache=isNotWorker,cachedir=cachJEC_dir)

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
        else:
            noSel = noSel.refine("withTrig", cut=makeMultiPrimaryDatasetTriggerSelection(sample, triggersPerPrimaryDataset))

        return tree,noSel,be,lumiArgs
        
    def definePlots(self, t, noSel, sample=None, sampleCfg=None): 
        # Some imports #
        from bamboo.analysisutils import forceDefine

        era = sampleCfg['era']
        # Get pile-up configs #
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
        isMC = self.isMC(sample)
        plots = []

        forceDefine(t._Muon.calcProd, noSel)

        #############################################################################
        ################################  Muons #####################################
        #############################################################################
        # Wp // 2016- 2017 -2018 : Muon_mediumId   // https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Muon_Isolation
        muonsByPt = op.sort(t.Muon, lambda mu : -mu.p4.Pt())
        muons = op.select(muonsByPt, lambda mu : op.AND(mu.p4.Pt() > 10., op.abs(mu.p4.Eta()) < 2.4, mu.tightId, mu.pfRelIso04_all<0.15))
      
        # Scalefactors #
        #if era=="2016":
        #    doubleMuTrigSF = get_scalefactor("dilepton", ("doubleMuLeg_HHMoriond17_2016"), systName="mumutrig")    
        #    muMediumIDSF = get_scalefactor("lepton", ("muon_{0}_{1}".format(era, sfTag), "id_medium"), combine="weight", systName="muid")
        #    muMediumISOSF = get_scalefactor("lepton", ("muon_{0}_{1}".format(era, sfTag), "iso_tight_id_medium"), combine="weight", systName="muiso")
        #    TrkIDSF = get_scalefactor("lepton", ("muon_{0}_{1}".format(era, sfTag), "highpt"), combine="weight")
        #    TrkISOSF = get_scalefactor("lepton", ("muon_{0}_{1}".format(era, sfTag), "isotrk_loose_idtrk_tightidandipcut"), combine="weight")
        #else:
        #    muMediumIDSF = get_scalefactor("lepton", ("muon_{0}_{1}".format(era, sfTag), "id_medium"), systName="muid")
        #    muMediumISOSF = get_scalefactor("lepton", ("muon_{0}_{1}".format(era, sfTag), "iso_tight_id_medium"), systName="muiso") 

        #############################################################################
        #############################  Electrons  ###################################
        #############################################################################
        #Wp  // 2016: Electron_cutBased_Sum16==3  -> medium     // 2017 -2018  : Electron_cutBased ==3   --> medium ( Fall17_V2)
        # asking for electrons to be in the Barrel region with dz<1mm & dxy< 0.5mm   //   Endcap region dz<2mm & dxy< 0.5mm 
        electronsByPt = op.sort(t.Electron, lambda ele : -ele.p4.Pt())
        electrons = op.select(electronsByPt, lambda ele : op.AND(ele.p4.Pt() > 15., op.abs(ele.p4.Eta()) < 2.5 , ele.cutBased>=3 )) # //cut-based ID Fall17 V2 the recomended one from POG for the FullRunII

        # Scalefactors #
        #elMediumIDSF = get_scalefactor("lepton", ("electron_{0}_{1}".format(era,sfTag), "id_medium"), systName="elid")
        #doubleEleTrigSF = get_scalefactor("dilepton", ("doubleEleLeg_HHMoriond17_2016"), systName="eleltrig")     

        #elemuTrigSF = get_scalefactor("dilepton", ("elemuLeg_HHMoriond17_2016"), systName="elmutrig")
        #mueleTrigSF = get_scalefactor("dilepton", ("mueleLeg_HHMoriond17_2016"), systName="mueltrig")

        OsElEl = op.combine(electrons, N=2, pred=lambda el1,el2 : op.AND(el1.charge != el2.charge , el1.p4.Pt() > 25, el2.p4.Pt() > 15 ))
        OsMuMu = op.combine(muons, N=2, pred=lambda mu1,mu2 : op.AND(mu1.charge != mu2.charge , mu1.p4.Pt() > 25, mu2.p4.Pt() > 15 ))
        OsElMu = op.combine((electrons, muons), pred=lambda el,mu : op.AND(el.charge != mu.charge ,el.p4.Pt() > 25, mu.p4.Pt() > 15 ))
        OsMuEl = op.combine((electrons, muons), pred=lambda el,mu : op.AND(el.charge != mu.charge ,el.p4.Pt() > 15, mu.p4.Pt() > 25 ))

        hasOsElEl = noSel.refine("hasOsElEl",cut=[ op.rng_len(OsElEl) >= 1 ])
        hasOsMuMu = noSel.refine("hasOsMuMu",cut=[ op.rng_len(OsMuMu) >= 1 ])
        hasOsElMu = noSel.refine("hasOsElMu",cut=[ op.rng_len(OsElMu) >= 1 ])
        hasOsMuEl = noSel.refine("hasOsMuEl",cut=[ op.rng_len(OsMuEl) >= 1 ])

        plots.append(Plot.make1D("ElEl_channel",op.rng_len(OsElEl),noSel,EquidistantBinning(10,0,10.),title='Number of dilepton events in ElEl channel',xTitle='N_{dilepton} (ElEl channel)'))
        plots.append(Plot.make1D("MuMu_channel",op.rng_len(OsMuMu),noSel,EquidistantBinning(10,0,10.),title='Number of dilepton events in MuMu channel',xTitle='N_{dilepton} (MuMu channel)'))
        plots.append(Plot.make1D("ElMu_channel",op.rng_len(OsElMu),noSel,EquidistantBinning(10,0,10.),title='Number of dilepton events in ElMu channel',xTitle='N_{dilepton} (ElMu channel)'))
        plots.append(Plot.make1D("MuEl_channel",op.rng_len(OsMuEl),noSel,EquidistantBinning(10,0,10.),title='Number of dilepton events in MuEl channel',xTitle='N_{dilepton} (MuEl channel)'))

        plots += makeDileptonPlots(self,
                                   sel         = hasOsElEl,
                                   dilepton    = OsElEl[0],
                                   suffix      = 'hasOsdilep',
                                   channel     = 'ElEl')
        plots += makeDileptonPlots(self,
                                   sel         = hasOsMuMu,
                                   dilepton    = OsMuMu[0],
                                   suffix      = 'hasOsdilep',
                                   channel     = 'MuMu')
        plots += makeDileptonPlots(self,
                                   sel         = hasOsElMu,
                                   dilepton    = OsElMu[0],
                                   suffix      = 'hasOsdilep',
                                   channel     = 'ElMu')
        plots += makeDileptonPlots(self,
                                   sel         = hasOsMuEl,
                                   dilepton    = OsMuEl[0],
                                   suffix      = 'hasOsdilep',
                                   channel     = 'MuEl')

   
        # Dilepton Z peak exclusion (charge already done in previous selection) # 
        lambda_mllLowerband = lambda dilep: op.in_range(12.,op.invariant_mass(dilep[0].p4, dilep[1].p4),80.)
        lambda_mllUpperband = lambda dilep: op.invariant_mass(dilep[0].p4, dilep[1].p4)>100.
        lambda_mllCut       = lambda dilep: op.OR(lambda_mllLowerband(dilep),lambda_mllUpperband(dilep))

        hasOsElElOutZ = hasOsElEl.refine("hasOsElElOutZ",cut=[lambda_mllCut(OsElEl[0])])
        hasOsMuMuOutZ = hasOsMuMu.refine("hasOsMuMuOutZ",cut=[lambda_mllCut(OsMuMu[0])])
        hasOsElMuOutZ = hasOsElMu.refine("hasOsElMuOutZ",cut=[lambda_mllCut(OsElMu[0])])
        hasOsMuElOutZ = hasOsMuEl.refine("hasOsMuElOutZ",cut=[lambda_mllCut(OsMuEl[0])])

        plots += makeDileptonPlots(self,
                                   sel         = hasOsElElOutZ,
                                   dilepton    = OsElEl[0],
                                   suffix      = 'hasOsdilep_OutZ',
                                   channel     = 'ElEl')
        plots += makeDileptonPlots(self,
                                   sel         = hasOsMuMuOutZ,
                                   dilepton    = OsMuMu[0],
                                   suffix      = 'hasOsdilep_OutZ',
                                   channel     = 'MuMu')
        plots += makeDileptonPlots(self,
                                   sel         = hasOsElMuOutZ,
                                   dilepton    = OsElMu[0],
                                   suffix      = 'hasOsdilep_OutZ',
                                   channel     = 'ElMu')
        plots += makeDileptonPlots(self,
                                   sel         = hasOsMuElOutZ,
                                   dilepton    = OsMuEl[0],
                                   suffix      = 'hasOsdilep_OutZ',
                                   channel     = 'MuEl')

        #############################################################################
        ################################  Jets  #####################################
        #############################################################################
        # select jets   // 2016 - 2017 - 2018   ( j.jetId &2) ->      tight jet ID
        jetsByPt = op.sort(t.Jet, lambda jet : -jet.p4.Pt())
        jetsSel = op.select(jetsByPt, lambda j : op.AND(j.p4.Pt() > 20., op.abs(j.p4.Eta())< 2.4, (j.jetId &2)))        # Jets = AK4 jets
        fatjetsByPt = op.sort(t.FatJet, lambda fatjet : -fatjet.p4.Pt())
        fatjetsSel = op.select(fatjetsByPt, lambda j : op.AND(j.p4.Pt() > 20., op.abs(j.p4.Eta())< 2.4, (j.jetId &2)))        # FatJets = AK8 jets
        # exclude from the jetsSel any jet that happens to include within its reconstruction cone a muon or an electron.
        jets = op.select(jetsSel, lambda j : op.AND(op.NOT(op.rng_any(electrons, lambda ele : op.deltaR(j.p4, ele.p4) < 0.3 )), op.NOT(op.rng_any(muons, lambda mu : op.deltaR(j.p4, mu.p4) < 0.3 ))))
        fatjets = op.select(fatjetsSel, lambda j : op.AND(op.NOT(op.rng_any(electrons, lambda ele : op.deltaR(j.p4, ele.p4) < 0.3 )), op.NOT(op.rng_any(muons, lambda mu : op.deltaR(j.p4, mu.p4) < 0.3 ))))
        # Boosted and resolved jets categories #
        if era == "2016": # Must check that subJet exists before looking at the btag
            lambda_boosted  = lambda fatjet : op.OR(op.AND(fatjet.subJet1._idx.result != -1,fatjet.subJet1.btagDeepB > 0.6321), op.AND(fatjet.subJet2._idx.result != -1,fatjet.subJet2.btagDeepB > 0.6321))
            lambda_resolved = lambda jet    : jet.btagDeepB > 0.6321
        elif era =="2017":
            lambda_boosted  = lambda fatjet : op.OR(op.AND(fatjet.subJet1._idx.result != -1,fatjet.subJet1.btagDeepB > 0.4941), op.AND(fatjet.subJet2._idx.result != -1,fatjet.subJet2.btagDeepB > 0.4941))
            lambda_resolved = lambda jet    : jet.btagDeepB > 0.4941
        elif era == "2018":
            lambda_boosted  = lambda fatjet : op.OR(op.AND(fatjet.subJet1._idx.result != -1,fatjet.subJet1.btagDeepB > 0.4184), op.AND(fatjet.subJet2._idx.result != -1,fatjet.subJet2.btagDeepB > 0.4184))
            lambda_resolved = lambda jet    : jet.btagDeepB > 0.4184

        # Select the bjets we want #
        bjetsResolved = op.select(jets, lambda_resolved)
        bjetsBoosted  = op.select(fatjets, lambda_boosted)
    
        # Define the boosted and Resolved (+exclusive) selections #
        hasBoostedJets = noSel.refine("hasBoostedJets",
                               cut=[op.rng_len(bjetsBoosted)>=1])
        hasNotBoostedJets = noSel.refine("hasNotBoostedJets",
                               cut=[op.rng_len(bjetsBoosted)==0])
        hasResolvedJets = noSel.refine("hasResolvedJets",
                               cut=[op.rng_len(jets)>=2,op.rng_len(bjetsResolved)>=1])
        hasNotResolvedJets = noSel.refine("hasNotResolvedJets",
                               cut=[op.OR(op.rng_len(jets)<=1,op.rng_len(bjetsResolved)==0)])
        hasBoostedAndResolvedJets = noSel.refine("hasBoostedAndResolvedJets", 
                               cut=[op.rng_len(bjetsBoosted)>=1,op.rng_len(jets)>=2,op.rng_len(bjetsResolved)>=1])
        hasNotBoostedAndResolvedJets = noSel.refine("hasNotBoostedAndResolvedJets", 
                               cut=[op.OR(op.rng_len(bjetsBoosted)==0,op.rng_len(jets)<=1,op.rng_len(bjetsResolved)==0)])
        hasExlusiveResolvedJets = noSel.refine("hasExlusiveResolved",
                               cut=[op.rng_len(jets)>=2, op.rng_len(bjetsResolved)>=1,op.rng_len(bjetsBoosted)==0])
        hasNotExlusiveResolvedJets = noSel.refine("hasNotExlusiveResolved",
                               cut=[op.OR(op.OR(op.rng_len(jets)<=1,op.rng_len(bjetsResolved)==0),op.AND(op.rng_len(bjetsBoosted)>=1,op.rng_len(jets)>=2,op.rng_len(bjetsResolved)>=1))])
        hasExlusiveBoostedJets = noSel.refine("hasExlusiveBoostedJets",
                               cut=[op.rng_len(bjetsBoosted)>=1,op.OR(op.rng_len(jets)<=1,op.rng_len(bjetsResolved)==0)])
        hasNotExlusiveBoostedJets = noSel.refine("hasNotExlusiveBoostedJets",
                               cut=[op.OR(op.rng_len(bjetsBoosted)==0,op.AND(op.rng_len(jets)>=2,op.rng_len(bjetsResolved)>=1))])

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
                                              hasExlusiveResolvedJets,
                                              EquidistantBinning(2,0.,2.),
                                              title='Passed Exclusive Resolved',
                                              xTitle='Passed Exclusive Resolved')
        FailedExclusiveResolved = Plot.make1D("FailedExclusiveResolved",
                                              op.c_int(0),
                                              hasNotExlusiveResolvedJets,
                                              EquidistantBinning(2,0.,2.),
                                              title='Failed Exclusive Resolved',
                                              xTitle='Failed Exclusive Resolved')
        plots.append(SummedPlot("ExclusiveResolvedCase",
                                [FailedExclusiveResolved,PassedExclusiveResolved],
                                xTitle="Exclusive Resolved selection"))

        # Passing Exclusive Boosted (Boosted AND NOT Resolved) #
        PassedExclusiveBoosted = Plot.make1D("PassedExclusiveBoosted",
                                             op.c_int(1),
                                             hasExlusiveBoostedJets,
                                             EquidistantBinning(2,0.,2.),
                                             title='Passed Exclusive Boosted',
                                             xTitle='Passed Exclusive Boosted')
        FailedExclusiveBoosted = Plot.make1D("FailedExclusiveBoosted",
                                             op.c_int(0),
                                             hasNotExlusiveBoostedJets,
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

        # Count number of boosted and resolved jets #
        plots.append(Plot.make1D("NBoostedJets",
                     op.rng_len(bjetsBoosted),
                     hasBoostedJets,
                     EquidistantBinning(5,0.,5.),
                     title='Number of boosted jets in boosted case',
                     xTitle='N boosted bjets'))
        plots.append(Plot.make1D("NResolvedJets",
                     op.rng_len(bjetsResolved),
                     hasExlusiveResolvedJets,
                     EquidistantBinning(5,0.,5.),
                     title='Number of resolved jets in exclusive resolved case',
                     xTitle='N resolved bjets'))

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

        # Plot jets quantities without the dilepton selections #
        plots += makeFatJetPlots(self,
                                 sel        = hasBoostedJets,
                                 fatjet     = bjetsBoosted[0],
                                 suffix     = "BoostedJets",
                                 channel    = "NoChannel")
        plots += makeJetsPlots(self,
                               sel          = hasResolvedJets,
                               jets         = bjetsResolved,
                               suffix       = "ResolvedJets",
                               channel      = "NoChannel")


        #############################################################################
        ##################### Jets + Dilepton combination ###########################
        #############################################################################
        # Combine dilepton and Resolved (Exclusive = NOT Boosted) selections #
        hasOsElElOutZResolvedJets = hasOsElElOutZ.refine("hasOsElElOutZResolvedJets", 
                                      cut=[op.rng_len(jets)>=2, op.rng_len(bjetsResolved)>=1,op.rng_len(bjetsBoosted)==0]) 
        hasOsMuMuOutZResolvedJets = hasOsMuMuOutZ.refine("hasOsMuMuOutZResolvedJets",  
                                      cut=[op.rng_len(jets)>=2, op.rng_len(bjetsResolved)>=1,op.rng_len(bjetsBoosted)==0]) 
        hasOsElMuOutZResolvedJets = hasOsElMuOutZ.refine("hasOsElMuOutZResolvedJets", 
                                      cut=[op.rng_len(jets)>=2, op.rng_len(bjetsResolved)>=1,op.rng_len(bjetsBoosted)==0])
        hasOsMuElOutZResolvedJets = hasOsMuElOutZ.refine("hasOsMuElOutZResolvedJets", 
                                      cut=[op.rng_len(jets)>=2, op.rng_len(bjetsResolved)>=1,op.rng_len(bjetsBoosted)==0])

        # Combine dilepton and Boosted selections #
        hasOsElElOutZBoostedJets = hasOsElElOutZ.refine("hasOsElElOutZBoostedJets", 
                                      cut=[op.rng_len(bjetsBoosted)>=1])
        hasOsMuMuOutZBoostedJets = hasOsMuMuOutZ.refine("hasOsMuMuOutZBoostedJets", 
                                      cut=[op.rng_len(bjetsBoosted)>=1])
        hasOsElMuOutZBoostedJets = hasOsElMuOutZ.refine("hasOsElMuOutZBoostedJets", 
                                      cut=[op.rng_len(bjetsBoosted)>=1])
        hasOsMuElOutZBoostedJets = hasOsMuElOutZ.refine("hasOsMuElOutZBoostedJets", 
                                      cut=[op.rng_len(bjetsBoosted)>=1])

        # Plot dilepton with OS, Z peak and Resolved jets selections #
        plots += makeDileptonPlots(self,
                                   sel         = hasOsElElOutZResolvedJets,
                                   dilepton    = OsElEl[0],
                                   suffix      = 'hasOsdilep_OutZ_ResolvedJets',
                                   channel     = 'ElEl')
        plots += makeDileptonPlots(self,
                                   sel         = hasOsMuMuOutZResolvedJets,
                                   dilepton    = OsMuMu[0],
                                   suffix      = 'hasOsdilep_OutZ_ResolvedJets',
                                   channel     = 'MuMu')
        plots += makeDileptonPlots(self,
                                   sel         = hasOsElMuOutZResolvedJets,
                                   dilepton    = OsElMu[0],
                                   suffix      = 'hasOsdilep_OutZ_ResolvedJets',
                                   channel     = 'ElMu')
        plots += makeDileptonPlots(self,
                                   sel         = hasOsMuElOutZResolvedJets,
                                   dilepton    = OsMuEl[0],
                                   suffix      = 'hasOsdilep_OutZ_ResolvedJets',
                                   channel     = 'MuEl')

        # Plot dilepton with OS dilepton, Z peak and Boosted jets selections #
        plots += makeDileptonPlots(self,
                                   sel         = hasOsElElOutZBoostedJets,
                                   dilepton    = OsElEl[0],
                                   suffix      = 'hasOsdilep_OutZ_BoostedJets',
                                   channel     = 'ElEl')
        plots += makeDileptonPlots(self,
                                   sel         = hasOsMuMuOutZBoostedJets,
                                   dilepton    = OsMuMu[0],
                                   suffix      = 'hasOsdilep_OutZ_BoostedJets',
                                   channel     = 'MuMu')
        plots += makeDileptonPlots(self,
                                   sel         = hasOsElMuOutZBoostedJets,
                                   dilepton    = OsElMu[0],
                                   suffix      = 'hasOsdilep_OutZ_BoostedJets',
                                   channel     = 'ElMu')
        plots += makeDileptonPlots(self,
                                   sel         = hasOsMuElOutZBoostedJets,
                                   dilepton    = OsMuEl[0],
                                   suffix      = 'hasOsdilep_OutZ_BoostedJets',
                                   channel     = 'MuEl')

    
        # Plotting the fatjet for OS dilepton, Z peak and Boosted jets selections #
        plots += makeFatJetPlots(self,
                                 sel       = hasOsElElOutZBoostedJets,
                                 fatjet    = bjetsBoosted[0],
                                 suffix    = "hasOsdilep_OutZ_BoostedJets",
                                 channel   = "ElEl")
        plots += makeFatJetPlots(self,
                                 sel       = hasOsMuMuOutZBoostedJets,
                                 fatjet    = bjetsBoosted[0],
                                 suffix    = "hasOsdilep_OutZ_BoostedJets",
                                 channel   = "MuMu")
        plots += makeFatJetPlots(self,
                                 sel       = hasOsElMuOutZBoostedJets,
                                 fatjet    = bjetsBoosted[0],
                                 suffix    = "hasOsdilep_OutZ_BoostedJets",
                                 channel   = "ElMu")
        plots += makeFatJetPlots(self,
                                 sel       = hasOsMuElOutZBoostedJets,
                                 fatjet    = bjetsBoosted[0],
                                 suffix    = "hasOsdilep_OutZ_BoostedJets",
                                 channel   = "MuEl")

        # Plotting the jets for OS dilepton, Z peak and Resolved jets selections #
        plots += makeJetsPlots(self,
                               sel         = hasOsElElOutZResolvedJets,
                               jets        = bjetsResolved,
                               suffix      = "hasOsdilep_OutZ_ResolvedJets",
                               channel     = "ElEl")
        plots += makeJetsPlots(self,
                               sel         = hasOsMuMuOutZResolvedJets,
                               jets        = bjetsResolved,
                               suffix      = "hasOsdilep_OutZ_ResolvedJets",
                               channel     = "MuMu")
        plots += makeJetsPlots(self,
                               sel         = hasOsElMuOutZResolvedJets,
                               jets        = bjetsResolved,
                               suffix      = "hasOsdilep_OutZ_ResolvedJets",
                               channel     = "ElMu")
        plots += makeJetsPlots(self,
                               sel         = hasOsMuElOutZResolvedJets,
                               jets        = bjetsResolved,
                               suffix      = "hasOsdilep_OutZ_ResolvedJets",
                               channel     = "MuEl")

       ## helper selection (OR) to make sure jet calculations are only done once
       #hasOSLL = noSel.refine("hasOSLL", cut=op.OR(*( hasOSLL_cmbRng(rng) for rng in osLLRng.values())))
       #forceDefine(t._Jet.calcProd, hasOSLL)
       #for varNm in t._Jet.available:
       #    forceDefine(t._Jet[varNm], hasOSLL)
        return plots

