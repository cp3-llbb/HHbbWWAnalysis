import sys

from bamboo.analysismodules import NanoAODHistoModule
from bamboo.analysisutils import makeMultiPrimaryDatasetTriggerSelection
from bamboo.scalefactors import binningVariables_nano

from bamboo import treefunctions as op
from bamboo import scalefactors
from bamboo.plots import Plot, EquidistantBinning, SummedPlot
from bamboo import treefunctions as op

class NanoHHTobbWW(NanoAODHistoModule):
    """ Example module: HH->bbW(->e/µ nu)W(->e/µ nu) histograms from NanoAOD """
    def __init__(self, args):
            super(NanoHHTobbWW, self).__init__(args)

    def prepareTree(self, tree, sample=None, sampleCfg=None, enableSystematics=None):
        # JEC's Recommendation for Full RunII: https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC
        # JER : -----------------------------: https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
        # Get base aguments #
        tree,noSel,be,lumiArgs = super(NanoHHTobbWW,self).prepareTree(tree, sample=sample, sampleCfg=sampleCfg)
        era = sampleCfg['era']
        #triggersPerPrimaryDataset = {}
        #from bamboo.analysisutils import configureJets ,configureRochesterCorrection

        # Check distributed option #
        #isNotWorker = (self.args.distributed != "worker") 

        # Rochester and JEC corrections (depends on era) #     
        #if era == "2016":
        #    configureRochesterCorrection(tree._Muon.calc,os.path.join(os.path.dirname(__file__), "data", "RoccoR2016.txt"))
        #    triggersPerPrimaryDataset = {
        #        "DoubleMuon" : [ tree.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL,
        #                         tree.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ,
        #                         tree.HLT.Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL,
        #                         tree.HLT.Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ ],
        #        "DoubleEG"   : [ tree.HLT.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ ],  # double electron (loosely isolated)
        #        "MuonEG"     : [ tree.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL ]
        #        }
        #    if self.isMC(sample) or "2016F" in sample or "2016G" in sample or "2016H" in sample:
        #        triggersPerPrimaryDataset["MuonEG"] += [ ## added from 2016F on
        #                tree.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,
        #                tree.HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ
        #                ]
        #    if "2016H" not in sample:
        #        triggersPerPrimaryDataset["MuonEG"] += [ ## removed for 2016H
        #                tree.HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL
        #                ]

        #    # Jet treatment #
        #    if self.isMC(sample):   # if MC -> needs spearing
        #        configureJets(tree, "Jet", "AK4PFchs",
        #            jec="Summer16_07Aug2017_V20_MC",
        #            smear="Summer16_25nsV1_MC",
        #            jesUncertaintySources=["Total"], mayWriteCache=isNotWorker)
        #    else:                   # If data -> extract info from config 
        #        if "2016B" in sample or "2016C" in sample or "2016D" in sample:
        #            configureJets(tree, "Jet", "AK4PFchs",
        #                jec="Summer16_07Aug2017BCD_V11_DATA", mayWriteCache=isNotWorker)
        #        elif "2016E" in sample or "2016F" in sample:
        #            configureJets(tree, "Jet", "AK4PFchs",
        #                jec="Summer16_07Aug2017EF_V11_DATA", mayWriteCache=isNotWorker)
        #        elif "2016G" in sample or "2016H" in sample:
        #            configureJets(tree, "Jet", "AK4PFchs",
        #                jec="Summer16_07Aug2017GH_V11_DATA", mayWriteCache=isNotWorker)

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
        #else:
        #    raise RuntimeError("Unknown era {0}".format(era))

        # Get weights #
        #if self.isMC(sample):
        #    noSel = noSel.refine("genWeight", weight=tree.genWeight, cut=op.OR(*chain.from_iterable(triggersPerPrimaryDataset.values())))
        #else:
        #    noSel = noSel.refine("withTrig", cut=makeMultiPrimaryDatasetTriggerSelection(sample, triggersPerPrimaryDataset))

        return tree,noSel,be,lumiArgs
        
    def makeDileptonPlots(self, sel, dilepton, suffix, channel):
        plots = []
        plots.append(Plot.make1D("%s_leadleptonPT_%s"%(channel,suffix), dilepton[0].p4.Pt(), sel, EquidistantBinning(100, 0., 150.), title="Transverse momentum of the leading lepton (channel %s)"%channel, xTitle= "P_{T}(leading lepton) [GeV]"))
        plots.append(Plot.make1D("%s_subleadleptonPT_%s"%(channel,suffix), dilepton[1].p4.Pt(), sel, EquidistantBinning(100, 0., 150.), title="Transverse momentum of the subleading lepton (channel %s)"%channel, xTitle= "P_{T}(leading lepton) [GeV]"))
        plots.append(Plot.make1D("%s_leadleptonEta_%s"%(channel,suffix), dilepton[0].p4.Eta(), sel, EquidistantBinning(20, -3., 3.), title="Pseudorapidity of the leading lepton (channel %s)"%channel, xTitle= "#eta (leading lepton) [GeV]"))
        plots.append(Plot.make1D("%s_subleadleptonEta_%s"%(channel,suffix), dilepton[1].p4.Eta(), sel, EquidistantBinning(20, -3., 3.), title="Pseudorapidity of the subleading lepton (channel %s)"%channel, xTitle= "#eta (subleading sublepton) [GeV]"))
        plots.append(Plot.make1D("%s_leadleptonPhi_%s"%(channel,suffix), dilepton[0].p4.Phi(), sel, EquidistantBinning(20, -3.2, 3.2), title="Azimutal angle of the leading lepton (channel %s)"%channel, xTitle= "#phi (leading lepton) [GeV]"))
        plots.append(Plot.make1D("%s_subleadleptonPhi_%s"%(channel,suffix), dilepton[1].p4.Phi(), sel, EquidistantBinning(20, -3.2, 3.2), title="Azimutal angle of the subleading lepton (channel %s)"%channel, xTitle= "#phi (subleading lepton) [GeV]"))
        plots.append(Plot.make1D("%s_leptonInvariantMass_%s"%(channel,suffix), op.invariant_mass(dilepton[0].p4,dilepton[1].p4), sel, EquidistantBinning(100, 0., 500.), title="Dilepton invariant mass (channel %s)"%channel, xTitle= "Invariant mass [GeV]"))

        return plots

    def makeJetsPlots(self, sel, jets, suffix, channel):
        plots = []

        # bjet 1 plots (always present by selection) #
        plots.append(Plot.make1D("%s_leadjetPT_%s"%(channel,suffix),jets[0].p4.Pt(),sel,EquidistantBinning(100,0,300.),title='Transverse momentum of the leading jet',xTitle="P_{T}(leading jet) [GeV]"))
        plots.append(Plot.make1D("%s_leadjetEta_%s"%(channel,suffix),jets[0].p4.Eta(),sel,EquidistantBinning(20,-3.,3.),title='Transverse momentum of the leading jet',xTitle="#eta(leading jet) [GeV]"))
        plots.append(Plot.make1D("%s_leadjetPhi_%s"%(channel,suffix),jets[0].p4.Phi(),sel,EquidistantBinning(20,-3.2,3.2),title='Azimutal angle of the leading jet',xTitle="#phi(leading jet) [GeV]"))

        # bjet 2 plots (not necessarily present) #
        hasAtLeastTwoBJets = sel.refine("hasAtLeastTwoBJets_%s_%s"%(channel,suffix),cut=[op.rng_len(jets)>=2]) 
        hasOnlyOneBJet = sel.refine("hasAtLeast2BJets_%s_%s"%(channel,suffix),cut=[op.rng_len(jets)==1]) 

        plot_hasTwoBJets_PT = Plot.make1D("%s_hasSubleadjetPT_%s"%(channel,suffix),jets[1].p4.Pt(),hasAtLeastTwoBJets,EquidistantBinning(100,-10.,300.),title='Transverse momentum of the subleading jet',xTitle="P_{T}(subleading jet) [GeV]")
        plot_notTwoBJets_PT = Plot.make1D("%s_notSubleadjetPT_%s"%(channel,suffix),op.c_float(-10.),hasOnlyOneBJet,EquidistantBinning(100,-10.,300.),title='Transverse momentum of the subleading jet',xTitle="P_{T}(subleading jet) [GeV]")
        plots.append(SummedPlot("%s_subleadjetPT_%s"%(channel,suffix),[plot_hasTwoBJets_PT,plot_notTwoBJets_PT],xTitle="P_{T}(subleading jet) [GeV]"))

        plot_hasTwoBJets_Eta = Plot.make1D("%s_hasSubleadjetEta_%s"%(channel,suffix),jets[1].p4.Eta(),hasAtLeastTwoBJets,EquidistantBinning(20,-3.,3.),title='Pseudorapidity of the subleading jet',xTitle="#eta(subleading jet) [GeV]")
        plot_notTwoBJets_Eta = Plot.make1D("%s_notSubleadjetEta_%s"%(channel,suffix),op.c_float(-3.),hasOnlyOneBJet,EquidistantBinning(20,-3.,3.),title='Pseudorapidity of the subleading jet',xTitle="#eta(subleading jet) [GeV]")
        plots.append(SummedPlot("%s_subleadjetEta_%s"%(channel,suffix),[plot_hasTwoBJets_Eta,plot_notTwoBJets_Eta],xTitle="#eta(subleading jet) [GeV]"))

        plot_hasTwoBJets_Phi = Plot.make1D("%s_hasSubleadjetPhi_%s"%(channel,suffix),jets[1].p4.Phi(),hasAtLeastTwoBJets,EquidistantBinning(20,-3.2,3.2),title='Azimutal angle of the subleading jet',xTitle="#phi(subleading jet) [GeV]")
        plot_notTwoBJets_Phi = Plot.make1D("%s_notSubleadjetPhi_%s"%(channel,suffix),op.c_float(-3.2),hasOnlyOneBJet,EquidistantBinning(20,-3.2,3.2),title='Azimutal angle of the subleading jet',xTitle="#phi(subleading jet) [GeV]")
        plots.append(SummedPlot("%s_subleadjetPhi_%s"%(channel,suffix),[plot_hasTwoBJets_Phi,plot_notTwoBJets_Phi],xTitle="#phi(subleading jet) [GeV]"))

        plot_hasTwoBJets_invMass = Plot.make1D("%s_hasTwoBJets_invMass_%s"%(channel,suffix),op.invariant_mass(jets[0].p4,jets[1].p4),hasAtLeastTwoBJets,EquidistantBinning(100,-10,500.),title="Dijet invariant mass (channel %s)"%channel, xTitle= "Invariant mass [GeV]")
        plot_notTwoBJets_invMass = Plot.make1D("%s_notTwoBJets_invMass_%s"%(channel,suffix),op.c_float(-10.),hasAtLeastTwoBJets,EquidistantBinning(100,-10,500.),title="Dijet invariant mass (channel %s)"%channel, xTitle= "Invariant mass [GeV]")
        plots.append(SummedPlot("%s_dijetInvariantMass_%s"%(channel,suffix),[plot_hasTwoBJets_invMass,plot_notTwoBJets_invMass],xTitle="Invariant mass [GeV]"))
        
        return plots

    def makeFatJetPlots(self, sel, fatjet, suffix, channel):
        plots = []
        lambda_subjet1 = lambda fatjet : fatjet.subJet1._idx.result != -1
        lambda_subjet2 = lambda fatjet : fatjet.subJet2._idx.result != -1
        lambda_notSubjet1 = lambda fatjet : fatjet.subJet1._idx.result == -1
        lambda_notSubjet2 = lambda fatjet : fatjet.subJet2._idx.result == -1
        has_subjet1 = sel.refine("has_subjet1_%s_%s"%(channel,suffix), cut=[lambda_subjet1(fatjet)])
        has_subjet2 = sel.refine("has_subjet2_%s_%s"%(channel,suffix), cut=[lambda_subjet2(fatjet)])
        has_notSubjet1 = sel.refine("has_notSubjet1_%s_%s"%(channel,suffix), cut=[lambda_notSubjet1(fatjet)])
        has_notSubjet2 = sel.refine("has_notSubjet2_%s_%s"%(channel,suffix), cut=[lambda_notSubjet2(fatjet)])
        has_bothSubjets = sel.refine("has_bothSubjets_%s_%s"%(channel,suffix), cut=[lambda_subjet1(fatjet),lambda_subjet2(fatjet)])
        has_notBothSubjets = sel.refine("has_notBothSubjets_%s_%s"%(channel,suffix), cut=[op.OR(lambda_notSubjet1(fatjet),lambda_notSubjet2(fatjet))])
        
        # PT Plot #
        plot_hasSubjet1_Pt = Plot.make1D("%s_fatjet_hasSubjet1PT_%s"%(channel,suffix),fatjet.subJet1.p4.Pt(),has_subjet1,EquidistantBinning(100,-10,500.),title='Transverse momentum of the first subjet',xTitle="P_{T} (first subjet) [GeV]")
        plot_notSubjet1_Pt = Plot.make1D("%s_fatjet_notSubjet1PT_%s"%(channel,suffix),op.c_float(-10.),has_notSubjet1,EquidistantBinning(100,-10,500.),title='Transverse momentum of the first subjet',xTitle="P_{T} (first subjet) [GeV]")
        plot_hasSubjet2_Pt = Plot.make1D("%s_fatjet_hasSubjet2PT_%s"%(channel,suffix),fatjet.subJet2.p4.Pt(),has_subjet2,EquidistantBinning(100,-10,500.),title='Transverse momentum of the second subjet',xTitle="P_{T} (second subjet) [GeV]")
        plot_notSubjet2_Pt = Plot.make1D("%s_fatjet_notSubjet2PT_%s"%(channel,suffix),op.c_float(-10.),has_notSubjet2,EquidistantBinning(100,-10,500.),title='Transverse momentum of the second subjet',xTitle="P_{T} (second subjet) [GeV]")

        plots.append(SummedPlot("%s_fatjet_subjet1PT_%s"%(channel,suffix),[plot_hasSubjet1_Pt,plot_notSubjet1_Pt],xTitle="P_{T} (first subjet) [GeV]"))
        plots.append(SummedPlot("%s_fatjet_subjet2PT_%s"%(channel,suffix),[plot_hasSubjet2_Pt,plot_notSubjet2_Pt],xTitle="P_{T} (second subjet) [GeV]"))

        # Eta plot #
        plot_hasSubjet1_Eta = Plot.make1D("%s_fatjet_hasSubjet1Eta_%s"%(channel,suffix),fatjet.subJet1.p4.Eta(),has_subjet1,EquidistantBinning(20,-3.,3.),title='Pseudorapidity of the first subjet',xTitle="#eta (first subjet) [GeV]")
        plot_notSubjet1_Eta = Plot.make1D("%s_fatjet_notSubjet1Eta_%s"%(channel,suffix),op.c_float(-3.),has_notSubjet1,EquidistantBinning(20,-3.,3.),title='Pseudorapidity of the first subjet',xTitle="#eta (first subjet) [GeV]")
        plot_hasSubjet2_Eta = Plot.make1D("%s_fatjet_hasSubjet2Eta_%s"%(channel,suffix),fatjet.subJet2.p4.Eta(),has_subjet2,EquidistantBinning(20,-3.,3.),title='Pseudorapidity of the second subjet',xTitle="#eta (second subjet) [GeV]")
        plot_notSubjet2_Eta = Plot.make1D("%s_fatjet_notSubjet2Eta_%s"%(channel,suffix),op.c_float(-3.),has_notSubjet2,EquidistantBinning(20,-3.,3.),title='Pseudorapidity of the second subjet',xTitle="#eta (second subjet) [GeV]")

        plots.append(SummedPlot("%s_fatjet_subjet1Eta_%s"%(channel,suffix),[plot_hasSubjet1_Eta,plot_notSubjet1_Eta],xTitle="#eta (first subjet) [GeV]"))
        plots.append(SummedPlot("%s_fatjet_subjet2Eta_%s"%(channel,suffix),[plot_hasSubjet2_Eta,plot_notSubjet2_Eta],xTitle="#eta (second subjet) [GeV]"))

        # Phi Plot #
        plot_hasSubjet1_Phi = Plot.make1D("%s_fatjet_hasSubjet1Phi_%s"%(channel,suffix),fatjet.subJet1.p4.Phi(),has_subjet1,EquidistantBinning(20,-3.2,3.2),title='Azimutal angle of the first subjet',xTitle="#phi (first subjet) [GeV]")
        plot_notSubjet1_Phi = Plot.make1D("%s_fatjet_notSubjet1Phi_%s"%(channel,suffix),op.c_float(-3.2),has_notSubjet1,EquidistantBinning(20,-3.2,3.2),title='Azimutal angle of the first subjet',xTitle="#phi (first subjet) [GeV]")
        plot_hasSubjet2_Phi = Plot.make1D("%s_fatjet_hasSubjet2Phi_%s"%(channel,suffix),fatjet.subJet2.p4.Phi(),has_subjet2,EquidistantBinning(20,-3.2,3.2),title='Azimutal angle of the second subjet',xTitle="#phi (second subjet) [GeV]")
        plot_notSubjet2_Phi = Plot.make1D("%s_fatjet_notSubjet2Phi_%s"%(channel,suffix),op.c_float(-3.2),has_notSubjet2,EquidistantBinning(20,-3.2,3.2),title='Azimutal angle of the second subjet',xTitle="#phi (second subjet) [GeV]")

        plots.append(SummedPlot("%s_fatjet_subjet1Phi_%s"%(channel,suffix),[plot_hasSubjet1_Phi,plot_notSubjet1_Phi],xTitle="#phi (first subjet) [GeV]"))
        plots.append(SummedPlot("%s_fatjet_subjet2Phi_%s"%(channel,suffix),[plot_hasSubjet2_Phi,plot_notSubjet2_Phi],xTitle="#phi (second subjet) [GeV]"))

        # Invariant mass #
        plot_hasBothSubjets_invMass = Plot.make1D("%s_fatjet_hasBothSubjets_invMass_%s"%(channel,suffix),op.invariant_mass(fatjet.subJet1.p4,fatjet.subJet2.p4), has_bothSubjets, EquidistantBinning(100, -10., 300.), title="Fatjet invariant mass (channel %s)"%channel, xTitle= "Invariant mass [GeV]")
        plot_hasNotBothSubjets_invMass = Plot.make1D("%s_fatjet_notBothSubjets_invMass_%s"%(channel,suffix),op.c_float(-10.), has_notBothSubjets, EquidistantBinning(100, -10., 300.), title="Fatjet invariant mass (channel %s)"%channel, xTitle= "Invariant mass [GeV]")

        plots.append(SummedPlot("%s_fatjet_invMass_%s"%(channel,suffix),[plot_hasBothSubjets_invMass,plot_hasNotBothSubjets_invMass],xTitle="Invariant mass [GeV]"))

        return plots
        

     

    def definePlots(self, t, noSel, sample=None, sampleCfg=None): 
        # Some imports #
        from bamboo.analysisutils import forceDefine

        era = sampleCfg['era']
        # Get pile-up configs #
        #puWeightsFile = None
        #if era == "2016":
        #    sfTag="94X"
        #    puWeightsFile = os.path.join(os.path.dirname(__file__), "data", "puweights2016.json")
        #elif era == "2017":
        #    sfTag="94X"     
        #    puWeightsFile = os.path.join(os.path.dirname(__file__), "data", "puweights2017.json")
        #elif era == "2018":
        #    sfTag="102X"
        #    puWeightsFile = os.path.join(os.path.dirname(__file__), "data", "puweights2018.json")
        #if self.isMC(sample) and puWeightsFile is not None:
        #    from bamboo.analysisutils import makePileupWeight
        #    noSel = noSel.refine("puWeight", weight=makePileupWeight(puWeightsFile, t.Pileup_nTrueInt, systName="pileup"))
        isMC = self.isMC(sample)
        plots = []

        #forceDefine(t._Muon.calcProd, noSel)

        #############################  Muons #####################################
        # Wp // 2016- 2017 -2018 : Muon_mediumId   // https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Muon_Isolation
        muonsByPt = op.sort(t.Muon, lambda mu : -mu.p4.Pt())
        muons = op.select(muonsByPt, lambda mu : op.AND(mu.p4.Pt() > 10., op.abs(mu.p4.Eta()) < 2.4, mu.tightId, mu.pfRelIso04_all<0.15))
      
        # Scalefactors #
        #if era=="2016":
        #    doubleMuTrigSF = get_scalefactor("dilepton", ("doubleMuLeg_HHMoriond17_2016"), systName="mumutrig")    
        #    muMediumIDSF = get_scalefactor("lepton", ("muon_{0}_{1}".format(era, sfTag), "id_medium"), combine="weight", systName="muid")
        #    muMediumISOSF = get_scalefactor("lepton", ("muon_{0}_{1}".format(era, sfTag), "iso_tight_id_medium"), combine="weight", systName="muiso")
        #    #TrkIDSF = get_scalefactor("lepton", ("muon_{0}_{1}".format(era, sfTag), "highpt"), combine="weight")
        #    #TrkISOSF = get_scalefactor("lepton", ("muon_{0}_{1}".format(era, sfTag), "isotrk_loose_idtrk_tightidandipcut"), combine="weight")
        #else:
        #    muMediumIDSF = get_scalefactor("lepton", ("muon_{0}_{1}".format(era, sfTag), "id_medium"), systName="muid")
        #    muMediumISOSF = get_scalefactor("lepton", ("muon_{0}_{1}".format(era, sfTag), "iso_tight_id_medium"), systName="muiso") 

        #############################  Electrons  #####################################
        #Wp  // 2016: Electron_cutBased_Sum16==3  -> medium     // 2017 -2018  : Electron_cutBased ==3   --> medium ( Fall17_V2)
        # asking for electrons to be in the Barrel region with dz<1mm & dxy< 0.5mm   //   Endcap region dz<2mm & dxy< 0.5mm 
        electronsByPt = op.sort(t.Electron, lambda ele : -ele.p4.Pt())
        electrons = op.select(electronsByPt, lambda ele : op.AND(ele.p4.Pt() > 15., op.abs(ele.p4.Eta()) < 2.5 , ele.cutBased>=3 )) # //cut-based ID Fall17 V2 the recomended one from POG for the FullRunII

        # Scalefactors #
        #elMediumIDSF = get_scalefactor("lepton", ("electron_{0}_{1}".format(era,sfTag), "id_medium"), systName="elid")
        #doubleEleTrigSF = get_scalefactor("dilepton", ("doubleEleLeg_HHMoriond17_2016"), systName="eleltrig")     

        #elemuTrigSF = get_scalefactor("dilepton", ("elemuLeg_HHMoriond17_2016"), systName="elmutrig")
        #mueleTrigSF = get_scalefactor("dilepton", ("mueleLeg_HHMoriond17_2016"), systName="mueltrig")

#------------- For debugging purpose; work on electronMVA to check if the e-mu & mu-e channel will screwd up in the same way as working only with Electron
        #working in the inner barrel  |eta|< 0.8
        #electronsIB = op.select(t.Electron, lambda ele : op.AND(ele.p4.Pt() > 10., op.abs(ele.p4.Eta()) < 0.8, ele.mvaFall17V2Iso_WPL))
        #OsElEl_IB = op.combine(electronsIB, N=2, pred=lambda ele1,ele2 : ele1.charge != ele2.charge)
        #firstoselel_IB = OsElEl_IB[0]
        #hasOsElEl_IB = noSel.refine("twoosElectronsIB", cut=[op.rng_len(OsElEl_IB) >= 1, firstoselel_IB[0].p4.Pt() > 10. , op.in_range(70, op.invariant_mass(firstoselel_IB[0].p4, firstoselel_IB[1].p4), 110.)], weight=([elMediumIDSF(firstoselel_IB[0]), elMediumIDSF(firstoselel_IB[1]), doubleEleTrigSF(firstoselel_IB) ]if isMC else None))
        ## working in the outer barrel  0.8< eta <1.44
        #electronsOB = op.select(t.Electron, lambda ele : op.AND(ele.p4.Pt() > 10., op.in_range(0.8,ele.p4.Eta(),1.44), ele.mvaFall17V2Iso_WPL))
        #OsElEl_OB = op.combine(electronsOB, N=2, pred=lambda ele1,ele2 : ele1.charge != ele2.charge)
        #firstoselel_OB = OsElEl_OB[0]
        #hasOsElEl_OB = noSel.refine("twoosElectronsOB", cut=[op.rng_len(OsElEl_OB) >= 1, firstoselel_OB[0].p4.Pt() > 10. , op.in_range(70, op.invariant_mass(firstoselel_OB[0].p4, firstoselel_OB[1].p4), 110.)], weight=([elMediumIDSF(firstoselel_OB[0]), elMediumIDSF(firstoselel_OB[1]), doubleEleTrigSF(firstoselel_OB) ]if isMC else None))
        ## working in the Endcap 1.44< eta < 1.57
        #electronsEC = op.select(t.Electron, lambda ele : op.AND(ele.p4.Pt() > 10., op.in_range(1.44,ele.p4.Eta(),1.57), ele.mvaFall17V2Iso_WPL))
        #OsElEl_EC = op.combine(electronsEC, N=2, pred=lambda ele1,ele2 : ele1.charge != ele2.charge)
        #firstoselel_EC = OsElEl_EC[0]
        #hasOsElEl_EC = noSel.refine("twoosElectronsEC", cut=[op.rng_len(OsElEl_EC) >= 1, firstoselel_EC[0].p4.Pt() > 10., op.in_range(70, op.invariant_mass(firstoselel_EC[0].p4, firstoselel_EC[1].p4), 110.)], weight=([elMediumIDSF(firstoselel_EC[0]), elMediumIDSF(firstoselel_EC[1]), doubleEleTrigSF(firstoselel_EC) ]if isMC else None))
        
        # Combining electron and muon into a pair #
        #OsElEl = op.combine(electrons, N=2, pred=lambda el1,el2 : op.AND(el1.charge != el2.charge , el1.p4.Pt() > 25, el2.p4.Pt() > 15 , el1.p4.Pt() > el2.p4.Pt()))
        #OsMuMu = op.combine(muons, N=2, pred=lambda mu1,mu2 : op.AND(mu1.charge != mu2.charge , mu1.p4.Pt() > 20, mu2.p4.Pt() > 10, mu1.p4.Pt() > mu2.p4.Pt() ))
        #OsElMu = op.combine((electrons, muons), pred=lambda el,mu : op.AND(el.charge != mu.charge ,el.p4.Pt() > 25, mu.p4.Pt() > 10, el.p4.Pt() > mu.p4.Pt() ))
        #OsMuEl = op.combine((electrons, muons), pred=lambda el,mu : op.AND(el.charge != mu.charge ,el.p4.Pt() > 15, mu.p4.Pt() > 25, el.p4.Pt() < mu.p4.Pt() ))

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

        plots += self.makeDileptonPlots(
                                            sel = hasOsElEl,
                                            dilepton = OsElEl[0],
                                            suffix = 'hasOsdilep',
                                            channel = 'ElEl'
                                       )
        plots += self.makeDileptonPlots(
                                            sel = hasOsMuMu,
                                            dilepton = OsMuMu[0],
                                            suffix = 'hasOsdilep',
                                            channel = 'MuMu'
                                       )
        plots += self.makeDileptonPlots(
                                            sel = hasOsElMu,
                                            dilepton = OsElMu[0],
                                            suffix = 'hasOsdilep',
                                            channel = 'ElMu'
                                       )
        plots += self.makeDileptonPlots(
                                            sel = hasOsMuEl,
                                            dilepton = OsMuEl[0],
                                            suffix = 'hasOsdilep',
                                            channel = 'MuEl'
                                       )

   
        # Dilepton Z peak exclusion (charge already done in previous selection) # 
        lambda_mllLowerband = lambda dilep: op.in_range(12.,op.invariant_mass(dilep[0].p4, dilep[1].p4),80.)
        lambda_mllUpperband = lambda dilep: op.invariant_mass(dilep[0].p4, dilep[1].p4)>100.
        lambda_mllCut       = lambda dilep: op.OR(lambda_mllLowerband(dilep),lambda_mllUpperband(dilep))

        hasOsElElOutZ = hasOsElEl.refine("hasOsElElOutZ",cut=[lambda_mllCut(OsElEl[0])])
        hasOsMuMuOutZ = hasOsMuMu.refine("hasOsMuMuOutZ",cut=[lambda_mllCut(OsMuMu[0])])
        hasOsElMuOutZ = hasOsElMu.refine("hasOsElMuOutZ",cut=[lambda_mllCut(OsElMu[0])])
        hasOsMuElOutZ = hasOsMuEl.refine("hasOsMuElOutZ",cut=[lambda_mllCut(OsMuEl[0])])

        plots += self.makeDileptonPlots(
                                            sel = hasOsElElOutZ,
                                            dilepton = OsElEl[0],
                                            suffix = 'hasOsdilep_OutZ',
                                            channel = 'ElEl'
                                       )
        plots += self.makeDileptonPlots(
                                            sel = hasOsMuMuOutZ,
                                            dilepton = OsMuMu[0],
                                            suffix = 'hasOsdilep_OutZ',
                                            channel = 'MuMu'
                                       )
        plots += self.makeDileptonPlots(
                                            sel = hasOsElMuOutZ,
                                            dilepton = OsElMu[0],
                                            suffix = 'hasOsdilep_OutZ',
                                            channel = 'ElMu'
                                       )
        plots += self.makeDileptonPlots(
                                            sel = hasOsMuElOutZ,
                                            dilepton = OsMuEl[0],
                                            suffix = 'hasOsdilep_OutZ',
                                            channel = 'MuEl'
                                       )

         #############################  Jets  #####################################
        # select jets   // 2016 - 2017 - 2018   ( j.jetId &2) ->      tight jet ID
        jetsByPt = op.sort(t.Jet, lambda jet : -jet.p4.Pt())
        jetsSel = op.select(jetsByPt, lambda j : op.AND(j.p4.Pt() > 20., op.abs(j.p4.Eta())< 2.4, (j.jetId &2)))        # Jets = AK4 jets
        fatjetsByPt = op.sort(t.FatJet, lambda fatjet : -fatjet.p4.Pt())
        fatjetsSel = op.select(fatjetsByPt, lambda j : op.AND(j.p4.Pt() > 20., op.abs(j.p4.Eta())< 2.4, (j.jetId &2)))        # FatJets = AK8 jets
        # exclude from the jetsSel any jet that happens to include within its reconstruction cone a muon or an electron.
        jets = op.select(jetsSel, lambda j : op.AND(op.NOT(op.rng_any(electrons, lambda ele : op.deltaR(j.p4, ele.p4) < 0.3 )), op.NOT(op.rng_any(muons, lambda mu : op.deltaR(j.p4, mu.p4) < 0.3 ))))
        fatjets = op.select(fatjetsSel, lambda j : op.AND(op.NOT(op.rng_any(electrons, lambda ele : op.deltaR(j.p4, ele.p4) < 0.3 )), op.NOT(op.rng_any(muons, lambda mu : op.deltaR(j.p4, mu.p4) < 0.3 ))))
        #N.B: DeepJet= DeepFlavour
        # ask for the jets to be a b-jets ## DeepCSV medium b-tag working point
        #if era == "2016":
        #    bjets = op.select(jets, lambda j : j.btagDeepB > 0.6321 ) 
        #    DeepB_discriVar = { "BTagDiscri": lambda j : j.btagDeepB }
        #    #deepBMediumSF = get_scalefactor("jet", ("btag_2016_94X", "DeepCSV_medium"), additionalVariables=DeepB_discriVar, systName="btagging2016")  
        #elif era =="2017":
        #    bjets = op.select(jets, lambda j : j.btagDeepB > 0.4941 ) 
        #    DeepB_discriVar = { "BTagDiscri": lambda j : j.btagDeepB }
        #    #deepBMediumSF = get_scalefactor("jet", ("btag_2017_94X", "DeepCSV_medium"), additionalVariables=DeepB_discriVar, systName="btagging2017")  
        #elif era == "2018":
        #    bjets = op.select(jets, lambda j : j.btagDeepB > 0.4184 )
        #    DeepB_discriVar = { "BTagDiscri": lambda j : j.btagDeepB }
        #    #deepBMediumSF = get_scalefactor("jet", ("btag_2018_102X", "DeepCSV_medium"), additionalVariables=DeepB_discriVar, systName="btagging2018")  
        #else:
        #    print("Could not find bjet eras")

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
        hasBoostedJets = noSel.refine("hasBoostedJets",cut=[op.rng_len(bjetsBoosted)>=1])
        hasNotBoostedJets = noSel.refine("hasNotBoostedJets",cut=[op.rng_len(bjetsBoosted)==0])
        hasResolvedJets = noSel.refine("hasResolvedJets",cut=[op.rng_len(jets)>=2,op.rng_len(bjetsResolved)>=1])
        hasNotResolvedJets = noSel.refine("hasNotResolvedJets",cut=[op.OR(op.rng_len(jets)<=1,op.rng_len(bjetsResolved)==0)])
        hasBoostedAndResolvedJets = noSel.refine("hasBoostedAndResolvedJets", cut=[op.rng_len(bjetsBoosted)>=1,op.rng_len(jets)>=2,op.rng_len(bjetsResolved)>=1])
        hasNotBoostedAndResolvedJets = noSel.refine("hasNotBoostedAndResolvedJets", cut=[op.OR(op.rng_len(bjetsBoosted)==0,op.rng_len(jets)<=1,op.rng_len(bjetsResolved)==0)])
        hasExlusiveResolvedJets = noSel.refine("hasExlusiveResolved",cut=[op.rng_len(jets)>=2, op.rng_len(bjetsResolved)>=1,op.rng_len(bjetsBoosted)==0])
        hasNotExlusiveResolvedJets = noSel.refine("hasNotExlusiveResolved",cut=[op.OR(op.OR(op.rng_len(jets)<=1,op.rng_len(bjetsResolved)==0),op.AND(op.rng_len(bjetsBoosted)>=1,op.rng_len(jets)>=2,op.rng_len(bjetsResolved)>=1))])
        hasExlusiveBoostedJets = noSel.refine("hasExlusiveBoostedJets",cut=[op.rng_len(bjetsBoosted)>=1,op.OR(op.rng_len(jets)<=1,op.rng_len(bjetsResolved)==0)])
        hasNotExlusiveBoostedJets = noSel.refine("hasNotExlusiveBoostedJets",cut=[op.OR(op.rng_len(bjetsBoosted)==0,op.AND(op.rng_len(jets)>=2,op.rng_len(bjetsResolved)>=1))])

        # Counting events from different selections for debugging #
        PassedBoosted = Plot.make1D("PassedBoosted",op.c_int(1),hasBoostedJets,EquidistantBinning(2,0.,2.),title='Passed Boosted',xTitle='Passed Boosted')
        FailedBoosted = Plot.make1D("FailedBoosted",op.c_int(0),hasNotBoostedJets,EquidistantBinning(2,0.,2.),title='Failed Boosted',xTitle='Failed Boosted')

        PassedResolved = Plot.make1D("PassedResolved",op.c_int(1),hasResolvedJets,EquidistantBinning(2,0.,2.),title='Passed Resolved',xTitle='Passed Resolved')
        FailedResolved = Plot.make1D("FailedResolved",op.c_int(0),hasNotResolvedJets,EquidistantBinning(2,0.,2.),title='Failed Resolved',xTitle='Failed Resolved')

        PassedExclusiveResolved = Plot.make1D("PassedExclusiveResolved",op.c_int(1),hasExlusiveResolvedJets,EquidistantBinning(2,0.,2.),title='Passed Exclusive Resolved',xTitle='Passed Exclusive Resolved')
        FailedExclusiveResolved = Plot.make1D("FailedExclusiveResolved",op.c_int(0),hasNotExlusiveResolvedJets,EquidistantBinning(2,0.,2.),title='Failed Exclusive Resolved',xTitle='Failed Exclusive Resolved')

        PassedExclusiveBoosted = Plot.make1D("PassedExclusiveBoosted",op.c_int(1),hasExlusiveBoostedJets,EquidistantBinning(2,0.,2.),title='Passed Exclusive Boosted',xTitle='Passed Exclusive Boosted')
        FailedExclusiveBoosted = Plot.make1D("FailedExclusiveBoosted",op.c_int(0),hasNotExlusiveBoostedJets,EquidistantBinning(2,0.,2.),title='Failed Exclusive Boosted',xTitle='Failed Exclusive Boosted')

        PassedBoth = Plot.make1D("PassedBoth",op.c_int(1),hasBoostedAndResolvedJets,EquidistantBinning(2,0.,2.),title='Passed Both Boosted and Resolved',xTitle='Passed Boosted and Resolved')
        FailedBoth = Plot.make1D("FailedBoth",op.c_int(0),hasNotBoostedAndResolvedJets,EquidistantBinning(2,0.,2.),title='Failed combination Boosted and Resolved',xTitle='Failed combination')

        plots.append(SummedPlot("BoostedCase",[FailedBoosted,PassedBoosted],xTitle="Boosted selection"))
        plots.append(SummedPlot("ResolvedCase",[FailedResolved,PassedResolved],xTitle="Resolved selection"))
        plots.append(SummedPlot("ExclusiveResolvedCase",[FailedExclusiveResolved,PassedExclusiveResolved],xTitle="Exclusive Resolved selection"))
        plots.append(SummedPlot("ExclusiveBoostedCase",[FailedExclusiveBoosted,PassedExclusiveBoosted],xTitle="Exclusive Boosted selection"))
        plots.append(SummedPlot("BoostedAndResolvedCase",[FailedBoth,PassedBoth],xTitle="Boosted and Resolved selection"))

        plots.append(Plot.make1D("NBoostedJets",op.rng_len(bjetsBoosted),hasBoostedJets,EquidistantBinning(5,0.,5.),title='Number of boosted jets in boosted case',xTitle='N boosted bjets'))
        plots.append(Plot.make1D("NResolvedJets",op.rng_len(bjetsResolved),hasExlusiveResolvedJets,EquidistantBinning(5,0.,5.),title='Number of resolved jets in exclusive resolved case',xTitle='N resolved bjets'))

        # Plot number of subjets in the boosted fatjets #
        lambda_noSubjet = lambda fatjet : op.AND(fatjet.subJet1._idx.result == -1, op.AND(fatjet.subJet2._idx.result == -1 ))
        lambda_oneSubjet = lambda fatjet : op.AND(fatjet.subJet1._idx.result != -1, op.AND(fatjet.subJet2._idx.result == -1 ))
        lambda_twoSubjet = lambda fatjet : op.AND(fatjet.subJet1._idx.result != -1, op.AND(fatjet.subJet2._idx.result != -1 ))
        hasNoSubjet = hasBoostedJets.refine("hasNoSubjet",cut=[lambda_noSubjet(bjetsBoosted[0])])
        hasOneSubjet = hasBoostedJets.refine("hasOneSubjet",cut=[lambda_oneSubjet(bjetsBoosted[0])])
        hasTwoSubjet = hasBoostedJets.refine("hasTwoSubjet",cut=[lambda_twoSubjet(bjetsBoosted[0])])
        plot_hasNoSubjet = Plot.make1D("plot_hasNoSubjet",op.c_int(0),hasNoSubjet,EquidistantBinning(3,0.,3.),title='Boosted jet without subjet')
        plot_hasOneSubjet = Plot.make1D("plot_hasOneSubjet",op.c_int(1),hasOneSubjet,EquidistantBinning(3,0.,3.),title='Boosted jet with one subjet')
        plot_hasTwoSubjet = Plot.make1D("plot_hasTwoSubjet",op.c_int(2),hasTwoSubjet,EquidistantBinning(3,0.,3.),title='Boosted jet with two subjets')
        plots.append(SummedPlot("NumberOfSubjets",[plot_hasNoSubjet,plot_hasOneSubjet,plot_hasTwoSubjet],xTitle="Number of subjets in boosted jet"))

        # Plot jets quantities without the dilepton selections #
        plots += self.makeFatJetPlots(
                                        sel = hasBoostedJets,
                                        fatjet = bjetsBoosted[0],
                                        suffix = "BoostedJets",
                                        channel = "NoChannel",
                                   )
        plots += self.makeJetsPlots(
                                        sel = hasResolvedJets,
                                        jets = bjetsResolved,
                                        suffix = "ResolvedJets",
                                        channel = "NoChannel",
                                   )


        # Combine dilepton and jet cuts #
        hasOsElElOutZResolvedJets = hasOsElElOutZ.refine("hasOsElElOutZResolvedJets", cut=[op.rng_len(jets)>=2, op.rng_len(bjetsResolved)>=1,op.rng_len(bjetsBoosted)==0])
        hasOsMuMuOutZResolvedJets = hasOsMuMuOutZ.refine("hasOsMuMuOutZResolvedJets", cut=[op.rng_len(jets)>=2, op.rng_len(bjetsResolved)>=1,op.rng_len(bjetsBoosted)==0])
        hasOsElMuOutZResolvedJets = hasOsElMuOutZ.refine("hasOsElMuOutZResolvedJets", cut=[op.rng_len(jets)>=2, op.rng_len(bjetsResolved)>=1,op.rng_len(bjetsBoosted)==0])
        hasOsMuElOutZResolvedJets = hasOsMuElOutZ.refine("hasOsMuElOutZResolvedJets", cut=[op.rng_len(jets)>=2, op.rng_len(bjetsResolved)>=1,op.rng_len(bjetsBoosted)==0])

        hasOsElElOutZBoostedJets = hasOsElElOutZ.refine("hasOsElElOutZBoostedJets", cut=[op.rng_len(bjetsBoosted)>=1])
        hasOsMuMuOutZBoostedJets = hasOsMuMuOutZ.refine("hasOsMuMuOutZBoostedJets", cut=[op.rng_len(bjetsBoosted)>=1])
        hasOsElMuOutZBoostedJets = hasOsElMuOutZ.refine("hasOsElMuOutZBoostedJets", cut=[op.rng_len(bjetsBoosted)>=1])
        hasOsMuElOutZBoostedJets = hasOsMuElOutZ.refine("hasOsMuElOutZBoostedJets", cut=[op.rng_len(bjetsBoosted)>=1])

        # Plot dilepton with OS, Z peak and jets selections #
        plots += self.makeDileptonPlots(
                                            sel = hasOsElElOutZResolvedJets,
                                            dilepton = OsElEl[0],
                                            suffix = 'hasOsdilep_OutZ_ResolvedJets',
                                            channel = 'ElEl'
                                       )
        plots += self.makeDileptonPlots(
                                            sel = hasOsMuMuOutZResolvedJets,
                                            dilepton = OsMuMu[0],
                                            suffix = 'hasOsdilep_OutZ_ResolvedJets',
                                            channel = 'MuMu'
                                       )
        plots += self.makeDileptonPlots(
                                            sel = hasOsElMuOutZResolvedJets,
                                            dilepton = OsElMu[0],
                                            suffix = 'hasOsdilep_OutZ_ResolvedJets',
                                            channel = 'ElMu'
                                       )
        plots += self.makeDileptonPlots(
                                            sel = hasOsMuElOutZResolvedJets,
                                            dilepton = OsMuEl[0],
                                            suffix = 'hasOsdilep_OutZ_ResolvedJets',
                                            channel = 'MuEl'
                                       )
        plots += self.makeDileptonPlots(
                                            sel = hasOsElElOutZBoostedJets,
                                            dilepton = OsElEl[0],
                                            suffix = 'hasOsdilep_OutZ_BoostedJets',
                                            channel = 'ElEl'
                                       )
        plots += self.makeDileptonPlots(
                                            sel = hasOsMuMuOutZBoostedJets,
                                            dilepton = OsMuMu[0],
                                            suffix = 'hasOsdilep_OutZ_BoostedJets',
                                            channel = 'MuMu'
                                       )
        plots += self.makeDileptonPlots(
                                            sel = hasOsElMuOutZBoostedJets,
                                            dilepton = OsElMu[0],
                                            suffix = 'hasOsdilep_OutZ_BoostedJets',
                                            channel = 'ElMu'
                                       )
        plots += self.makeDileptonPlots(
                                            sel = hasOsMuElOutZBoostedJets,
                                            dilepton = OsMuEl[0],
                                            suffix = 'hasOsdilep_OutZ_BoostedJets',
                                            channel = 'MuEl'
                                       )

    
        # Plotting the jets for OS, Z peak and jets selections #
        plots += self.makeFatJetPlots(
                                        sel = hasOsElElOutZBoostedJets,
                                        fatjet = bjetsBoosted[0],
                                        suffix = "hasOsdilep_OutZ_BoostedJets",
                                        channel = "ElEl",
                                   )
        plots += self.makeFatJetPlots(
                                        sel = hasOsMuMuOutZBoostedJets,
                                        fatjet = bjetsBoosted[0],
                                        suffix = "hasOsdilep_OutZ_BoostedJets",
                                        channel = "MuMu",
                                   )
        plots += self.makeFatJetPlots(
                                        sel = hasOsElMuOutZBoostedJets,
                                        fatjet = bjetsBoosted[0],
                                        suffix = "hasOsdilep_OutZ_BoostedJets",
                                        channel = "ElMu",
                                   )
        plots += self.makeFatJetPlots(
                                        sel = hasOsMuElOutZBoostedJets,
                                        fatjet = bjetsBoosted[0],
                                        suffix = "hasOsdilep_OutZ_BoostedJets",
                                        channel = "MuEl",
                                   )
        plots += self.makeJetsPlots(
                                        sel = hasOsElElOutZResolvedJets,
                                        jets = bjetsResolved,
                                        suffix = "hasOsdilep_OutZ_ResolvedJets",
                                        channel = "ElEl",
                                   )
        plots += self.makeJetsPlots(
                                        sel = hasOsMuMuOutZResolvedJets,
                                        jets = bjetsResolved,
                                        suffix = "hasOsdilep_OutZ_ResolvedJets",
                                        channel = "MuMu",
                                   )
        plots += self.makeJetsPlots(
                                        sel = hasOsElMuOutZResolvedJets,
                                        jets = bjetsResolved,
                                        suffix = "hasOsdilep_OutZ_ResolvedJets",
                                        channel = "ElMu",
                                   )
        plots += self.makeJetsPlots(
                                        sel = hasOsMuElOutZResolvedJets,
                                        jets = bjetsResolved,
                                        suffix = "hasOsdilep_OutZ_ResolvedJets",
                                        channel = "MuEl",
                                   )

                #hasOsDilepOutZ_cmbRng = lambda cmbRng : op.AND(op.rng_len(cmbRng) > 0, cmbRng[0][0].p4.Pt() > 25.) # The leading pT for the µµ channel should be above 20 Gev !

        ## helper selection (OR) to make sure jet calculations are only done once
        #hasOSLL = noSel.refine("hasOSLL", cut=op.OR(*( hasOSLL_cmbRng(rng) for rng in osLLRng.values())))
        #forceDefine(t._Jet.calcProd, hasOSLL)
        #for varNm in t._Jet.available:
        #    forceDefine(t._Jet[varNm], hasOSLL)

        #llSFs = {
        #    "MuMu" : (lambda ll : [ muMediumIDSF(ll[0]), muMediumIDSF(ll[1]), muMediumISOSF(ll[0]), muMediumISOSF(ll[1]), doubleMuTrigSF(ll) ]),#,TrkIDSF(ll), TrkISOSF(ll)
        #    "ElMu" : (lambda ll : [ elMediumIDSF(ll[0]), muMediumIDSF(ll[1]), muMediumISOSF(ll[1]), elemuTrigSF(ll) ]),
        #    "MuEl" : (lambda ll : [ muMediumIDSF(ll[0]), muMediumISOSF(ll[0]), elMediumIDSF(ll[1]), mueleTrigSF(ll) ]),
        #    "ElEl" : (lambda ll : [ elMediumIDSF(ll[0]), elMediumIDSF(ll[1]), doubleEleTrigSF(ll) ])
        #    }
        #categories = dict((catN, (catLLRng[0], hasOSLL.refine("hasOS{0}".format(catN), cut=hasOSLL_cmbRng(catLLRng), weight=(llSFs[catN](catLLRng[0]) if isMC else None)))) for catN, catLLRng in osLLRng.items())
        #for catN, (dilepton, catSel) in categories.items():

        #    plots.append(Plot.make1D("{0}_leadleptonPT".format(catN), dilepton[0].p4.Pt(), catSel, EquidistantBinning(100, 0., 600.), title="Transverse momentum of the leading lepton", xTitle= "pT(leading lepton)(GeV)"))
        #    plots.append(Plot.make1D("{0}_SubleadleptonPT".format(catN), dilepton[1].p4.Pt(), catSel, EquidistantBinning(100, 0., 600.), title="Transverse momentum of the Subleading lepton", xTitle= "pT(Sub-leading lepton)(GeV)"))
        #    plots.append(Plot.make1D("{0}_leadleptonETA".format(catN), dilepton[0].p4.eta(), catSel, EquidistantBinning(10, -3, 3), title="Pseudo-rapidity of the leading lepton", xTitle= "Eta(leading lepton)"))
        #    plots.append(Plot.make1D("{0}_subleadleptonETA".format(catN), dilepton[1].p4.eta(), catSel, EquidistantBinning(10, -3, 3), title="Pseudo-rapidity of the sub-leading lepton", xTitle= "Eta(Sub-leading lepton)"))
        #    plots.append(Plot.make1D("{0}_mll".format(catN), op.invariant_mass(dilepton[0].p4, dilepton[1].p4), catSel, EquidistantBinning(100, 70., 110.), title=" dilepton invariant mass", xTitle= "mll(GeV)"))
        #    plots.append(Plot.make1D("{0}_llpT".format(catN), (dilepton[0].p4.Pt() + dilepton[1].p4.Pt()), catSel, EquidistantBinning(100,0.,600.),title= "dilepton transverse momentum" , xTitle= "dilepton pT (GeV)"))
        #    plots.append(Plot.make1D("{0}_nVX".format(catN), t.PV.npvs, catSel, EquidistantBinning(10, 0., 60.), title="Distrubtion of the number of the reconstructed vertices", xTitle="nPVX"))
        #    plots.append(Plot.make2D("{0}_Electron_dzdxy".format(catN), (dilepton[0].dz ,dilepton[0].dxy),catSel, (EquidistantBinning(10, 0., 2.),EquidistantBinning(10, 0., 2.)) ,title="Electron in Barrel/EndCAP region" ))

        #    TwoJetsTwoLeptons=catSel.refine("twoJet{0}Sel".format(catN), cut=[ op.rng_len(jets) > 1 ]) 

        #    plots.append(Plot.make1D("{0}_leadJetPT".format(catN), jets[0].p4.Pt(), TwoJetsTwoLeptons, EquidistantBinning(100, 0., 600.), title="Transverse momentum of the leading jet PT", xTitle= "pT(leading Jet)(GeV)"))
        #    plots.append(Plot.make1D("{0}_subleadJetPT".format(catN), jets[1].p4.Pt(), TwoJetsTwoLeptons,EquidistantBinning(100, 0., 600.), title="Transverse momentum of the sub-leading jet PT", xTitle= "pT(Sub-leading Jet)(GeV)"))
        #    plots.append(Plot.make1D("{0}_leadJetETA".format(catN), jets[0].p4.eta(), TwoJetsTwoLeptons, EquidistantBinning(10, -3, 3), title="Pseudo-rapidity of the leading jet", xTitle="Eta(leading Jet"))
        #    plots.append(Plot.make1D("{0}_subleadJetETA".format(catN), jets[1].p4.eta(), TwoJetsTwoLeptons, EquidistantBinning(10, -3, 3), title="Pseudo-rapidity of the sub-leading jet", xTitle="Eta(Sub-leading Jet"))
        #    plots.append(Plot.make1D("{0}_nJets".format(catN), op.rng_len(jets), TwoJetsTwoLeptons, EquidistantBinning(5, 2., 6.), title="Number of jets", xTitle= "nbr Jets"))
        #    plots.append(Plot.make1D("{0}_jjpT".format(catN), (jets[0].p4.Pt() + jets[1].p4.Pt()), TwoJetsTwoLeptons, EquidistantBinning(100,0.,600.),title= "dijet transverse momentum" , xTitle= "dijet pT (GeV)"))
        #    plots.append(Plot.make1D("{0}_mjj".format(catN),op.invariant_mass(jets[0].p4, jets[1].p4) , TwoJetsTwoLeptons, EquidistantBinning(100, 0., 800.), title="mjj", xTitle= "mjj(GeV)"))
        #    plots.append(Plot.make1D("{0}_mlljj".format(catN), (dilepton[0].p4 +dilepton[1].p4+jets[0].p4+jets[1].p4).M(), TwoJetsTwoLeptons, EquidistantBinning(100, 0., 800.), title="mlljj", xTitle="mlljj(GeV)"))
        #    plots.append(Plot.make2D("{0}_mlljjvsmjj".format(catN), (op.invariant_mass(jets[0].p4, jets[1].p4),(dilepton[0].p4 +dilepton[1].p4+jets[0].p4+jets[1].p4).M()), TwoJetsTwoLeptons, (EquidistantBinning(1000, 0., 1000.), EquidistantBinning(1000, 0., 1000.)), title="mlljj vs mjj invariant mass"))

        #    # asking for bjets -----------------------           
        #    #TwoLeptonsTwoBjets = TwoJetsTwoLeptons.refine("TwoElTwoBjets{0}".format(catN), cut=[ op.rng_len(bJets) > 1 ],weight=([ deepBMediumSF(bJets[0]), deepBMediumSF(bJets[1]) ]if isMC else None))
        #    TwoLeptonsTwoBjets = catSel.refine("TwoElTwoBjets{0}".format(catN), cut=[ op.rng_len(bJets) > 1 ],weight=([ deepBMediumSF(bJets[0]), deepBMediumSF(bJets[1]) ]if isMC else None))
        #    plots.append(Plot.make1D("{0}_lead_BJetPT".format(catN), jets[0].p4.Pt(), TwoLeptonsTwoBjets, EquidistantBinning(100, 0., 600.), title="Transverse momentum of the leading bjet PT", xTitle= "pT(leading b-Jet)(GeV)"))
        #    plots.append(Plot.make1D("{0}_sublead_BJetPT".format(catN), jets[1].p4.Pt(), TwoLeptonsTwoBjets,EquidistantBinning(100, 0., 600.), title="Transverse momentum of the sub-leading bjet PT", xTitle= "pT(Sub-leading b-Jet)(GeV)"))
        #    plots.append(Plot.make1D("{0}_lead_BJetETA".format(catN), jets[0].p4.eta(), TwoLeptonsTwoBjets, EquidistantBinning(10, -3, 3), title="Pseudo-rapidity of the leading bjet", xTitle="Eta(leading b-Jet"))
        #    plots.append(Plot.make1D("{0}_sublead_BJetETA".format(catN), jets[1].p4.eta(), TwoLeptonsTwoBjets, EquidistantBinning(10, -3, 3), title="Pseudo-rapidity of the sub-leading bjet", xTitle="Eta(Sub-leading b-Jet"))
        #    plots.append(Plot.make1D("{0}_nBJets".format(catN), op.rng_len(jets), TwoLeptonsTwoBjets, EquidistantBinning(5, 2., 6.), title="Number of bjets", xTitle= "nbr b-Jets"))
        #    plots.append(Plot.make1D("{0}_twoBtaggedjetspT".format(catN), (jets[0].p4.Pt() + jets[1].p4.Pt()), TwoLeptonsTwoBjets, EquidistantBinning(100,0.,600.),title= "di-bjet transverse momentum" , xTitle= "di-bjet pT (GeV)"))
        #    plots.append(Plot.make1D("{0}_mjj_btagged".format(catN),op.invariant_mass(jets[0].p4, jets[1].p4) , TwoLeptonsTwoBjets, EquidistantBinning(100, 0., 800.), title="mass of two b-tagged jets", xTitle= "mjj(GeV)"))
        #    plots.append(Plot.make1D("{0}_mlljj_btagged".format(catN), (dilepton[0].p4 +dilepton[1].p4+jets[0].p4+jets[1].p4).M(),TwoLeptonsTwoBjets, EquidistantBinning(100, 0., 800.), title="invariant mass of 2 leptons two b-tagged jets", xTitle="mlljj(GeV)"))
        #    plots.append(Plot.make2D("{0}_mlljjvsmjj_btagged".format(catN), (op.invariant_mass(jets[0].p4, jets[1].p4),(dilepton[0].p4 +dilepton[1].p4+jets[0].p4+jets[1].p4).M()),TwoLeptonsTwoBjets, (EquidistantBinning(1000, 0., 1000.), EquidistantBinning(1000, 0., 1000.)), title="mlljj vs mjj invariant mass"))
        #    plots.append(Plot.make1D("{0}_mll_btagged".format(catN), op.invariant_mass(dilepton[0].p4, dilepton[1].p4), TwoLeptonsTwoBjets, EquidistantBinning(100, 70., 110.), title=" dilepton invariant mass", xTitle= "mll(GeV)"))

        return plots

