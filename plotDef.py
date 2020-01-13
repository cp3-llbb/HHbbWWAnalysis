from bamboo.plots import Plot, EquidistantBinning, SummedPlot
from bamboo import treefunctions as op

##########################  DILEPTON PLOT #################################
def makeDileptonPlots(self, sel, dilepton, suffix, channel):
    """
    Make dilepton basic plots
    sel         = refine selection 
    dilepton    = dilepton object (ElEl, MuEl, ElMu, MuMu)
    suffix      = string identifying the selecton 
    channel     = string identifying the channel of the dilepton (can be "NoChannel")
    """
    plots = []
    # PT plot #
    plots.append(Plot.make1D("%s_leadleptonPT_%s"%(channel,suffix), 
                             dilepton[0].p4.Pt(), 
                             sel, 
                             EquidistantBinning(100, 0., 150.), 
                             title="Transverse momentum of the leading lepton (channel %s)"%channel, 
                             xTitle= "P_{T}(leading lepton) [GeV]"))
    plots.append(Plot.make1D("%s_subleadleptonPT_%s"%(channel,suffix), 
                             dilepton[1].p4.Pt(), 
                             sel, 
                             EquidistantBinning(100, 0., 150.), 
                             title="Transverse momentum of the subleading lepton (channel %s)"%channel, 
                             xTitle= "P_{T}(leading lepton) [GeV]"))

    # Eta plot #
    plots.append(Plot.make1D("%s_leadleptonEta_%s"%(channel,suffix), 
                             dilepton[0].p4.Eta(), 
                             sel, 
                             EquidistantBinning(20, -3., 3.), 
                             title="Pseudorapidity of the leading lepton (channel %s)"%channel, 
                             xTitle= "#eta (leading lepton) [GeV]"))
    plots.append(Plot.make1D("%s_subleadleptonEta_%s"%(channel,suffix), 
                             dilepton[1].p4.Eta(), 
                             sel, 
                             EquidistantBinning(20, -3., 3.), 
                             title="Pseudorapidity of the subleading lepton (channel %s)"%channel, 
                             xTitle= "#eta (subleading sublepton) [GeV]"))

    # Phi plot #
    plots.append(Plot.make1D("%s_leadleptonPhi_%s"%(channel,suffix), 
                             dilepton[0].p4.Phi(), 
                             sel, 
                             EquidistantBinning(20, -3.2, 3.2), 
                             title="Azimutal angle of the leading lepton (channel %s)"%channel, 
                             xTitle= "#phi (leading lepton) [GeV]"))
    plots.append(Plot.make1D("%s_subleadleptonPhi_%s"%(channel,suffix), 
                             dilepton[1].p4.Phi(), 
                             sel, 
                             EquidistantBinning(20, -3.2, 3.2), 
                             title="Azimutal angle of the subleading lepton (channel %s)"%channel, 
                             xTitle= "#phi (subleading lepton) [GeV]"))

    # InvMass plot #
    plots.append(Plot.make1D("%s_leptonInvariantMass_%s"%(channel,suffix), 
                             op.invariant_mass(dilepton[0].p4,dilepton[1].p4), 
                             sel, 
                             EquidistantBinning(100, 0., 500.), 
                             title="Dilepton invariant mass (channel %s)"%channel, xTitle= "Invariant mass [GeV]"))

    return plots

##########################  JETS (AK4) PLOT #################################
def makeJetsPlots(self, sel, jets, suffix, channel):
    """
    Make basic plots
    sel         = refine selection 
    jets        = bjet container with len(jets)>=1 (we require at least one bjet) 
    suffix      = string identifying the selecton 
    channel     = string identifying the channel of the dilepton (can be "NoChannel")
    """
 
    plots = []

    # bjet 1 plots (always present by selection) #
    plots.append(Plot.make1D("%s_leadjetPT_%s"%(channel,suffix),
                             jets[0].p4.Pt(),
                             sel,
                             EquidistantBinning(100,0,300.),
                             title='Transverse momentum of the leading jet',
                             xTitle="P_{T}(leading jet) [GeV]"))
    plots.append(Plot.make1D("%s_leadjetEta_%s"%(channel,suffix),
                             jets[0].p4.Eta(),
                             sel,
                             EquidistantBinning(20,-3.,3.),
                             title='Transverse momentum of the leading jet',
                             xTitle="#eta(leading jet) [GeV]"))
    plots.append(Plot.make1D("%s_leadjetPhi_%s"%(channel,suffix),
                             jets[0].p4.Phi(),
                             sel,
                             EquidistantBinning(20,-3.2,3.2),
                             title='Azimutal angle of the leading jet',
                             xTitle="#phi(leading jet) [GeV]"))

    # bjet 2 plots (not necessarily present) #
    hasAtLeastTwoBJets = sel.refine("hasAtLeastTwoBJets_%s_%s"%(channel,suffix),
                                    cut=[op.rng_len(jets)>=2]) 
    hasOnlyOneBJet = sel.refine("hasAtLeast2BJets_%s_%s"%(channel,suffix), 
                                cut=[op.rng_len(jets)==1]) 
        # Fill a specific bin apart so that the information is not lost

    # PT plot #
    plot_hasTwoBJets_PT = Plot.make1D("%s_hasSubleadjetPT_%s"%(channel,suffix),
                                      jets[1].p4.Pt(),
                                      hasAtLeastTwoBJets,
                                      EquidistantBinning(100,-10.,300.),
                                      title='Transverse momentum of the subleading jet',
                                      xTitle="P_{T}(subleading jet) [GeV]")
    plot_notTwoBJets_PT = Plot.make1D("%s_notSubleadjetPT_%s"%(channel,suffix),
                                      op.c_float(-10.),
                                      hasOnlyOneBJet,
                                      EquidistantBinning(100,-10.,300.),
                                      title='Transverse momentum of the subleading jet',
                                      xTitle="P_{T}(subleading jet) [GeV]")
    plots.append(SummedPlot("%s_subleadjetPT_%s"%(channel,suffix),
                            [plot_hasTwoBJets_PT,plot_notTwoBJets_PT],
                            xTitle="P_{T}(subleading jet) [GeV]"))

    # Eta plot #
    plot_hasTwoBJets_Eta = Plot.make1D("%s_hasSubleadjetEta_%s"%(channel,suffix),
                                       jets[1].p4.Eta(),
                                       hasAtLeastTwoBJets,
                                       EquidistantBinning(20,-3.,3.),
                                       title='Pseudorapidity of the subleading jet',
                                       xTitle="#eta(subleading jet) [GeV]")
    plot_notTwoBJets_Eta = Plot.make1D("%s_notSubleadjetEta_%s"%(channel,suffix),
                                       op.c_float(-3.),
                                       hasOnlyOneBJet,
                                       EquidistantBinning(20,-3.,3.),
                                       title='Pseudorapidity of the subleading jet',
                                       xTitle="#eta(subleading jet) [GeV]")
    plots.append(SummedPlot("%s_subleadjetEta_%s"%(channel,suffix),
                            [plot_hasTwoBJets_Eta,plot_notTwoBJets_Eta],
                            xTitle="#eta(subleading jet) [GeV]"))

    # Phi plot #
    plot_hasTwoBJets_Phi = Plot.make1D("%s_hasSubleadjetPhi_%s"%(channel,suffix),
                                       jets[1].p4.Phi(),
                                       hasAtLeastTwoBJets,
                                       EquidistantBinning(20,-3.2,3.2),
                                       title='Azimutal angle of the subleading jet',
                                       xTitle="#phi(subleading jet) [GeV]")
    plot_notTwoBJets_Phi = Plot.make1D("%s_notSubleadjetPhi_%s"%(channel,suffix),
                                       op.c_float(-3.2),
                                       hasOnlyOneBJet,
                                       EquidistantBinning(20,-3.2,3.2),
                                       title='Azimutal angle of the subleading jet',
                                       xTitle="#phi(subleading jet) [GeV]")
    plots.append(SummedPlot("%s_subleadjetPhi_%s"%(channel,suffix),
                            [plot_hasTwoBJets_Phi,plot_notTwoBJets_Phi],
                            xTitle="#phi(subleading jet) [GeV]"))

    # InvMass plot # 
    plot_hasTwoBJets_invMass = Plot.make1D("%s_hasTwoBJets_invMass_%s"%(channel,suffix),
                                           op.invariant_mass(jets[0].p4,jets[1].p4),
                                           hasAtLeastTwoBJets,
                                           EquidistantBinning(100,-10,500.),
                                           title="Dijet invariant mass (channel %s)"%channel, 
                                           xTitle= "Invariant mass [GeV]")
    plot_notTwoBJets_invMass = Plot.make1D("%s_notTwoBJets_invMass_%s"%(channel,suffix),
                                           op.c_float(-10.),
                                           hasAtLeastTwoBJets,
                                           EquidistantBinning(100,-10,500.),
                                           title="Dijet invariant mass (channel %s)"%channel, 
                                           xTitle= "Invariant mass [GeV]")
    plots.append(SummedPlot("%s_dijetInvariantMass_%s"%(channel,suffix),
                            [plot_hasTwoBJets_invMass,plot_notTwoBJets_invMass],
                            xTitle="Invariant mass [GeV]"))
    
    return plots

##########################  FATJETS (AK8) PLOT #################################
def makeFatJetPlots(self, sel, fatjet, suffix, channel):
    """
    Make fatjet basic plots
    sel         = refine selection 
    fatjet      = fatjet object (AK8 jet) (can contain 0,1 or 2 subjets), where at least on subjet is btagged
    suffix      = string identifying the selecton 
    channel     = string identifying the channel of the dilepton (can be "NoChannel")
    """
 
    plots = []
    # Refine to check what and how many subjets were present #
    lambda_subjet1 = lambda fatjet : fatjet.subJet1._idx.result != -1
    lambda_subjet2 = lambda fatjet : fatjet.subJet2._idx.result != -1
    lambda_notSubjet1 = lambda fatjet : fatjet.subJet1._idx.result == -1
    lambda_notSubjet2 = lambda fatjet : fatjet.subJet2._idx.result == -1
    has_subjet1 = sel.refine("has_subjet1_%s_%s"%(channel,suffix), 
                             cut=[lambda_subjet1(fatjet)])
    has_subjet2 = sel.refine("has_subjet2_%s_%s"%(channel,suffix), 
                             cut=[lambda_subjet2(fatjet)])
    has_notSubjet1 = sel.refine("has_notSubjet1_%s_%s"%(channel,suffix), 
                             cut=[lambda_notSubjet1(fatjet)])
    has_notSubjet2 = sel.refine("has_notSubjet2_%s_%s"%(channel,suffix), 
                             cut=[lambda_notSubjet2(fatjet)])
    has_bothSubjets = sel.refine("has_bothSubjets_%s_%s"%(channel,suffix), 
                             cut=[lambda_subjet1(fatjet),lambda_subjet2(fatjet)])
    has_notBothSubjets = sel.refine("has_notBothSubjets_%s_%s"%(channel,suffix), 
                             cut=[op.OR(lambda_notSubjet1(fatjet),lambda_notSubjet2(fatjet))])
    
    # PT Plot #
    plot_hasSubjet1_Pt = Plot.make1D("%s_fatjet_hasSubjet1PT_%s"%(channel,suffix),
                                     fatjet.subJet1.p4.Pt(),
                                     has_subjet1,
                                     EquidistantBinning(100,-10,500.),
                                     title='Transverse momentum of the first subjet',
                                     xTitle="P_{T} (first subjet) [GeV]")
    plot_notSubjet1_Pt = Plot.make1D("%s_fatjet_notSubjet1PT_%s"%(channel,suffix),
                                     op.c_float(-10.),
                                     has_notSubjet1,
                                     EquidistantBinning(100,-10,500.),
                                     title='Transverse momentum of the first subjet',
                                     xTitle="P_{T} (first subjet) [GeV]")
    plots.append(SummedPlot("%s_fatjet_subjet1PT_%s"%(channel,suffix),
                            [plot_hasSubjet1_Pt,plot_notSubjet1_Pt],
                            xTitle="P_{T} (first subjet) [GeV]"))

    plot_hasSubjet2_Pt = Plot.make1D("%s_fatjet_hasSubjet2PT_%s"%(channel,suffix),
                                     fatjet.subJet2.p4.Pt(),
                                     has_subjet2,
                                     EquidistantBinning(100,-10,500.),
                                     title='Transverse momentum of the second subjet',
                                     xTitle="P_{T} (second subjet) [GeV]")
    plot_notSubjet2_Pt = Plot.make1D("%s_fatjet_notSubjet2PT_%s"%(channel,suffix),
                                     op.c_float(-10.),
                                     has_notSubjet2,
                                     EquidistantBinning(100,-10,500.),
                                     title='Transverse momentum of the second subjet',
                                     xTitle="P_{T} (second subjet) [GeV]")
    plots.append(SummedPlot("%s_fatjet_subjet2PT_%s"%(channel,suffix),
                            [plot_hasSubjet2_Pt,plot_notSubjet2_Pt],
                            xTitle="P_{T} (second subjet) [GeV]"))

    # Eta plot #
    plot_hasSubjet1_Eta = Plot.make1D("%s_fatjet_hasSubjet1Eta_%s"%(channel,suffix),
                                      fatjet.subJet1.p4.Eta(),
                                      has_subjet1,
                                      EquidistantBinning(20,-3.,3.),
                                      title='Pseudorapidity of the first subjet',
                                      xTitle="#eta (first subjet) [GeV]")
    plot_notSubjet1_Eta = Plot.make1D("%s_fatjet_notSubjet1Eta_%s"%(channel,suffix),
                                      op.c_float(-3.),
                                      has_notSubjet1,
                                      EquidistantBinning(20,-3.,3.),
                                      title='Pseudorapidity of the first subjet',
                                      xTitle="#eta (first subjet) [GeV]")
    plots.append(SummedPlot("%s_fatjet_subjet1Eta_%s"%(channel,suffix),
                                      [plot_hasSubjet1_Eta,plot_notSubjet1_Eta],
                                      xTitle="#eta (first subjet) [GeV]"))

    plot_hasSubjet2_Eta = Plot.make1D("%s_fatjet_hasSubjet2Eta_%s"%(channel,suffix),
                                      fatjet.subJet2.p4.Eta(),
                                      has_subjet2,
                                      EquidistantBinning(20,-3.,3.),
                                      title='Pseudorapidity of the second subjet',
                                      xTitle="#eta (second subjet) [GeV]")
    plot_notSubjet2_Eta = Plot.make1D("%s_fatjet_notSubjet2Eta_%s"%(channel,suffix),
                                      op.c_float(-3.),
                                      has_notSubjet2,
                                      EquidistantBinning(20,-3.,3.),
                                      title='Pseudorapidity of the second subjet',
                                      xTitle="#eta (second subjet) [GeV]")
    plots.append(SummedPlot("%s_fatjet_subjet2Eta_%s"%(channel,suffix),
                            [plot_hasSubjet2_Eta,plot_notSubjet2_Eta],
                            xTitle="#eta (second subjet) [GeV]"))

    # Phi Plot #
    plot_hasSubjet1_Phi = Plot.make1D("%s_fatjet_hasSubjet1Phi_%s"%(channel,suffix),
                                      fatjet.subJet1.p4.Phi(),
                                      has_subjet1,
                                      EquidistantBinning(20,-3.2,3.2),
                                      title='Azimutal angle of the first subjet',
                                      xTitle="#phi (first subjet) [GeV]")
    plot_notSubjet1_Phi = Plot.make1D("%s_fatjet_notSubjet1Phi_%s"%(channel,suffix),
                                      op.c_float(-3.2),
                                      has_notSubjet1,
                                      EquidistantBinning(20,-3.2,3.2),
                                      title='Azimutal angle of the first subjet',
                                      xTitle="#phi (first subjet) [GeV]")
    plots.append(SummedPlot("%s_fatjet_subjet1Phi_%s"%(channel,suffix),
                            [plot_hasSubjet1_Phi,plot_notSubjet1_Phi],
                            xTitle="#phi (first subjet) [GeV]"))

    plot_hasSubjet2_Phi = Plot.make1D("%s_fatjet_hasSubjet2Phi_%s"%(channel,suffix),
                                      fatjet.subJet2.p4.Phi(),
                                      has_subjet2,
                                      EquidistantBinning(20,-3.2,3.2),
                                      title='Azimutal angle of the second subjet',
                                      xTitle="#phi (second subjet) [GeV]")
    plot_notSubjet2_Phi = Plot.make1D("%s_fatjet_notSubjet2Phi_%s"%(channel,suffix),
                                      op.c_float(-3.2),
                                      has_notSubjet2,
                                      EquidistantBinning(20,-3.2,3.2),
                                      title='Azimutal angle of the second subjet',
                                      xTitle="#phi (second subjet) [GeV]")
    plots.append(SummedPlot("%s_fatjet_subjet2Phi_%s"%(channel,suffix),
                            [plot_hasSubjet2_Phi,plot_notSubjet2_Phi],
                            xTitle="#phi (second subjet) [GeV]"))

    # Invariant mass #
    plot_hasBothSubjets_invMass = Plot.make1D("%s_fatjet_hasBothSubjets_invMass_%s"%(channel,suffix),
                                              op.invariant_mass(fatjet.subJet1.p4,fatjet.subJet2.p4), 
                                              has_bothSubjets, EquidistantBinning(100, -10., 300.), 
                                              title="Fatjet invariant mass (channel %s)"%channel, 
                                              xTitle= "Invariant mass [GeV]")
    plot_hasNotBothSubjets_invMass = Plot.make1D("%s_fatjet_notBothSubjets_invMass_%s"%(channel,suffix),
                                                 op.c_float(-10.), 
                                                 has_notBothSubjets, 
                                                 EquidistantBinning(100, -10., 300.), 
                                                 title="Fatjet invariant mass (channel %s)"%channel, 
                                                 xTitle= "Invariant mass [GeV]")

    plots.append(SummedPlot("%s_fatjet_invMass_%s"%(channel,suffix),
                            [plot_hasBothSubjets_invMass,plot_hasNotBothSubjets_invMass],
                            xTitle="Invariant mass [GeV]"))

    return plots

