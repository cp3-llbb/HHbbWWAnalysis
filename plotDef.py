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
    plots.append(Plot.make1D("%s_%s_leadlepton_pt"%(channel,suffix), 
                             dilepton[0].p4.Pt(), 
                             sel, 
                             EquidistantBinning(100, 0., 300.), 
                             title="Transverse momentum of the leading lepton (channel %s)"%channel, 
                             xTitle= "P_{T}(leading lepton) [GeV]"))
    plots.append(Plot.make1D("%s_%s_subleadlepton_pt"%(channel,suffix), 
                             dilepton[1].p4.Pt(), 
                             sel, 
                             EquidistantBinning(100, 0., 300.), 
                             title="Transverse momentum of the subleading lepton (channel %s)"%channel, 
                             xTitle= "P_{T}(leading lepton) [GeV]"))

    # Eta plot #
    plots.append(Plot.make1D("%s_%s_leadlepton_eta"%(channel,suffix), 
                             dilepton[0].p4.Eta(), 
                             sel, 
                             EquidistantBinning(20, -3., 3.), 
                             title="Pseudorapidity of the leading lepton (channel %s)"%channel, 
                             xTitle= "#eta (leading lepton) [GeV]"))
    plots.append(Plot.make1D("%s_%s_subleadlepton_eta"%(channel,suffix), 
                             dilepton[1].p4.Eta(), 
                             sel, 
                             EquidistantBinning(20, -3., 3.), 
                             title="Pseudorapidity of the subleading lepton (channel %s)"%channel, 
                             xTitle= "#eta (subleading sublepton) [GeV]"))

    # Phi plot #
    plots.append(Plot.make1D("%s_%s_leadlepton_phi"%(channel,suffix), 
                             dilepton[0].p4.Phi(), 
                             sel, 
                             EquidistantBinning(20, -3.2, 3.2), 
                             title="Azimutal angle of the leading lepton (channel %s)"%channel, 
                             xTitle= "#phi (leading lepton) [GeV]"))
    plots.append(Plot.make1D("%s_%s_subleadlepton_phi"%(channel,suffix), 
                             dilepton[1].p4.Phi(), 
                             sel, 
                             EquidistantBinning(20, -3.2, 3.2), 
                             title="Azimutal angle of the subleading lepton (channel %s)"%channel, 
                             xTitle= "#phi (subleading lepton) [GeV]"))

    # InvMass plot #
    plots.append(Plot.make1D("%s_%s_dilepton_invariantMass"%(channel,suffix), 
                             op.invariant_mass(dilepton[0].p4,dilepton[1].p4), 
                             sel, 
                             EquidistantBinning(100, 0., 400.), 
                             title="Dilepton invariant mass (channel %s)"%channel, 
                             xTitle= "Invariant mass [GeV]"))

    return plots

########################## DELTA PLOTS #################################
def makeDeltaRPlots(self,sel,cont1,cont2,suffix,channel,isMC):
    plots = [] 
    mixedCont = op.combine((cont1,cont2))
    plots.append(Plot.make1D("%s_%s_DeltaR"%(channel,suffix), 
                             op.map(mixedCont,lambda c : op.deltaR(c[0].p4,c[1].p4)),
                             sel, 
                             EquidistantBinning(100,0.,5.), 
                             title="#Delta R %s"%channel, 
                             xTitle= "#Delta R"))

    mixedCont_inDeltaR_0p4 = op.combine((cont1,cont2),pred=lambda c1,c2 : op.deltaR(c1.p4,c2.p4)<0.4)
    mixedCont_outDeltaR_0p4 = op.combine((cont1,cont2),pred=lambda c1,c2 : op.deltaR(c1.p4,c2.p4)>=0.4)
    if isMC: # We don't have gen level info on data
        plots.append(Plot.make1D("%s_%s_ID_DeltaRCut_Below0p4"%(channel,suffix), 
                                 op.map(mixedCont_inDeltaR_0p4,lambda c : c[0].genPartFlav), # Print gen flavour of first cont (aka lepton)
                                 sel, 
                                 EquidistantBinning(22,0.,22.), 
                                 title="%s genFlavour ID : #Delta R < 0.4"%channel, 
                                 xTitle= "ID"))
        plots.append(Plot.make1D("%s_%s_ID_DeltaRCut_Above0p4"%(channel,suffix), 
                                 op.map(mixedCont_outDeltaR_0p4,lambda c : c[0].genPartFlav), # Print gen flavour of first cont (aka lepton)
                                 sel, 
                                 EquidistantBinning(22,0.,22.), 
                                 title="%s genFlavour ID : #Delta R > 0.4"%channel, 
                                 xTitle= "#ID"))
 
        # Electron/Muon  genPartFlav :    1       -> prompt electron
        #                                 3,4,5   -> light quark, c quark, b quark
        #                                 15      -> tau 
        #                                 22      -> photon conversion (only electron)



    return plots

########################## JETS PLOTS #################################
def makeJetsPlots(self,sel,bjets,lightjets,alljets,suffix,channel):
    plots = []

    # Selection for bjets and mixed categories #
    bjetsCategory = sel.refine(suffix+"bjetsCategory",cut=[op.rng_len(bjets) >= 2])
    mixedCategory = sel.refine(suffix+"mixedCategory",cut=[op.rng_len(bjets) == 1])

    ## Inclusive plots (len(alljets)>=2 by definition ): leading = highest Pt jet, subleading : second highest Pt jet #
    plots.extend(makeSeparateJetsPlots(self          = self,
                                       sel           = sel,
                                       leadjet       = alljets[0],
                                       subleadjet    = alljets[1],
                                       suffix        = suffix,
                                       channel       = channel,
                                       plot_type     = "inclusive"))

    ## Bjets plots (if len(bjets)>=2) : leading = highest Pt bjet, subleading : second highest Pt bjet #
    # Must check len(bjets)>=2
    plots.extend(makeSeparateJetsPlots(self          = self,
                                       sel           = bjetsCategory,
                                       leadjet       = bjets[0],
                                       subleadjet    = bjets[1],
                                       suffix        = suffix,
                                       channel       = channel,
                                       plot_type     = "bjets"))

    ## Mixed plots (if len(bjets)==1, >0 by definition) : leading = bjet, subleading = highest Pt non-bjet
    # Must check len(bjets)==1
    plots.extend(makeSeparateJetsPlots(self          = self,
                                       sel           = mixedCategory,
                                       leadjet       = bjets[0],
                                       subleadjet    = lightjets[0],
                                       suffix        = suffix,
                                       channel       = channel,
                                       plot_type     = "mixed"))
    return plots

##########################  JETS SEPARATE PLOTS #################################
def makeSeparateJetsPlots(self, sel, leadjet, subleadjet, suffix, channel, plot_type):
    """
    Make basic plots
    sel         = refine selection 
    leadjet     = Depends on scenario (plot_type) : highest-Pt jet ("inclusive"), highest-Pt bjet ("bjets"), only bjet ("mixed")
    subleadjet  = Depends on scenario (plot_type) : second-highest-Pt jet ("inclusive"), second-highest-Pt bjet ("bjets"), highest-Pt lightjet ("mixed")
    jets        = jet container (content depends on plot_type) 
    suffix      = string identifying the selection 
    channel     = string identifying the channel of the dilepton (can be "NoChannel")
    plot_type   = defines the scenario : "inclusive", "bjets", "mixed"
    """
 
    plots = []
    if plot_type == "inclusive":
        lead_base_name      = "%s_%s_inclusive_leadjet_{var}"%(channel,suffix)
        lead_base_title     = "leading jet"
        sublead_base_name   = "%s_%s_inclusive_subleadjet_{var}"%(channel,suffix)
        sublead_base_title  = "subleading jet"
    elif plot_type == "bjets":
        lead_base_name      = "%s_%s_bjets_leadjet_{var}"%(channel,suffix)
        lead_base_title     = "leading bjet"
        sublead_base_name   = "%s_%s_bjets_subleadjet_{var}"%(channel,suffix)
        sublead_base_title  = "subleading bjet"
    elif plot_type == "mixed":
        lead_base_name      = "%s_%s_mixed_bjet_{var}"%(channel,suffix)
        lead_base_title     = "bjet"
        sublead_base_name   = "%s_%s_mixed_lightjet_{var}"%(channel,suffix)
        sublead_base_title  = "lightjet"
    else:
        print ("[ERROR] Jet plot type no understood")
        sys.exit(1)

    # leadjet plots #
    plots.append(Plot.make1D(lead_base_name.format(var="pt"),
                             leadjet.p4.pt(),
                             sel,
                             EquidistantBinning(100,0,300.),
                             title='Transverse momentum of the %s'%lead_base_title,
                             xTitle="P_{T}(%s) [GeV]"%lead_base_title,))
    plots.append(Plot.make1D(lead_base_name.format(var="eta"),
                             leadjet.p4.eta(),
                             sel,
                             EquidistantBinning(20,-3.,3.),
                             title='Pseudorapidity of the %s'%lead_base_title,
                             xTitle="#eta(%s)"%lead_base_title))
    plots.append(Plot.make1D(lead_base_name.format(var="phi"),
                             leadjet.p4.phi(),
                             sel,
                             EquidistantBinning(20,-3.2,3.2),
                             title='Azimutal angle of the %s'%lead_base_title,
                             xTitle="#phi(%s)"%lead_base_title))
    # subleadjet plots #
    plots.append(Plot.make1D(sublead_base_name.format(var="pt"),
                             subleadjet.p4.pt(),
                             sel,
                             EquidistantBinning(100,0,300.),
                             title='Transverse momentum of the %s'%sublead_base_title,
                             xTitle="P_{T}(%s) [GeV]"%sublead_base_title,))
    plots.append(Plot.make1D(sublead_base_name.format(var="eta"),
                             subleadjet.p4.eta(),
                             sel,
                             EquidistantBinning(20,-3.,3.),
                             title='Pseudorapidity of the %s'%sublead_base_title,
                             xTitle="#eta(%s)"%sublead_base_title))
    plots.append(Plot.make1D(sublead_base_name.format(var="phi"),
                             subleadjet.p4.phi(),
                             sel,
                             EquidistantBinning(20,-3.2,3.2),
                             title='Azimutal angle of the %s'%sublead_base_title,
                             xTitle="#phi(%s)"%sublead_base_title))
    # invariant mass plot #
    plots.append(Plot.make1D("%s_%s_%s_dijet_invariantMass"%(channel,suffix,plot_type),
                             op.invariant_mass(leadjet.p4,subleadjet.p4),
                             sel,
                             EquidistantBinning(100, 0., 800.), 
                             title="Dijet invariant mass (channel %s)"%channel, 
                             xTitle= "Invariant mass [GeV]"))

    return plots

##########################  JETS (AK4) PLOT #################################
def makeFatJetPlots(self, sel, fatjets, suffix, channel):
    """
    Make fatjet subjet basic plots
    sel         = refine selection 
    fatjets     = fatjet container (AK8 jet) 
    suffix      = string identifying the selecton 
    channel     = string identifying the channel of the dilepton (can be "NoChannel")
    """
 
    plots = []

    # fatjet plots (always present by selection) #
    plots.append(Plot.make1D("%s_%s_fatjet_pt"%(channel,suffix),
                             fatjets[0].p4.pt(),
                             sel,
                             EquidistantBinning(100,150,600.),
                             title='Transverse momentum of the fatjet',
                             xTitle="P_{T}(fatjet) [GeV]"))
    plots.append(Plot.make1D("%s_%s_fatjet_eta"%(channel,suffix),
                             fatjets[0].p4.eta(),
                             sel,
                             EquidistantBinning(20,-3.,3.),
                             title='Pseudorapidity of the fatjet',
                             xTitle="#eta(fatjet)"))
    plots.append(Plot.make1D("%s_%s_fatjet_phi"%(channel,suffix),
                             fatjets[0].p4.phi(),
                             sel,
                             EquidistantBinning(20,-3.2,3.2),
                             title='Azimutal angle of the fatjet',
                             xTitle="#phi(fatjet)"))
    plots.append(Plot.make1D("%s_%s_fatjet_mass"%(channel,suffix),
                             fatjets[0].mass,
                             sel,
                             EquidistantBinning(100,0.,200.),
                             title='Invariant mass of the fatjet',
                             xTitle="M(fatjet)"))
    plots.append(Plot.make1D("%s_%s_fatjet_softdropmass"%(channel,suffix),
                             fatjets[0].msoftdrop,
                             sel,
                             EquidistantBinning(100,0.,200.),
                             title='Soft Drop mass of the fatjet',
                             xTitle="M_{Soft Drop}(fatjet) [GeV]"))
    return plots

##########################  FATJETS SUBJETS (AK8) PLOT #################################
def makeFatJetSubJetPlots(self, sel, fatjet, suffix, channel): # Mostly debugging 
    """
    Make fatjet subjet basic plots
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

