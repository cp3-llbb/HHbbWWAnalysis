import math
from bamboo.plots import Plot, EquidistantBinning, SummedPlot
from bamboo import treefunctions as op

# TODO : remove the self

######################## Channel title #################################
def channelTitleLabel(channel):
    if (channel == "ElEl"):
        channel = "e^{+}e^{-}"
    elif (channel == "MuMu"):
        channel ="#mu^{+}#mu^{-}"
    elif (channel == "ElMu"):
        channel = "e^{#pm} #mu^{#mp}"
    channelLabel = {'labels': [{'text': channel, 'position': [0.23, 0.87], 'size': 36}]}
    return channelLabel 

##########################  YIELD PLOT #################################
def makeYieldPlot(self, sel, name, title, order):
    """
    Make Yield plot and use it alos in the latex yield table
    sel     = refine selection
    name    = name of the PDF to be produced
    title   = title that will be used in the LateX yield table
    order   = int that gives the entry order in yield table (user must ensure the order)
    """
    plot = Plot.make1D("Yield_"+name,   
                       op.c_int(1),
                       sel,
                       EquidistantBinning(1, 0., 1.),
                       title = title + " Yield",
                       xTitle = title + " Yield",
                       plotopts = {"for-yields":True, "yields-title":title, 'yields-table-order':order})
    return plot

##########################  MET PLOT #################################
def makeMETPlots(self, sel, met, suffix, channel):
    """
    Make MET basic plots
    sel         = refine selection 
    met         = MET object
    suffix      = string identifying the selecton 
    channel     = string identifying the channel of the dilepton (can be "NoChannel")
    """

    plots = []

    channelLabel = channelTitleLabel(channel)

    # PT plot #
    plots.append(Plot.make1D("%s_%s_met_pt"%(channel,suffix), 
                             met.pt, 
                             sel, 
                             EquidistantBinning(50, 0., 300.), 
                             title="Transverse momentum of the MET (channel %s)"%channel, 
                             xTitle= "P_{T}(MET) [GeV]",
                             plotopts = channelLabel))
    # Phi plot #
    plots.append(Plot.make1D("%s_%s_met_phi"%(channel,suffix), 
                             met.phi, 
                             sel, 
                             EquidistantBinning(20, -3.2, 3.2), 
                             title="Azimutal angle of the MET (channel %s)"%channel, 
                             xTitle= "#phi (MET)",
                             plotopts = channelLabel))

    return plots

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

    channelLabel = channelTitleLabel(channel)

    # PT plot #
    plots.append(Plot.make1D("%s_%s_firstlepton_pt"%(channel,suffix), 
                             dilepton[0].p4.Pt(), 
                             sel, 
                             EquidistantBinning(50,0.,300.),
                             title="Transverse momentum of the first lepton (channel %s)"%channel, 
                             xTitle= "P_{T} (first lepton) [GeV]",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_secondlepton_pt"%(channel,suffix), 
                             dilepton[1].p4.Pt(), 
                             sel, 
                             EquidistantBinning(50,0.,300.),
                             title="Transverse momentum of the second lepton (channel %s)"%channel, 
                             xTitle= "P_{T} (second lepton) [GeV]",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_dilepton_pt"%(channel,suffix), 
                             (dilepton[0].p4+dilepton[1].p4).Pt(), 
                             sel, 
                             EquidistantBinning(50,0.,300.),
                             title="Transverse momentum of the dilepton (channel %s)"%channel, 
                             xTitle= "P_{T} (dilepton) [GeV]",
                             plotopts = channelLabel))

    # Eta plot #
    plots.append(Plot.make1D("%s_%s_firstlepton_eta"%(channel,suffix), 
                             dilepton[0].p4.Eta(), 
                             sel, 
                             EquidistantBinning(20, -3., 3.), 
                             title="Pseudorapidity of the first lepton (channel %s)"%channel, 
                             xTitle= "#eta (first lepton)",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_secondlepton_eta"%(channel,suffix), 
                             dilepton[1].p4.Eta(), 
                             sel, 
                             EquidistantBinning(20, -3., 3.), 
                             title="Pseudorapidity of the second lepton (channel %s)"%channel, 
                             xTitle= "#eta (second lepton)",
                             plotopts = channelLabel))

    # Phi plot #
    plots.append(Plot.make1D("%s_%s_firstlepton_phi"%(channel,suffix), 
                             dilepton[0].p4.Phi(), 
                             sel, 
                             EquidistantBinning(20, -3.2, 3.2), 
                             title="Azimutal angle of the first lepton (channel %s)"%channel, 
                             xTitle= "#phi (first lepton)",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_secondlepton_phi"%(channel,suffix), 
                             dilepton[1].p4.Phi(), 
                             sel, 
                             EquidistantBinning(20, -3.2, 3.2), 
                             title="Azimutal angle of the second lepton (channel %s)"%channel, 
                             xTitle= "#phi (second lepton)",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_dilepton_deltaPhi"%(channel,suffix), 
                             op.abs(op.deltaPhi(dilepton[0].p4,dilepton[1].p4)),
                             sel, 
                             EquidistantBinning(20,0, 3.2), 
                             title="Azimutal angle difference of the dilepton (channel %s)"%channel, 
                             xTitle= "|#Delta \phi (dilepton)|",
                             plotopts = channelLabel))

    # DeltaR plot #
    plots.append(Plot.make1D("%s_%s_dilepton_deltaR"%(channel,suffix), 
                             op.deltaR(dilepton[0].p4,dilepton[1].p4), 
                             sel, 
                             EquidistantBinning(50, 0., 5.), 
                             title="Dilepton Delta R (channel %s)"%channel, 
                             xTitle= "#Delta R (dilepton)",
                             plotopts = channelLabel))
    # InvMass plot #
    plots.append(Plot.make1D("%s_%s_dilepton_invariantMass"%(channel,suffix), 
                             op.invariant_mass(dilepton[0].p4,dilepton[1].p4), 
                             sel, 
                             EquidistantBinning(50, 0., 500.), 
                             title="Dilepton invariant mass (channel %s)"%channel, 
                             xTitle= "Invariant mass (dilepton) [GeV]",
                             plotopts = channelLabel))

    return plots

########################## DELTA PLOTS #################################
def makeDeltaRPlots(self,sel,cont1,cont2,suffix,channel,isMC):
    plots = [] 

    channelLabel = channelTitleLabel(channel)

    mixedCont = op.combine((cont1,cont2))
    plots.append(Plot.make1D("%s_%s_DeltaR"%(channel,suffix), 
                             op.map(mixedCont,lambda c : op.deltaR(c[0].p4,c[1].p4)),
                             sel, 
                             EquidistantBinning(100,0.,5.), 
                             title="#Delta R %s"%channel, 
                             xTitle= "#Delta R",
                             plotopts = channelLabel))

#    mixedCont_inDeltaR_0p4 = op.combine((cont1,cont2),pred=lambda c1,c2 : op.deltaR(c1.p4,c2.p4)<0.4)
#    mixedCont_outDeltaR_0p4 = op.combine((cont1,cont2),pred=lambda c1,c2 : op.deltaR(c1.p4,c2.p4)>=0.4)
#    if isMC: # We don't have gen level info on data
#        plots.append(Plot.make1D("%s_%s_ID_DeltaRCut_Below0p4"%(channel,suffix), 
#                                 op.map(mixedCont_inDeltaR_0p4,lambda c : c[0].genPartFlav), # Print gen flavour of first cont (aka lepton)
#                                 sel, 
#                                 EquidistantBinning(22,0.,22.), 
#                                 title="%s genFlavour ID : #Delta R < 0.4"%channel, 
#                                 xTitle= "ID"))
#        plots.append(Plot.make1D("%s_%s_ID_DeltaRCut_Above0p4"%(channel,suffix), 
#                                 op.map(mixedCont_outDeltaR_0p4,lambda c : c[0].genPartFlav), # Print gen flavour of first cont (aka lepton)
#                                 sel, 
#                                 EquidistantBinning(22,0.,22.), 
#                                 title="%s genFlavour ID : #Delta R > 0.4"%channel, 
#                                 xTitle= "#ID"))
 
        # Electron/Muon  genPartFlav :    1       -> prompt electron
        #                                 3,4,5   -> light quark, c quark, b quark
        #                                 15      -> tau 
        #                                 22      -> photon conversion (only electron)
        # script : https://github.com/cms-nanoAOD/cmssw/blob/6f76b77b43c852475cff452a97203441a9c962d1/PhysicsTools/NanoAOD/plugins/CandMCMatchTableProducer.cc#L99-L138



    return plots

########################## JETS PLOTS #################################
def makeJetsPlots(self,sel,bjets,lightjets,alljets,suffix,channel):
    plots = []

    channelLabel = channelTitleLabel(channel)

    # Selection for bjets and mixed categories #
    bjetsCategory = sel.refine(suffix+"bjetsCategory",cut=[op.rng_len(bjets) >= 2]) # TODO : might want to veto more than 2 bjets 
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

    # Number of Resolved jets #
    plots.append(Plot.make1D("%s_%s_resolvedjets_N"%(channel,suffix),                                   
                             op.rng_len(bjets),
                             sel,
                             EquidistantBinning(5,0.,5.),                                              
                             title='Number of Resolved jets (%s channel)'%channel,                     
                             xTitle='N Resolved jets',
                             plotopts = channelLabel))
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

    channelLabel = channelTitleLabel(channel)

    if plot_type == "inclusive":
        lead_base_name      = "%s_%sinclusive_leadjet_{var}"%(channel,suffix)
        lead_base_title     = "leading jet"
        sublead_base_name   = "%s_%sinclusive_subleadjet_{var}"%(channel,suffix)
        sublead_base_title  = "subleading jet"
    elif plot_type == "bjets":
        lead_base_name      = "%s_%sbjets_leadjet_{var}"%(channel,suffix)
        lead_base_title     = "leading bjet"
        sublead_base_name   = "%s_%sbjets_subleadjet_{var}"%(channel,suffix)
        sublead_base_title  = "subleading bjet"
    elif plot_type == "mixed":
        lead_base_name      = "%s_%smixed_bjet_{var}"%(channel,suffix)
        lead_base_title     = "bjet"
        sublead_base_name   = "%s_%smixed_lightjet_{var}"%(channel,suffix)
        sublead_base_title  = "lightjet"
    else:
        print ("[ERROR] Jet plot type no understood")
        sys.exit(1)

    # leadjet plots #
    plots.append(Plot.make1D(lead_base_name.format(var="pt"),
                             leadjet.p4.pt(),
                             sel,
                             EquidistantBinning(50,0.,300.),
                             title='Transverse momentum of the %s'%lead_base_title,
                             xTitle="P_{T}(%s) [GeV]"%lead_base_title,
                             plotopts = channelLabel))
    plots.append(Plot.make1D(lead_base_name.format(var="eta"),
                             leadjet.p4.eta(),
                             sel,
                             EquidistantBinning(20,-3.,3.),
                             title='Pseudorapidity of the %s'%lead_base_title,
                             xTitle="#eta(%s)"%lead_base_title,
                             plotopts = channelLabel))
    plots.append(Plot.make1D(lead_base_name.format(var="phi"),
                             leadjet.p4.phi(),
                             sel,
                             EquidistantBinning(20,-3.2,3.2),
                             title='Azimutal angle of the %s'%lead_base_title,
                             xTitle="#phi(%s)"%lead_base_title,
                             plotopts = channelLabel))
    # subleadjet plots #
    plots.append(Plot.make1D(sublead_base_name.format(var="pt"),
                             subleadjet.p4.pt(),
                             sel,
                             EquidistantBinning(50,0.,300.),
                             title='Transverse momentum of the %s'%sublead_base_title,
                             xTitle="P_{T}(%s) [GeV]"%sublead_base_title,
                             plotopts = channelLabel))
    plots.append(Plot.make1D(sublead_base_name.format(var="eta"),
                             subleadjet.p4.eta(),
                             sel,
                             EquidistantBinning(20,-3.,3.),
                             title='Pseudorapidity of the %s'%sublead_base_title,
                             xTitle="#eta(%s)"%sublead_base_title,
                             plotopts = channelLabel))
    plots.append(Plot.make1D(sublead_base_name.format(var="phi"),
                             subleadjet.p4.phi(),
                             sel,
                             EquidistantBinning(20,-3.2,3.2),
                             title='Azimutal angle of the %s'%sublead_base_title,
                             xTitle="#phi(%s)"%sublead_base_title,
                             plotopts = channelLabel))
    # Dijet Pt plot #
    plots.append(Plot.make1D("%s_%s%s_dijet_pt"%(channel,suffix,plot_type), 
                             (leadjet.p4+subleadjet.p4).Pt(), 
                             sel, 
                             EquidistantBinning(50,0.,300.),
                             title="Transverse momentum of the dijet (channel %s)"%channel, 
                             xTitle= "P_{T}(dijet) [GeV]",
                             plotopts = channelLabel))
    # DeltaPhi plot #
    plots.append(Plot.make1D("%s_%s%s_dijet_deltaPhi"%(channel,suffix,plot_type), 
                             op.abs(op.deltaPhi(leadjet.p4,subleadjet.p4)), 
                             sel, 
                             EquidistantBinning(20,0, 3.2),
                             title="Azimutal angle difference of the dijet (channel %s)"%channel,
                             xTitle= "| #Delta \phi (dijet)|",
                             plotopts = channelLabel))

    # DeltaR plot #
    plots.append(Plot.make1D("%s_%s%s_dijet_deltaR"%(channel,suffix,plot_type), 
                             op.deltaR(leadjet.p4,subleadjet.p4), 
                             sel, 
                             EquidistantBinning(50, 0., 5.), 
                             title="Dilepton Delta R (channel %s)"%channel, 
                             xTitle= "#Delta R (dijet)",
                             plotopts = channelLabel))
    # invariant mass plot #
    plots.append(Plot.make1D("%s_%s%s_dijet_invariantMass"%(channel,suffix,plot_type),
                             op.invariant_mass(leadjet.p4,subleadjet.p4),
                             sel,
                             EquidistantBinning(50, 0., 500.), 
                             title="Dijet invariant mass (channel %s)"%channel, 
                             xTitle= "Invariant mass (dijet) [GeV]",
                             plotopts = channelLabel))

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

    channelLabel = channelTitleLabel(channel)

    # fatjet plots (always present by selection) #
    plots.append(Plot.make1D("%s_%s_fatjet_pt"%(channel,suffix),
                             fatjets[0].p4.pt(),
                             sel,
                             EquidistantBinning(50,200,600.),
                             title='Transverse momentum of the fatjet',
                             xTitle="P_{T}(fatjet) [GeV]",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_fatjet_eta"%(channel,suffix),
                             fatjets[0].p4.eta(),
                             sel,
                             EquidistantBinning(10,-3.,3.),
                             title='Pseudorapidity of the fatjet',
                             xTitle="#eta(fatjet)",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_fatjet_phi"%(channel,suffix),
                             fatjets[0].p4.phi(),
                             sel,
                             EquidistantBinning(10,-3.2,3.2),
                             title='Azimutal angle of the fatjet',
                             xTitle="#phi(fatjet)",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_fatjet_mass"%(channel,suffix),
                             fatjets[0].mass,
                             sel,
                             EquidistantBinning(50,0.,200.),
                             title='Invariant mass of the fatjet',
                             xTitle="M(fatjet)",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_fatjet_softdropmass"%(channel,suffix),
                             fatjets[0].msoftdrop,
                             sel,
                             EquidistantBinning(50,0.,200.),
                             title='Soft Drop mass of the fatjet',
                             xTitle="M_{Soft Drop}(fatjet) [GeV]",
                             plotopts = channelLabel))

    # Number of Boosted jets #
    plots.append(Plot.make1D("%s_%s_boostedjets_N"%(channel,suffix),                                   
                             op.rng_len(fatjets),
                             sel,
                             EquidistantBinning(5,0.,5.),                                              
                             title='Number of Boosted jets (%s channel)'%channel,                     
                             xTitle='N Boosted jets',   
                             plotopts = channelLabel))

    return plots

#########################  High-level quantities ################################
def makeHighLevelQuantities(self,sel,dilepton,met,jets,resolvedjets,lightjets,boostedjets,suffix,channel):
    plots = []

    # Categories spltting #
    # Check if boosted category : make selection on number of boosted jets #
    isBoosted = sel.refine(suffix+"highlevelBoosted",cut=[op.rng_len(boostedjets)>=1]) 
        # If has boosted jet : goes into boosted (even if resolved jet present) -> inclusive
    isResolved = sel.refine(suffix+"highlevelResolved",cut=[op.rng_len(boostedjets)==0]) 
        # If not boosted : must be resolved -> exclusive
    # Resolved subcategories 
    # isResolved -> inclusive : pass the jet leading and sybleading
        # bjets : pass the resolved jets 
    isResolved_bjets = isResolved.refine(suffix+"highlevelResolvedBjets",
                                         cut=[op.rng_len(resolvedjets)>=2])
        # Mixed : pass leading bjet and leading lightjet
    isResolved_mixed = isResolved.refine(suffix+"highlevelResolvedMixed",
                                         cut=[op.rng_len(resolvedjets)==1])

    # Pass to function #
    # Boosted #
    plots.extend(makeSeparateHighLevelQuantities(self       = self,
                                                 sel        = isBoosted,
                                                 met        = met,
                                                 l1         = dilepton[0],
                                                 l2         = dilepton[1],
                                                 j1         = boostedjets[0].subJet1,
                                                 j2         = boostedjets[0].subJet2,
                                                 suffix     = suffix,
                                                 channel    = channel))
    # Resolved inclusive #
    plots.extend(makeSeparateHighLevelQuantities(self       = self,
                                                 sel        = isResolved,
                                                 met        = met,
                                                 l1         = dilepton[0],
                                                 l2         = dilepton[1],
                                                 j1         = jets[0],
                                                 j2         = jets[1],
                                                 suffix     = suffix+"inclusive",
                                                 channel    = channel))
    # Resolved bjets #
    plots.extend(makeSeparateHighLevelQuantities(self       = self,
                                                 sel        = isResolved_bjets,
                                                 met        = met,
                                                 l1         = dilepton[0],
                                                 l2         = dilepton[1],
                                                 j1         = resolvedjets[0],
                                                 j2         = resolvedjets[1],
                                                 suffix     = suffix+"bjets",
                                                 channel    = channel))
    # Resolved bjets #
    plots.extend(makeSeparateHighLevelQuantities(self       = self,
                                                 sel        = isResolved_mixed,
                                                 met        = met,
                                                 l1         = dilepton[0],
                                                 l2         = dilepton[1],
                                                 j1         = resolvedjets[0],
                                                 j2         = lightjets[0],
                                                 suffix     = suffix+"mixed",
                                                 channel    = channel))

    return plots

def makeSeparateHighLevelQuantities(self,sel,met,l1,l2,j1,j2,suffix,channel):
    plots = []

    channelLabel = channelTitleLabel(channel)

    # Useful lambdas #
    ll_p4 = lambda l1,l2 : l1.p4+l2.p4
    lljj_p4 = lambda l1,l2,j1,j2 : l1.p4+l2.p4+j1.p4+j2.p4
 
    # dilepton-MET plots #
    plots.append(Plot.make1D("%s_%s_highlevelvariable_DilepMETdeltaPhi"%(channel,suffix),
                             op.abs(ll_p4(l1,l2).Phi()-met.phi),
                             sel,
                             EquidistantBinning(20,0.,3.2),
                             title='Azimutal angle between dilepton and MET (%s channel)'%channel,
                             xTitle="|#Delta \phi (ll,MET)|",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_highlevelvariable_DilepMETpt"%(channel,suffix),
                             op.sqrt(op.pow(met.pt*op.cos(met.phi)+ll_p4(l1,l2).Px(),2)+op.pow(met.pt*op.sin(met.phi)+ll_p4(l1,l2).Py(),2)),
                             sel,
                             EquidistantBinning(50,0.,500.),
                             title='Transverse momentum of dilepton and MET (%s channel)'%channel,
                             xTitle="P_{T}(ll,MET) [GeV]",
                             plotopts = channelLabel))

    # Transverse mass plots #
    # m_T = sqrt(2*P_T^{l1l2}*P_T^{miss} * (1-cos(\Delta\Phi(ll,P_T^{miss}))))
    mTll = lambda l1,l2,met : op.sqrt(2*ll_p4(l1,l2).Pt()*met.pt*(1-op.cos(ll_p4(l1,l2).Phi()-met.phi)))
    mTlljj = lambda l1,l2,j1,j2,met : op.sqrt(2*lljj_p4(l1,l2,j1,j2).Pt()*met.pt*(1-op.cos(lljj_p4(l1,l2,j1,j2).Phi()-met.phi)))
        # mT for dilepton #
    plots.append(Plot.make1D("%s_%s_highlevelvariable_mTll"%(channel,suffix),
                             mTll(l1,l2,met),
                             sel,
                             EquidistantBinning(50,0.,1000.),
                             title='Transverse mass of dilepton and MET (%s channel)'%channel,
                             xTitle="M_{T}(ll,MET) [GeV]",
                             plotopts = channelLabel))
        # mT for lljj #
    plots.append(Plot.make1D("%s_%s_highlevelvariable_mTlljj"%(channel,suffix),
                             mTlljj(l1,l2,j1,j2,met),
                             sel,
                             EquidistantBinning(50,0.,1000.),
                             title='Transverse mass of dilepton+dijet and MET (%s channel)'%channel,
                             xTitle="M_{T}(lljj,MET) [GeV]",
                             plotopts = channelLabel))


    # Scalar magnitude sum #
    HT2 = lambda met,l1,l2,j1,j2 : op.sqrt(op.pow(met.pt*op.cos(met.phi)+l1.p4.Px()+l2.p4.Px(),2)+op.pow(met.pt*op.sin(met.phi)+l1.p4.Py()+l2.p4.Py(),2)) + op.abs((j1.p4+j2.p4).Pt())
    HT2R = lambda met,l1,l2,j1,j2 : HT2(met,l1,l2,j1,j2)/(met.pt+l1.p4.Pt()+l2.p4.Pt()+j1.p4.Pt()+j2.p4.Pt())

    plots.append(Plot.make1D("%s_%s_highlevelvariable_HT2"%(channel,suffix),
                             HT2(met,l1,l2,j1,j2),
                             sel,
                             EquidistantBinning(50,0.,1000.),
                             title='Di-Higgs magnitude (%s channel)'%channel,
                             xTitle="H_{T2} [GeV]",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_highlevelvariable_HT2R"%(channel,suffix),
                             HT2R(met,l1,l2,j1,j2),
                             sel,
                             EquidistantBinning(50,0.,1.),
                             title='Di-Higgs magnitude ratio (%s channel)'%channel,
                             xTitle="H_{T2}^{R} [GeV]",
                             plotopts = channelLabel))

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

