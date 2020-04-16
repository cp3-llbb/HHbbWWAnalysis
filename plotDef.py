import math
from bamboo.plots import Plot, EquidistantBinning, SummedPlot
from bamboo import treefunctions as op

# TODO : remove the self

########################   Channel title   #############################
def channelTitleLabel(channel):
    if (channel == "ElEl"):
        channel = "e^{+}e^{-}"
    elif (channel == "MuMu"):
        channel ="#mu^{+}#mu^{-}"
    elif (channel == "ElMu"):
        channel = "e^{#pm} #mu^{#mp}"
    channelLabel = {'labels': [{'text': channel, 'position': [0.23, 0.87], 'size': 36}]}
    return channelLabel 

######################  Plot Objects Number  ###########################
def objectsNumberPlot(channel,objName,suffix,objCont,sel,Nmax,xTitle):
    channelLabel = channelTitleLabel(channel)
    return Plot.make1D(channel+"_"+suffix+"_"+objName+"_N",
                       op.rng_len(objCont),
                       sel,
                       EquidistantBinning(Nmax,0.,Nmax),
                       xTitle = xTitle,
                       plotopts = channelLabel)

#######################   Channel Plot   ###############################
def channelPlot(sel,DilepElEl,DilepMuMu,DilepElMu,suffix,channel):
    """
    Plots the number of events in each channel
    """
    plots = []
    plots.append(objectsNumberPlot(channel = channel,
                                   suffix  = suffix,
                                   objName = "DilepElEl",
                                   objCont = DilepElEl,
                                   sel     = sel,
                                   Nmax    = 5,
                                   xTitle  = "N(e^{+}e^{-})"))
    plots.append(objectsNumberPlot(channel = channel,
                                   suffix  = suffix,
                                   objName = "DilepMuMu",
                                   objCont = DilepMuMu,
                                   sel     = sel,
                                   Nmax    = 5,
                                   xTitle  = "N(#mu^{+}#mu^{-})"))
    plots.append(objectsNumberPlot(channel = channel,
                                   suffix  = suffix,
                                   objName = "DilepElMu",
                                   objCont = DilepElMu,
                                   sel     = sel,
                                   Nmax    = 5,
                                   xTitle  = "N(e^{#pm} #mu^{#mp})"))
    return plots

##########################  YIELD PLOT #################################
class makeYieldPlots:
    def __init__(self):
        self.calls = 0
        self.plots = []
    def addYield(self, sel, name, title):
        """
        Make Yield plot and use it also in the latex yield table
        sel     = refine selection
        name    = name of the PDF to be produced
        title   = title that will be used in the LateX yield table
        """
        self.plots.append(Plot.make1D("Yield_"+name,   
                                      op.c_int(0),
                                      sel,
                                      EquidistantBinning(1, 0., 1.),
                                      title = title + " Yield",
                                      xTitle = title + " Yield",
                                      plotopts = {"for-yields":True, "yields-title":title, 'yields-table-order':self.calls}))
        self.calls += 1

    def returnPlots(self):
        return self.plots

######################   TRIGGER PLOTS  ##############################
def triggerPlots(sel,triggerDict,suffix,channel):
    plots = []
    channelLabel = channelTitleLabel(channel)
    
    # Single Electron #
    plots.append(Plot.make1D("%s_%s_trigger_singleElectron"%(channel,suffix), 
                             op.OR(*triggerDict["SingleElectron"]),
                             sel, 
                             EquidistantBinning(2, 0., 2.), 
                             xTitle= "SingleElectron trigger",
                             plotopts = channelLabel))
    # Single Muon #
    plots.append(Plot.make1D("%s_%s_trigger_singleMuon"%(channel,suffix), 
                             op.OR(*triggerDict['SingleMuon']), 
                             sel, 
                             EquidistantBinning(2, 0., 2.), 
                             xTitle= "SingleMuon trigger",
                             plotopts = channelLabel))
    # Double Electron #
    plots.append(Plot.make1D("%s_%s_trigger_doubleElectron"%(channel,suffix), 
                             op.OR(*triggerDict['DoubleEGamma']), 
                             sel, 
                             EquidistantBinning(2, 0., 2.), 
                             xTitle= "DoubleEGamma trigger",
                             plotopts = channelLabel))
    # Double Muon #
    plots.append(Plot.make1D("%s_%s_trigger_doubleMuon"%(channel,suffix), 
                             op.OR(*triggerDict['DoubleMuon']), 
                             sel, 
                             EquidistantBinning(2, 0., 2.), 
                             xTitle= "DoubleMuon trigger",
                             plotopts = channelLabel))
    # Muon-Electron #
    plots.append(Plot.make1D("%s_%s_trigger_MuonElectron"%(channel,suffix), 
                             op.OR(*triggerDict['MuonEG']), 
                             sel, 
                             EquidistantBinning(2, 0., 2.), 
                             xTitle= "MuonEG trigger",
                             plotopts = channelLabel))

    return plots
##########################  MET PLOT #################################
def makeMETPlots(sel, met, suffix, channel):
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
def makeDileptonPlots(sel, dilepton, suffix, channel, isMC=False):
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
                             dilepton[0].pt, 
                             sel, 
                             EquidistantBinning(50,0.,300.),
                             title="Transverse momentum of the first lepton (channel %s)"%channel, 
                             xTitle= "P_{T} (first lepton) [GeV]",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_secondlepton_pt"%(channel,suffix), 
                             dilepton[1].pt, 
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
                             dilepton[0].eta, 
                             sel, 
                             EquidistantBinning(20, -3., 3.), 
                             title="Pseudorapidity of the first lepton (channel %s)"%channel, 
                             xTitle= "#eta (first lepton)",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_secondlepton_eta"%(channel,suffix), 
                             dilepton[1].eta, 
                             sel, 
                             EquidistantBinning(20, -3., 3.), 
                             title="Pseudorapidity of the second lepton (channel %s)"%channel, 
                             xTitle= "#eta (second lepton)",
                             plotopts = channelLabel))

    # Phi plot #
    plots.append(Plot.make1D("%s_%s_firstlepton_phi"%(channel,suffix), 
                             dilepton[0].phi, 
                             sel, 
                             EquidistantBinning(20, -3.2, 3.2), 
                             title="Azimutal angle of the first lepton (channel %s)"%channel, 
                             xTitle= "#phi (first lepton)",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_secondlepton_phi"%(channel,suffix), 
                             dilepton[1].phi, 
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

    # GenPartFlav (if isMC) #
#    if isMC:
#        plots.append(Plot.make1D("%s_%s_firstlepton_genPartFlav"%(channel,suffix), 
#                                 dilepton[0].genPartFlav,
#                                 sel, 
#                                 EquidistantBinning(22, 0., 22.), 
#                                 title="Flavour of genParticle (channel %s)"%channel, 
#                                 xTitle= "GenParticle flavour (first lepton)",
#                                 plotopts = channelLabel))
#        plots.append(Plot.make1D("%s_%s_secondlepton_genPartFlav"%(channel,suffix), 
#                                 dilepton[1].genPartFlav,
#                                 sel, 
#                                 EquidistantBinning(22, 0., 22.), 
#                                 title="Flavour of genParticle (channel %s)"%channel, 
#                                 xTitle= "GenParticle flavour (second lepton)",
#                                 plotopts = channelLabel))
#
    return plots

######################  DILEPTON CHECK PLOTS #############################
def makeDileptonCheckPlots(self,sel, dilepton, suffix, channel):
    """
    Make dilepton basic plots
    sel         = refine selection 
    dilepton    = dilepton object (ElEl, ElMu, MuMu)
    suffix      = string identifying the selecton 
    channel     = string identifying the channel of the dilepton (can be "NoChannel")
    """
    plots = []

    channelLabel = channelTitleLabel(channel)

    # Lepton MVA #
    plots.append(Plot.make1D("%s_%s_firstlepton_checkMVAttH"%(channel,suffix), 
                             dilepton[0].mvaTTH, 
                             sel, 
                             EquidistantBinning(50,-1.,1.),
                             xTitle= "ttH Lepton MVA (first lepton)",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_secondlepton_checkMVAttH"%(channel,suffix), 
                             dilepton[1].mvaTTH, 
                             sel, 
                             EquidistantBinning(50,-1.,1.),
                             xTitle= "ttH Lepton MVA (second lepton)",
                             plotopts = channelLabel))

    # JetDeepJet (needs to check if associated jet) #
    firstHasAssociatedJet = sel.refine(suffix+"firstHasAssociatedJet",cut=[self.lambda_hasAssociatedJet(dilepton[0])])
    firstHasNotAssociatedJet = sel.refine(suffix+"firstHasNotAssociatedJet",cut=[op.NOT(self.lambda_hasAssociatedJet(dilepton[0]))])
    secondHasAssociatedJet = sel.refine(suffix+"secondHasAssociatedJet",cut=[self.lambda_hasAssociatedJet(dilepton[1])])
    secondHasNotAssociatedJet = sel.refine(suffix+"secondHasNotAssociatedJet",cut=[op.NOT(self.lambda_hasAssociatedJet(dilepton[1]))])

    firstDeepJet_has = Plot.make1D("%s_%s_firstlepton_checkhasJetDeepFlavB"%(channel,suffix), 
                                 dilepton[0].jet.btagDeepFlavB, 
                                 firstHasAssociatedJet, 
                                 EquidistantBinning(50,-1,1.),
                                 xTitle= "jetDeepJet (first lepton)",
                                 plotopts = channelLabel)
    firstDeepJet_hasNot = Plot.make1D("%s_%s_firstlepton_checkhasNotJetDeepFlavB"%(channel,suffix),
                                 op.c_float(-1.),
                                 firstHasNotAssociatedJet, 
                                 EquidistantBinning(50,-1.,1.),
                                 xTitle= "jetDeepJet (first lepton)",
                                 plotopts = channelLabel)
    plots.append(SummedPlot("%s_%s_firstlepton_checkjetDeepJet"%(channel,suffix),
                            [firstDeepJet_has,firstDeepJet_hasNot],
                            xTitle="jetDeepJet (first lepton)"))
    secondDeepJet_has = Plot.make1D("%s_%s_secondlepton_checkhasJetDeepFlavB"%(channel,suffix), 
                                 dilepton[1].jet.btagDeepFlavB, 
                                 secondHasAssociatedJet, 
                                 EquidistantBinning(50,0.,1.),
                                 xTitle= "jetDeepJet (second lepton)",
                                 plotopts = channelLabel)
    secondDeepJet_hasNot = Plot.make1D("%s_%s_secondlepton_checkhasNotJetDeepFlavB"%(channel,suffix),
                                 op.c_float(-1.),
                                 secondHasNotAssociatedJet, 
                                 EquidistantBinning(50,0.,1.),
                                 xTitle= "jetDeepJet (second lepton)",
                                 plotopts = channelLabel)
    plots.append(SummedPlot("%s_%s_secondlepton_checkjetDeepJet"%(channel,suffix),
                            [secondDeepJet_has,secondDeepJet_hasNot],
                            xTitle="jetDeepJet (second lepton)",
                            plotopts = channelLabel))


    return plots

########################## DELTA PLOTS #################################
def makeDeltaRPlots(sel,cont1,cont2,suffix,channel,isMC):
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

########################## BJETS PLOTS #################################
def makeAk4BJetsPlots(sel,bjets,lightjets,alljets,suffix,channel):
    plots = []

    channelLabel = channelTitleLabel(channel)

    # Selection for bjets and mixed categories #
    bjetsCategory = sel.refine(channel+suffix+"bjetsCategory",cut=[op.rng_len(bjets) >= 2]) # TODO : might want to veto more than 2 bjets 
    mixedCategory = sel.refine(channel+suffix+"mixedCategory",cut=[op.rng_len(bjets) == 1])

    ## Inclusive plots (len(alljets)>=2 by definition ): leading = highest Pt jet, subleading : second highest Pt jet #
    plots.extend(makeSeparateJetsPlots(sel           = sel,
                                       leadjet       = alljets[0],
                                       subleadjet    = alljets[1],
                                       suffix        = suffix,
                                       channel       = channel,
                                       tag     = "inclusive"))

    ## Bjets plots (if len(bjets)>=2) : leading = highest Pt bjet, subleading : second highest Pt bjet #
    # Must check len(bjets)>=2
    plots.extend(makeSeparateJetsPlots(sel           = bjetsCategory,
                                       leadjet       = bjets[0],
                                       subleadjet    = bjets[1],
                                       suffix        = suffix,
                                       channel       = channel,
                                       tag     = "bjets"))

    ## Mixed plots (if len(bjets)==1, >0 by definition) : leading = bjet, subleading = highest Pt non-bjet
    # Must check len(bjets)==1
    plots.extend(makeSeparateJetsPlots(sel           = mixedCategory,
                                       leadjet       = bjets[0],
                                       subleadjet    = lightjets[0],
                                       suffix        = suffix,
                                       channel       = channel,
                                       tag     = "mixed"))

    return plots

##########################  JETS SEPARATE PLOTS #################################
def makeAk4JetsPlots(sel, leadjet, subleadjet, suffix, channel, lead_is_b=False, sublead_is_b=False):
    """
    Make basic plots
    sel         = refine selection 
    leadjet     = Depends on scenario
    subleadjet  = Depends on scenario
    suffix      = string identifying the selection 
    channel     = string identifying the channel of the dilepton (can be "NoChannel")
    lead_is_b   = boolean for the titles
    sublead_is_b= boolean for the titles
    """
 
    plots = []

    channelLabel = channelTitleLabel(channel)

    if lead_is_b:
        lead_base_name      = "%s_%s_leadbjet_{var}"%(channel,suffix)
        lead_base_title     = "leading bjet"
    else:
        lead_base_name      = "%s_%s_leadjet_{var}"%(channel,suffix)
        lead_base_title     = "leading jet"
    if sublead_is_b:
        sublead_base_name   = "%s_%s_subleadbjet_{var}"%(channel,suffix)
        sublead_base_title  = "subleading bjet"
    else:
        sublead_base_name   = "%s_%s_subleadjet_{var}"%(channel,suffix)
        sublead_base_title  = "subleading jet"

    # leadjet plots #
    plots.append(Plot.make1D(lead_base_name.format(var="pt"),
                             leadjet.p4.Pt(),
                             sel,
                             EquidistantBinning(50,0.,300.),
                             title='Transverse momentum of the %s'%lead_base_title,
                             xTitle="P_{T}(%s) [GeV]"%lead_base_title,
                             plotopts = channelLabel))
    plots.append(Plot.make1D(lead_base_name.format(var="eta"),
                             leadjet.eta,
                             sel,
                             EquidistantBinning(20,-3.,3.),
                             title='Pseudorapidity of the %s'%lead_base_title,
                             xTitle="#eta(%s)"%lead_base_title,
                             plotopts = channelLabel))
    plots.append(Plot.make1D(lead_base_name.format(var="phi"),
                             leadjet.phi,
                             sel,
                             EquidistantBinning(20,-3.2,3.2),
                             title='Azimutal angle of the %s'%lead_base_title,
                             xTitle="#phi(%s)"%lead_base_title,
                             plotopts = channelLabel))
    # subleadjet plots #
    plots.append(Plot.make1D(sublead_base_name.format(var="pt"),
                             subleadjet.pt,
                             sel,
                             EquidistantBinning(50,0.,300.),
                             title='Transverse momentum of the %s'%sublead_base_title,
                             xTitle="P_{T}(%s) [GeV]"%sublead_base_title,
                             plotopts = channelLabel))
    plots.append(Plot.make1D(sublead_base_name.format(var="eta"),
                             subleadjet.eta,
                             sel,
                             EquidistantBinning(20,-3.,3.),
                             title='Pseudorapidity of the %s'%sublead_base_title,
                             xTitle="#eta(%s)"%sublead_base_title,
                             plotopts = channelLabel))
    plots.append(Plot.make1D(sublead_base_name.format(var="phi"),
                             subleadjet.phi,
                             sel,
                             EquidistantBinning(20,-3.2,3.2),
                             title='Azimutal angle of the %s'%sublead_base_title,
                             xTitle="#phi(%s)"%sublead_base_title,
                             plotopts = channelLabel))
    # Dijet Pt plot #
    plots.append(Plot.make1D("%s_%s_dijet_pt"%(channel,suffix), 
                             (leadjet.p4+subleadjet.p4).Pt(), 
                             sel, 
                             EquidistantBinning(50,0.,300.),
                             title="Transverse momentum of the dijet (channel %s)"%channel, 
                             xTitle= "P_{T}(dijet) [GeV]",
                             plotopts = channelLabel))
    # DeltaPhi plot #
    plots.append(Plot.make1D("%s_%s_dijet_deltaPhi"%(channel,suffix), 
                             op.abs(op.deltaPhi(leadjet.p4,subleadjet.p4)), 
                             sel, 
                             EquidistantBinning(20,0, 3.2),
                             title="Azimutal angle difference of the dijet (channel %s)"%channel,
                             xTitle= "| #Delta \phi (dijet)|",
                             plotopts = channelLabel))

    # DeltaR plot #
    plots.append(Plot.make1D("%s_%s_dijet_deltaR"%(channel,suffix), 
                             op.deltaR(leadjet.p4,subleadjet.p4), 
                             sel, 
                             EquidistantBinning(50, 0., 5.), 
                             title="Dilepton Delta R (channel %s)"%channel, 
                             xTitle= "#Delta R (dijet)",
                             plotopts = channelLabel))
    # invariant mass plot #
    plots.append(Plot.make1D("%s_%s_dijet_invariantMass"%(channel,suffix),
                             op.invariant_mass(leadjet.p4,subleadjet.p4),
                             sel,
                             EquidistantBinning(50, 0., 500.), 
                             title="Dijet invariant mass (channel %s)"%channel, 
                             xTitle= "Invariant mass (dijet) [GeV]",
                             plotopts = channelLabel))

    return plots

##########################  JETS (AK4) PLOT #################################
def makeAk8JetsPlots(sel, fatjet, suffix, channel):
    """
    Make fatjet subjet basic plots
    sel         = refine selection 
    fatjet      = fatjet (AK8 jet) 
    suffix      = string identifying the selecton 
    channel     = string identifying the channel of the dilepton (can be "NoChannel")
    """
 
    plots = []

    channelLabel = channelTitleLabel(channel)

    # fatjet plots (always present by selection) #
    plots.append(Plot.make1D("%s_%s_fatjet_pt"%(channel,suffix),
                             fatjet.p4.pt(),
                             sel,
                             EquidistantBinning(50,200,600.),
                             title='Transverse momentum of the fatjet',
                             xTitle="P_{T}(fatjet) [GeV]",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_fatjet_eta"%(channel,suffix),
                             fatjet.p4.eta(),
                             sel,
                             EquidistantBinning(10,-3.,3.),
                             title='Pseudorapidity of the fatjet',
                             xTitle="#eta(fatjet)",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_fatjet_phi"%(channel,suffix),
                             fatjet.p4.phi(),
                             sel,
                             EquidistantBinning(10,-3.2,3.2),
                             title='Azimutal angle of the fatjet',
                             xTitle="#phi(fatjet)",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_fatjet_mass"%(channel,suffix),
                             fatjet.mass,
                             sel,
                             EquidistantBinning(50,0.,200.),
                             title='Invariant mass of the fatjet',
                             xTitle="M(fatjet)",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_fatjet_softdropmass"%(channel,suffix),
                             fatjet.msoftdrop,
                             sel,
                             EquidistantBinning(50,0.,200.),
                             title='Soft Drop mass of the fatjet',
                             xTitle="M_{Soft Drop}(fatjet) [GeV]",
                             plotopts = channelLabel))

    return plots

#########################  High-level quantities ################################
def makeHighLevelQuantities(sel,met,l1,l2,j1,j2,suffix,channel):
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

    # Invariant mass plots #
    plots.append(Plot.make1D("%s_%s_highlevelvariable_DilepJetInvariantMass"%(channel,suffix),
                             op.invariant_mass(lljj_p4(l1,l2,j1,j2)),    
                             sel,
                             EquidistantBinning(50,0.,200.),
                             title='Dilepton-jets invariant mass (%s channel)'%channel,
                             xTitle="M(lljj)",
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

