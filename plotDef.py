import math
from bamboo.plots import Plot, EquidistantBinning, SummedPlot
from bamboo import treefunctions as op
from bamboo.plots import CutFlowReport

from highlevelLambdas import *

# TODO : remove the self

########################   Channel title   #############################
def DoubleLeptonChannelTitleLabel(channel):
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
    channelLabel = DoubleLeptonChannelTitleLabel(channel)
    return Plot.make1D(channel+"_"+suffix+"_"+objName+"_N",
                       op.rng_len(objCont),
                       sel,
                       EquidistantBinning(Nmax,0.,Nmax),
                       xTitle = xTitle,
                       plotopts = channelLabel)

#######################   Channel Plot   ###############################
def doubleLeptonChannelPlot(sel,DilepElEl,DilepMuMu,DilepElMu,suffix,channel):
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
    def __init__(self,unitYield=False):
        self.calls = 0
        self.plots = []
        self.unitYield = unitYield
    def addYield(self, sel, name, title):
        """
        Make Yield plot and use it also in the latex yield table
        sel     = refine selection
        name    = name of the PDF to be produced
        title   = title that will be used in the LateX yield table
        """
        if self.unitYield:
            sel = sel.refine("Yield_"+name,weight=op.c_float(1.)/sel.weight)
        self.plots.append(Plot.make1D("Yield_"+name,   
                                      op.c_int(0),
                                      sel,
                                      EquidistantBinning(1, 0., 1.),
                                      title = "Yield_"+name,
                                      xTitle = "Yield_"+name,
                                      plotopts = {"for-yields":True, "yields-title":title, 'yields-table-order':self.calls}))
        self.calls += 1

    def returnPlots(self):
        return self.plots 

######################   TRIGGER PLOTS  ##############################
def triggerPlots(sel,triggerDict,suffix,channel):
    plots = []
    channelLabel = DoubleLeptonChannelTitleLabel(channel)
    
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

    channelLabel = DoubleLeptonChannelTitleLabel(channel)

    # PT plot #
    plots.append(Plot.make1D("%s_%s_met_pt"%(channel,suffix), 
                             met.pt, 
                             sel, 
                             EquidistantBinning(40, 0., 200.), 
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
def makeDileptonPlots(sel, dilepton, suffix, channel, is_MC=False):
    """
    Make dilepton basic plots
    sel         = refine selection 
    dilepton    = dilepton object (ElEl, MuEl, ElMu, MuMu)
    suffix      = string identifying the selecton 
    channel     = string identifying the channel of the dilepton (can be "NoChannel")
    """
    plots = []

    channelLabel = DoubleLeptonChannelTitleLabel(channel)

    # PT plot #
    plots.append(Plot.make1D("%s_%s_firstlepton_pt"%(channel,suffix), 
                             dilepton[0].pt, 
                             sel, 
                             EquidistantBinning(60,0.,300.),
                             title="Transverse momentum of the first lepton (channel %s)"%channel, 
                             xTitle= "P_{T} (first lepton) [GeV]",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_secondlepton_pt"%(channel,suffix), 
                             dilepton[1].pt, 
                             sel, 
                             EquidistantBinning(40,0.,200.),
                             title="Transverse momentum of the second lepton (channel %s)"%channel, 
                             xTitle= "P_{T} (second lepton) [GeV]",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_dilepton_pt"%(channel,suffix), 
                             (dilepton[0].p4+dilepton[1].p4).Pt(), 
                             sel, 
                             EquidistantBinning(60,0.,300.),
                             title="Transverse momentum of the dilepton (channel %s)"%channel, 
                             xTitle= "P_{T} (dilepton) [GeV]",
                             plotopts = channelLabel))

    # PT-eta plots #
    plots.append(Plot.make2D("%s_%s_firstlepton_ptVSeta"%(channel,suffix), 
                             [dilepton[0].eta, dilepton[0].pt],
                             sel, 
                             [EquidistantBinning(22, -3., 3.),EquidistantBinning(60,0.,300.)],
                             xTitle= "P_{T} (first lepton) [GeV]",
                             yTitle= "#eta (second lepton)",
                             plotopts = channelLabel))

    # Eta plot #
    plots.append(Plot.make1D("%s_%s_firstlepton_eta"%(channel,suffix), 
                             dilepton[0].eta, 
                             sel, 
                             EquidistantBinning(22, -3., 3.), 
                             title="Pseudorapidity of the first lepton (channel %s)"%channel, 
                             xTitle= "#eta (first lepton)",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_secondlepton_eta"%(channel,suffix), 
                             dilepton[1].eta, 
                             sel, 
                             EquidistantBinning(22, -3., 3.), 
                             title="Pseudorapidity of the second lepton (channel %s)"%channel, 
                             xTitle= "#eta (second lepton)",
                             plotopts = channelLabel))

    # Phi plot #
    plots.append(Plot.make1D("%s_%s_firstlepton_phi"%(channel,suffix), 
                             dilepton[0].phi, 
                             sel, 
                             EquidistantBinning(16, -3.2, 3.2), 
                             title="Azimutal angle of the first lepton (channel %s)"%channel, 
                             xTitle= "#phi (first lepton)",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_secondlepton_phi"%(channel,suffix), 
                             dilepton[1].phi, 
                             sel, 
                             EquidistantBinning(16, -3.2, 3.2), 
                             title="Azimutal angle of the second lepton (channel %s)"%channel, 
                             xTitle= "#phi (second lepton)",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_dilepton_deltaPhi"%(channel,suffix), 
                             op.abs(op.deltaPhi(dilepton[0].p4,dilepton[1].p4)),
                             sel, 
                             EquidistantBinning(32,0, 3.2), 
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
    plots.append(Plot.make1D("%s_%s_firstlepton_genPartFlav"%(channel,suffix), 
                             dilepton[0].genPartFlav if is_MC else op.c_int(-1),
                             sel, 
                             EquidistantBinning(23, -1., 22.), 
                             title="Flavour of genParticle (channel %s)"%channel, 
                             xTitle= "GenParticle flavour (first lepton)",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_secondlepton_genPartFlav"%(channel,suffix), 
                             dilepton[1].genPartFlav if is_MC else op.c_int(-1),
                             sel, 
                             EquidistantBinning(23, -1., 22.), 
                             title="Flavour of genParticle (channel %s)"%channel, 
                             xTitle= "GenParticle flavour (second lepton)",
                             plotopts = channelLabel))

    return plots

########################## DELTA PLOTS #################################
def makeDeltaRPlots(sel,cont1,cont2,suffix,channel,isMC):
    plots = [] 

    channelLabel = DoubleLeptonChannelTitleLabel(channel)

    mixedCont = op.combine((cont1,cont2))
    plots.append(Plot.make1D("%s_%s_DeltaR"%(channel,suffix), 
                             op.map(mixedCont,lambda c : op.deltaR(c[0].p4,c[1].p4)),
                             sel, 
                             EquidistantBinning(100,0.,5.), 
                             title="#Delta R %s"%channel, 
                             xTitle= "#Delta R",
                             plotopts = channelLabel))

    return plots

##########################  JETS SEPARATE PLOTS #################################
def makeTwoAk4JetsPlots(sel, leadjet, subleadjet, suffix, channel, lead_is_b=False, sublead_is_b=False,is_MC=False):
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

    channelLabel = DoubleLeptonChannelTitleLabel(channel)

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
                             leadjet.pt,
                             sel,
                             EquidistantBinning(40,0.,200.),
                             title='Transverse momentum of the %s'%lead_base_title,
                             xTitle="P_{T}(%s) [GeV]"%lead_base_title,
                             plotopts = channelLabel))
    plots.append(Plot.make1D(lead_base_name.format(var="eta"),
                             leadjet.eta,
                             sel,
                             EquidistantBinning(22,-3.,3.),
                             title='Pseudorapidity of the %s'%lead_base_title,
                             xTitle="#eta(%s)"%lead_base_title,
                             plotopts = channelLabel))
    plots.append(Plot.make1D(lead_base_name.format(var="phi"),
                             leadjet.phi,
                             sel,
                             EquidistantBinning(16,-3.2,3.2),
                             title='Azimutal angle of the %s'%lead_base_title,
                             xTitle="#phi(%s)"%lead_base_title,
                             plotopts = channelLabel))

    plots.append(Plot.make1D(lead_base_name.format(var="hadronflavour"),
                             op.abs(leadjet.hadronFlavour) if is_MC else op.c_int(-1),
                             sel,
                             EquidistantBinning(23,-1.,22.),
                             title='Hadron flavour of the %s'%lead_base_title,
                             xTitle="Hadron flavour (%s) [GeV]"%lead_base_title,
                             plotopts = channelLabel))
    plots.append(Plot.make1D(lead_base_name.format(var="partonflavour"),
                             op.abs(leadjet.partonFlavour) if is_MC else op.c_int(-1),
                             sel,
                             EquidistantBinning(23,-1.,22.),
                             title='Parton flavour of the %s'%lead_base_title,
                             xTitle="Parton flavour (%s) [GeV]"%lead_base_title,
                             plotopts = channelLabel))


    # subleadjet plots #
    plots.append(Plot.make1D(sublead_base_name.format(var="pt"),
                             subleadjet.pt,
                             sel,
                             EquidistantBinning(40,0.,200.),
                             title='Transverse momentum of the %s'%sublead_base_title,
                             xTitle="P_{T}(%s) [GeV]"%sublead_base_title,
                             plotopts = channelLabel))
    plots.append(Plot.make1D(sublead_base_name.format(var="eta"),
                             subleadjet.eta,
                             sel,
                             EquidistantBinning(22,-3.,3.),
                             title='Pseudorapidity of the %s'%sublead_base_title,
                             xTitle="#eta(%s)"%sublead_base_title,
                             plotopts = channelLabel))
    plots.append(Plot.make1D(sublead_base_name.format(var="phi"),
                             subleadjet.phi,
                             sel,
                             EquidistantBinning(16,-3.2,3.2),
                             title='Azimutal angle of the %s'%sublead_base_title,
                             xTitle="#phi(%s)"%sublead_base_title,
                             plotopts = channelLabel))
    plots.append(Plot.make1D(sublead_base_name.format(var="hadronflavour"),
                             op.abs(subleadjet.hadronFlavour) if is_MC else op.c_int(-1),
                             sel,
                             EquidistantBinning(23,-1.,22.),
                             title='Hadron flavour of the %s'%sublead_base_title,
                             xTitle="Hadron flavour (%s) [GeV]"%sublead_base_title,
                             plotopts = channelLabel))
    plots.append(Plot.make1D(sublead_base_name.format(var="partonflavour"),
                             op.abs(subleadjet.partonFlavour) if is_MC else op.c_int(-1),
                             sel,
                             EquidistantBinning(23,-1.,22.),
                             title='Parton flavour of the %s'%sublead_base_title,
                             xTitle="Parton flavour (%s) [GeV]"%sublead_base_title,
                             plotopts = channelLabel))

    # Dijet Pt plot #
    plots.append(Plot.make1D("%s_%s_dijet_pt"%(channel,suffix), 
                             (leadjet.p4+subleadjet.p4).Pt(), 
                             sel, 
                             EquidistantBinning(80,0.,400.),
                             title="Transverse momentum of the dijet (channel %s)"%channel, 
                             xTitle= "P_{T}(dijet) [GeV]",
                             plotopts = channelLabel))
    # DeltaPhi plot #
    plots.append(Plot.make1D("%s_%s_dijet_deltaPhi"%(channel,suffix), 
                             op.abs(op.deltaPhi(leadjet.p4,subleadjet.p4)), 
                             sel, 
                             EquidistantBinning(32,0.,3.2),
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

    channelLabel = DoubleLeptonChannelTitleLabel(channel)

    # fatjet plots (always present by selection) #
    plots.append(Plot.make1D("%s_%s_fatjet_pt"%(channel,suffix),
                             fatjet.p4.pt(),
                             sel,
                             EquidistantBinning(60,0.,600.),
                             title='Transverse momentum of the fatjet',
                             xTitle="P_{T}(fatjet) [GeV]",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_fatjet_eta"%(channel,suffix),
                             fatjet.p4.eta(),
                             sel,
                             EquidistantBinning(12,-2.4,2.4),
                             title='Pseudorapidity of the fatjet',
                             xTitle="#eta(fatjet)",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_fatjet_phi"%(channel,suffix),
                             fatjet.p4.phi(),
                             sel,
                             EquidistantBinning(8,-3.2,3.2),
                             title='Azimutal angle of the fatjet',
                             xTitle="#phi(fatjet)",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_fatjet_mass"%(channel,suffix),
                             fatjet.mass,
                             sel,
                             EquidistantBinning(40,0.,200.),
                             title='Invariant mass of the fatjet',
                             xTitle="M(fatjet)",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_fatjet_softdropmass"%(channel,suffix),
                             fatjet.msoftdrop,
                             sel,
                             EquidistantBinning(40,0.,200.),
                             title='Soft Drop mass of the fatjet',
                             xTitle="M_{Soft Drop}(fatjet) [GeV]",
                             plotopts = channelLabel))

    return plots

#########################  High-level quantities ################################
def makeHighLevelQuantities(sel,met,l1,l2,j1,j2,suffix,channel):
    plots = []

    channelLabel = DoubleLeptonChannelTitleLabel(channel)

    # Lepton PT versus Jet Pt #
    plots.append(Plot.make2D("%s_%s_highlevelvariable_firstLeptonPtVSLeadjetPt"%(channel,suffix),
                             [l1.pt,j1.pt],
                             sel,
                             [EquidistantBinning(60,0.,300.),EquidistantBinning(60,0.,300.)],
                             xTitle="First lepton P_{T}",
                             yTitle="Leading jet P_{T}",
                             plotopts = channelLabel))

    # dilepton-MET plots #
    plots.append(Plot.make1D("%s_%s_highlevelvariable_DilepMETdeltaPhi"%(channel,suffix),
                             op.abs(DilepMET_deltaPhi(l1,l2,met)),
                             sel,
                             EquidistantBinning(32,0.,3.2),
                             title='Azimutal angle between dilepton and MET (%s channel)'%channel,
                             xTitle="|#Delta \phi (ll,MET)|",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_highlevelvariable_DilepMETpt"%(channel,suffix),
                             DilepMET_Pt(l1,l2,met),
                             sel,
                             EquidistantBinning(50,0.,500.),
                             title='Transverse momentum of dilepton and MET (%s channel)'%channel,
                             xTitle="P_{T}(ll,MET) [GeV]",
                             plotopts = channelLabel))

    # Invariant mass plots #
    plots.append(Plot.make1D("%s_%s_highlevelvariable_DilepJetInvariantMass"%(channel,suffix),
                             M_lljj(l1,l2,j1,j2),
                             sel,
                             EquidistantBinning(50,0.,1000.),
                             title='Dilepton-jets invariant mass (%s channel)'%channel,
                             xTitle="M_{lljj}",
                             plotopts = channelLabel))


    # Transverse mass plots #
        # mT for dilepton #
    plots.append(Plot.make1D("%s_%s_highlevelvariable_MTll"%(channel,suffix),
                             MT_ll(l1,l2,met),
                             sel,
                             EquidistantBinning(50,0.,1000.),
                             title='Transverse mass of dilepton and MET (%s channel)'%channel,
                             xTitle="M_{T}(ll,MET) [GeV]",
                             plotopts = channelLabel))
        # mT for lljj #
    plots.append(Plot.make1D("%s_%s_highlevelvariable_MTlljj"%(channel,suffix),
                             MT_lljj(l1,l2,j1,j2,met),
                             sel,
                             EquidistantBinning(50,0.,1000.),
                             title='Transverse mass of dilepton+dijet and MET (%s channel)'%channel,
                             xTitle="M_{T}(lljj,MET) [GeV]",
                             plotopts = channelLabel))


    # Scalar magnitude sum #

    plots.append(Plot.make1D("%s_%s_highlevelvariable_HT2"%(channel,suffix),
                             HT2(l1,l2,j1,j2,met),
                             sel,
                             EquidistantBinning(60,0.,600.),
                             title='Di-Higgs magnitude (%s channel)'%channel,
                             xTitle="H_{T2} [GeV]",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_highlevelvariable_HT2R"%(channel,suffix),
                             HT2R(l1,l2,j1,j2,met),
                             sel,
                             EquidistantBinning(50,0.,1.),
                             title='Di-Higgs magnitude ratio (%s channel)'%channel,
                             xTitle="H_{T2}^{R} [GeV]",
                             plotopts = channelLabel))

    return plots 

