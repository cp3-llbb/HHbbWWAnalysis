import math
from bamboo.plots import Plot, EquidistantBinning, SummedPlot, VariableBinning
from bamboo import treefunctions as op
from bamboo.plots import CutFlowReport
from highlevelLambdas import highlevelLambdas

########################   Channel title   #############################
def SingleLeptonChannelTitleLabel(channel):
    if (channel == "El"):
        channel = "e^{+,-}"
    elif (channel == "Mu"):
        channel ="#mu^{+,-}"
    channelLabel = {'labels': [{'text': channel, 'position': [0.23, 0.87], 'size': 36}]}
    return channelLabel 

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
def singleLeptonTriggerPlots(sel,triggerDict,suffix,channel):
    plots = []
    channelLabel = SingleLeptonChannelTitleLabel(channel)
    
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
    return plots

def doubleLeptonTriggerPlots(sel,triggerDict,suffix,channel):
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


##########################  LHE PLOT #################################
def makeLHEPlots(sel, LHE):
    plots = []
    plots.append(Plot.make1D("LHE_Njets",
                             LHE.Njets,
                             sel,
                             EquidistantBinning(5,0.,5),
                             xTitle='LHE Njets'))
    plots.append(Plot.make1D("LHE_HT",
                             LHE.HT,
                             sel,
                             VariableBinning([70,100,200,400,600,800,1200,2500]),
                             xTitle='LHE HT'))
    plots.append(Plot.make1D("LHE_HT_finer",
                             LHE.HT,
                             sel,
                             EquidistantBinning(40, 0., 400.), 
                             xTitle='LHE HT'))
    plots.append(Plot.make2D("LHE_HTVSNjets",
                             [LHE.HT,LHE.Njets],
                             sel,
                             [VariableBinning([0,70,100,200,400,600,800,1200,2500]),EquidistantBinning(5,0.,5)],
                             yTitle='LHE Njets',
                             xTitle='LHE HT'))
    plots.append(Plot.make2D("LHE_HTVSNjets_finer",
                             [LHE.HT,LHE.Njets],
                             sel,
                             [EquidistantBinning(40, 0., 400.),EquidistantBinning(5,0.,5)],
                             yTitle='LHE Njets',
                             xTitle='LHE HT'))
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

##########################  LEPTON PLOT #################################
def makeSinleptonPlots(sel, lep, suffix, channel, is_MC=False):
    """
    Make dilepton basic plots
    sel         = refine selection 
    sinlepton   = sinlepton object (El, Mu)
    suffix      = string identifying the selecton 
    channel     = string identifying the channel of the sinlepton (can be "NoChannel")
    """
    plots = []

    channelLabel = SingleLeptonChannelTitleLabel(channel)

    # PT plot #
    plots.append(Plot.make1D("%s_%s_lepton_pt"%(channel,suffix), 
                             lep.pt, 
                             sel, 
                             EquidistantBinning(60,0.,300.),
                             title="Transverse momentum of the lepton (channel %s)"%channel, 
                             xTitle= "P_{T} (lepton) [GeV]",
                             plotopts = channelLabel))

    # Eta plot #
    plots.append(Plot.make1D("%s_%s_lepton_eta"%(channel,suffix), 
                             lep.eta, 
                             sel, 
                             EquidistantBinning(22, -3., 3.), 
                             title="Pseudorapidity of the lepton (channel %s)"%channel, 
                             xTitle= "#eta (lepton)",
                             plotopts = channelLabel))

    # PT-eta plots #
    plots.append(Plot.make2D("%s_%s_lepton_ptVSeta"%(channel,suffix), 
                             [lep.pt, lep.eta],
                             sel, 
                             [EquidistantBinning(60,0.,300.),EquidistantBinning(22, -3., 3.)],
                             xTitle= "P_{T} (lepton) [GeV]",
                             yTitle= "#eta (lepton)",
                             plotopts = channelLabel))

    # Phi plot #
    plots.append(Plot.make1D("%s_%s_lepton_phi"%(channel,suffix), 
                             lep.phi, 
                             sel, 
                             EquidistantBinning(20, -3.2, 3.2), 
                             title="Azimutal angle of the lepton (channel %s)"%channel, 
                             xTitle= "#phi (lepton)",
                             plotopts = channelLabel))

    # GenPartFlav (if isMC) #
    plots.append(Plot.make1D("%s_%s_lepton_genPartFlav"%(channel,suffix), 
                             lep.genPartFlav if is_MC else op.c_int(-1),
                             sel, 
                             EquidistantBinning(23, -1., 22.), 
                             title="Flavour of genParticle (channel %s)"%channel, 
                             xTitle= "GenParticle flavour (lepton)",
                             plotopts = channelLabel))

    return plots


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
    # ttH mva plot #
    plots.append(Plot.make1D("%s_%s_firstlepton_ttHmva"%(channel,suffix), 
                             dilepton[0].mvaTTH, 
                             sel, 
                             EquidistantBinning(50,0.,1.),
                             title="ttH MVA of the first lepton (channel %s)"%channel, 
                             xTitle= "MVA_{ttH} (first lepton) [GeV]",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_secondlepton_ttHmva"%(channel,suffix), 
                             dilepton[1].mvaTTH, 
                             sel, 
                             EquidistantBinning(50,0.,1.),
                             title="ttH MVA of the second lepton (channel %s)"%channel, 
                             xTitle= "MVA_{ttH} (second lepton) [GeV]",
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
def makeAk4JetsPlots (sel,j1,j2,j3,j4,channel,suffix,nJet,nbJet,HLL,is_MC=False):

    #print ("===============>>>Ak4Jet Plots__channel:%s__sel:%s"%(channel,suffix))
    plots = []

    channelLabel = SingleLeptonChannelTitleLabel(channel)
    if nbJet == 0:
        j1_name = "leadAk4"
        j2_name = "subLeadAk4"
        j3_name = "thirdAk4"
        j4_name = "fourthAk4"
    else:
        j1_name = "leadAk4B"
        j2_name = "leadAk4"if nbJet == 1 else "subLeadAk4B"
        j3_name = "subLeadAk4" if nbJet == 1 else "leadAk4"
        j4_name = "thirdAk4" if nbJet == 1 else "subLeadAk4"

    
    # leadjet plots #
    plots.append(Plot.make1D("%s_%s_pT_%s"%(channel,suffix,j1_name),
                             #bJetCorrPt(j1),
                             HLL.bJetCorrP4(j1).Pt(),
                             sel,
                             EquidistantBinning(40,0.,200.),
                             title='Transverse momentum of the %s'%j1_name,
                             xTitle="P_{T}(%s) [GeV]"%j1_name,
                             plotopts = channelLabel))
    
    plots.append(Plot.make1D("%s_%s_eta_%s"%(channel,suffix,j1_name),
                             j1.eta,
                             sel,
                             EquidistantBinning(22,-3.,3.),
                             title='Pseudorapidity of the %s'%j1_name,
                             xTitle="#eta(%s)"%j1_name,
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_phi_%s"%(channel,suffix,j1_name),
                             j1.phi,
                             sel,
                             EquidistantBinning(20,-3.2,3.2),
                             title='Azimutal angle of the %s'%j1_name,
                             xTitle="#phi(%s)"%j1_name,
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_hadronFlv_%s"%(channel,suffix,j1_name),
                             op.abs(j1.hadronFlavour) if is_MC else op.c_int(-1),
                             sel,
                             EquidistantBinning(23,-1.,22.),
                             title='Hadron flavour of the %s'%j1_name,
                             xTitle="Hadron flavour (%s) [GeV]"%j1_name,
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_partonFlv_%s"%(channel,suffix,j1_name),
                             op.abs(j1.partonFlavour) if is_MC else op.c_int(-1),
                             sel,
                             EquidistantBinning(23,-1.,22.),
                             title='Parton flavour of the %s'%j1_name,
                             xTitle="Parton flavour (%s) [GeV]"%j1_name,
                             plotopts = channelLabel))



    # subleadjet plots #
    plots.append(Plot.make1D("%s_%s_pT_%s"%(channel,suffix,j2_name),
                             HLL.bJetCorrP4(j2).Pt(),
                             sel,
                             EquidistantBinning(40,0.,200.),
                             title='Transverse momentum of the %s'%j2_name,
                             xTitle="P_{T}(%s) [GeV]"%j2_name,
                             plotopts = channelLabel))
        
    plots.append(Plot.make1D("%s_%s_eta_%s"%(channel,suffix,j2_name),
                             j2.eta,
                             sel,
                             EquidistantBinning(22,-3.,3.),
                             title='Pseudorapidity of the %s'%j2_name,
                             xTitle="#eta(%s)"%j2_name,
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_phi_%s"%(channel,suffix,j2_name),
                             j2.phi,
                             sel,
                             EquidistantBinning(20,-3.2,3.2),
                             title='Azimutal angle of the %s'%j2_name,
                             xTitle="#phi(%s)"%j2_name,
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_hadronFlv_%s"%(channel,suffix,j2_name),
                             op.abs(j2.hadronFlavour) if is_MC else op.c_int(-1),
                             sel,
                             EquidistantBinning(23,-1.,22.),
                             title='Hadron flavour of the %s'%j2_name,
                             xTitle="Hadron flavour (%s) [GeV]"%j2_name,
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_partonFlv_%s"%(channel,suffix,j2_name),
                             op.abs(j2.partonFlavour) if is_MC else op.c_int(-1),
                             sel,
                             EquidistantBinning(23,-1.,22.),
                             title='Parton flavour of the %s'%j2_name,
                             xTitle="Parton flavour (%s) [GeV]"%j2_name,
                             plotopts = channelLabel))
    

    # thirdjet plots #
    plots.append(Plot.make1D("%s_%s_pT_%s"%(channel,suffix,j3_name),
                             j3.pt,
                             sel,
                             EquidistantBinning(40,0.,200.),
                             title='Transverse momentum of the %s'%j3_name,
                             xTitle="P_{T}(%s) [GeV]"%j3_name,
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_eta_%s"%(channel,suffix,j3_name),
                             j3.eta,
                             sel,
                             EquidistantBinning(22,-3.,3.),
                             title='Pseudorapidity of the %s'%j3_name,
                             xTitle="#eta(%s)"%j3_name,
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_phi_%s"%(channel,suffix,j3_name),
                             j3.phi,
                             sel,
                             EquidistantBinning(20,-3.2,3.2),
                             title='Azimutal angle of the %s'%j3_name,
                             xTitle="#phi(%s)"%j3_name,
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_hadronFlv_%s"%(channel,suffix,j3_name),
                             op.abs(j3.hadronFlavour) if is_MC else op.c_int(-1),
                             sel,
                             EquidistantBinning(23,-1.,22.),
                             title='Hadron flavour of the %s'%j3_name,
                             xTitle="Hadron flavour (%s) [GeV]"%j3_name,
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_partonFlv_%s"%(channel,suffix,j3_name),
                             op.abs(j3.partonFlavour) if is_MC else op.c_int(-1),
                             sel,
                             EquidistantBinning(23,-1.,22.),
                             title='Parton flavour of the %s'%j3_name,
                             xTitle="Parton flavour (%s) [GeV]"%j3_name,
                             plotopts = channelLabel))

    # Dijet Pt plot #
    plots.append(Plot.make1D("%s_%s_DiJetPT_%s_%s"%(channel,suffix,j1_name,j2_name), 
                             #(j1.p4+j2.p4).Pt(), 
                             (HLL.bJetCorrP4(j1)+HLL.bJetCorrP4(j2)).Pt(),
                             sel, 
                             EquidistantBinning(80,0.,400.),
                             title="Transverse momentum of the dijet_[%s_%s] (channel %s)"%(j1_name,j2_name,channel), 
                             xTitle= "DiJetPT [%s_%s]"%(j1_name, j2_name),
                             plotopts = channelLabel))
        
    # DeltaPhi plot #
    plots.append(Plot.make1D("%s_%s_DeltaPhi_%s_%s"%(channel,suffix,j1_name,j2_name), 
                             op.abs(op.deltaPhi(j1.p4,j2.p4)), 
                             sel, 
                             EquidistantBinning(20,0, 3.2),
                             title="Azimutal angle difference of the dijet_[%s_%s] (channel %s)"%(j1_name,j2_name,channel),
                             xTitle= "DeltaPhi [%s_%s]"%(j1_name, j2_name),
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_DeltaPhi_%s_%s"%(channel,suffix,j1_name,j3_name), 
                             op.abs(op.deltaPhi(j1.p4,j3.p4)), 
                             sel, 
                             EquidistantBinning(20,0, 3.2),
                             title="Azimutal angle difference of the dijet_[%s_%s] (channel %s)"%(j1_name,j3_name,channel),
                             xTitle= "DeltaPhi [%s_%s]"%(j1_name, j3_name),
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_DeltaPhi_%s_%s"%(channel,suffix,j2_name,j3_name), 
                             op.abs(op.deltaPhi(j2.p4,j3.p4)), 
                             sel, 
                             EquidistantBinning(20,0, 3.2),
                             title="Azimutal angle difference of the dijet_[%s_%s] (channel %s)"%(j2_name,j3_name,channel),
                             xTitle= "DeltaPhi [%s_%s]"%(j2_name, j3_name),
                             plotopts = channelLabel))

    # DeltaR plot #
    plots.append(Plot.make1D("%s_%s_DeltaR_%s_%s"%(channel,suffix,j1_name,j2_name), 
                             op.deltaR(j1.p4,j2.p4), 
                             sel, 
                             EquidistantBinning(50,0, 5.0),
                             title="DeltaR of the dijet_[%s_%s] (channel %s)"%(j1_name,j2_name,channel),
                             xTitle= "DeltaR [%s_%s]"%(j1_name, j2_name),
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_DeltaR_%s_%s"%(channel,suffix,j1_name,j3_name), 
                             op.deltaR(j1.p4,j3.p4), 
                             sel, 
                             EquidistantBinning(50,0, 5.0),
                             title="DeltaR of the dijet_[%s_%s] (channel %s)"%(j1_name,j3_name,channel),
                             xTitle= "DeltaR [%s_%s]"%(j1_name, j3_name),
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_DeltaR_%s_%s"%(channel,suffix,j2_name,j3_name), 
                             op.deltaR(j2.p4,j3.p4), 
                             sel, 
                             EquidistantBinning(50,0, 5.0),
                             title="DeltaR of the dijet_[%s_%s] (channel %s)"%(j2_name,j3_name,channel),
                             xTitle= "DeltaR [%s_%s]"%(j2_name, j3_name),
                             plotopts = channelLabel))
    
    # minDR
    plots.append(Plot.make1D("%s_%s_minDeltaR_DiJet"%(channel,suffix), 
                             HLL.MinDiJetDRLoose(j1,j2,j3), 
                             sel, 
                             EquidistantBinning(50,0, 5.0),
                             title="MinDeltaR between dijet (channel %s)"%channel,
                             xTitle= "MinDeltaR [Di-Jet]",
                             plotopts = channelLabel))
    
    # invariant mass plot #
    plots.append(Plot.make1D("%s_%s_InvariantMass_%s_%s"%(channel,suffix,j1_name,j2_name),
                             #op.invariant_mass(j1.p4,j2.p4),
                             op.invariant_mass(HLL.bJetCorrP4(j1),HLL.bJetCorrP4(j2)),
                             sel,
                             EquidistantBinning(50, 0., 500.), 
                             title="Invariant Mass of the dijet_[%s_%s] (channel %s)"%(j1_name,j2_name,channel), 
                             xTitle= "InvariantMass [%s_%s]"%(j1_name, j2_name),
                             plotopts = channelLabel))
        
    if nJet == 4:
        # fourthjet plots #
        plots.append(Plot.make1D("%s_%s_pT_%s"%(channel,suffix,j4_name),
                                 j4.pt,
                                 sel,
                                 EquidistantBinning(40,0.,200.),
                                 title='Transverse momentum of the %s'%j4_name,
                                 xTitle="P_{T}(%s) [GeV]"%j4_name,
                                 plotopts = channelLabel))
        plots.append(Plot.make1D("%s_%s_eta_%s"%(channel,suffix,j4_name),
                                 j4.eta,
                                 sel,
                                 EquidistantBinning(22,-3.,3.),
                                 title='Pseudorapidity of the %s'%j4_name,
                                 xTitle="#eta(%s)"%j4_name,
                                 plotopts = channelLabel))
        plots.append(Plot.make1D("%s_%s_phi_%s"%(channel,suffix,j4_name),
                                 j4.phi,
                                 sel,
                                 EquidistantBinning(20,-3.2,3.2),
                                 title='Azimutal angle of the %s'%j4_name,
                                 xTitle="#phi(%s)"%j4_name,
                                 plotopts = channelLabel))
        plots.append(Plot.make1D("%s_%s_hadronFlv_%s"%(channel,suffix,j4_name),
                                 op.abs(j4.hadronFlavour) if is_MC else op.c_int(-1),
                                 sel,
                                 EquidistantBinning(23,-1.,22.),
                                 title='Hadron flavour of the %s'%j4_name,
                                 xTitle="Hadron flavour (%s) [GeV]"%j4_name,
                                 plotopts = channelLabel))
        plots.append(Plot.make1D("%s_%s_partonFlv_%s"%(channel,suffix,j4_name),
                                 op.abs(j4.partonFlavour) if is_MC else op.c_int(-1),
                                 sel,
                                 EquidistantBinning(23,-1.,22.),
                                 title='Parton flavour of the %s'%j4_name,
                                 xTitle="Parton flavour (%s) [GeV]"%j4_name,
                                 plotopts = channelLabel))
        # DiJet Pt Plot #
        plots.append(Plot.make1D("%s_%s_DiJetPT_%s_%s"%(channel,suffix,j3_name,j4_name), 
                                 (j3.p4+j4.p4).Pt(), 
                                 sel, 
                                 EquidistantBinning(80,0.,400.),
                                 title="Transverse momentum of the dijet_[%s_%s] (channel %s)"%(j3_name,j4_name,channel), 
                                 xTitle= "DiJetPT [%s_%s]"%(j3_name, j4_name),
                                 plotopts = channelLabel))
        
        # DeltaPhi plot #
        plots.append(Plot.make1D("%s_%s_DeltaPhi_%s_%s"%(channel,suffix,j1_name,j4_name), 
                                 op.abs(op.deltaPhi(j1.p4,j4.p4)), 
                                 sel, 
                                 EquidistantBinning(20,0, 3.2),
                                 title="Azimutal angle difference of the dijet_[%s_%s] (channel %s)"%(j1_name,j4_name,channel),
                                 xTitle= "DeltaPhi [%s_%s]"%(j1_name, j4_name),
                                 plotopts = channelLabel))
        plots.append(Plot.make1D("%s_%s_DeltaPhi_%s_%s"%(channel,suffix,j3_name,j4_name), 
                                 op.abs(op.deltaPhi(j3.p4,j4.p4)), 
                                 sel, 
                                 EquidistantBinning(20,0, 3.2),
                                 title="Azimutal angle difference of the dijet_[%s_%s] (channel %s)"%(j3_name,j4_name,channel),
                                 xTitle= "DeltaPhi [%s_%s]"%(j3_name, j4_name),
                                 plotopts = channelLabel))
        plots.append(Plot.make1D("%s_%s_DeltaPhi_%s_%s"%(channel,suffix,'Hbb','Wjj'), 
                                 op.abs(op.deltaPhi((j1.p4+j2.p4),(j3.p4+j4.p4))), 
                                 sel, 
                                 EquidistantBinning(20,0, 3.2),
                                 title="Azimutal angle difference_[%s_%s] (channel %s)"%('Hbb','Wjj',channel),
                                 xTitle= "DeltaPhi [%s_%s]"%('Hbb', 'Wjj'),
                                 plotopts = channelLabel))
    
        # DeltaR plots #
        plots.append(Plot.make1D("%s_%s_DeltaR_%s_%s"%(channel,suffix,j1_name,j4_name), 
                                 op.deltaR(j1.p4,j4.p4), 
                                 sel, 
                                 EquidistantBinning(50,0, 5.0),
                                 title="DeltaR of the dijet_[%s_%s] (channel %s)"%(j1_name,j4_name,channel),
                                 xTitle= "DeltaR [%s_%s]"%(j1_name, j4_name),
                                 plotopts = channelLabel))
        plots.append(Plot.make1D("%s_%s_DeltaR_%s_%s"%(channel,suffix,j3_name,j4_name), 
                                 op.deltaR(j3.p4,j4.p4), 
                                 sel, 
                                 EquidistantBinning(20,0, 3.2),
                                 title="DeltaR of the dijet_[%s_%s] (channel %s)"%(j3_name,j4_name,channel),
                                 xTitle= "DeltaR [%s_%s]"%(j3_name, j4_name),
                                 plotopts = channelLabel))        
        plots.append(Plot.make1D("%s_%s_DeltaR_%s_%s"%(channel,suffix,'Hbb','Wjj'), 
                                 op.deltaR((j1.p4+j2.p4),(j3.p4+j4.p4)), 
                                 sel, 
                                 EquidistantBinning(50,0, 5.0),
                                 title="DeltaR_[%s_%s] (channel %s)"%('Hbb','Wjj',channel),
                                 xTitle= "DeltaR [%s_%s]"%('Hbb', 'Wjj'),
                                 plotopts = channelLabel))

        # invariant mass plot #
        plots.append(Plot.make1D("%s_%s_InvariantMass_%s_%s"%(channel,suffix,j3_name,j4_name),
                                 op.invariant_mass(j3.p4,j4.p4),
                                 sel,
                                 EquidistantBinning(50, 0., 500.), 
                                 title="Invariant Mass of the dijet_[%s_%s] (channel %s)"%(j3_name,j4_name,channel), 
                                 xTitle= "InvariantMass [%s_%s]"%(j3_name, j4_name),
                                 plotopts = channelLabel))

    
    return plots


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
def makeSingleLeptonAk8JetsPlots(sel,j1,j2,j3,suffix,channel,has1fat=False,has2fat=False,has1fat1slim=False,has1fat2slim=False):
    """
    Make fatjet subjet basic plots
    sel         = refine selection 

    fatjet      = fatjet (AK8 jet) 
    suffix      = string identifying the selecton 
    channel     = string identifying the channel of the dilepton (can be "NoChannel")
    if j1 j2 are present in selection:
      if both are fat, keep the bjet in the 1st pos
    """
    #print ("===============>>>Ak8Jet plots__channel:%s__sel:%s"%(channel,suffix))

    plots = []

    channelLabel = SingleLeptonChannelTitleLabel(channel)

    #j1_name = "Ak8Jet" if has1fat else "Ak8bJet"
    #j2_name = "Ak8Jet" if has2fat else "Ak4Jet"
    j1_name = "Ak8BJet"
    j2_name = "Ak8LSJet" if has2fat else "LeadAk4NonBJet"
    j3_name = "SubLeadAk4NonBJet"

    # 1st-jet : FatJet #
    plots.append(Plot.make1D("%s_%s_pT_%s"%(channel,suffix,j1_name),
                             j1.p4.pt(),
                             sel,
                             EquidistantBinning(80,100.,500.),
                             title='Transverse momentum of the %s'%j1_name,
                             xTitle="P_{T} %s[GeV]"%j1_name,
                             plotopts = channelLabel))
    
    plots.append(Plot.make1D("%s_%s_SubJet1_pT_%s"%(channel,suffix,j1_name),
                             j1.subJet1.p4.pt(),
                             sel,
                             EquidistantBinning(40,0.,200.),
                             title='Transverse momentum of the subJet1 of %s'%j1_name,
                             xTitle="subJet1_P_{T} %s[GeV]"%j1_name,
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_SubJet2_pT_%s"%(channel,suffix,j1_name),
                             j1.subJet2.p4.pt(),
                             sel,
                             EquidistantBinning(40,0.,200.),
                             title='Transverse momentum of the subJet2 of %s'%j1_name,
                             xTitle="subJet2_P_{T} %s[GeV]"%j1_name,
                             plotopts = channelLabel))
    
    plots.append(Plot.make1D("%s_%s_eta_%s"%(channel,suffix,j1_name),
                             j1.p4.eta(),
                             sel,
                             EquidistantBinning(22,-3.,3.),
                             title='Pseudorapidity of the %s'%j1_name,
                             xTitle="#eta %s"%j1_name,
                             plotopts = channelLabel))

    plots.append(Plot.make1D("%s_%s_SubJet1_eta_%s"%(channel,suffix,j1_name),
                             j1.subJet1.p4.eta(),
                             sel,
                             EquidistantBinning(22,-3.,3.),
                             title='Pseudorapidity of the subJet1 of %s'%j1_name,
                             xTitle="subJet1_#eta %s[GeV]"%j1_name,
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_SubJet2_eta_%s"%(channel,suffix,j1_name),
                             j1.subJet2.p4.eta(),
                             sel,
                             EquidistantBinning(22,-3.,3.),
                             title='Pseudorapidity of the subJet2 of %s'%j1_name,
                             xTitle="subJet2_#eta %s[GeV]"%j1_name,
                             plotopts = channelLabel))

    plots.append(Plot.make1D("%s_%s_phi_%s"%(channel,suffix,j1_name),
                             j1.p4.phi(),
                             sel,
                             EquidistantBinning(20,-3.2,3.2),
                             title='Azimutal angle of the %s'%j1_name,
                             xTitle="#phi %s"%j1_name,
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_mass_%s"%(channel,suffix,j1_name),
                             j1.mass,
                             sel,
                             EquidistantBinning(40,0.,200.),
                             title='mass of the %s'%j1_name,
                             xTitle="M %s"%j1_name,
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_massSD_%s"%(channel,suffix,j1_name),
                             j1.msoftdrop,
                             sel,
                             EquidistantBinning(40,0.,200.),
                             title='Soft Drop mass of the %s'%j1_name,
                             xTitle="M_{Soft Drop} %s[GeV]"%j1_name,
                             plotopts = channelLabel))
    '''
    plots.append(Plot.make1D("%s_%s_btagDDBvL_%s"%(channel,suffix,j1_name),
                             j1.btagDDBvL,
                             sel,
                             EquidistantBinning(10,0.,1.),
                             title='btagDDBvL of %s'%j1_name,
                             xTitle="btagDDBvL %s[GeV]"%j1_name,
                             plotopts = channelLabel))        
    plots.append(Plot.make1D("%s_%s_btagDDBvL_noMD_%s"%(channel,suffix,j1_name),
                             j1.btagDDBvL_noMD,
                             sel,
                             EquidistantBinning(10,0.,1.),
                             title='btagDDBvL_noMD of %s'%j1_name,
                             xTitle="btagDDBvL_noMD %s[GeV]"%j1_name,
                             plotopts = channelLabel))        
    plots.append(Plot.make1D("%s_%s_btagDDCvB_%s"%(channel,suffix,j1_name),
                             j1.btagDDCvB,
                             sel,
                             EquidistantBinning(10,0.,1.),
                             title='btagDDCvB of %s'%j1_name,
                             xTitle="btagDDCvB %s[GeV]"%j1_name,
                             plotopts = channelLabel))        
    plots.append(Plot.make1D("%s_%s_btagDDCvB_noMD_%s"%(channel,suffix,j1_name),
                             j1.btagDDCvB_noMD,
                             sel,
                             EquidistantBinning(10,0.,1.),
                             title='btagDDCvB_noMD of %s'%j1_name,
                             xTitle="btagDDCvB_noMD %s[GeV]"%j1_name,
                             plotopts = channelLabel))        
    plots.append(Plot.make1D("%s_%s_btagDDCvL_%s"%(channel,suffix,j1_name),
                             j1.btagDDCvL,
                             sel,
                             EquidistantBinning(10,0.,1.),
                             title='btagDDCvL of %s'%j1_name,
                             xTitle="btagDDCvL %s[GeV]"%j1_name,
                             plotopts = channelLabel))        
    plots.append(Plot.make1D("%s_%s_btagDDCvL_noMD_%s"%(channel,suffix,j1_name),
                             j1.btagDDCvL_noMD,
                             sel,
                             EquidistantBinning(10,0.,1.),
                             title='btagDDCvL_noMD of %s'%j1_name,
                             xTitle="btagDDCvL_noMD %s[GeV]"%j1_name,
                             plotopts = channelLabel))        
    plots.append(Plot.make1D("%s_%s_btagDeepB_%s"%(channel,suffix,j1_name),
                             j1.btagDeepB,
                             sel,
                             EquidistantBinning(10,0.,1.),
                             title='btagDeepB of %s'%j1_name,
                             xTitle="btagDeepB %s[GeV]"%j1_name,
                             plotopts = channelLabel))        

    '''
    if has2fat or has1fat1slim or has1fat2slim:
        # 2nd-jet : Fat or Slim Jet #
        plots.append(Plot.make1D("%s_%s_pT_%s"%(channel,suffix,j2_name),
                                 j2.p4.pt(),
                                 sel,
                                 EquidistantBinning(40,0.,200.),
                                 title='Transverse momentum of the %s'%j2_name,
                                 xTitle="P_{T} %s[GeV]"%j2_name,
                                 plotopts = channelLabel))
        plots.append(Plot.make1D("%s_%s_eta_%s"%(channel,suffix,j2_name),
                                 j2.p4.eta(),
                                 sel,
                                 EquidistantBinning(22,-3.,3.),
                                 title='Pseudorapidity of the %s'%j2_name,
                                 xTitle="#eta %s"%j2_name,
                                 plotopts = channelLabel))
        plots.append(Plot.make1D("%s_%s_phi_%s"%(channel,suffix,j2_name),
                                 j2.p4.phi(),
                                 sel,
                                 EquidistantBinning(20,-3.2,3.2),
                                 title='Azimutal angle of the %s'%j2_name,
                                 xTitle="#phi %s"%j2_name,
                                 plotopts = channelLabel))
        plots.append(Plot.make1D("%s_%s_mass_%s"%(channel,suffix,j2_name),
                                 j2.mass,
                                 sel,
                                 EquidistantBinning(40,0.,200.),
                                 title='mass of the %s'%j2_name,
                                 xTitle="M %s"%j2_name,
                                 plotopts = channelLabel))
        if has2fat:
            plots.append(Plot.make1D("%s_%s_massSD_%s"%(channel,suffix,j2_name),
                                     j2.msoftdrop,
                                     sel,
                                     EquidistantBinning(40,0.,200.),
                                     title='Soft Drop mass of the %s'%j2_name,
                                     xTitle="M_{Soft Drop} %s[GeV]"%j2_name,
                                     plotopts = channelLabel))        

        # Di-Jet Plots #
        # pT Sum #
        plots.append(Plot.make1D("%s_%s_DijetPT_%s_%s"%(channel,suffix,j1_name,j2_name), 
                                 (j1.p4+j2.p4).Pt(), 
                                 sel, 
                                 EquidistantBinning(80,0.,400.),
                                 title="Transverse momentum of the dijet_%s_%s (channel %s)"%(j1_name,j2_name,channel), 
                                 xTitle= "DiJetPT [%s_%s] [GeV]"%(j1_name,j2_name),
                                 plotopts = channelLabel))

        # DeltaR #
        plots.append(Plot.make1D("%s_%s_DeltaR_%s_%s"%(channel,suffix,j1_name,j2_name), 
                                 op.deltaR(j1.p4,j2.p4), 
                                 sel, 
                                 EquidistantBinning(50, 0., 5.), 
                                 title="Dijet DeltaR_%s_%s (channel %s)"%(j1_name,j2_name,channel), 
                                 xTitle= "DeltaR [%s_%s]"%(j1_name,j2_name),
                                 plotopts = channelLabel))
        # DeltaPhi #
        plots.append(Plot.make1D("%s_%s_DeltaPhi_%s_%s"%(channel,suffix,j1_name,j2_name), 
                                 op.abs(op.deltaPhi(j1.p4,j2.p4)), 
                                 sel, 
                                 EquidistantBinning(20, 0., 3.2), 
                                 title="Dijet DeltaPhi_%s_%s (channel %s)"%(j1_name,j2_name,channel), 
                                 xTitle= "DeltaPhi [%s_%s]"%(j1_name,j2_name),
                                 plotopts = channelLabel))

        # Invariant Mass #
        plots.append(Plot.make1D("%s_%s_InvariantMass_%s_%s"%(channel,suffix,j1_name,j2_name),
                                 op.invariant_mass(j1.p4,j2.p4),
                                 sel,
                                 EquidistantBinning(50, 0., 500.), 
                                 title="Dijet invariant mass_%s_%s (channel %s)"%(j1_name,j2_name,channel), 
                                 xTitle= "InvariantMass [%s_%s] [GeV]"%(j1_name,j2_name),
                                 plotopts = channelLabel))

        if has1fat2slim:
            # Di-Jet Plots #
            # pT Sum #
            plots.append(Plot.make1D("%s_%s_DijetPT_%s_%s"%(channel,suffix,j2_name,j3_name), 
                                     (j2.p4+j3.p4).Pt(), 
                                     sel, 
                                     EquidistantBinning(80,0.,400.),
                                     title="Transverse momentum of the dijet_%s_%s (channel %s)"%(j2_name,j3_name,channel), 
                                     xTitle= "DiJetPT [%s_%s] [GeV]"%(j2_name,j3_name),
                                     plotopts = channelLabel))
            
            # DeltaR #
            plots.append(Plot.make1D("%s_%s_DeltaR_%s_%s"%(channel,suffix,j2_name,j3_name), 
                                     op.deltaR(j2.p4,j3.p4), 
                                     sel, 
                                     EquidistantBinning(50, 0., 5.), 
                                     title="Dijet DeltaR_%s_%s (channel %s)"%(j2_name,j3_name,channel), 
                                     xTitle= "DeltaR [%s_%s]"%(j2_name,j3_name),
                                     plotopts = channelLabel))
            plots.append(Plot.make1D("%s_%s_DeltaR_%s_%s"%(channel,suffix,'Hbb','Wjj'), 
                                     op.deltaR(j1.p4,(j2.p4+j3.p4)), 
                                     sel, 
                                     EquidistantBinning(50, 0., 5.), 
                                     title="DeltaR_%s_%s (channel %s)"%('Hbb','Wjj',channel), 
                                     xTitle= "DeltaR [%s_%s]"%('Hbb','Wjj'),
                                     plotopts = channelLabel))
            # DeltaPhi #
            plots.append(Plot.make1D("%s_%s_DeltaPhi_%s_%s"%(channel,suffix,j2_name,j3_name), 
                                     op.abs(op.deltaPhi(j2.p4,j3.p4)), 
                                     sel, 
                                     EquidistantBinning(20, 0., 3.2), 
                                     title="Dijet DeltaPhi_%s_%s (channel %s)"%(j2_name,j3_name,channel), 
                                     xTitle= "DeltaPhi [%s_%s]"%(j2_name,j3_name),
                                     plotopts = channelLabel))
            plots.append(Plot.make1D("%s_%s_DeltaPhi_%s_%s"%(channel,suffix,'Hbb','Wjj'), 
                                     op.abs(op.deltaPhi(j1.p4,(j2.p4+j3.p4))), 
                                     sel, 
                                     EquidistantBinning(20, 0., 3.2), 
                                     title="DeltaPhi_%s_%s (channel %s)"%('Hbb','Wjj',channel), 
                                     xTitle= "DeltaPhi [%s_%s]"%('Hbb','Wjj'),
                                     plotopts = channelLabel))

            # Invariant Mass #
            plots.append(Plot.make1D("%s_%s_InvariantMass_%s_%s"%(channel,suffix,j2_name,j3_name),
                                     op.invariant_mass(j2.p4,j3.p4),
                                     sel,
                                     EquidistantBinning(50, 0., 500.), 
                                     title="Dijet invariant mass_%s_%s (channel %s)"%(j2_name,j3_name,channel), 
                                     xTitle= "InvariantMass [%s_%s] [GeV]"%(j2_name,j3_name),
                                     plotopts = channelLabel))
            

    return plots


def makeDoubleLeptonAk8JetsPlots(sel, fatjet, suffix, channel):
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
    plots.append(Plot.make1D("%s_%s_fatjet_btagDeepB"%(channel,suffix),
                             fatjet.btagDeepB,
                             sel,
                             EquidistantBinning(40,0.,1.),
                             title='btagDeepB score of the fatjet',
                             xTitle="btagDeepB",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_fatjet_btagHbb"%(channel,suffix),
                             fatjet.btagHbb,
                             sel,
                             EquidistantBinning(40,0.,1.),
                             title='btagHbb score of the fatjet',
                             xTitle="btagHbb",
                             plotopts = channelLabel))

    return plots

#########################  High-level quantities ################################
def makeHighLevelPlotsResolved(sel,met,lep,j1,j2,j3,j4,channel,suffix,nJet,nbJet,HLL):

    #print ("===============>>>HighLevel Plots for Resolved__channel:%s__sel:%s"%(channel,suffix))
    plots = []

    channelLabel = SingleLeptonChannelTitleLabel(channel)

    if nbJet == 0:
        j1_name = "leadAk4"
        j2_name = "subLeadAk4"
        j3_name = "thirdAk4"
        j4_name = "fourthAk4"
    else:
        j1_name = "leadAk4B"
        j2_name = "leadAk4"if nbJet == 1 else "subLeadAk4B"
        j3_name = "subLeadAk4" if nbJet == 1 else "leadAk4"
        j4_name = "thirdAk4" if nbJet == 1 else "subLeadAk4"

    # mT for Single lepton #
    plots.append(Plot.make1D("%s_%s_highlevelvariable_MT"%(channel,suffix),
                             HLL.MT(lep,met),
                             sel,
                             EquidistantBinning(50,0.,1000.),
                             title='Transverse mass of lepton and MET (%s channel)'%channel,
                             xTitle="M_{T}(lep,MET) [GeV]",
                             plotopts = channelLabel))
    # lepton-MET plots #
    plots.append(Plot.make1D("%s_%s_highlevelvariable_SinglepMET_Pt"%(channel,suffix),
                             op.abs(HLL.SinglepMet_Pt(lep, met)),
                             sel,
                             EquidistantBinning(20,0.,3.2),
                             title='Reultant pT of lepton and MET (%s channel)'%channel,
                             xTitle="|Pt (lep,MET)|",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_highlevelvariable_SinglepMETdeltaPhi"%(channel,suffix),
                             op.abs(HLL.SinglepMet_dPhi(lep, met)),
                             sel,
                             EquidistantBinning(20,0.,3.2),
                             title='Azimutal angle between lepton and MET (%s channel)'%channel,
                             xTitle="|#Delta \phi (lep,MET)|",
                             plotopts = channelLabel))
    
    # jet-MET plots
    plots.append(Plot.make1D("%s_%s_highlevelvariable_%s_METdeltaPhi"%(channel,suffix,j1_name),
                             op.abs(HLL.SinglepMet_dPhi(j1, met)),
                             sel,
                             EquidistantBinning(20,0.,3.2),
                             title='Azimutal angle between %s and MET (%s channel)'%(j1_name,channel),
                             xTitle="|#Delta \phi (%s,MET)|"%j1_name,
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_highlevelvariable_%s_METdeltaPhi"%(channel,suffix,j2_name),
                             op.abs(HLL.SinglepMet_dPhi(j2, met)),
                             sel,
                             EquidistantBinning(20,0.,3.2),
                             title='Azimutal angle between %s and MET (%s channel)'%(j2_name,channel),
                             xTitle="|#Delta \phi (%s,MET)|"%j2_name,
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_highlevelvariable_%s_METdeltaPhi"%(channel,suffix,j3_name),
                             op.abs(HLL.SinglepMet_dPhi(j3, met)),
                             sel,
                             EquidistantBinning(20,0.,3.2),
                             title='Azimutal angle between %s and MET (%s channel)'%(j3_name,channel),
                             xTitle="|#Delta \phi (%s,MET)|"%j3_name,
                             plotopts = channelLabel))
    

    # lepton_Jet DeltaR plots #
    plots.append(Plot.make1D("%s_%s_highlevelvariable_DeltaR_lep_%s"%(channel,suffix,j1_name),
                             op.deltaR(lep.p4,j1.p4),
                             sel,
                             EquidistantBinning(50,0.,5.),
                             title="DeltaR_%s_%s (channel %s)"%(channel,j1_name,channel),
                             xTitle="DeltaR [%s_%s]"%(channel,j1_name),
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_highlevelvariable_DeltaR_lep_%s"%(channel,suffix,j2_name),
                             op.deltaR(lep.p4,j2.p4),
                             sel,
                             EquidistantBinning(50,0.,5.),
                             title="DeltaR_%s_%s (channel %s)"%(channel,j2_name,channel),
                             xTitle="DeltaR [%s_%s]"%(channel,j2_name),
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_highlevelvariable_DeltaR_lep_%s"%(channel,suffix,j3_name),
                             op.deltaR(lep.p4,j3.p4),
                             sel,
                             EquidistantBinning(50,0.,5.),
                             title="DeltaR_%s_%s (channel %s)"%(channel,j3_name,channel),
                             xTitle="DeltaR [%s_%s]"%(channel,j3_name),
                             plotopts = channelLabel))
    
    # min-DR
    plots.append(Plot.make1D("%s_%s_highlevelvariable_minDeltaR_lep_jet"%(channel,suffix),
                             HLL.MinDR_lep3j(lep,j1,j2,j3),
                             sel,
                             EquidistantBinning(50,0.,5.),
                             title="minDeltaR_lep_jet (channel %s)"%channel,
                             xTitle="minDeltaR [%s_Jet]"%channel,
                             plotopts = channelLabel))
    
    # lepton_Jet DeltaPhi plots #
    plots.append(Plot.make1D("%s_%s_highlevelvariable_DeltaPhi_lep_%s"%(channel,suffix,j1_name),
                             op.abs(op.deltaPhi(lep.p4,j1.p4)),
                             sel,
                             EquidistantBinning(20,0.,3.2),
                             title="DeltaPhi_%s_%s (channel %s)"%(channel,j1_name,channel),
                             xTitle="DeltaPhi [%s_%s]"%(channel,j1_name),
                             plotopts=channelLabel))
    plots.append(Plot.make1D("%s_%s_highlevelvariable_DeltaPhi_lep_%s"%(channel,suffix,j2_name),
                             op.abs(op.deltaPhi(lep.p4,j2.p4)),
                             sel,
                             EquidistantBinning(20,0.,3.2),
                             title="DeltaPhi_%s_%s (channel %s)"%(channel,j2_name,channel),
                             xTitle="DeltaPhi [%s_%s]"%(channel,j2_name),
                             plotopts=channelLabel))
    plots.append(Plot.make1D("%s_%s_highlevelvariable_DeltaPhi_lep_%s"%(channel,suffix,j3_name),
                             op.abs(op.deltaPhi(lep.p4,j3.p4)),
                             sel,
                             EquidistantBinning(20,0.,3.2),
                             title="DeltaPhi_%s_%s (channel %s)"%(channel,j3_name,channel),
                             xTitle="DeltaPhi [%s_%s]"%(channel,j3_name),
                             plotopts=channelLabel))
    plots.append(Plot.make1D("%s_%s_highlevelvariable_MT_lep_%s"%(channel,suffix,j3_name),
                             HLL.MT_W1W2_lj(lep,j3,met),
                             sel,
                             EquidistantBinning(50,0.,1000.),
                             title='Transverse mass of lepton+jet and MET (%s channel)'%channel,
                             xTitle="M_{T}(lepj,MET) [GeV]",
                             plotopts = channelLabel))
    # Scalar magnitude sum #
    plots.append(Plot.make1D("%s_%s_highlevelvariable_HT2_l3jmet"%(channel,suffix),
                             HLL.HT2_l3jmet(lep,j1,j2,j3,met),
                             sel,
                             EquidistantBinning(60,0.,600.),
                             title='Di-Higgs magnitude_l3jmet (%s channel)'%channel,
                             xTitle="H_{T2} [GeV]",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_highlevelvariable_HT2R_l3jmet"%(channel,suffix),
                             HLL.HT2R_l3jmet(lep,j1,j2,j3,met),
                             sel,
                             EquidistantBinning(40,0.,1.),
                             title='Di-Higgs magnitude_Ratio_l3jmet (%s channel)'%channel,
                             xTitle="H_{T2}_Ratio [GeV]",
                             plotopts = channelLabel))
    if nJet == 4:
        plots.append(Plot.make1D("%s_%s_highlevelvariable_%s_METdeltaPhi"%(channel,suffix,j4_name),
                                 op.abs(HLL.SinglepMet_dPhi(j4, met)),
                                 sel,
                                 EquidistantBinning(20,0.,3.2),
                                 title='Azimutal angle between %s and MET (%s channel)'%(j4_name,channel),
                                 xTitle="|#Delta \phi (%s,MET)|"%j4_name,
                                 plotopts = channelLabel))

        plots.append(Plot.make1D("%s_%s_highlevelvariable_DeltaR_lep_%s"%(channel,suffix,j4_name),
                                 op.deltaR(lep.p4,j4.p4),
                                 sel,
                                 EquidistantBinning(50,0.,5.),
                                 title="DeltaR_%s_%s (channel %s)"%(channel,j4_name,channel),
                                 xTitle="DeltaR [%s_%s]"%(channel,j4_name),
                                 plotopts=channelLabel))
        plots.append(Plot.make1D("%s_%s_highlevelvariable_DeltaPhi_lep_%s"%(channel,suffix,j4_name),
                                 op.abs(op.deltaPhi(lep.p4,j4.p4)),
                                 sel,
                                 EquidistantBinning(20,0.,3.2),
                                 title="DeltaPhi_%s_%s (channel %s)"%(channel,j4_name,channel),
                                 xTitle="DeltaPhi [%s_%s]"%(channel,j4_name),
                                 plotopts=channelLabel))
        plots.append(Plot.make1D("%s_%s_highlevelvariable_HT2_l4jmet"%(channel,suffix),
                                 HLL.HT2_l4jmet(lep,j1,j2,j3,j4,met),
                                 sel,
                                 EquidistantBinning(60,0.,600.),
                                 title='Di-Higgs magnitude_l4jmet (%s channel)'%channel,
                                 xTitle="H_{T2} [GeV]",
                                 plotopts = channelLabel))
        plots.append(Plot.make1D("%s_%s_highlevelvariable_HT2R_l4jmet"%(channel,suffix),
                                 HLL.HT2R_l4jmet(lep,j1,j2,j3,j4,met),
                                 sel,
                                 EquidistantBinning(40,0.,1.),
                                 title='Di-Higgs magnitude_Ratio_l4jmet (%s channel)'%channel,
                                 xTitle="H_{T2}_Ratio [GeV]",
                                 plotopts = channelLabel))
        # Transverse mass plots #
        if nbJet == 2:
            plots.append(Plot.make1D("%s_%s_highlevelvariable_MT_lep_%s_%s"%(channel,suffix,j3_name,j4_name),
                                     HLL.MT_W1W2_ljj(lep,j3,j4,met),
                                     sel,
                                     EquidistantBinning(50,0.,1000.),
                                     title='Transverse mass of lepton+dijet and MET (%s channel)'%channel,
                                     xTitle="M_{T}(lepjj,MET) [GeV]",
                                     plotopts = channelLabel))

    return plots 

# Make HighLevel plots for semi boosted and boosted categories
def makeHighLevelPlotsBoosted(sel,met,lep,j1,j2,channel,suffix,HLL,bothAreFat=False):
    #print ("===============>>>HighLevel Plots for Boosted__channel:%s__sel:%s"%(channel,suffix))
    if bothAreFat:
        j1_name   = "Ak8BJet"
        j2_name   = "Ak8Jet"
        
    else:
        j1_name   = "Ak8BJet"
        j2_name   = "Ak4NonbJet"
    
    plots = []
    channelLabel = SingleLeptonChannelTitleLabel(channel)

    # lepton_Jet DeltaR plots #
    plots.append(Plot.make1D("%s_%s_highlevelvariable_DeltaR_%s_%s"%(channel,suffix,channel,j1_name),
                             op.deltaR(lep.p4,j1.p4),
                             sel,
                             EquidistantBinning(50,0.,5.),
                             title="DeltaR_(%s_%s)"%(channel, j1_name),
                             xTitle="DeltaR [%s_%s]"%(channel, j1_name),
                             plotopts=channelLabel))
    plots.append(Plot.make1D("%s_%s_highlevelvariable_DeltaR_%s_%s"%(channel,suffix,channel,j2_name),
                             op.deltaR(lep.p4,j2.p4),
                             sel,
                             EquidistantBinning(50,0.,5.),
                             title="DeltaR_(%s_%s)"%(channel,j1_name),
                             xTitle="DeltaR [%s_%s]"%(channel, j2_name),
                             plotopts=channelLabel))

    # lepton_Jet DeltaPhi plots #
    plots.append(Plot.make1D("%s_%s_highlevelvariable_DeltaPhi_%s_%s"%(channel,suffix,channel,j1_name),
                             op.abs(op.deltaPhi(lep.p4,j1.p4)),
                             sel,
                             EquidistantBinning(20,0.,3.2),
                             title="DeltaPhi_(%s_%s)"%(channel, j1_name),
                             xTitle="DeltaPhi [%s_%s]"%(channel, j1_name),
                             plotopts=channelLabel))
    plots.append(Plot.make1D("%s_%s_highlevelvariable_DeltaPhi_%s_%s"%(channel,suffix,channel,j2_name),
                             op.abs(op.deltaPhi(lep.p4,j2.p4)),
                             sel,
                             EquidistantBinning(20,0.,3.2),
                             title="DeltaPhi(%s %s)"%(channel, j2_name),
                             xTitle="DeltaPhi [%s_%s]"%(channel, j2_name),
                             plotopts=channelLabel))
    
    # lepton-MET plots #
    plots.append(Plot.make1D("%s_%s_highlevelvariable_SinglepMET_Pt"%(channel,suffix),
                             op.abs(HLL.SinglepMet_Pt(lep, met)),
                             sel,
                             EquidistantBinning(20,0.,3.2),
                             title='Reultant pT of lepton and MET (%s channel)'%channel,
                             xTitle="|Pt (lep,MET)|",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_highlevelvariable_SinglepMETdeltaPhi"%(channel,suffix),
                             op.abs(HLL.SinglepMet_dPhi(lep, met)),
                             sel,
                             EquidistantBinning(20,0.,3.2),
                             title='Azimutal angle between lepton and MET (%s channel)'%channel,
                             xTitle="|#Delta \phi (lep,MET)|",
                             plotopts = channelLabel))
    # jet-MET plots
    plots.append(Plot.make1D("%s_%s_highlevelvariable_%s_METdeltaPhi"%(channel,suffix,j1_name),
                             op.abs(HLL.SinglepMet_dPhi(j1, met)),
                             sel,
                             EquidistantBinning(20,0.,3.2),
                             title='Azimutal angle between %s and MET (%s channel)'%(j1_name,channel),
                             xTitle="|#Delta \phi (%s,MET)|"%j1_name,
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_highlevelvariable_%s_METdeltaPhi"%(channel,suffix,j2_name),
                             op.abs(HLL.SinglepMet_dPhi(j2, met)),
                             sel,
                             EquidistantBinning(20,0.,3.2),
                             title='Azimutal angle between %s and MET (%s channel)'%(j2_name,channel),
                             xTitle="|#Delta \phi (%s,MET)|"%j2_name,
                             plotopts = channelLabel))


    # Transverse mass plots #
    # mT for Single lepton #
    plots.append(Plot.make1D("%s_%s_highlevelvariable_MT"%(channel,suffix),
                             HLL.MT(lep,met),
                             sel,
                             EquidistantBinning(50,0.,1000.),
                             title='Transverse mass of lepton and MET (%s channel)'%channel,
                             xTitle="M_{T}(lep,MET) [GeV]",
                             plotopts = channelLabel))
    
    return plots

def makeBtagDDplots(sel,jet,channel,suffix):

    plots = []
    channelLabel = SingleLeptonChannelTitleLabel(channel)
    jname = 'ak8Jet'
    
    plots.append(Plot.make1D("%s_%s_btagDDBvL_%s"%(channel,suffix,jname),
                             jet.btagDDBvL,
                             sel,
                             EquidistantBinning(20,0.,1.),
                             title='btagDDBvL of %s'%jname,
                             xTitle="btagDDBvL %s[GeV]"%jname,
                             plotopts = channelLabel))        
    plots.append(Plot.make1D("%s_%s_btagDDBvL_noMD_%s"%(channel,suffix,jname),
                             jet.btagDDBvL_noMD,
                             sel,
                             EquidistantBinning(20,0.,1.),
                             title='btagDDBvL_noMD of %s'%jname,
                             xTitle="btagDDBvL_noMD %s[GeV]"%jname,
                             plotopts = channelLabel))        
    plots.append(Plot.make1D("%s_%s_btagDDCvB_%s"%(channel,suffix,jname),
                             jet.btagDDCvB,
                             sel,
                             EquidistantBinning(20,0.,1.),
                             title='btagDDCvB of %s'%jname,
                             xTitle="btagDDCvB %s[GeV]"%jname,
                             plotopts = channelLabel))        
    plots.append(Plot.make1D("%s_%s_btagDDCvB_noMD_%s"%(channel,suffix,jname),
                             jet.btagDDCvB_noMD,
                             sel,
                             EquidistantBinning(20,0.,1.),
                             title='btagDDCvB_noMD of %s'%jname,
                             xTitle="btagDDCvB_noMD %s[GeV]"%jname,
                             plotopts = channelLabel))        
    plots.append(Plot.make1D("%s_%s_btagDDCvL_%s"%(channel,suffix,jname),
                             jet.btagDDCvL,
                             sel,
                             EquidistantBinning(20,0.,1.),
                             title='btagDDCvL of %s'%jname,
                             xTitle="btagDDCvL %s[GeV]"%jname,
                             plotopts = channelLabel))        
    plots.append(Plot.make1D("%s_%s_btagDDCvL_noMD_%s"%(channel,suffix,jname),
                             jet.btagDDCvL_noMD,
                             sel,
                             EquidistantBinning(20,0.,1.),
                             title='btagDDCvL_noMD of %s'%jname,
                             xTitle="btagDDCvL_noMD %s[GeV]"%jname,
                             plotopts = channelLabel))        
    plots.append(Plot.make1D("%s_%s_btagDeepB_%s"%(channel,suffix,jname),
                             jet.btagDeepB,
                             sel,
                             EquidistantBinning(20,0.,1.),
                             title='btagDeepB of %s'%jname,
                             xTitle="btagDeepB %s[GeV]"%jname,
                             plotopts = channelLabel))        

    return plots


def makeDoubleLeptonHighLevelQuantities (sel,met,l1,l2,j1,j2,suffix,channel,HLL):
    plots = []

    channelLabel = DoubleLeptonChannelTitleLabel(channel)

    # Bjets corr mass #
#    plots.append(Plot.make2D("%s_%s_highlevelvariable_bbRegcorrMass"%(channel,suffix),
#                             op.invariant_mass(HLL.bJetCorrP4(j1), HLL.bJetCorrP4(j2)),
#                             sel,
#                             EquidistantBinning(60,0.,300.),
#                             xTitle="m_{bb}^{Regcorr}",
#                             plotopts = channelLabel))


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
                             op.abs(HLL.DilepMET_deltaPhi(l1,l2,met)),
                             sel,
                             EquidistantBinning(32,0.,3.2),
                             title='Azimutal angle between dilepton and MET (%s channel)'%channel,
                             xTitle="|#Delta \phi (ll,MET)|",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_highlevelvariable_DilepMETpt"%(channel,suffix),
                             HLL.DilepMET_Pt(l1,l2,met),
                             sel,
                             EquidistantBinning(50,0.,500.),
                             title='Transverse momentum of dilepton and MET (%s channel)'%channel,
                             xTitle="P_{T}(ll,MET) [GeV]",
                             plotopts = channelLabel))

    # Invariant mass plots #
    plots.append(Plot.make1D("%s_%s_highlevelvariable_DilepJetInvariantMass"%(channel,suffix),
                             HLL.M_lljj(l1,l2,j1,j2),
                             sel,
                             EquidistantBinning(50,0.,1000.),
                             title='Dilepton-jets invariant mass (%s channel)'%channel,
                             xTitle="M_{lljj}",
                             plotopts = channelLabel))


    # Transverse mass plots #
        # mT for dilepton #
    plots.append(Plot.make1D("%s_%s_highlevelvariable_MTll"%(channel,suffix),
                             HLL.MT_ll(l1,l2,met),
                             sel,
                             EquidistantBinning(50,0.,1000.),
                             title='Transverse mass of dilepton and MET (%s channel)'%channel,
                             xTitle="M_{T}(ll,MET) [GeV]",
                             plotopts = channelLabel))
        # mT for lljj #
    plots.append(Plot.make1D("%s_%s_highlevelvariable_MTlljj"%(channel,suffix),
                             HLL.MT_lljj(l1,l2,j1,j2,met),
                             sel,
                             EquidistantBinning(50,0.,1000.),
                             title='Transverse mass of dilepton+dijet and MET (%s channel)'%channel,
                             xTitle="M_{T}(lljj,MET) [GeV]",
                             plotopts = channelLabel))
    # M_HH #
    plots.append(Plot.make1D("%s_%s_highlevelvariable_MHH"%(channel,suffix),
                             HLL.M_HH(l1,l2,j1,j2,met),
                             sel,
                             EquidistantBinning(100,0.,1500.),
                             title='HH invariant mass (%s channel)'%channel,
                             xTitle="M_{HH}(lljj,MET) [GeV]",
                             plotopts = channelLabel))


    # Scalar magnitude sum #

    plots.append(Plot.make1D("%s_%s_highlevelvariable_HT2"%(channel,suffix),
                             HLL.HT2(l1,l2,j1,j2,met),
                             sel,
                             EquidistantBinning(60,0.,600.),
                             title='Di-Higgs magnitude (%s channel)'%channel,
                             xTitle="H_{T2} [GeV]",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_highlevelvariable_HT2R"%(channel,suffix),
                             HLL.HT2R(l1,l2,j1,j2,met),
                             sel,
                             EquidistantBinning(50,0.,1.),
                             title='Di-Higgs magnitude ratio (%s channel)'%channel,
                             xTitle="H_{T2}^{R} [GeV]",
                             plotopts = channelLabel))

    return plots 

#########################  Machine Learning  ################################
def makeSingleLeptonMachineLearningPlotsBDT(sel,fakeLepColl,lep,met,jets,bJets,lJets,j1,j2,j3,j4,suffix,channel,model_even,model_odd,event,HLL):
    plots = []
    channelLabel = SingleLeptonChannelTitleLabel(channel)

    #output_odd = model_even(op.c_float(0.),op.c_float(0.),op.c_float(0.),op.c_float(0.),op.c_float(0.),op.c_float(0.),op.c_float(0.),op.c_float(0.),op.c_float(0.),op.c_float(0.),op.c_float(0.),op.c_float(0.),op.c_float(0.),op.c_float(0.),op.c_float(0.),op.c_float(0.),op.c_float(0.),op.c_float(0.),op.c_float(0.),op.c_float(0.),op.c_float(0.))
    #output_even = model_odd(op.c_float(0.),op.c_float(0.),op.c_float(0.),op.c_float(0.),op.c_float(0.),op.c_float(0.),op.c_float(0.),op.c_float(0.),op.c_float(0.),op.c_float(0.),op.c_float(0.),op.c_float(0.),op.c_float(0.),op.c_float(0.),op.c_float(0.),op.c_float(0.),op.c_float(0.),op.c_float(0.),op.c_float(0.),op.c_float(0.),op.c_float(0.))
    inputs = [op.c_float(0.)]*21 
    # TODO : 
    output = op.switch(event%2,model_odd(*inputs),model_even(*inputs))

    '''
    output_odd = model_even(
        mindr_lep1_jet(lep,jets),                             # mindr_lep1_jet
        op.invariant_mass(bJetCorrP4(j1), bJetCorrP4(j2)),    # m_Hbb_regCorr
        HHP4_simple_met(bJetCorrP4(j1)+bJetCorrP4(j2), j3.p4, j4.p4, lep.p4, met.p4).M(),   # mHH_simple_met
        Wlep_met_simple(lep.p4, met.p4).M(),           # mWlep_met_simple
        HWW_met_simple(j3.p4,j4.p4,lep.p4,met.p4).M(), # mWW_simple_met
        Wjj_simple(j3.p4, j4.p4).M(),                  # mWjj_simple
        comp_cosThetaS(bJetCorrP4(j1),bJetCorrP4(j2)), # cosThetaS_Hbb
        comp_cosThetaS(j3.p4, j4.p4),                  # cosThetaS_Wjj_simple
        comp_cosThetaS(Wjj_simple(j3.p4,j4.p4),Wlep_met_simple(lep.p4,met.p4)),                      # cosThetaS_WW_simple_met
        comp_cosThetaS(bJetCorrP4(j1)+bJetCorrP4(j2), HWW_met_simple(j3.p4,j4.p4,lep.p4,met.p4)), # cosThetaS_HH_simple_met
        op.rng_len(jets),         # nJet
        op.rng_len(bJets),        # nBJetMedium
        op.deltaR(j1.p4, lep.p4), # dR_b1lep
        op.deltaR(j2.p4, lep.p4), # dR_b2lep
        op.switch(op.abs(lep.pdgId)==11, lambda_conept_electron(lep), lambda_conept_muon(lep)),          # lep_conePt
        #lepConePt(lep),           # lep_conePt
        #op.c_float(1.0),
        bJetCorrP4(j1).Pt(),      # selJet1_Hbb_pT
        bJetCorrP4(j2).Pt(),      # selJet2_Hbb_pT 
        MET_LD(met,jets,fakeLepColl), # met_LD
        HT(jets),                     # HT
        op.min(mT2(j3,lep,met), mT2(j4,lep,met)),       # mT_top_3particle
        MT(lep,met)                   # mT_W
    )
    output_even = model_odd(
        mindr_lep1_jet(lep,jets),                             # mindr_lep1_jet
        op.invariant_mass(bJetCorrP4(j1), bJetCorrP4(j2)),    # m_Hbb_regCorr
        HHP4_simple_met(bJetCorrP4(j1)+bJetCorrP4(j2), j3.p4, j4.p4, lep.p4, met.p4).M(),   # mHH_simple_met
        Wlep_met_simple(lep.p4, met.p4).M(),           # mWlep_met_simple
        HWW_met_simple(j3.p4,j4.p4,lep.p4,met.p4).M(), # mWW_simple_met
        Wjj_simple(j3.p4, j4.p4).M(),                  # mWjj_simple
        comp_cosThetaS(bJetCorrP4(j1),bJetCorrP4(j2)), # cosThetaS_Hbb
        comp_cosThetaS(j3.p4, j4.p4),                  # cosThetaS_Wjj_simple
        comp_cosThetaS(Wjj_simple(j3.p4,j4.p4),Wlep_met_simple(lep.p4,met.p4)),                      # cosThetaS_WW_simple_met
        comp_cosThetaS(bJetCorrP4(j1)+bJetCorrP4(j2), HWW_met_simple(j3.p4,j4.p4,lep.p4,met.p4)), # cosThetaS_HH_simple_met
        op.rng_len(jets),         # nJet
        op.rng_len(bJets),        # nBJetMedium
        op.deltaR(j1.p4, lep.p4), # dR_b1lep
        op.deltaR(j2.p4, lep.p4), # dR_b2lep
        op.switch(op.abs(lep.pdgId)==11, lambda_conept_electron(lep), lambda_conept_muon(lep)),          # lep_conePt
        #lepConePt(lep),           # lep_conePt
        #op.c_float(1.0),
        bJetCorrP4(j1).Pt(),      # selJet1_Hbb_pT
        bJetCorrP4(j2).Pt(),      # selJet2_Hbb_pT
        MET_LD(met,jets,fakeLepColl), # met_LD 
        HT(jets),                 # HT
        op.min(mT2(j3,lep,met), mT2(j4,lep,met)),       # mT_top_3particle
        MT(lep,met)               # mT_W
    )
    '''
    for i in range(10):
        plots.append(Plot.make1D("%s_%s_BDTOutput_%d"%(channel,suffix,i),
                                 output[i],
                                 sel,
                                 EquidistantBinning(50,0.,1.),
                                 xTitle = 'BDT output %d'%i,
                                 plotopts = channelLabel))
    
    return plots

'''
def makeSingleLeptonMachineLearningPlotsBDTmissReco(sel,lep,met,jets,bJets,lJets,j1,j2,j3,j4,suffix,channel,model):
    plots = []
    channelLabel = SingleLeptonChannelTitleLabel(channel)
    output = model(op.c_float(0.),
                   MT(lep, met),
                   op.sort(jets, lambda j : op.deltaR(lep.p4, j.p4))[0],
                   op.deltaR(j1.p4, lep.p4),
                   op.deltaR(j2.p4, lep.p4),
                   op.invariant_mass(bJetCorrP4(j1),bJetCorrP4(j2)),
                   bJetCorrP4(j1).Pt(),
                   bJetCorrP4(j2).Pt(),
                   op.deltaR(lep.p4, j3.p4),
                   op.rng_len(bJets),
                   lep.pt,
                   met.pt,
                   op.rng_sum(jets, lambda j : j.pt))

    for i,node in enumerate(['HH', 'H', 'DY', 'ST', 'ttbar', 'ttVX', 'VVV', 'Rare']):
        plots.append(Plot.make1D("%s_%s_BDTOutput_%s"%(channel,suffix,node),
                                 output[i],
                                 sel,
                                 EquidistantBinning(50,0.,1.),
                                 xTitle = 'BDT output %s'%node,
                                 plotopts = channelLabel))

    return plots
    '''
def makeDoubleLeptonMachineLearningInputPlots(sel,suffix,channel,inputs):
    plots = []
    channelLabel = DoubleLeptonChannelTitleLabel(channel)
    for (varname,vartitle,binning),variable in inputs.items():
        plots.append(Plot.make1D("%s_%s_DNNInput_%s"%(channel,suffix,varname),
                                 variable,
                                 sel,
                                 EquidistantBinning(*binning),
                                 xTitle = '%s'%vartitle,
                                 plotopts = channelLabel))
    return plots
    
def makeDoubleLeptonMachineLearningOutputPlots(selObjNodesDict,output,nodes,channel):
    plots = []

    channelLabel = DoubleLeptonChannelTitleLabel(channel)
    for i,node in enumerate(nodes):
        suffix = selObjNodesDict[node].selName
        sel = selObjNodesDict[node].sel
        plots.append(Plot.make1D("%s_%s_DNNOutput_%s"%(channel,suffix,node),
                                 output[i],
                                 sel,
                                 EquidistantBinning(50,0.,1.),
                                 xTitle = 'DNN output %s'%node,
                                 plotopts = channelLabel))
        plots.append(Plot.make1D("%s_%s_DNNOutputRebin_%s"%(channel,suffix,node),
                                 output[i],
                                 sel,
                                 EquidistantBinning(1,0.,1.),
                                 xTitle = 'DNN output %s'%node,
                                 plotopts = channelLabel))


    return plots    

