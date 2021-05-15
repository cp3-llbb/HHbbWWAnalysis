import math
from bamboo.plots import Plot, EquidistantBinning, SummedPlot, VariableBinning
from bamboo import treefunctions as op
from copy import copy
from bamboo.plots import CutFlowReport
from highlevelLambdas import highlevelLambdas
#from JPA import bJetCorrPT
from JPA import *
#from bamboo.root import loadHeader
#loadHeader("/home/ucl/cp3/gsaha/bamboodev/HHbbWWAnalysis/jpa.h")

########################   Channel title   #############################
def SingleLeptonChannelTitleLabel(channel):
    if (channel == "El"):
        channel = "e^{#pm}"
    elif (channel == "Mu"):
        channel ="#mu^{#pm}"
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
    #plots.append(Plot.make1D("%s_%s_trigger_singleElectron"%(channel,suffix), 
    #                         op.OR(*triggerDict["SingleElectron"]),
    #                         sel, 
    #                         EquidistantBinning(2, 0., 2.), 
    #                         xTitle= "SingleElectron trigger",
    #                         plotopts = channelLabel))
    # Single Muon #
    #plots.append(Plot.make1D("%s_%s_trigger_singleMuon"%(channel,suffix), 
    #                         op.OR(*triggerDict['SingleMuon']), 
    #                         sel, 
    #                         EquidistantBinning(2, 0., 2.), 
    #                         xTitle= "SingleMuon trigger",
    #                         plotopts = channelLabel))
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
    #plots.append(Plot.make2D("%s_%s_lepton_ptVSeta"%(channel,suffix), 
    #                         [lep.pt, lep.eta],
    #                         sel, 
    #                         [EquidistantBinning(60,0.,300.),EquidistantBinning(22, -3., 3.)],
    #                         xTitle= "P_{T} (lepton) [GeV]",
    #                         yTitle= "#eta (lepton)",
    #                         plotopts = channelLabel))
    # Phi plot #
    plots.append(Plot.make1D("%s_%s_lepton_phi"%(channel,suffix), 
                             lep.phi, 
                             sel, 
                             EquidistantBinning(20, -3.2, 3.2), 
                             title="Azimutal angle of the lepton (channel %s)"%channel, 
                             xTitle= "#phi (lepton)",
                             plotopts = channelLabel))

    # GenPartFlav (if isMC) #
    #plots.append(Plot.make1D("%s_%s_lepton_genPartFlav"%(channel,suffix), 
    #                         lep.genPartFlav if is_MC else op.c_int(-1),
    #                         sel, 
    #                         EquidistantBinning(23, -1., 22.), 
    #                         title="Flavour of genParticle (channel %s)"%channel, 
    #                         xTitle= "GenParticle flavour (lepton)",
    #                         plotopts = channelLabel))

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
    #plots.append(Plot.make2D("%s_%s_firstlepton_ptVSeta"%(channel,suffix), 
    #                         [dilepton[0].eta, dilepton[0].pt],
    #                         sel, 
    #                         [EquidistantBinning(22, -3., 3.),EquidistantBinning(60,0.,300.)],
    #                         xTitle= "P_{T} (first lepton) [GeV]",
    #                         yTitle= "#eta (second lepton)",
    #                         plotopts = channelLabel))

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
def makeAk4JetsPlots (sel,jet1,jet2,jet3,jet4,channel,suffix,nJet,nbJet,HLL,is_MC=False):

    #print ("===============>>>Ak4Jet Plots__channel:%s__sel:%s"%(channel,suffix))
    plots = []
    j1 = jet1
    j2 = jet2
    j3 = jet3
    j4 = jet4
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
                             j1.pt,
                             sel,
                             EquidistantBinning(40,0.,200.),
                             title='Transverse momentum of the %s'%j1_name,
                             xTitle="P_{T}(%s) [GeV]"%j1_name,
                             plotopts = channelLabel))

    plots.append(Plot.make1D("%s_%s_btagDeepFlavB_%s"%(channel,suffix,j1_name),
                             j1.btagDeepFlavB,
                             sel,
                             EquidistantBinning(20,0.,1.),
                             title='btagDeepFlavB of the %s'%j1_name,
                             xTitle="btagDeepFlavB_%s"%j1_name,
                             plotopts = channelLabel))

    #plots.append(Plot.make1D("%s_%s_RegCorrPt_%s"%(channel,suffix,j1_name),
    #                         bJetCorrPT(j1),
    #                         sel,
    #                         EquidistantBinning(40,0.,200.),
    #                         title='Corr p_{T} of the %s'%j1_name,
    #                         xTitle="RegCorrP_{T}(%s) [GeV]"%j1_name,
    #                         plotopts = channelLabel))

    plots.append(Plot.make1D("%s_%s_eta_%s"%(channel,suffix,j1_name),
                             j1.eta,
                             sel,
                             EquidistantBinning(22,-3.,3.),
                             title='Pseudorapidity of the %s'%j1_name,
                             xTitle="#eta(%s)"%j1_name,
                             plotopts = channelLabel))
    #plots.append(Plot.make1D("%s_%s_phi_%s"%(channel,suffix,j1_name),
    #                         j1.phi,
    #                         sel,
    #                         EquidistantBinning(20,-3.2,3.2),
    #                         title='Azimutal angle of the %s'%j1_name,
    #                         xTitle="#phi(%s)"%j1_name,
    #                         plotopts = channelLabel))

    #plots.append(Plot.make1D("%s_%s_hadronFlv_%s"%(channel,suffix,j1_name),
    #                         op.abs(j1.hadronFlavour) if is_MC else op.c_int(-1),
    #                         sel,
    #                         EquidistantBinning(23,-1.,22.),
    #                         title='Hadron flavour of the %s'%j1_name,
    #                         xTitle="Hadron flavour (%s) [GeV]"%j1_name,
    #                         plotopts = channelLabel))
    #plots.append(Plot.make1D("%s_%s_partonFlv_%s"%(channel,suffix,j1_name),
    #                         op.abs(j1.partonFlavour) if is_MC else op.c_int(-1),
    #                         sel,
    #                         EquidistantBinning(23,-1.,22.),
    #                         title='Parton flavour of the %s'%j1_name,
    #                         xTitle="Parton flavour (%s) [GeV]"%j1_name,
    #                         plotopts = channelLabel))

    if nJet > 1 :
        plots.append(Plot.make1D("%s_%s_btagDeepFlavB_%s"%(channel,suffix,j2_name),
                                 j2.btagDeepFlavB,
                                 sel,
                                 EquidistantBinning(20,0.,1.),
                                 title='btagDeepFlavB of the %s'%j2_name,
                                 xTitle="btagDeepFlavB_%s"%j2_name,
                                 plotopts = channelLabel))
        #plots.append(Plot.make1D("%s_%s_qgDiscr_%s"%(channel,suffix,j2_name),
        #                         j2.qgl,
        #                         sel,
        #                         EquidistantBinning(20,0.,1.),
        #                         title='qgDiscr of the %s'%j2_name,
        #                         xTitle="qgDiscr_%s"%j2_name,
        #                         plotopts = channelLabel))
        
        #plots.append(Plot.make1D("%s_%s_Hbb_massReg_%s_%s"%(channel,suffix,j1_name,j2_name),
        #                         (HLL.bJetCorrP4(j1) + HLL.bJetCorrP4(j2)).M(),
        #                         sel,
        #                         EquidistantBinning(50, 0., 500.), 
        #                         title="Hbb_massReg_[%s_%s] (channel %s)"%(j1_name,j2_name,channel), 
        #                         xTitle= "Hbb_massReg",
        #                         plotopts = channelLabel))
        
        # subleadjet plots #
        plots.append(Plot.make1D("%s_%s_pT_%s"%(channel,suffix,j2_name),
                                 j2.pt,
                                 sel,
                                 EquidistantBinning(40,0.,200.),
                                 title='Transverse momentum of the %s'%j2_name,
                                 xTitle="P_{T}(%s) [GeV]"%j2_name,
                                 plotopts = channelLabel))

        #plots.append(Plot.make1D("%s_%s_RegCorrPt_%s"%(channel,suffix,j2_name),
        #                         HLL.bJetCorrP4(j2).Pt(),
        #                         sel,
        #                         EquidistantBinning(40,0.,200.),
        #                         title='RegCorr Transverse momentum of the %s'%j2_name,
        #                         xTitle="RegCorr P_{T}(%s) [GeV]"%j2_name,
        #                         plotopts = channelLabel))
        
        plots.append(Plot.make1D("%s_%s_eta_%s"%(channel,suffix,j2_name),
                                 j2.eta,
                                 sel,
                                 EquidistantBinning(22,-3.,3.),
                                 title='Pseudorapidity of the %s'%j2_name,
                                 xTitle="#eta(%s)"%j2_name,
                                 plotopts = channelLabel))
        #plots.append(Plot.make1D("%s_%s_phi_%s"%(channel,suffix,j2_name),
        #                         j2.phi,
        #                         sel,
        #                         EquidistantBinning(20,-3.2,3.2),
        #                         title='Azimutal angle of the %s'%j2_name,
        #                         xTitle="#phi(%s)"%j2_name,
        #                         plotopts = channelLabel))

        #plots.append(Plot.make1D("%s_%s_hadronFlv_%s"%(channel,suffix,j2_name),
        #                         op.abs(j2.hadronFlavour) if is_MC else op.c_int(-1),
        #                         sel,
        #                         EquidistantBinning(23,-1.,22.),
        #                         title='Hadron flavour of the %s'%j2_name,
        #                         xTitle="Hadron flavour (%s) [GeV]"%j2_name,
        #                         plotopts = channelLabel))
        #plots.append(Plot.make1D("%s_%s_partonFlv_%s"%(channel,suffix,j2_name),
        #                         op.abs(j2.partonFlavour) if is_MC else op.c_int(-1),
        #                         sel,
        #                         EquidistantBinning(23,-1.,22.),
        #                         title='Parton flavour of the %s'%j2_name,
        #                         xTitle="Parton flavour (%s) [GeV]"%j2_name,
        #                         plotopts = channelLabel))

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
        #plots.append(Plot.make1D("%s_%s_DeltaPhi_%s_%s"%(channel,suffix,j1_name,j2_name), 
        #                         op.abs(op.deltaPhi(j1.p4,j2.p4)), 
        #                         sel, 
        #                         EquidistantBinning(20,0, 3.2),
        #                         title="Azimutal angle difference of the dijet_[%s_%s] (channel %s)"%(j1_name,j2_name,channel),
        #                         xTitle= "DeltaPhi [%s_%s]"%(j1_name, j2_name),
        #                         plotopts = channelLabel))
        
        # DeltaR plot #
        #plots.append(Plot.make1D("%s_%s_DeltaR_%s_%s"%(channel,suffix,j1_name,j2_name), 
        #                         op.deltaR(j1.p4,j2.p4), 
        #                         sel, 
        #                         EquidistantBinning(50,0, 5.0),
        #                         title="DeltaR of the dijet_[%s_%s] (channel %s)"%(j1_name,j2_name,channel),
        #                         xTitle= "DeltaR [%s_%s]"%(j1_name, j2_name),
        #                         plotopts = channelLabel))

        # invariant mass plot #
        #plots.append(Plot.make1D("%s_%s_InvariantMass_%s_%s"%(channel,suffix,j1_name,j2_name),
        #                         #op.invariant_mass(j1.p4,j2.p4),
        #                         op.invariant_mass(HLL.bJetCorrP4(j1),HLL.bJetCorrP4(j2)),
        #                         sel,
        #                         EquidistantBinning(50, 0., 500.), 
        #                         title="Invariant Mass of the dijet_[%s_%s] (channel %s)"%(j1_name,j2_name,channel), 
        #                         xTitle= "InvariantMass [%s_%s]"%(j1_name, j2_name),
        #                         plotopts = channelLabel))
        
        if nJet > 2 :
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
            plots.append(Plot.make1D("%s_%s_btagDeepFlavB_%s"%(channel,suffix,j3_name),
                                     j3.btagDeepFlavB,
                                     sel,
                                     EquidistantBinning(20,0.,1.),
                                     title='btagDeepFlavB of the %s'%j3_name,
                                     xTitle="btagDeepFlavB_%s"%j3_name,
                                     plotopts = channelLabel))
            
            #plots.append(Plot.make1D("%s_%s_phi_%s"%(channel,suffix,j3_name),
            #                         j3.phi,
            #                         sel,
            #                         EquidistantBinning(20,-3.2,3.2),
            #                         title='Azimutal angle of the %s'%j3_name,
            #                         xTitle="#phi(%s)"%j3_name,
            #                         plotopts = channelLabel))
            #plots.append(Plot.make1D("%s_%s_hadronFlv_%s"%(channel,suffix,j3_name),
            #                         op.abs(j3.hadronFlavour) if is_MC else op.c_int(-1),
            #                         sel,
            #                         EquidistantBinning(23,-1.,22.),
            #                         title='Hadron flavour of the %s'%j3_name,
            #                         xTitle="Hadron flavour (%s) [GeV]"%j3_name,
            #                         plotopts = channelLabel))
            #plots.append(Plot.make1D("%s_%s_partonFlv_%s"%(channel,suffix,j3_name),
            #                         op.abs(j3.partonFlavour) if is_MC else op.c_int(-1),
            #                         sel,
            #                         EquidistantBinning(23,-1.,22.),
            #                         title='Parton flavour of the %s'%j3_name,
            #                         xTitle="Parton flavour (%s) [GeV]"%j3_name,
            #                         plotopts = channelLabel))
            #plots.append(Plot.make1D("%s_%s_DeltaPhi_%s_%s"%(channel,suffix,j1_name,j3_name), 
            #                         op.abs(op.deltaPhi(j1.p4,j3.p4)), 
            #                         sel, 
            #                         EquidistantBinning(20,0, 3.2),
            #                         title="Azimutal angle difference of the dijet_[%s_%s] (channel %s)"%(j1_name,j3_name,channel),
            #                         xTitle= "DeltaPhi [%s_%s]"%(j1_name, j3_name),
            #                         plotopts = channelLabel))
            #plots.append(Plot.make1D("%s_%s_DeltaPhi_%s_%s"%(channel,suffix,j2_name,j3_name), 
            #                         op.abs(op.deltaPhi(j2.p4,j3.p4)), 
            #                         sel, 
            #                         EquidistantBinning(20,0, 3.2),
            #                         title="Azimutal angle difference of the dijet_[%s_%s] (channel %s)"%(j2_name,j3_name,channel),
            #                         xTitle= "DeltaPhi [%s_%s]"%(j2_name, j3_name),
            #                         plotopts = channelLabel))
            #plots.append(Plot.make1D("%s_%s_DeltaR_%s_%s"%(channel,suffix,j1_name,j3_name), 
            #                         op.deltaR(j1.p4,j3.p4), 
            #                         sel, 
            #                         EquidistantBinning(50,0, 5.0),
            #                         title="DeltaR of the dijet_[%s_%s] (channel %s)"%(j1_name,j3_name,channel),
            #                         xTitle= "DeltaR [%s_%s]"%(j1_name, j3_name),
            #                         plotopts = channelLabel))
            #plots.append(Plot.make1D("%s_%s_DeltaR_%s_%s"%(channel,suffix,j2_name,j3_name), 
            #                         op.deltaR(j2.p4,j3.p4), 
            #                         sel, 
            #                         EquidistantBinning(50,0, 5.0),
            #                         title="DeltaR of the dijet_[%s_%s] (channel %s)"%(j2_name,j3_name,channel),
            #                         xTitle= "DeltaR [%s_%s]"%(j2_name, j3_name),
            #                         plotopts = channelLabel))
            '''
            # minDR
            plots.append(Plot.make1D("%s_%s_minDeltaR_DiJet_3Jets"%(channel,suffix), 
                                     HLL.MinDiJetDRLoose(j1,j2,j3), 
                                     sel, 
                                     EquidistantBinning(50,0, 5.0),
                                     title="MinDeltaR between dijet (channel %s)"%channel,
                                     xTitle= "MinDeltaR [Di-Jet]",
                                     plotopts = channelLabel))
        
            '''
    if nJet > 3 :
        # fourthjet plots #
        #plots.append(Plot.make1D("%s_%s_maxWjetbtagCSV"%(channel,suffix), 
        #                         op.max(j3.btagDeepFlavB, j4.btagDeepFlavB), 
        #                         sel, 
        #                         EquidistantBinning(20, 0, 1.0),
        #                         title=" maxWjetbtagCSV (channel %s)"%channel,
        #                         xTitle= "max-Wjet_btagCSV",
        #                         plotopts = channelLabel))

        #plots.append(Plot.make1D("%s_%s_RegCorrPt_%s"%(channel,suffix,j4_name),
        #                         bJetCorrPT(j4),
        #                         #j4.pt*j4.bRegCorr,
        #                         sel,
        #                         EquidistantBinning(40,0.,200.),
        #                         title='wjet2_ptReg of the %s'%j4_name,
        #                         xTitle="RegCorrP_{T}(%s) [GeV]"%j4_name,
        #                         plotopts = channelLabel))

        #plots.append(Plot.make1D("%s_%s_qgDiscr_%s"%(channel,suffix,j4_name),
        #                         j4.qgl,
        #                         sel,
        #                         EquidistantBinning(20,0.,1.),
        #                         title='qgDiscr of the %s'%j4_name,
        #                         xTitle="qgDiscr_%s"%j4_name,
        #                         plotopts = channelLabel))
        
        #plots.append(Plot.make1D("%s_%s_pT_%s"%(channel,suffix,j4_name),
        #                         j4.pt,
        #                         sel,
        #                         EquidistantBinning(40,0.,200.),
        #                         title='Transverse momentum of the %s'%j4_name,
        #                         xTitle="P_{T}(%s) [GeV]"%j4_name,
        #                         plotopts = channelLabel))
        plots.append(Plot.make1D("%s_%s_eta_%s"%(channel,suffix,j4_name),
                                 j4.eta,
                                 sel,
                                 EquidistantBinning(22,-3.,3.),
                                 title='Pseudorapidity of the %s'%j4_name,
                                 xTitle="#eta(%s)"%j4_name,
                                 plotopts = channelLabel))
        plots.append(Plot.make1D("%s_%s_btagDeepFlavB_%s"%(channel,suffix,j4_name),
                                 j4.btagDeepFlavB,
                                 sel,
                                 EquidistantBinning(20,0.,1.),
                                 title='btagDeepFlavB of the %s'%j4_name,
                                 xTitle="btagDeepFlavB_%s"%j4_name,
                                 plotopts = channelLabel))
        #plots.append(Plot.make1D("%s_%s_phi_%s"%(channel,suffix,j4_name),
        #                         j4.phi,
        #                         sel,
        #                         EquidistantBinning(20,-3.2,3.2),
        #                         title='Azimutal angle of the %s'%j4_name,
        #                         xTitle="#phi(%s)"%j4_name,
        #                         plotopts = channelLabel))
        #plots.append(Plot.make1D("%s_%s_hadronFlv_%s"%(channel,suffix,j4_name),
        #                         op.abs(j4.hadronFlavour) if is_MC else op.c_int(-1),
        #                         sel,
        #                         EquidistantBinning(23,-1.,22.),
        #                         title='Hadron flavour of the %s'%j4_name,
        #                         xTitle="Hadron flavour (%s) [GeV]"%j4_name,
        #                         plotopts = channelLabel))
        #plots.append(Plot.make1D("%s_%s_partonFlv_%s"%(channel,suffix,j4_name),
        #                         op.abs(j4.partonFlavour) if is_MC else op.c_int(-1),
        #                         sel,
        #                         EquidistantBinning(23,-1.,22.),
        #                         title='Parton flavour of the %s'%j4_name,
        #                         xTitle="Parton flavour (%s) [GeV]"%j4_name,
        #                         plotopts = channelLabel))
        # DiJet Pt Plot #
        #plots.append(Plot.make1D("%s_%s_DiJetPT_%s_%s"%(channel,suffix,j3_name,j4_name), 
        #                         (j3.p4+j4.p4).Pt(), 
        #                         sel, 
        #                         EquidistantBinning(80,0.,400.),
        #                         title="Transverse momentum of the dijet_[%s_%s] (channel %s)"%(j3_name,j4_name,channel), 
        #                         xTitle= "DiJetPT [%s_%s]"%(j3_name, j4_name),
        #                         plotopts = channelLabel))
        
        # DeltaPhi plot #
        #plots.append(Plot.make1D("%s_%s_DeltaPhi_%s_%s"%(channel,suffix,j1_name,j4_name), 
        #                         op.abs(op.deltaPhi(j1.p4,j4.p4)), 
        #                         sel, 
        #                         EquidistantBinning(20,0, 3.2),
        #                         title="Azimutal angle difference of the dijet_[%s_%s] (channel %s)"%(j1_name,j4_name,channel),
        #                         xTitle= "DeltaPhi [%s_%s]"%(j1_name, j4_name),
        #                         plotopts = channelLabel))
        #plots.append(Plot.make1D("%s_%s_DeltaPhi_%s_%s"%(channel,suffix,j3_name,j4_name), 
        #                         op.abs(op.deltaPhi(j3.p4,j4.p4)), 
        #                         sel, 
        #                         EquidistantBinning(20,0, 3.2),
        #                         title="Azimutal angle difference of the dijet_[%s_%s] (channel %s)"%(j3_name,j4_name,channel),
        #                         xTitle= "DeltaPhi [%s_%s]"%(j3_name, j4_name),
        #                         plotopts = channelLabel))
        #plots.append(Plot.make1D("%s_%s_DeltaPhi_%s_%s"%(channel,suffix,'Hbb','Wjj'), 
        #                         op.abs(op.deltaPhi((j1.p4+j2.p4),(j3.p4+j4.p4))), 
        #                         sel, 
        #                         EquidistantBinning(20,0, 3.2),
        #                         title="Azimutal angle difference_[%s_%s] (channel %s)"%('Hbb','Wjj',channel),
        #                         xTitle= "DeltaPhi [%s_%s]"%('Hbb', 'Wjj'),
        #                         plotopts = channelLabel))
    
        # DeltaR plots #
        #plots.append(Plot.make1D("%s_%s_DeltaR_%s_%s"%(channel,suffix,j1_name,j4_name), 
        #                         op.deltaR(j1.p4,j4.p4), 
        #                         sel, 
        #                         EquidistantBinning(50,0, 5.0),
        #                         title="DeltaR of the dijet_[%s_%s] (channel %s)"%(j1_name,j4_name,channel),
        #                         xTitle= "DeltaR [%s_%s]"%(j1_name, j4_name),
        #                         plotopts = channelLabel))
        #plots.append(Plot.make1D("%s_%s_DeltaR_%s_%s"%(channel,suffix,'Hbb','Wjj'), 
        #                         op.deltaR((j1.p4+j2.p4),(j3.p4+j4.p4)), 
        #                         sel, 
        #                         EquidistantBinning(50,0, 5.0),
        #                         title="DeltaR_[%s_%s] (channel %s)"%('Hbb','Wjj',channel),
        #                         xTitle= "DeltaR [%s_%s]"%('Hbb', 'Wjj'),
        #                         plotopts = channelLabel))

        # invariant mass plot #
        #plots.append(Plot.make1D("%s_%s_InvariantMass_%s_%s"%(channel,suffix,j3_name,j4_name),
        #                         op.invariant_mass(j3.p4,j4.p4),
        #                         sel,
        #                         EquidistantBinning(50, 0., 500.), 
        #                         title="Invariant Mass of the dijet_[%s_%s] (channel %s)"%(j3_name,j4_name,channel), 
        #                         xTitle= "InvariantMass [%s_%s]"%(j3_name, j4_name),
        #                         plotopts = channelLabel))


        #plots.append(Plot.make1D("%s_%s_HadW_mass_%s_%s"%(channel,suffix,j3_name,j4_name),
        #                         HLL.Wjj_simple(j3.p4, j4.p4).M(),
        #                         sel,
        #                         EquidistantBinning(25, 0., 250.), 
        #                         title="Invariant Mass of the dijet_[%s_%s] (channel %s)"%(j3_name,j4_name,channel), 
        #                         xTitle= "HadW_mass",
        #                         plotopts = channelLabel))
        #plots.append(Plot.make1D("%s_%s_maxDeltaR_%s_%s"%(channel,suffix,'HadW','bJet'), 
        #                         op.max(HLL.dR_HadW_bjet(j1.p4, j3.p4, j4.p4),   
        #                                HLL.dR_HadW_bjet(j2.p4, j3.p4, j4.p4)),
        #                         sel, 
        #                         EquidistantBinning(50,0, 5.0),
        #                         title="maxDeltaR_[%s_%s] (channel %s)"%('HadW','bJet',channel),
        #                         xTitle= "maxDeltaR [%s_%s]"%('HadW', 'bJet'),
        #                         plotopts = channelLabel))

        #plots.append(Plot.make1D("%s_%s_DeltaR_%s_%s"%(channel,suffix,j3_name,j4_name), 
        #                         op.deltaR(j3.p4,j4.p4), 
        #                         sel, 
        #                         EquidistantBinning(50,0, 5.0),
        #                         title="DeltaR of the dijet_[%s_%s] (channel %s)"%(j3_name,j4_name,channel),
        #                         xTitle= "DeltaR [%s_%s]"%(j3_name, j4_name),
        #                         plotopts = channelLabel))        

        #plots.append(Plot.make1D("%s_%s_mTop1"%(channel,suffix), 
        #                         (HLL.bJetCorrP4(j1) + HLL.Wjj_simple(j3.p4, j4.p4)).M(), 
        #                         sel, 
        #                         EquidistantBinning(50,0, 500.0),
        #                         title="mTop1 (channel %s)"%channel,
        #                         xTitle= "mTop1",
        #                         plotopts = channelLabel))        
        
        #plots.append(Plot.make1D("%s_%s_nBJetMedium"%(channel,suffix), 
        #                         op.c_int(nbJet),
        #                         sel, 
        #                         EquidistantBinning(10,-0.5, 9.5),
        #                         title="nBJetMedium (channel %s)"%channel,
        #                         xTitle= "nBJetMedium",
        #                         plotopts = channelLabel))        

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
def makeSingleLeptonAk8JetsPlots(sel,j1,j2,j3,suffix,channel,nMedBJets,HLL,has1fat=False,has2fat=False,has1fat1slim=False,has1fat2slim=False):
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

    j1_name = "Ak8BJet"
    j2_name = "Ak8LSJet" if has2fat else "LeadAk4NonBJet"
    j3_name = "SubLeadAk4NonBJet"

    #plots.append(Plot.make1D("%s_%s_Hbb_Pt"%(channel,suffix),
    #                         (j1.subJet1.p4 + j1.subJet2.p4).Pt(),
    #                         sel,
    #                         EquidistantBinning(50,0.,500.),
    #                         title='Hbb_Pt',
    #                         xTitle="Hbb_P_{T} [GeV]",
    #                         plotopts = channelLabel))

    
    # 1st-jet : FatJet #
    #plots.append(Plot.make1D("%s_%s_nMedSubBjets_%s"%(channel,suffix,j1_name),
    #                         nMedBJets,
    #                         sel,
    #                         EquidistantBinning(6,-0.5,5.5),
    #                         title='nMediumBtagged SubJets of %s'%j1_name,
    #                         xTitle="nSubJets_MediumBTagged%s"%j1_name,
    #                         plotopts = channelLabel))

    plots.append(Plot.make1D("%s_%s_pT_%s"%(channel,suffix,j1_name),
                             j1.p4.pt(),
                             sel,
                             EquidistantBinning(80,100.,500.),
                             title='Transverse momentum of the %s'%j1_name,
                             xTitle="P_{T} %s[GeV]"%j1_name,
                             plotopts = channelLabel))
    
    #plots.append(Plot.make1D("%s_%s_SubJet1_pT_%s"%(channel,suffix,j1_name),
    #                         j1.subJet1.p4.pt(),
    #                         sel,
    #                         EquidistantBinning(40,0.,200.),
    #                         title='Transverse momentum of the subJet1 of %s'%j1_name,
    #                         xTitle="subJet1_P_{T} %s[GeV]"%j1_name,
    #                         plotopts = channelLabel))
    #plots.append(Plot.make1D("%s_%s_SubJet2_pT_%s"%(channel,suffix,j1_name),
    #                         j1.subJet2.p4.pt(),
    #                         sel,
    #                         EquidistantBinning(40,0.,200.),
    #                         title='Transverse momentum of the subJet2 of %s'%j1_name,
    #                         xTitle="subJet2_P_{T} %s[GeV]"%j1_name,
    #                         plotopts = channelLabel))
    
    plots.append(Plot.make1D("%s_%s_eta_%s"%(channel,suffix,j1_name),
                             j1.p4.eta(),
                             sel,
                             EquidistantBinning(22,-3.,3.),
                             title='Pseudorapidity of the %s'%j1_name,
                             xTitle="#eta %s"%j1_name,
                             plotopts = channelLabel))

    #plots.append(Plot.make1D("%s_%s_SubJet1_eta_%s"%(channel,suffix,j1_name),
    #                         j1.subJet1.p4.eta(),
    #                         sel,
    #                         EquidistantBinning(22,-3.,3.),
    #                         title='Pseudorapidity of the subJet1 of %s'%j1_name,
    #                         xTitle="subJet1_#eta %s[GeV]"%j1_name,
    #                         plotopts = channelLabel))
    #plots.append(Plot.make1D("%s_%s_SubJet2_eta_%s"%(channel,suffix,j1_name),
    #                         j1.subJet2.p4.eta(),
    #                         sel,
    #                         EquidistantBinning(22,-3.,3.),
    #                         title='Pseudorapidity of the subJet2 of %s'%j1_name,
    #                         xTitle="subJet2_#eta %s[GeV]"%j1_name,
    #                         plotopts = channelLabel))

    #plots.append(Plot.make1D("%s_%s_phi_%s"%(channel,suffix,j1_name),
    #                         j1.p4.phi(),
    #                         sel,
    #                         EquidistantBinning(20,-3.2,3.2),
    #                         title='Azimutal angle of the %s'%j1_name,
    #                         xTitle="#phi %s"%j1_name,
    #                         plotopts = channelLabel))
    #plots.append(Plot.make1D("%s_%s_mass_%s"%(channel,suffix,j1_name),
    #                         j1.mass,
    #                         sel,
    #                         EquidistantBinning(40,0.,200.),
    #                         title='mass of the %s'%j1_name,
    #                         xTitle="M %s"%j1_name,
    #                         plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_massSD_%s"%(channel,suffix,j1_name),
                             j1.msoftdrop,
                             sel,
                             EquidistantBinning(40,0.,200.),
                             title='Soft Drop mass of the %s'%j1_name,
                             xTitle="M_{Soft Drop} %s[GeV]"%j1_name,
                             plotopts = channelLabel))

    #plots.append(Plot.make1D("%s_%s_btagDDBvL_%s"%(channel,suffix,j1_name),
    #                         j1.btagDDBvL,
    #                         sel,
    #                         EquidistantBinning(10,0.,1.),
    #                         title='btagDDBvL of %s'%j1_name,
    #                         xTitle="btagDDBvL %s[GeV]"%j1_name,
    #                         plotopts = channelLabel))        
    #plots.append(Plot.make1D("%s_%s_btagDDBvL_noMD_%s"%(channel,suffix,j1_name),
    #                         j1.btagDDBvL_noMD,
    #                         sel,
    #                         EquidistantBinning(10,0.,1.),
    #                         title='btagDDBvL_noMD of %s'%j1_name,
    #                         xTitle="btagDDBvL_noMD %s[GeV]"%j1_name,
    #                         plotopts = channelLabel))        
    #plots.append(Plot.make1D("%s_%s_btagDDCvB_%s"%(channel,suffix,j1_name),
    #                         j1.btagDDCvB,
    #                         sel,
    #                         EquidistantBinning(10,0.,1.),
    #                         title='btagDDCvB of %s'%j1_name,
    #                         xTitle="btagDDCvB %s[GeV]"%j1_name,
    #                         plotopts = channelLabel))        
    #plots.append(Plot.make1D("%s_%s_btagDDCvB_noMD_%s"%(channel,suffix,j1_name),
    #                         j1.btagDDCvB_noMD,
    #                         sel,
    #                         EquidistantBinning(10,0.,1.),
    #                         title='btagDDCvB_noMD of %s'%j1_name,
    #                         xTitle="btagDDCvB_noMD %s[GeV]"%j1_name,
    #                         plotopts = channelLabel))        
    #plots.append(Plot.make1D("%s_%s_btagDDCvL_%s"%(channel,suffix,j1_name),
    #                         j1.btagDDCvL,
    #                         sel,
    #                         EquidistantBinning(10,0.,1.),
    #                         title='btagDDCvL of %s'%j1_name,
    #                         xTitle="btagDDCvL %s[GeV]"%j1_name,
    #                         plotopts = channelLabel))        
    #plots.append(Plot.make1D("%s_%s_btagDDCvL_noMD_%s"%(channel,suffix,j1_name),
    #                         j1.btagDDCvL_noMD,
    #                         sel,
    #                         EquidistantBinning(10,0.,1.),
    #                         title='btagDDCvL_noMD of %s'%j1_name,
    #                         xTitle="btagDDCvL_noMD %s[GeV]"%j1_name,
    #                         plotopts = channelLabel))        
    #plots.append(Plot.make1D("%s_%s_btagDeepB_%s"%(channel,suffix,j1_name),
    #                         j1.btagDeepB,
    #                         sel,
    #                         EquidistantBinning(10,0.,1.),
    #                         title='btagDeepB of %s'%j1_name,
    #                         xTitle="btagDeepB %s[GeV]"%j1_name,
    #                         plotopts = channelLabel))        

    if has2fat or has1fat1slim or has1fat2slim:
        # 2nd-jet : Fat or Slim Jet #
        #plots.append(Plot.make1D("%s_%s_btagDeepB_%s"%(channel,suffix,j2_name),
        #                         j2.btagDeepB,
        #                         sel,
        #                         EquidistantBinning(20,0.,1.),
        #                         title='btagDeepB of %s'%j2_name,
        #                         xTitle="btagDeepB %s[GeV]"%j2_name,
        #                         plotopts = channelLabel))        
        #plots.append(Plot.make1D("%s_%s_qgDiscr_%s"%(channel,suffix,j2_name),
        #                         j2.qgl,
        #                         sel,
        #                         EquidistantBinning(20,0.,1.),
        #                         title='qgDiscr of the %s'%j2_name,
        #                         xTitle="qgDiscr_%s"%j2_name,
        #                         plotopts = channelLabel))
        
        #plots.append(Plot.make1D("%s_%s_pT_%s"%(channel,suffix,j2_name),
        #                         j2.p4.pt(),
        #                         sel,
        #                         EquidistantBinning(40,0.,200.),
        #                         title='Transverse momentum of the %s'%j2_name,
        #                         xTitle="P_{T} %s[GeV]"%j2_name,
        #                         plotopts = channelLabel))
        #plots.append(Plot.make1D("%s_%s_eta_%s"%(channel,suffix,j2_name),
        #                          j2.p4.eta(),
        #                          sel,
        #                          EquidistantBinning(22,-3.,3.),
        #                          title='Pseudorapidity of the %s'%j2_name,
        #                          xTitle="#eta %s"%j2_name,
        #                          plotopts = channelLabel))
        #plots.append(Plot.make1D("%s_%s_phi_%s"%(channel,suffix,j2_name),
        #                         j2.p4.phi(),
        #                         sel,
        #                         EquidistantBinning(20,-3.2,3.2),
        #                         title='Azimutal angle of the %s'%j2_name,
        #                         xTitle="#phi %s"%j2_name,
        #                         plotopts = channelLabel))
        #plots.append(Plot.make1D("%s_%s_mass_%s"%(channel,suffix,j2_name),
        #                         j2.mass,
        #                         sel,
        #                         EquidistantBinning(40,0.,200.),
        #                         title='mass of the %s'%j2_name,
        #                         xTitle="M %s"%j2_name,
        #                         plotopts = channelLabel))
        if has2fat:
            plots.append(Plot.make1D("%s_%s_massSD_%s"%(channel,suffix,j2_name),
                                     j2.msoftdrop,
                                     sel,
                                     EquidistantBinning(40,0.,200.),
                                     title='Soft Drop mass of the %s'%j2_name,
                                     xTitle="M_{Soft Drop} %s[GeV]"%j2_name,
                                     plotopts = channelLabel))        

        #plots.append(Plot.make1D("%s_%s_DijetPT_%s_%s"%(channel,suffix,j1_name,j2_name), 
        #                         (j1.p4+j2.p4).Pt(), 
        #                         sel, 
        #                         EquidistantBinning(80,0.,400.),
        #                         title="Transverse momentum of the dijet_%s_%s (channel %s)"%(j1_name,j2_name,channel), 
        #                         xTitle= "DiJetPT [%s_%s] [GeV]"%(j1_name,j2_name),
        #                         plotopts = channelLabel))#

        # DeltaR #
        #plots.append(Plot.make1D("%s_%s_DeltaR_%s_%s"%(channel,suffix,j1_name,j2_name), 
        #                         op.deltaR(j1.p4,j2.p4), 
        #                         sel, 
        #                         EquidistantBinning(50, 0., 5.), 
        #                         title="Dijet DeltaR_%s_%s (channel %s)"%(j1_name,j2_name,channel), 
        #                         xTitle= "DeltaR [%s_%s]"%(j1_name,j2_name),
        #                         plotopts = channelLabel))
        ## DeltaPhi #
        #plots.append(Plot.make1D("%s_%s_DeltaPhi_%s_%s"%(channel,suffix,j1_name,j2_name), 
        #                         op.abs(op.deltaPhi(j1.p4,j2.p4)), 
        #                         sel, 
        #                         EquidistantBinning(20, 0., 3.2), 
        #                         title="Dijet DeltaPhi_%s_%s (channel %s)"%(j1_name,j2_name,channel), 
        #                         xTitle= "DeltaPhi [%s_%s]"%(j1_name,j2_name),
        #                         plotopts = channelLabel))

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
            #plots.append(Plot.make1D("%s_%s_btagDeepB_%s"%(channel,suffix,j3_name),
            #                         j3.btagDeepB,
            #                         sel,
            #                         EquidistantBinning(20,0.,1.),
            #                         title='btagDeepB of %s'%j3_name,
            #                         xTitle="btagDeepB %s[GeV]"%j3_name,
            #                         plotopts = channelLabel))        
            #plots.append(Plot.make1D("%s_%s_qgDiscr_%s"%(channel,suffix,j3_name),
            #                         j3.qgl,
            #                         sel,
            #                         EquidistantBinning(20,0.,1.),
            #                         title='qgDiscr of the %s'%j3_name,
            #                         xTitle="qgDiscr_%s"%j3_name,
            #                         plotopts = channelLabel))
            #plots.append(Plot.make1D("%s_%s_HadW_mass"%(channel,suffix),
            #                         HLL.Wjj_simple(j2.p4, j3.p4).M(), 
            #                         sel,
            #                         EquidistantBinning(25, 0., 250.), 
            #                         title="HadW_mass (channel %s)"%channel, 
            #                         xTitle= "HadW_mass [GeV]",
            #                         plotopts = channelLabel))
            #plots.append(Plot.make1D("%s_%s_HadW_cosTheta"%(channel,suffix), 
            #                         HLL.comp_cosThetaS(j2.p4, j3.p4),
            #                         sel, 
            #                         EquidistantBinning(20,0.,1.),
            #                         title="HadW_cosTheta (channel %s)"%channel, 
            #                         xTitle= "HadW_cosTheta",
            #                         plotopts = channelLabel))
            plots.append(Plot.make1D("%s_%s_HadW_cosTheta_pDavid"%(channel,suffix), 
                                     op.extMethod("HHbbWWJPA::cosThetaS", returnType="float")(j2.p4, j3.p4),
                                     sel, 
                                     EquidistantBinning(20,0.,1.),
                                     title="HadW_cosTheta (channel %s)"%channel, 
                                     xTitle= "HadW_cosTheta",
                                     plotopts = channelLabel))
            #plots.append(Plot.make1D("%s_%s_min_dR_HadW_bjet"%(channel,suffix), 
            #                         op.min(HLL.dR_HadW_bjet(j1.subJet1.p4, j2.p4, j3.p4), 
            #                                HLL.dR_HadW_bjet(j1.subJet2.p4, j2.p4, j3.p4)),
            #                         sel, 
            #                         EquidistantBinning(50,0.,5.),
            #                         title="min_dR_HadW_bjet (channel %s)"%channel, 
            #                         xTitle= "min_dR_HadW_bjet",
            #                         plotopts = channelLabel))
            #plots.append(Plot.make1D("%s_%s_dR_wjet1wjet2"%(channel,suffix), 
            #                         op.deltaR(j2.p4, j3.p4),
            #                         sel, 
            #                         EquidistantBinning(50,0.,5.),
            #                         title="dR_wjet1wjet2 (channel %s)"%channel, 
            #                         xTitle= "dR_wjet1wjet2",
            #                         plotopts = channelLabel))

            #plots.append(Plot.make1D("%s_%s_DijetPT_%s_%s"%(channel,suffix,j2_name,j3_name), 
            #                         (j2.p4+j3.p4).Pt(), 
            #                         sel, 
            #                         EquidistantBinning(80,0.,400.),
            #                         title="Transverse momentum of the dijet_%s_%s (channel %s)"%(j2_name,j3_name,channel), 
            #                         xTitle= "DiJetPT [%s_%s] [GeV]"%(j2_name,j3_name),
            #                         plotopts = channelLabel))
            
            # DeltaR #
            #plots.append(Plot.make1D("%s_%s_DeltaR_%s_%s"%(channel,suffix,j2_name,j3_name), 
            #                         op.deltaR(j2.p4,j3.p4), 
            #                         sel, 
            #                         EquidistantBinning(50, 0., 5.), 
            #                         title="Dijet DeltaR_%s_%s (channel %s)"%(j2_name,j3_name,channel), 
            #                         xTitle= "DeltaR [%s_%s]"%(j2_name,j3_name),
            #                         plotopts = channelLabel))
            #plots.append(Plot.make1D("%s_%s_DeltaR_%s_%s"%(channel,suffix,'Hbb','Wjj'), 
            #                         op.deltaR(j1.p4,(j2.p4+j3.p4)), 
            #                         sel, 
            #                         EquidistantBinning(50, 0., 5.), 
            #                         title="DeltaR_%s_%s (channel %s)"%('Hbb','Wjj',channel), 
            #                         xTitle= "DeltaR [%s_%s]"%('Hbb','Wjj'),
            #                         plotopts = channelLabel))
            ## DeltaPhi #
            #plots.append(Plot.make1D("%s_%s_DeltaPhi_%s_%s"%(channel,suffix,j2_name,j3_name), 
            #                         op.abs(op.deltaPhi(j2.p4,j3.p4)), 
            #                         sel, 
            #                         EquidistantBinning(20, 0., 3.2), 
            #                         title="Dijet DeltaPhi_%s_%s (channel %s)"%(j2_name,j3_name,channel), 
            #                         xTitle= "DeltaPhi [%s_%s]"%(j2_name,j3_name),
            #                         plotopts = channelLabel))
            #plots.append(Plot.make1D("%s_%s_DeltaPhi_%s_%s"%(channel,suffix,'Hbb','Wjj'), 
            #                         op.abs(op.deltaPhi(j1.p4,(j2.p4+j3.p4))), 
            #                         sel, 
            #                         EquidistantBinning(20, 0., 3.2), 
            #                         title="DeltaPhi_%s_%s (channel %s)"%('Hbb','Wjj',channel), 
            #                         xTitle= "DeltaPhi [%s_%s]"%('Hbb','Wjj'),
            #                         plotopts = channelLabel))

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
                             EquidistantBinning(25,0.,250.),
                             title='Transverse mass of lepton and MET (%s channel)'%channel,
                             xTitle="M_{T}(lep,MET) [GeV]",
                             plotopts = channelLabel))
    # lepton-MET plots #
    #plots.append(Plot.make1D("%s_%s_highlevelvariable_SinglepMET_Pt"%(channel,suffix),
    #                         op.abs(HLL.SinglepMet_Pt(lep, met)),
    #                         sel,
    #                         EquidistantBinning(20,0.,3.2),
    #                         title='Reultant pT of lepton and MET (%s channel)'%channel,
    #                         xTitle="|Pt (lep,MET)|",
    #                         plotopts = channelLabel))
    #plots.append(Plot.make1D("%s_%s_highlevelvariable_SinglepMETdeltaPhi"%(channel,suffix),
    #                         op.abs(HLL.SinglepMet_dPhi(lep, met)),
    #                         sel,
    #                         EquidistantBinning(20,0.,3.2),
    #                         title='Azimutal angle between lepton and MET (%s channel)'%channel,
    #                         xTitle="|#Delta \phi (lep,MET)|",
    #                         plotopts = channelLabel))
    
    # jet-MET plots
    #plots.append(Plot.make1D("%s_%s_highlevelvariable_%s_METdeltaPhi"%(channel,suffix,j1_name),
    #                         op.abs(HLL.SinglepMet_dPhi(j1, met)),
    #                         sel,
    #                         EquidistantBinning(20,0.,3.2),
    #                         title='Azimutal angle between %s and MET (%s channel)'%(j1_name,channel),
    #                         xTitle="|#Delta \phi (%s,MET)|"%j1_name,
    #                         plotopts = channelLabel))

    # lepton_Jet DeltaR plots #
    #plots.append(Plot.make1D("%s_%s_highlevelvariable_DeltaR_lep_%s"%(channel,suffix,j1_name),
    #                         op.deltaR(lep.p4,j1.p4),
    #                         sel,
    #                         EquidistantBinning(50,0.,5.),
    #                         title="DeltaR_%s_%s (channel %s)"%(channel,j1_name,channel),
    #                         xTitle="DeltaR [%s_%s]"%(channel,j1_name),
    #                         plotopts = channelLabel))

    # lepton_Jet DeltaPhi plots #
    #plots.append(Plot.make1D("%s_%s_highlevelvariable_DeltaPhi_lep_%s"%(channel,suffix,j1_name),
    #                         op.abs(op.deltaPhi(lep.p4,j1.p4)),
    #                         sel,
    #                         EquidistantBinning(20,0.,3.2),
    #                         title="DeltaPhi_%s_%s (channel %s)"%(channel,j1_name,channel),
    #                         xTitle="DeltaPhi [%s_%s]"%(channel,j1_name),
    #                         plotopts=channelLabel))

    if nJet > 1 :
        '''
        plots.append(Plot.make1D("%s_%s_highlevelvariable_%s_METdeltaPhi"%(channel,suffix,j2_name),
                                 op.abs(HLL.SinglepMet_dPhi(j2, met)),
                                 sel,
                                 EquidistantBinning(20,0.,3.2),
                                 title='Azimutal angle between %s and MET (%s channel)'%(j2_name,channel),
                                 xTitle="|#Delta \phi (%s,MET)|"%j2_name,
                                 plotopts = channelLabel))

        plots.append(Plot.make1D("%s_%s_highlevelvariable_DeltaR_lep_%s"%(channel,suffix,j2_name),
                                 op.deltaR(lep.p4,j2.p4),
                                 sel,
                                 EquidistantBinning(50,0.,5.),
                                 title="DeltaR_%s_%s (channel %s)"%(channel,j2_name,channel),
                                 xTitle="DeltaR [%s_%s]"%(channel,j2_name),
                                 plotopts = channelLabel))
        
        plots.append(Plot.make1D("%s_%s_highlevelvariable_DeltaPhi_lep_%s"%(channel,suffix,j2_name),
                                 op.abs(op.deltaPhi(lep.p4,j2.p4)),
                                 sel,
                                 EquidistantBinning(20,0.,3.2),
                                 title="DeltaPhi_%s_%s (channel %s)"%(channel,j2_name,channel),
                                 xTitle="DeltaPhi [%s_%s]"%(channel,j2_name),
                                 plotopts=channelLabel))
        '''
    if nJet > 2 :
        '''
        plots.append(Plot.make1D("%s_%s_highlevelvariable_%s_METdeltaPhi"%(channel,suffix,j3_name),
                                 op.abs(HLL.SinglepMet_dPhi(j3, met)),
                                 sel,
                                 EquidistantBinning(20,0.,3.2),
                                 title='Azimutal angle between %s and MET (%s channel)'%(j3_name,channel),
                                 xTitle="|#Delta \phi (%s,MET)|"%j3_name,
                                 plotopts = channelLabel))
    

        plots.append(Plot.make1D("%s_%s_highlevelvariable_DeltaR_lep_%s"%(channel,suffix,j3_name),
                                 op.deltaR(lep.p4,j3.p4),
                                 sel,
                                 EquidistantBinning(50,0.,5.),
                                 title="DeltaR_%s_%s (channel %s)"%(channel,j3_name,channel),
                                 xTitle="DeltaR [%s_%s]"%(channel,j3_name),
                                 plotopts = channelLabel))
        
        # min-DR
        plots.append(Plot.make1D("%s_%s_highlevelvariable_minDeltaR_lep_jets123"%(channel,suffix),
                                 HLL.MinDR_lep3j(lep,j1,j2,j3),
                                 sel,
                                 EquidistantBinning(50,0.,5.),
                                 title="minDeltaR_lep_jet (channel %s)"%channel,
                                 xTitle="minDeltaR [%s_Jet]"%channel,
                                 plotopts = channelLabel))

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
        '''
    
        # Scalar magnitude sum #
        #plots.append(Plot.make1D("%s_%s_highlevelvariable_HT2_l3jmet"%(channel,suffix),
        #                         HLL.HT2_l3jmet(lep,j1,j2,j3,met),
        #                         sel,
        #                         EquidistantBinning(60,0.,600.),
        #                         title='Di-Higgs magnitude_l3jmet (%s channel)'%channel,
        #                         xTitle="H_{T2} [GeV]",
        #                         plotopts = channelLabel))
        #plots.append(Plot.make1D("%s_%s_highlevelvariable_HT2R_l3jmet"%(channel,suffix),
        #                         HLL.HT2R_l3jmet(lep,j1,j2,j3,met),
        #                         sel,
        #                         EquidistantBinning(40,0.,1.),
        #                         title='Di-Higgs magnitude_Ratio_l3jmet (%s channel)'%channel,
        #                         xTitle="H_{T2}_Ratio [GeV]",
        #                         plotopts = channelLabel))
    if nJet > 3:
        '''
        plots.append(Plot.make1D("%s_%s_DeltaR_lep_%s"%(channel,suffix,'HadW'),
                                 HLL.Wjj_simple(j3.p4, j4.p4).M(),
                                 sel,
                                 EquidistantBinning(50,0.,5.),
                                 title="DeltaR_%s_%s (channel %s)"%(channel,j4_name,channel),
                                 xTitle="dR_HadW_%s"%channel,
                                 plotopts=channelLabel))

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
        '''
        plots.append(Plot.make1D("%s_%s_highlevelvariable_HT2_l4jmet"%(channel,suffix),
                                 HLL.HT2_l4jmet(lep,j1,j2,j3,j4,met),
                                 sel,
                                 EquidistantBinning(60,0.,600.),
                                 title='Di-Higgs magnitude_l4jmet (%s channel)'%channel,
                                 xTitle="H_{T2} [GeV]",
                                 plotopts = channelLabel))
        #plots.append(Plot.make1D("%s_%s_highlevelvariable_HT2R_l4jmet"%(channel,suffix),
        #                         HLL.HT2R_l4jmet(lep,j1,j2,j3,j4,met),
        #                         sel,
        #                         EquidistantBinning(40,0.,1.),
        #                         title='Di-Higgs magnitude_Ratio_l4jmet (%s channel)'%channel,
        #                         xTitle="H_{T2}_Ratio [GeV]",
        #                         plotopts = channelLabel))

        '''
        # cosThetaS variables
        plots.append(Plot.make1D("%s_%s_highlevelvariable_cosThetaS_Wjj_simple"%(channel,suffix),
                                 HLL.comp_cosThetaS(j3.p4, j4.p4),
                                 sel,
                                 EquidistantBinning(50,0.,1.),
                                 title='cosThetaS_Wjj_simple (%s channel)'%channel,
                                 xTitle="cosThetaS_Wjj_simple",
                                 plotopts = channelLabel))

        plots.append(Plot.make1D("%s_%s_highlevelvariable_cosThetaS_WW_simple_met"%(channel,suffix),
                                 HLL.comp_cosThetaS(HLL.Wjj_simple(j3.p4,j4.p4),HLL.Wlep_met_simple(lep.p4,met.p4)),
                                 sel,
                                 EquidistantBinning(50,0.,1.),
                                 title='cosThetaS_WW_simple_met (%s channel)'%channel,
                                 xTitle="cosThetaS_WW_simple_met",
                                 plotopts = channelLabel))
        #####################
        plots.append(Plot.make1D("%s_%s_highlevelvariable_cosThetaS_Wjj_simple_pDavid"%(channel,suffix),
                                 op.extMethod("HHbbWWJPA::cosThetaS", returnType="float")(j3.p4, j4.p4),
                                 sel,
                                 EquidistantBinning(50,0.,1.),
                                 title='cosThetaS_Wjj_simple (%s channel)'%channel,
                                 xTitle="cosThetaS_Wjj_simple",
                                 plotopts = channelLabel))

        plots.append(Plot.make1D("%s_%s_highlevelvariable_cosThetaS_WW_simple_met_pDavid"%(channel,suffix),
                                 op.extMethod("HHbbWWJPA::cosThetaS", returnType="float")(HLL.Wjj_simple(j3.p4,j4.p4), HLL.Wlep_met_simple(lep.p4,met.p4)),
                                 sel,
                                 EquidistantBinning(50,0.,1.),
                                 title='cosThetaS_WW_simple_met (%s channel)'%channel,
                                 xTitle="cosThetaS_WW_simple_met",
                                 plotopts = channelLabel))
        #####################
        # Transverse mass plots #
        #if nbJet == 2:
        plots.append(Plot.make1D("%s_%s_highlevelvariable_MT_lep_%s_%s"%(channel,suffix,j3_name,j4_name),
                                 HLL.MT_W1W2_ljj(lep,j3,j4,met),
                                 sel,
                                 EquidistantBinning(50,0.,1000.),
                                 title='Transverse mass of lepton+dijet and MET (%s channel)'%channel,
                                 xTitle="M_{T}(lepjj,MET) [GeV]",
                                 plotopts = channelLabel))
        
        # BDT variables
        plots.append(Plot.make1D("%s_%s_highlevelvariable_mHH_simple_met"%(channel,suffix),
                                 HLL.HHP4_simple_met(HLL.bJetCorrP4(j1)+HLL.bJetCorrP4(j2), j3.p4, j4.p4, lep.p4, met.p4).M(),
                                 sel,
                                 EquidistantBinning(50,0.,1000.),
                                 title='mHH_simple_met (%s channel)'%channel,
                                 xTitle="mHH_simple_met [GeV]",
                                 plotopts = channelLabel))
        
        plots.append(Plot.make1D("%s_%s_highlevelvariable_mT_top_3particle"%(channel,suffix),
                                 op.min(HLL.mT2(HLL.bJetCorrP4(j1),lep.p4,met.p4), HLL.mT2(HLL.bJetCorrP4(j2),lep.p4,met.p4)),
                                 sel,
                                 EquidistantBinning(50,0.,1000.),
                                 title='mT_top_3particle (%s channel)'%channel,
                                 xTitle="mT_top_3particle [GeV]",
                                 plotopts = channelLabel))
        
        # cosThetaS variables
        plots.append(Plot.make1D("%s_%s_highlevelvariable_cosThetaS_Hbb"%(channel,suffix),
                                 HLL.comp_cosThetaS(HLL.bJetCorrP4(j1),HLL.bJetCorrP4(j2)),
                                 sel,
                                 EquidistantBinning(50,0.,1.),
                                 title='cosThetaS_Hbb (%s channel)'%channel,
                                 xTitle="cosThetaS_Hbb",
                                 plotopts = channelLabel))
        
        plots.append(Plot.make1D("%s_%s_highlevelvariable_cosThetaS_HH_simple_met"%(channel,suffix),
                                 HLL.comp_cosThetaS(HLL.bJetCorrP4(j1)+HLL.bJetCorrP4(j2), HLL.HWW_met_simple(j3.p4,j4.p4,lep.p4,met.p4)),
                                 sel,
                                 EquidistantBinning(50,0.,1.),
                                 title='cosThetaS_HH_simple_met (%s channel)'%channel,
                                 xTitle="cosThetaS_HH_simple_met",
                                 plotopts = channelLabel))
        #####################
        plots.append(Plot.make1D("%s_%s_highlevelvariable_cosThetaS_Hbb_pDavid"%(channel,suffix),
                                 op.extMethod("HHbbWWJPA::cosThetaS", returnType="float")(j1.p4, j2.p4),
                                 sel,
                                 EquidistantBinning(50,0.,1.),
                                 title='cosThetaS_Hbb (%s channel)'%channel,
                                 xTitle="cosThetaS_Hbb",
                                 plotopts = channelLabel))
        
        plots.append(Plot.make1D("%s_%s_highlevelvariable_cosThetaS_HH_simple_met_pDavid"%(channel,suffix),
                                 op.extMethod("HHbbWWJPA::cosThetaS", returnType="float")(HLL.bJetCorrP4(j1)+HLL.bJetCorrP4(j2), HLL.HWW_met_simple(j3.p4,j4.p4,lep.p4,met.p4)),
                                 sel,
                                 EquidistantBinning(50,0.,1.),
                                 title='cosThetaS_HH_simple_met (%s channel)'%channel,
                                 xTitle="cosThetaS_HH_simple_met",
                                 plotopts = channelLabel))
        #####################
        '''
    return plots 

# Make HighLevel plots for semi boosted and boosted categories
def makeHighLevelPlotsBoosted(sel,met,lep,j1,j2,j3,channel,suffix,HLL,has1fat2slim,bothAreFat):
    #print ("===============>>>HighLevel Plots for Boosted__channel:%s__sel:%s"%(channel,suffix))
    if bothAreFat:
        j1_name   = "Ak8BJet"
        j2_name   = "Ak8Jet"
        
    else:
        j1_name   = "Ak8BJet"
        j2_name   = "Ak4LeadingNonbJet"
        j3_name   = "Ak4SubLeadingNonbJet"
    
    plots = []
    channelLabel = SingleLeptonChannelTitleLabel(channel)
    '''
    # dR between lepton and fat-bJet
    plots.append(Plot.make1D("%s_%s_highlevelvariable_DeltaR_%s_%s"%(channel,suffix,channel,j1_name),
                             op.deltaR(lep.p4,j1.p4),
                             sel,
                             EquidistantBinning(50,0.,5.),
                             title="DeltaR_(%s_%s)"%(channel, j1_name),
                             xTitle="DeltaR [%s_%s]"%(channel, j1_name),
                             plotopts=channelLabel))
    # dR between lepton and fat-bSubJet1
    plots.append(Plot.make1D("%s_%s_highlevelvariable_DeltaR_%s_%s_SubJet1"%(channel,suffix,channel,j1_name),
                             op.deltaR(lep.p4,j1.subJet1.p4),
                             sel,
                             EquidistantBinning(50,0.,5.),
                             title="DeltaR_(%s_%s)_SubJet1"%(channel, j1_name),
                             xTitle="DeltaR [%s_%s]_SubJet1"%(channel, j1_name),
                             plotopts=channelLabel))
    # dR between lepton and fat-bSubJet2
    plots.append(Plot.make1D("%s_%s_highlevelvariable_DeltaR_%s_%s_SubJet2"%(channel,suffix,channel,j1_name),
                             op.deltaR(lep.p4,j1.subJet2.p4),
                             sel,
                             EquidistantBinning(50,0.,5.),
                             title="DeltaR_(%s_%s)_SubJet1"%(channel, j1_name),
                             xTitle="DeltaR [%s_%s]_SubJet1"%(channel, j1_name),
                             plotopts=channelLabel))

    # dPhi between lepton and fat-bJet
    plots.append(Plot.make1D("%s_%s_highlevelvariable_DeltaPhi_%s_%s"%(channel,suffix,channel,j1_name),
                             op.abs(op.deltaPhi(lep.p4,j1.p4)),
                             sel,
                             EquidistantBinning(20,0.,3.2),
                             title="DeltaPhi_(%s_%s)"%(channel, j1_name),
                             xTitle="DeltaPhi [%s_%s]"%(channel, j1_name),
                             plotopts=channelLabel))
    # dPhi between lepton and fat-b-SubJet1
    plots.append(Plot.make1D("%s_%s_highlevelvariable_DeltaPhi_%s_%s_SubJet1"%(channel,suffix,channel,j1_name),
                             op.abs(op.deltaPhi(lep.p4,j1.subJet1.p4)),
                             sel,
                             EquidistantBinning(20,0.,3.2),
                             title="DeltaPhi_(%s_%s)_SubJet1"%(channel, j1_name),
                             xTitle="DeltaPhi [%s_%s]_SubJet1"%(channel, j1_name),
                             plotopts=channelLabel))
    # dPhi between lepton and fat-b-SubJet2
    plots.append(Plot.make1D("%s_%s_highlevelvariable_DeltaPhi_%s_%s_SubJet2"%(channel,suffix,channel,j1_name),
                             op.abs(op.deltaPhi(lep.p4,j1.subJet2.p4)),
                             sel,
                             EquidistantBinning(20,0.,3.2),
                             title="DeltaPhi_(%s_%s)_SubJet2"%(channel, j1_name),
                             xTitle="DeltaPhi [%s_%s]_SubJet2"%(channel, j1_name),
                             plotopts=channelLabel))

    # Resultant pT of Lepton and MET
    plots.append(Plot.make1D("%s_%s_highlevelvariable_SinglepMET_Pt"%(channel,suffix),
                             op.abs(HLL.SinglepMet_Pt(lep, met)),
                             sel,
                             EquidistantBinning(20,0.,3.2),
                             title='Resultant pT of lepton and MET (%s channel)'%channel,
                             xTitle="|Pt (lep,MET)|",
                             plotopts = channelLabel))
    # DPhi between Lepton and MET
    plots.append(Plot.make1D("%s_%s_highlevelvariable_SinglepMETdeltaPhi"%(channel,suffix),
                             op.abs(HLL.SinglepMet_dPhi(lep, met)),
                             sel,
                             EquidistantBinning(20,0.,3.2),
                             title='Azimutal angle between lepton and MET (%s channel)'%channel,
                             xTitle="|#Delta \phi (lep,MET)|",
                             plotopts = channelLabel))
    # DPhi between fat-bJet and MET
    plots.append(Plot.make1D("%s_%s_highlevelvariable_%s_METdeltaPhi"%(channel,suffix,j1_name),
                             op.abs(HLL.SinglepMet_dPhi(j1, met)),
                             sel,
                             EquidistantBinning(20,0.,3.2),
                             title='Azimutal angle between %s and MET (%s channel)'%(j1_name,channel),
                             xTitle="|#Delta \phi (%s,MET)|"%j1_name,
                             plotopts = channelLabel))
    # DPhi between fat-b_SubJet1 and MET 
    plots.append(Plot.make1D("%s_%s_highlevelvariable_%s_SubJet1_METdeltaPhi"%(channel,suffix,j1_name),
                             op.abs(HLL.SinglepMet_dPhi(j1.subJet1, met)),
                             sel,
                             EquidistantBinning(20,0.,3.2),
                             title='Azimutal angle between %s SubJet1 and MET (%s channel)'%(j1_name,channel),
                             xTitle="|#Delta \phi (%s SubJet1, MET)|"%j1_name,
                             plotopts = channelLabel))
    # DPhi between fat-b_SubJet2 and MET 
    plots.append(Plot.make1D("%s_%s_highlevelvariable_%s_SubJet2_METdeltaPhi"%(channel,suffix,j1_name),
                             op.abs(HLL.SinglepMet_dPhi(j1.subJet2, met)),
                             sel,
                             EquidistantBinning(20,0.,3.2),
                             title='Azimutal angle between %s SubJet2 and MET (%s channel)'%(j1_name,channel),
                             xTitle="|#Delta \phi (%s SubJet2, MET)|"%j1_name,
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
 
    # cosThetaS
    plots.append(Plot.make1D("%s_%s_highlevelvariable_cosThetaS_Hbb"%(channel,suffix),
                             HLL.comp_cosThetaS(j1.subJet1.p4, j1.subJet2.p4),
                             sel,
                             EquidistantBinning(50,0.,1.),
                             title='cosThetaS_Hbb (%s channel)'%channel,
                             xTitle="cosThetaS_Hbb",
                             plotopts = channelLabel))
    # cosThetaS
    plots.append(Plot.make1D("%s_%s_highlevelvariable_cosThetaS_Hbb_pDavid"%(channel,suffix),
                             op.extMethod("HHbbWWJPA::cosThetaS", returnType="float")(j1.subJet1.p4, j1.subJet2.p4),
                             sel,
                             EquidistantBinning(50,0.,1.),
                             title='cosThetaS_Hbb (%s channel)'%channel,
                             xTitle="cosThetaS_Hbb",
                             plotopts = channelLabel))
    
   

    if bothAreFat:
        # dR between lepton and the 2nd Jet
        plots.append(Plot.make1D("%s_%s_highlevelvariable_DeltaR_%s_%s"%(channel,suffix,channel,j2_name),
                                 op.deltaR(lep.p4,j2.p4),
                                 sel,
                                 EquidistantBinning(50,0.,5.),
                                 title="DeltaR_(%s_%s)"%(channel,j1_name),
                                 xTitle="DeltaR [%s_%s]"%(channel, j2_name),
                                 plotopts=channelLabel))

        # dPhi between lepton and 2nd jet
        plots.append(Plot.make1D("%s_%s_highlevelvariable_DeltaPhi_%s_%s"%(channel,suffix,channel,j2_name),
                                 op.abs(op.deltaPhi(lep.p4,j2.p4)),
                                 sel,
                                 EquidistantBinning(20,0.,3.2),
                                 title="DeltaPhi(%s %s)"%(channel, j2_name),
                                 xTitle="DeltaPhi [%s_%s]"%(channel, j2_name),
                                 plotopts=channelLabel))
    
        # DPhi between 2ndJet and MET 
        plots.append(Plot.make1D("%s_%s_highlevelvariable_%s_METdeltaPhi"%(channel,suffix,j2_name),
                                 op.abs(HLL.SinglepMet_dPhi(j2, met)),
                                 sel,
                                 EquidistantBinning(20,0.,3.2),
                                 title='Azimutal angle between %s and MET (%s channel)'%(j2_name,channel),
                                 xTitle="|#Delta \phi (%s,MET)|"%j2_name,
                                 plotopts = channelLabel))

    if has1fat2slim:
        # cosThetaS variables
        plots.append(Plot.make1D("%s_%s_highlevelvariable_cosThetaS_Wjj_simple"%(channel,suffix),
                                 HLL.comp_cosThetaS(j2.p4, j3.p4),
                                 sel,
                                 EquidistantBinning(50,0.,1.),
                                 title='cosThetaS_Wjj_simple (%s channel)'%channel,
                                 xTitle="cosThetaS_Wjj_simple",
                                 plotopts = channelLabel))

        plots.append(Plot.make1D("%s_%s_highlevelvariable_cosThetaS_WW_simple_met"%(channel,suffix),
                                 HLL.comp_cosThetaS(HLL.Wjj_simple(j2.p4,j3.p4),HLL.Wlep_met_simple(lep.p4,met.p4)),
                                 sel,
                                 EquidistantBinning(50,0.,1.),
                                 title='cosThetaS_WW_simple_met (%s channel)'%channel,
                                 xTitle="cosThetaS_WW_simple_met",
                                 plotopts = channelLabel))

        plots.append(Plot.make1D("%s_%s_highlevelvariable_cosThetaS_HH_simple_met"%(channel,suffix),
                                 HLL.comp_cosThetaS(j1.subJet1.p4 + j1.subJet2.p4, HLL.HWW_met_simple(j2.p4,j3.p4,lep.p4,met.p4)),
                                 sel,
                                 EquidistantBinning(50,0.,1.),
                                 title='cosThetaS_HH_simple_met (%s channel)'%channel,
                                 xTitle="cosThetaS_HH_simple_met",
                                 plotopts = channelLabel))

        plots.append(Plot.make1D("%s_%s_Hww_mass"%(channel,suffix), 
                                 HLL.HWW_simple(j2.p4, j3.p4, lep.p4, met).M(),
                                 sel, 
                                 EquidistantBinning(50,0.,500.),
                                 title="Hww_mass (channel %s)"%channel, 
                                 xTitle= "Hww_mass",
                                 plotopts = channelLabel))
        plots.append(Plot.make1D("%s_%s_dR_HadW_lep"%(channel,suffix), 
                                 op.deltaR(HLL.Wjj_simple(j2.p4, j3.p4), lep.p4),
                                 sel, 
                                 EquidistantBinning(50,0.,5.),
                                 title="dR_HadW_lep (channel %s)"%channel, 
                                 xTitle= "dR_HadW_lep",
                                 plotopts = channelLabel))
    '''
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
    #    plots.append(Plot.make2D("%s_%s_highlevelvariable_firstLeptonPtVSLeadjetPt"%(channel,suffix),
    #                             [l1.pt,j1.pt],
    #                             sel,
    #                             [EquidistantBinning(60,0.,300.),EquidistantBinning(60,0.,300.)],
    #                             xTitle="First lepton P_{T}",
    #                             yTitle="Leading jet P_{T}",
    #                             plotopts = channelLabel))
    #
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


####################################  Machine Learning [BDT implementation :: Single Lepton] #####################################
def makeSingleLeptonMachineLearningPlotsBDTfullRecoResolved(sel,fakeLepColl,lep,met,jets,bJets,lJets,j1,j2,j3,j4,suffix,channel,model_even,model_odd,nBins,event,HLL):
    plots = []
    channelLabel = SingleLeptonChannelTitleLabel(channel)

    #inputs = [op.c_float(0.)]*21
    inputs = [HLL.mindr_lep1_jet(lep,jets),                                                                  # mindr_lep1_jet
              op.invariant_mass(HLL.bJetCorrP4(j1), HLL.bJetCorrP4(j2)),                                     # m_Hbb_regCorr 
              HLL.HHP4_simple_met(HLL.bJetCorrP4(j1)+HLL.bJetCorrP4(j2), j3.p4, j4.p4, lep.p4, met.p4).M(),  # mHH_simple_met 
              HLL.Wlep_met_simple(lep.p4, met.p4).M(),                                                       # mWlep_met_simple
              HLL.HWW_met_simple(j3.p4,j4.p4,lep.p4,met.p4).M(),                                             # mWW_simple_met
              HLL.Wjj_simple(j3.p4, j4.p4).M(),                                                              # mWjj_simple
              HLL.comp_cosThetaS(HLL.bJetCorrP4(j1),HLL.bJetCorrP4(j2)),                                     # cosThetaS_Hbb                   
              HLL.comp_cosThetaS(j3.p4, j4.p4),                                                                          # cosThetaS_Wjj_simple     
              HLL.comp_cosThetaS(HLL.Wjj_simple(j3.p4,j4.p4),HLL.Wlep_met_simple(lep.p4,met.p4)),                        # cosThetaS_WW_simple_met  
              HLL.comp_cosThetaS(HLL.bJetCorrP4(j1)+HLL.bJetCorrP4(j2), HLL.HWW_met_simple(j3.p4,j4.p4,lep.p4,met.p4)),  # cosThetaS_HH_simple_met     
              op.rng_len(jets),                                                                       # nJet
              op.rng_len(bJets),                                                                      # nBJetMedium
              op.deltaR(j1.p4, lep.p4),                                                               # dR_b1lep
              op.deltaR(j2.p4, lep.p4),                                                               # dR_b2lep
              HLL.lambdaConePt(lep),                                                                  # lep_conePt
              HLL.bJetCorrP4(j1).Pt(),                                                                # selJet1_Hbb_pT
              HLL.bJetCorrP4(j2).Pt(),                                                                # selJet2_Hbb_pT 
              HLL.MET_LD(met,jets,fakeLepColl),                                                       # met_LD
              HLL.HTfull(fakeLepColl,HLL.bJetCorrP4(j1),HLL.bJetCorrP4(j2),j3.p4,j4.p4),              # HT
              op.min(HLL.mT2(HLL.bJetCorrP4(j1),lep.p4 ,met.p4), HLL.mT2(HLL.bJetCorrP4(j2), lep.p4, met.p4)), # mT_top_3particle
              HLL.MT(lep,met)                                                                         # mT_W
          ]

    output = op.switch(event%2,model_odd(*inputs),model_even(*inputs))

    plots.append(Plot.make1D("%s_%s_nClasses"%(channel,suffix),
                             op.rng_len(output),
                             sel,
                             EquidistantBinning(10,-0.5,9.5),
                             xTitle = 'nClasses (nBDT outputs)',
                             plotopts = channelLabel))

    plots.append(Plot.make1D("%s_%s_BDTOutput"%(channel,suffix),
                             output[0],
                             sel,
                             EquidistantBinning(nBins,0.,1.),
                             xTitle = 'BDT Response',
                             plotopts = channelLabel))

    return plots


def makeSingleLeptonMachineLearningPlotsBDTfullRecoBoosted(sel,fakeLepColl,lep,met,jets,bJets,lJets,j1,j2,j3,j4,suffix,channel,model_even,model_odd,nBins,nMedBJets,event,HLL):
    plots = []
    channelLabel = SingleLeptonChannelTitleLabel(channel)

    #inputs = [op.c_float(0.)]*21
    inputs = [HLL.mindr_lep1_jet(lep,jets),                                                                  # mindr_lep1_jet
              op.invariant_mass(j1.p4, j2.p4),                                                               # m_Hbb_regCorr 
              HLL.HHP4_simple_met(j1.p4+j2.p4, j3.p4, j4.p4, lep.p4, met.p4).M(),                            # mHH_simple_met 
              HLL.Wlep_met_simple(lep.p4, met.p4).M(),                                                       # mWlep_met_simple
              HLL.HWW_met_simple(j3.p4,j4.p4,lep.p4,met.p4).M(),                                             # mWW_simple_met
              HLL.Wjj_simple(j3.p4, j4.p4).M(),                                                              # mWjj_simple
              HLL.comp_cosThetaS(j1.p4, j2.p4),                                                              # cosThetaS_Hbb                   
              HLL.comp_cosThetaS(j3.p4, j4.p4),                                                                          # cosThetaS_Wjj_simple     
              HLL.comp_cosThetaS(HLL.Wjj_simple(j3.p4,j4.p4),HLL.Wlep_met_simple(lep.p4,met.p4)),                        # cosThetaS_WW_simple_met  
              HLL.comp_cosThetaS(j1.p4+j2.p4, HLL.HWW_met_simple(j3.p4,j4.p4,lep.p4,met.p4)),         # cosThetaS_HH_simple_met     
              op.rng_len(jets),                                                                       # nJet
              #op.rng_len(bJets),                                                                      # nBJetMedium
              nMedBJets,
              op.deltaR(j1.p4, lep.p4),                                                               # dR_b1lep
              op.deltaR(j2.p4, lep.p4),                                                               # dR_b2lep
              HLL.lambdaConePt(lep),                                                                  # lep_conePt
              j1.pt,                                                                                  # selJet1_Hbb_pT
              j2.pt,                                                                                  # selJet2_Hbb_pT 
              HLL.MET_LD(met,jets,fakeLepColl),                                                       # met_LD
              HLL.HTfull(fakeLepColl,j1.p4,j2.p4,j3.p4,j4.p4),                                        # HT
              op.min(HLL.mT2(j1.p4, lep.p4 ,met.p4), HLL.mT2(j2.p4, lep.p4, met.p4)),                 # mT_top_3particle
              HLL.MT(lep,met)                                                                         # mT_W
          ]

    output = op.switch(event%2,model_odd(*inputs),model_even(*inputs))

    plots.append(Plot.make1D("%s_%s_nClasses"%(channel,suffix),
                             op.rng_len(output),
                             sel,
                             EquidistantBinning(10,-0.5,9.5),
                             xTitle = 'nClasses (nBDT outputs)',
                             plotopts = channelLabel))

    plots.append(Plot.make1D("%s_%s_BDTOutput"%(channel,suffix),
                             output[0],
                             sel,
                             EquidistantBinning(nBins,0.,1.),
                             xTitle = 'BDT Response',
                             plotopts = channelLabel))

    return plots

def makeSingleLeptonMachineLearningPlotsBDTmissRecoResolved(sel,fakeLepColl,lep,met,jets,bJets,lJets,j1,j2,j3,j4,suffix,channel,model_even,model_odd,nBins,event,HLL):
    plots = []
    channelLabel = SingleLeptonChannelTitleLabel(channel)

    #inputs = [op.c_float(0.)]*13
    inputs = [op.min(HLL.mT2(HLL.bJetCorrP4(j1),lep.p4 ,met.p4), HLL.mT2(HLL.bJetCorrP4(j2), lep.p4, met.p4)),# mT_top_3particle
              HLL.MT(lep,met),                                                                        # mT_W
              HLL.mindr_lep1_jet(lep,jets),                                                           # mindr_lep1_jet
              op.deltaR(j1.p4, lep.p4),                                                               # dR_b1lep
              op.deltaR(j2.p4, lep.p4),                                                               # dR_b2lep             
              op.invariant_mass(HLL.bJetCorrP4(j1), HLL.bJetCorrP4(j2)),                              # m_Hbb_regCorr 
              HLL.bJetCorrP4(j1).Pt(),                                                                # selJet1_Hbb_pT
              HLL.bJetCorrP4(j2).Pt(),                                                                # selJet2_Hbb_pT 
              op.deltaR(j3.p4, HLL.Wlep_simple(j1.p4,j2.p4,lep.p4,met)),                              # dr_Wj1_lep_simple
              op.rng_len(bJets),                                                                      # nBJetMedium
              HLL.lambdaConePt(lep),                                                                  # lep_conePt
              HLL.MET_LD(met,jets,fakeLepColl),                                                       # met_LD
              HLL.HTmiss(fakeLepColl,HLL.bJetCorrP4(j1),HLL.bJetCorrP4(j2),j3.p4)                     # HT
          ]

    output = op.switch(event%2,model_odd(*inputs),model_even(*inputs))

    plots.append(Plot.make1D("%s_%s_nClasses"%(channel,suffix),
                                 op.rng_len(output),
                                 sel,
                                 EquidistantBinning(10,-0.5,9.5),
                                 xTitle = 'nClasses (nBDT outputs)',
                                 plotopts = channelLabel))

    for i in range(1):
        plots.append(Plot.make1D("%s_%s_BDTOutput_%d"%(channel,suffix,i),
                                 output[i],
                                 sel,
                                 EquidistantBinning(nBins,0.,1.),
                                 xTitle = 'BDT output %d'%i,
                                 plotopts = channelLabel))
        
    return plots


def makeSingleLeptonMachineLearningPlotsBDTmissRecoBoosted(sel,fakeLepColl,lep,met,jets,bJets,lJets,j1,j2,j3,j4,suffix,channel,model_even,model_odd,nBins,nMedBJets,event,HLL):
    plots = []
    channelLabel = SingleLeptonChannelTitleLabel(channel)

    #inputs = [op.c_float(0.)]*13
    inputs = [op.min(HLL.mT2(j1.p4, lep.p4 ,met.p4), HLL.mT2(j2.p4, lep.p4, met.p4)),                 # mT_top_3particle
              HLL.MT(lep,met),                                                                        # mT_W
              HLL.mindr_lep1_jet(lep,jets),                                                           # mindr_lep1_jet
              op.deltaR(j1.p4, lep.p4),                                                               # dR_b1lep
              op.deltaR(j2.p4, lep.p4),                                                               # dR_b2lep             
              op.invariant_mass(j1.p4, j2.p4),                                                        # m_Hbb_regCorr 
              j1.pt,                                                                                  # selJet1_Hbb_pT
              j2.pt,                                                                                  # selJet2_Hbb_pT 
              op.deltaR(j3.p4, HLL.Wlep_simple(j1.p4,j2.p4,lep.p4,met)),                              # dr_Wj1_lep_simple
              nMedBJets,
              #op.rng_len(bJets),                                                                      # nBJetMedium
              HLL.lambdaConePt(lep),                                                                  # lep_conePt
              HLL.MET_LD(met,jets,fakeLepColl),                                                       # met_LD
              HLL.HTmiss(fakeLepColl,j1.p4,j2.p4,j3.p4)                                               # HT
          ]

    output = op.switch(event%2,model_odd(*inputs),model_even(*inputs))

    plots.append(Plot.make1D("%s_%s_nClasses"%(channel,suffix),
                                 op.rng_len(output),
                                 sel,
                                 EquidistantBinning(10,-0.5,9.5),
                                 xTitle = 'nClasses (nBDT outputs)',
                                 plotopts = channelLabel))

    for i in range(1):
        plots.append(Plot.make1D("%s_%s_BDTOutput_%d"%(channel,suffix,i),
                                 output[i],
                                 sel,
                                 EquidistantBinning(nBins,0.,1.),
                                 xTitle = 'BDT output %d'%i,
                                 plotopts = channelLabel))
        
    return plots

####################################  Machine Learning [BDT implementation :: Single Lepton] #####################################
def plotSingleLeptonBDTResponse(channel,sel,suffix,nBins,output):
    plots = []

    channelLabel = SingleLeptonChannelTitleLabel(channel)
    plots.append(Plot.make1D("%s_%s_Output"%(channel,suffix),
                             output,
                             sel,
                             EquidistantBinning(nBins,0.,1.),
                             xTitle = 'BDT Response',
                             plotopts = channelLabel))
    
    return plots    


####################################  Machine Learning [DNN implementation :: Single Lepton] #####################################
def plotSingleLeptonDNNResponse(selObjNodesDict,output,nodes,channel):
    plots = []

    channelLabel = SingleLeptonChannelTitleLabel(channel)

    for i,node in enumerate(nodes):
        suffix = selObjNodesDict[node].selName
        sel = selObjNodesDict[node].sel
        plots.append(Plot.make1D("%s_%s_DNNOutput_%s"%(channel,suffix,node),
                                 output[i],
                                 sel,
                                 EquidistantBinning(50,0.,1.),
                                 xTitle = 'DNN output %s'%node,
                                 plotopts = channelLabel))
        
    return plots    



####################################  Machine Learning [DNN implementation :: Double Lepton] #####################################
def makeDoubleLeptonMachineLearningPlots(selObjNodesDict,output,nodes,channel):
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
        
    return plots    



def makeJPAtestPlots(sel, JPAscorePerCat, JPAjetIdxsPerCat, JPAnodeList, channel):
    plots = []
    channelLabel = SingleLeptonChannelTitleLabel(channel)
    suffix = sel.selName
    plots.append(Plot.make1D("%s_%s_nBDTMaxScores"%(channel,suffix),
                             op.rng_len(JPAscorePerCat),
                             sel.sel,
                             EquidistantBinning(10, 0.0, 10.0),
                             xTitle = 'BDT output nNodes',
                             plotopts = channelLabel))
    

    return plots



def plotBDTout(channel, out, jets, sel, suffix):
    plots = []
    channelLabel = SingleLeptonChannelTitleLabel(channel)

    plots.append(Plot.make1D("%s_%s_BDTout"%(channel,suffix),
                             out,
                             sel,
                             EquidistantBinning(100, -1.0, 1.0),
                             xTitle = 'BDT output',
                             plotopts = channelLabel))
    
    return plots


def makeTestPlot(channel, var, sel, suffix):
    print(type(var))
    plots = []
    channelLabel = SingleLeptonChannelTitleLabel(channel)

    plots.append(Plot.make1D("%s_%s_BDT-MaxOut"%(channel,suffix),
                             var,
                             sel,
                             EquidistantBinning(10, -0.5, 9.5),
                             xTitle = 'BDT MaxIntVal',
                             plotopts = channelLabel))

    return plots


def makeDoubleLeptonSelectedResolvedVariables(sel,l1,l2,b1,b2,jets,met,suffix,channel,HLL):
    plots = []
    channelLabel = DoubleLeptonChannelTitleLabel(channel)
    if channel == "ElEl":
        l1conept = lambda l1 : HLL.electron_conept(l1)
        l2conept = lambda l2 : HLL.electron_conept(l2)
    elif channel == "MuMu":
        l1conept = lambda l1 : HLL.muon_conept(l1)
        l2conept = lambda l2 : HLL.muon_conept(l2)
    elif channel == "ElMu":
        l1conept = lambda l1 : HLL.electron_conept(l1)
        l2conept = lambda l2 : HLL.muon_conept(l2)
    else:
        raise RuntimeError('Could not find correct channel %s'%channel)

    plots.append(Plot.make1D("%s_%s_l1_Pt"%(channel,suffix),
                             op.switch(l1conept(l1) >= l2conept(l2) , l1.pt , l2.pt),
                             sel,
                             EquidistantBinning(60,0.,300.),
                             xTitle = "Lead lepton P_{T} [GeV]",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_b1_Pt"%(channel,suffix),
                             HLL.getCorrBp4(b1).Pt(),
                             sel,
                             EquidistantBinning(60,0.,300.),
                             xTitle = "Lead bjet P_{T} [GeV]",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_met_Pt"%(channel,suffix),
                             met.pt,
                             sel,
                             EquidistantBinning(80,0.,400.),
                             xTitle = "MET P_{T} [GeV]",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_combined_mll"%(channel,suffix),
                             op.invariant_mass(l1.p4,l2.p4),
                             sel,
                             EquidistantBinning(80,0.,400.),
                             xTitle = "M_{ll} [GeV]",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_combined_DRll"%(channel,suffix),
                             op.deltaR(l1.p4,l2.p4),
                             sel,
                             EquidistantBinning(50,0.,5.),
                             xTitle = "#Delta R_{ll}",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_combined_DRbb"%(channel,suffix),
                             op.deltaR(HLL.getCorrBp4(b1),HLL.getCorrBp4(b2)),
                             sel,
                             EquidistantBinning(50,0.,5.),
                             xTitle = "#Delta R_{bb}",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_combined_DRllbb"%(channel,suffix),
                             op.deltaR((l1.p4+l2.p4),(HLL.getCorrBp4(b1)+HLL.getCorrBp4(b2))),
                             sel,
                             EquidistantBinning(50,0.,5.),
                             xTitle = "#Delta R_{ll,bb}",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_combined_mHH"%(channel,suffix),
                             op.invariant_mass(l1.p4,l2.p4,met.p4,HLL.getCorrBp4(b1),HLL.getCorrBp4(b2)),
                             sel,
                             EquidistantBinning(100,0.,1000.),
                             xTitle = "M_{HH} [GeV]",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_combined_HT"%(channel,suffix),
                             op.rng_sum(jets, lambda j : j.pt),
                             sel,
                             EquidistantBinning(1000,0.,1000.),
                             xTitle = "H_{T} [GeV]",
                             plotopts = channelLabel))

    return plots

def makeDoubleLeptonSelectedBoostedVariables(sel,l1,l2,B,jets,met,suffix,channel,HLL):
    plots = []
    channelLabel = DoubleLeptonChannelTitleLabel(channel)
    if channel == "ElEl":
        l1conept = lambda l1 : HLL.electron_conept(l1)
        l2conept = lambda l2 : HLL.electron_conept(l2)
    elif channel == "MuMu":
        l1conept = lambda l1 : HLL.muon_conept(l1)
        l2conept = lambda l2 : HLL.muon_conept(l2)
    elif channel == "ElMu":
        l1conept = lambda l1 : HLL.electron_conept(l1)
        l2conept = lambda l2 : HLL.muon_conept(l2)
    else:
        raise RuntimeError('Could not find correct channel %s'%channel)

    plots.append(Plot.make1D("%s_%s_l1_Pt"%(channel,suffix),
                             op.switch(l1conept(l1) >= l2conept(l2) , l1.pt , l2.pt),
                             sel,
                             EquidistantBinning(50,0.,500.),
                             xTitle = "Lead lepton P_{T} [GeV]",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_B_Pt"%(channel,suffix),
                             B.pt,
                             sel,
                             EquidistantBinning(50,0.,1000.),
                             xTitle = "B P_{T} [GeV]",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_B_M"%(channel,suffix),
                             B.msoftdrop,
                             sel,
                             EquidistantBinning(30,0.,300.),
                             xTitle = "B M_{soft drop} [GeV]",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_met_Pt"%(channel,suffix),
                             met.pt,
                             sel,
                             EquidistantBinning(50,0.,500.),
                             xTitle = "MET P_{T} [GeV]",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_combined_mll"%(channel,suffix),
                             op.invariant_mass(l1.p4,l2.p4),
                             sel,
                             EquidistantBinning(40,0.,400.),
                             xTitle = "M_{ll} [GeV]",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_combined_DRll"%(channel,suffix),
                             op.deltaR(l1.p4,l2.p4),
                             sel,
                             EquidistantBinning(25,0.,5.),
                             xTitle = "#Delta R_{ll}",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_combined_DRllB"%(channel,suffix),
                             op.deltaR((l1.p4+l2.p4),B.p4),
                             sel,
                             EquidistantBinning(25,0.,5.),
                             xTitle = "#Delta R_{ll,B}",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_combined_mHH"%(channel,suffix),
                             op.invariant_mass(l1.p4,l2.p4,met.p4,B.p4),
                             sel,
                             EquidistantBinning(50,0.,1500.),
                             xTitle = "M_{HH} [GeV]",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_combined_HT"%(channel,suffix),
                             op.rng_sum(jets, lambda j : j.pt),
                             sel,
                             EquidistantBinning(100,0.,1000.),
                             xTitle = "H_{T} [GeV]",
                             plotopts = channelLabel))

    return plots

def makeDoubleLeptonHMEPlots(sel,suffix,channel,HME):
    plots = []
    channelLabel = DoubleLeptonChannelTitleLabel(channel)
    plots.append(Plot.make1D("%s_%s_HME"%(channel,suffix),
                             HME,
                             sel,
                             EquidistantBinning(200,0.,2000.),
                             xTitle = "HME [GeV]",
                             plotopts = channelLabel))

    return plots




def makeDoubleLeptonMachineLearningInputPlots(sel,suffix,channel,inputs):
    plots = []
    channelLabel = DoubleLeptonChannelTitleLabel(channel)
    for (varname,vartitle,binning),variable in inputs.items():
        plots.append(Plot.make1D("%s_%s_DNNInput_%s"%(channel,suffix,varname),
                                 variable,
                                 sel,
                                 EquidistantBinning(*binning),
                                 xTitle = vartitle,
                                 plotopts = channelLabel))
    return plots
    
def makeDoubleLeptonMachineLearningExclusiveOutputPlots(selObjNodesDict,output,nodes,channel):
    plots = []

    channelLabel = DoubleLeptonChannelTitleLabel(channel)
    plotopts = {**channelLabel}
    for i,node in enumerate(selObjNodesDict.keys()):
        if node not in nodes:
            continue
        suffix = selObjNodesDict[node].selName
        sel = selObjNodesDict[node].sel
        if node in ["GGF","VBF"]:
            plotopts.update({'blinded-range':[0.5,1]})
        plots.append(Plot.make1D("%s_%s_DNNOutput_%s"%(channel,suffix,node),
                                 output[i],
                                 sel,
                                 EquidistantBinning(400,0.,1.),
                                 #EquidistantBinning(10,0.,1.),
                                 xTitle = 'DNN output %s'%node,
                                 plotopts = plotopts))

    return plots    

def makeDoubleLeptonMachineLearningExclusiveOutputPlotsWithHME(selObjNodesDict,output,nodes,channel,HME):
    plots = []

    channelLabel = DoubleLeptonChannelTitleLabel(channel)
    plotopts = {**channelLabel}
    for i,node in enumerate(selObjNodesDict.keys()):
        if node not in nodes:
            continue
        suffix = selObjNodesDict[node].selName
        sel = selObjNodesDict[node].sel
        plots.append(Plot.make2D("%s_%s_DNNOutputvsHME_%s"%(channel,suffix,node),
                                 [output[i],HME],
                                 sel,
                                 [EquidistantBinning(400,0.,1.),EquidistantBinning(400,0,2000.)],
                                 xTitle = 'DNN output %s'%node,
                                 yTitle = 'HME',
                                 plotopts = plotopts))

    return plots    


def makeDoubleLeptonMachineLearningInclusiveOutputPlots(selObjNodesDict,output,nodes,channel):
    plots = []

    channelLabel = DoubleLeptonChannelTitleLabel(channel)
    plotopts = {**channelLabel}
    for i,node in enumerate(selObjNodesDict.keys()):
        if node not in nodes:
            continue
        suffix = selObjNodesDict[node].selName
        sel = selObjNodesDict[node].sel
        for i,plotNode in enumerate(nodes):
            if selectNode in ["GGF","VBF"]:
                plotopts.update({'blinded-range':[0.5,1]})
            plots.append(Plot.make1D("%s_%s_DNNOutput_%s"%(channel,suffix,plotNode),
                                     output[i],
                                     sel,
                                     EquidistantBinning(400,0.,1.),
                                     xTitle = 'DNN output %s'%plotNode,
                                     plotopts = plotopts))

    #channelLabel = DoubleLeptonChannelTitleLabel(channel)
    channelLabel = SingleLeptonChannelTitleLabel(channel)
    for i,node in enumerate(nodes):
        plotots = {**channelLabel}
        suffix = selObjNodesDict[node].selName
        sel = selObjNodesDict[node].sel
        if node in ["GGF","VBF"]:
            plotots.update({'blinded-range':[0.5,1]})
            # just to plot the VBF and GGF nodes
            # slide the below to the left if you want make plots for all nodes
            plots.append(Plot.make1D("%s_%s_DNNOutput_%s"%(channel,suffix,node),
                                     output[i],
                                     sel,
                                     EquidistantBinning(400,0.,1.),
                                     #EquidistantBinning(10,0.,1.),
                                     xTitle = 'DNN output %s'%node,
                                     plotopts = plotots))
    return plots    
