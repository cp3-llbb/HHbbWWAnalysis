import math
from bamboo.plots import Plot, EquidistantBinning, SummedPlot
from bamboo import treefunctions as op

from highlevelLambdasSL import *
from bamboo.plots import CutFlowReport

# TODO : remove the self

########################   Channel title   #############################
def channelTitleLabel(channel):
    if (channel == "El"):
        channel = "e^{+,-}"
    elif (channel == "Mu"):
        channel ="#mu^{+,-}"
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
def channelPlot(sel,SinlepEl,SinlepMu,suffix,channel):
    """
    Plots the number of events in each channel
    """
    plots = []
    plots.append(objectsNumberPlot(channel = channel,
                                   suffix  = suffix,
                                   objName = "SinlepEl",
                                   objCont = SinlepEl,
                                   sel     = sel,
                                   Nmax    = 5,
                                   xTitle  = "N(e^{+}/e^{-})"))
    plots.append(objectsNumberPlot(channel = channel,
                                   suffix  = suffix,
                                   objName = "SinlepMu",
                                   objCont = SinlepMu,
                                   sel     = sel,
                                   Nmax    = 5,
                                   xTitle  = "N(#mu^{+}/#mu^{-})"))
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
def makeSinleptonPlots(sel, sinlepton, suffix, channel, is_MC=False):
    """
    Make dilepton basic plots
    sel         = refine selection 
    sinlepton   = sinlepton object (El, Mu)
    suffix      = string identifying the selecton 
    channel     = string identifying the channel of the sinlepton (can be "NoChannel")
    """
    plots = []

    channelLabel = channelTitleLabel(channel)

    # PT plot #
    plots.append(Plot.make1D("%s_%s_lepton_pt"%(channel,suffix), 
                             sinlepton.pt, 
                             sel, 
                             EquidistantBinning(60,0.,300.),
                             title="Transverse momentum of the lepton (channel %s)"%channel, 
                             xTitle= "P_{T} (lepton) [GeV]",
                             plotopts = channelLabel))

    # Eta plot #
    plots.append(Plot.make1D("%s_%s_lepton_eta"%(channel,suffix), 
                             sinlepton.eta, 
                             sel, 
                             EquidistantBinning(22, -3., 3.), 
                             title="Pseudorapidity of the lepton (channel %s)"%channel, 
                             xTitle= "#eta (lepton)",
                             plotopts = channelLabel))

    # PT-eta plots #
    plots.append(Plot.make2D("%s_%s_lepton_ptVSeta"%(channel,suffix), 
                             [sinlepton.pt, sinlepton.eta],
                             sel, 
                             [EquidistantBinning(60,0.,300.),EquidistantBinning(22, -3., 3.)],
                             xTitle= "P_{T} (lepton) [GeV]",
                             yTitle= "#eta (lepton)",
                             plotopts = channelLabel))

    # Phi plot #
    plots.append(Plot.make1D("%s_%s_lepton_phi"%(channel,suffix), 
                             sinlepton.phi, 
                             sel, 
                             EquidistantBinning(20, -3.2, 3.2), 
                             title="Azimutal angle of the lepton (channel %s)"%channel, 
                             xTitle= "#phi (lepton)",
                             plotopts = channelLabel))

    # GenPartFlav (if isMC) #
    plots.append(Plot.make1D("%s_%s_lepton_genPartFlav"%(channel,suffix), 
                             sinlepton.genPartFlav if is_MC else op.c_int(-1),
                             sel, 

                             EquidistantBinning(23, -1., 22.), 
                             title="Flavour of genParticle (channel %s)"%channel, 
                             xTitle= "GenParticle flavour (lepton)",
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

##########################  JETS SEPARATE PLOTS #################################
def makeAk4JetsPlots (sel,j1,j2,j3,j4,channel,suffix,nJet,nbJet,is_MC=False):

    #print ("===============>>>Ak4Jet Plots__channel:%s__sel:%s"%(channel,suffix))
    plots = []

    channelLabel = channelTitleLabel(channel)
    '''
    #    if nbJet == 0:
    j1_name = "leadingAk4Jet"
    j2_name = "subleadingAk4Jet"
    j3_name = "thirdAk4Jet"
    j4_name = "fourthAk4Jet"

    elif nbJet == 1:
        j1_name = "leadingAk4-bJet"
        j2_name = "leadingAk4Jet"
        j3_name = "subleadingAk4Jet"
        j4_name = "thirdAk4Jet"

    elif nbJet == 2:
        j1_name = "leadingAk4-bJet"
        j2_name = "subleadingAk4-bJet"
        j3_name = "leadingAk4Jet"
        j4_name = "subleadingAk4Jet"
            
    else: raise RuntimeError ("Error in makeAk4JetsPlots args!!!")
    '''
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

    '''
    print ("nJet :%s , nbJet :%s"%(nJet,nbJet))
    print (type(j1),"__%s"%j1_name)
    print (type(j2),"__%s"%j2_name)
    print (type(j3),"__%s"%j3_name)
    print (type(j4),"__%s"%j4_name)
    '''

    # leadjet plots #
    plots.append(Plot.make1D("%s_%s_pT_%s"%(channel,suffix,j1_name),
                             j1.pt,
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
                             j2.pt,
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
                             (j1.p4+j2.p4).Pt(), 
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
                             MinDiJetDRLoose(j1,j2,j3), 
                             sel, 
                             EquidistantBinning(50,0, 5.0),
                             title="MinDeltaR between dijet (channel %s)"%channel,
                             xTitle= "MinDeltaR [Di-Jet]",
                             plotopts = channelLabel))
    
    # invariant mass plot #
    plots.append(Plot.make1D("%s_%s_InvariantMass_%s_%s"%(channel,suffix,j1_name,j2_name),
                             op.invariant_mass(j1.p4,j2.p4),
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
        # invariant mass plot #
        plots.append(Plot.make1D("%s_%s_InvariantMass_%s_%s"%(channel,suffix,j3_name,j4_name),
                                 op.invariant_mass(j3.p4,j4.p4),
                                 sel,
                                 EquidistantBinning(50, 0., 500.), 
                                 title="Invariant Mass of the dijet_[%s_%s] (channel %s)"%(j3_name,j4_name,channel), 
                                 xTitle= "InvariantMass [%s_%s]"%(j3_name, j4_name),
                                 plotopts = channelLabel))

    
    return plots

##########################  JETS (AK8) PLOT #################################
def makeAk8JetsPlots(sel,j1,j2,suffix,channel,has1fat=False,has2fat=False,has1fat1slim=False):
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

    channelLabel = channelTitleLabel(channel)

    j1_name = "Ak8Jet" if has1fat else "Ak8bJet"
    j2_name = "Ak8Jet" if has2fat else "Ak4Jet"

    # 1st-jet : FatJet #
    plots.append(Plot.make1D("%s_%s_pT_%s"%(channel,suffix,j1_name),
                             j1.p4.pt(),
                             sel,
                             EquidistantBinning(40,0.,200.),
                             title='Transverse momentum of the %s'%j1_name,
                             xTitle="P_{T} %s[GeV]"%j1_name,
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_eta_%s"%(channel,suffix,j1_name),
                             j1.p4.eta(),
                             sel,
                             EquidistantBinning(22,-3.,3.),
                             title='Pseudorapidity of the %s'%j1_name,
                             xTitle="#eta %s"%j1_name,
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

    if has2fat or has1fat1slim:
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

    return plots

#########################  High-level quantities ################################
def makeHighLevelPlotsResolved(sel,met,lep,j1,j2,j3,j4,channel,suffix,nJet,nbJet):

    #print ("===============>>>HighLevel Plots for Resolved__channel:%s__sel:%s"%(channel,suffix))
    plots = []

    channelLabel = channelTitleLabel(channel)
    '''
    #if nbJet == 0:
    j1_name = "leadingAk4Jet"
    j2_name = "subleadingAk4Jet"
    j3_name = "thirdAk4Jet"
    j4_name = "fourthAk4Jet"
    
    elif nbJet == 1:
        j1_name = "leadingAk4-bJet"
        j2_name = "leadingAk4Jet"
        j3_name = "subleadingAk4Jet"
        j4_name = "thirdAk4Jet"
    elif nbJet == 2:
        j1_name = "leadingAk4-bJet"
        j2_name = "subleadingAk4-bJet"
        j3_name = "leadingAk4Jet"
        j4_name = "subleadingAk4Jet"
    else: raise RuntimeError ("Please correct the values of nbJet!!!")
    '''

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
                             MT(lep,met),
                             sel,
                             EquidistantBinning(50,0.,1000.),
                             title='Transverse mass of lepton and MET (%s channel)'%channel,
                             xTitle="M_{T}(lep,MET) [GeV]",
                             plotopts = channelLabel))
    # lepton-MET plots #
    plots.append(Plot.make1D("%s_%s_highlevelvariable_SinglepMET_Pt"%(channel,suffix),
                             op.abs(SinglepMet_Pt(lep, met)),
                             sel,
                             EquidistantBinning(20,0.,3.2),
                             title='Reultant pT of lepton and MET (%s channel)'%channel,
                             xTitle="|Pt (lep,MET)|",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_highlevelvariable_SinglepMETdeltaPhi"%(channel,suffix),
                             op.abs(SinglepMet_dPhi(lep, met)),
                             sel,
                             EquidistantBinning(20,0.,3.2),
                             title='Azimutal angle between lepton and MET (%s channel)'%channel,
                             xTitle="|#Delta \phi (lep,MET)|",
                             plotopts = channelLabel))
    
    # jet-MET plots
    plots.append(Plot.make1D("%s_%s_highlevelvariable_%s_METdeltaPhi"%(channel,suffix,j1_name),
                             op.abs(SinglepMet_dPhi(j1, met)),
                             sel,
                             EquidistantBinning(20,0.,3.2),
                             title='Azimutal angle between %s and MET (%s channel)'%(j1_name,channel),
                             xTitle="|#Delta \phi (%s,MET)|"%j1_name,
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_highlevelvariable_%s_METdeltaPhi"%(channel,suffix,j2_name),
                             op.abs(SinglepMet_dPhi(j2, met)),
                             sel,
                             EquidistantBinning(20,0.,3.2),
                             title='Azimutal angle between %s and MET (%s channel)'%(j2_name,channel),
                             xTitle="|#Delta \phi (%s,MET)|"%j2_name,
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_highlevelvariable_%s_METdeltaPhi"%(channel,suffix,j3_name),
                             op.abs(SinglepMet_dPhi(j3, met)),
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
                             MinDR_lep3j(lep,j1,j2,j3),
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
                             MT_W1W2_lj(lep,j3,met),
                             sel,
                             EquidistantBinning(50,0.,1000.),
                             title='Transverse mass of lepton+jet and MET (%s channel)'%channel,
                             xTitle="M_{T}(lepj,MET) [GeV]",
                             plotopts = channelLabel))
    # Scalar magnitude sum #
    plots.append(Plot.make1D("%s_%s_highlevelvariable_HT2_l3jmet"%(channel,suffix),
                             HT2_l3jmet(lep,j1,j2,j3,met),
                             sel,
                             EquidistantBinning(60,0.,600.),
                             title='Di-Higgs magnitude_l3jmet (%s channel)'%channel,
                             xTitle="H_{T2} [GeV]",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_highlevelvariable_HT2R_l3jmet"%(channel,suffix),
                             HT2R_l3jmet(lep,j1,j2,j3,met),
                             sel,
                             EquidistantBinning(40,0.,1.),
                             title='Di-Higgs magnitude_Ratio_l3jmet (%s channel)'%channel,
                             xTitle="H_{T2}_Ratio [GeV]",
                             plotopts = channelLabel))
    if nJet == 4:
        plots.append(Plot.make1D("%s_%s_highlevelvariable_%s_METdeltaPhi"%(channel,suffix,j4_name),
                                 op.abs(SinglepMet_dPhi(j4, met)),
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
                                 HT2_l4jmet(lep,j1,j2,j3,j4,met),
                                 sel,
                                 EquidistantBinning(60,0.,600.),
                                 title='Di-Higgs magnitude_l4jmet (%s channel)'%channel,
                                 xTitle="H_{T2} [GeV]",
                                 plotopts = channelLabel))
        plots.append(Plot.make1D("%s_%s_highlevelvariable_HT2R_l4jmet"%(channel,suffix),
                                 HT2R_l4jmet(lep,j1,j2,j3,j4,met),
                                 sel,
                                 EquidistantBinning(40,0.,1.),
                                 title='Di-Higgs magnitude_Ratio_l4jmet (%s channel)'%channel,
                                 xTitle="H_{T2}_Ratio [GeV]",
                                 plotopts = channelLabel))
        # Transverse mass plots #
        if nbJet == 2:
            plots.append(Plot.make1D("%s_%s_highlevelvariable_MT_lep_%s_%s"%(channel,suffix,j3_name,j4_name),
                                     MT_W1W2_ljj(lep,j3,j4,met),
                                     sel,
                                     EquidistantBinning(50,0.,1000.),
                                     title='Transverse mass of lepton+dijet and MET (%s channel)'%channel,
                                     xTitle="M_{T}(lepjj,MET) [GeV]",
                                     plotopts = channelLabel))

    return plots 

# Make HighLevel plots for semi boosted and boosted categories
def makeHighLevelPlotsBoosted(sel,met,lep,j1,j2,channel,suffix,bothAreFat=False):
    #print ("===============>>>HighLevel Plots for Boosted__channel:%s__sel:%s"%(channel,suffix))
    if bothAreFat:
        j1_name   = "Ak8BJet"
        j2_name   = "Ak8Jet"
        
    else:
        j1_name   = "Ak8BJet"
        j2_name   = "Ak4Jet"
    
    plots = []
    channelLabel = channelTitleLabel(channel)

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
                             op.abs(SinglepMet_Pt(lep, met)),
                             sel,
                             EquidistantBinning(20,0.,3.2),
                             title='Reultant pT of lepton and MET (%s channel)'%channel,
                             xTitle="|Pt (lep,MET)|",
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_highlevelvariable_SinglepMETdeltaPhi"%(channel,suffix),
                             op.abs(SinglepMet_dPhi(lep, met)),
                             sel,
                             EquidistantBinning(20,0.,3.2),
                             title='Azimutal angle between lepton and MET (%s channel)'%channel,
                             xTitle="|#Delta \phi (lep,MET)|",
                             plotopts = channelLabel))
    # jet-MET plots
    plots.append(Plot.make1D("%s_%s_highlevelvariable_%s_METdeltaPhi"%(channel,suffix,j1_name),
                             op.abs(SinglepMet_dPhi(j1, met)),
                             sel,
                             EquidistantBinning(20,0.,3.2),
                             title='Azimutal angle between %s and MET (%s channel)'%(j1_name,channel),
                             xTitle="|#Delta \phi (%s,MET)|"%j1_name,
                             plotopts = channelLabel))
    plots.append(Plot.make1D("%s_%s_highlevelvariable_%s_METdeltaPhi"%(channel,suffix,j2_name),
                             op.abs(SinglepMet_dPhi(j2, met)),
                             sel,
                             EquidistantBinning(20,0.,3.2),
                             title='Azimutal angle between %s and MET (%s channel)'%(j2_name,channel),
                             xTitle="|#Delta \phi (%s,MET)|"%j2_name,
                             plotopts = channelLabel))


    # Transverse mass plots #
    # mT for Single lepton #
    plots.append(Plot.make1D("%s_%s_highlevelvariable_MT"%(channel,suffix),
                             MT(lep,met),
                             sel,
                             EquidistantBinning(50,0.,1000.),
                             title='Transverse mass of lepton and MET (%s channel)'%channel,
                             xTitle="M_{T}(lep,MET) [GeV]",
                             plotopts = channelLabel))
    
    return plots

