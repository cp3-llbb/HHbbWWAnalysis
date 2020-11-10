import os
import sys
from copy import copy, deepcopy
from bamboo import treefunctions as op
from bamboo.plots import SelectionWithDataDriven

#===============================================================================================#
#                                  SelectionObject class                                        #
#===============================================================================================#
class SelectionObject:
    """
        Helper class to hold a selection including :
        sel         : refine object
        selName     : name of the selection used in the refine
        yieldTitle  : Title for the yield table in LateX
        NB : copy() might be useful to be able to branch out the selection
    """
    def __init__(self,sel,selName,yieldTitle):
        self.sel        = sel
        self.selName    = selName
        self.yieldTitle = yieldTitle
    
    def refine(self,cut=None,weight=None,autoSyst=True):
        """
        Overload the refine function to avoid repeating the selName arg 
        """
        self.sel = self.sel.refine(name    = self.selName,
                                   cut     = cut,
                                   weight  = weight,
                                   autoSyst= autoSyst)

    def create(self, ddSuffix, cut=None, weight=None, ddCut=None, ddWeight=None, autoSyst=True, enable=True):
        """
        Overload the create function for data-driven estimation to avoid repeating the sel and selName arg 
        """
        self.sel = SelectionWithDataDriven.create(parent   = self.sel,
                                                  name     = self.selName,
                                                  ddSuffix = ddSuffix,
                                                  cut      = cut,
                                                  weight   = weight,
                                                  ddCut    = ddCut,
                                                  ddWeight = ddWeight,
                                                  autoSyst = autoSyst,
                                                  enable   = enable)

    def makeYield(self,yieldObject):
        """
        Record an entry in the yield table 
        """
        yieldObject.addYield(sel    = self.sel,
                             name   = self.selName,
                             title  = self.yieldTitle)

#===============================================================================================#
#                                        Selections                                             #
#===============================================================================================#
def makeSingleLeptonSelection(self,baseSel,plot_yield=False,use_dd=True): 
    """
    Produces the requested lepton selection (encapsulated in SelectionObject class objects)
    Will produce a dict :
        - key = level required (Preselected, Fakeable, Tight and/or FakeExtrapolation)
        - value = list of SelectionObject class objects per channel [ElEl,MuMu,ElMu] 
    We start by leptons so no need to pass selObject
    Code organized such that selections are not repeated and hopefully optimzed the RooDataFrame
    """
    # Select level #
    sel_level = []   # What lepton type selection to perform (e.g. Tight needs Fakeable and Preselected)
    lepton_level = ["Preselected","Fakeable","Tight","FakeExtrapolation"]
    # What lepton type selection to save for plotting
    save_level = [arg for (arg,boolean) in self.args.__dict__.items() if arg in lepton_level and boolean]  
        # Make list of requested lepton selection in the arguments
    if len(save_level) == 0: save_level = lepton_level # If nothing asked, all will be done

    if "Tight" in save_level:
        sel_level.extend(["Preselected","Fakeable","Tight"]) # Tight needs Fakeable and Preselected
    if "FakeExtrapolation" in save_level:
        sel_level.extend(["Preselected","Fakeable","FakeExtrapolation"]) # FakeExtrapolation needs Fakeable and Preselected
    if "Fakeable" in save_level:
        sel_level.extend(["Preselected","Fakeable"]) # Fakeable needs Preselected
    if "Preselected" in save_level:
        sel_level.extend(["Preselected"]) 
            # TODO : find cleaner-> if args.Tight And args.FakeExtrapolation : will be redundancies (not a problem though)
    #--- Lambdas ---#
    # Fakeable lambdas #
    lambdaLowPtCut  = lambda lepColl: lepColl[0].pt > 25 # leading above 15 GeV
    lambdaHighPtCut = lambda lepColl: lepColl[0].pt > 32 # leading above 25 GeV

    #Z-nominal mass
    Zmass = 91.1876
    # Mll cut lambdas #
    lambda_lowMllCut = lambda dileptons: op.NOT(op.rng_any(dileptons, lambda dilep : op.invariant_mass(dilep[0].p4, dilep[1].p4) < 12.))
    # If any dilepton preselected below 12GeV, returns False
    lambda_outZ      = lambda dileptons: op.NOT(op.rng_any(dileptons, lambda dilep : op.in_range(Zmass - 10., op.invariant_mass(dilep[0].p4, dilep[1].p4), Zmass + 10.)))
    # If any dilepton preselected within Z peak, returns False
    lambda_inZ       = lambda dilep : op.in_range(Zmass - 10., op.invariant_mass(dilep[0].p4, dilep[1].p4), Zmass + 10.)


    if self.is_MC:
        # Loose SF #
        ElLooseSF = lambda lepColl : [op.defineOnFirstUse(self.ttH_singleElectron_trigSF(lepColl[0]))] + self.lambda_ElectronLooseSF(lepColl[0])
        MuLooseSF = lambda lepColl : [op.defineOnFirstUse(self.ttH_singleMuon_trigSF(lepColl[0]))] + self.lambda_MuonLooseSF(lepColl[0])
        # Tight SF #
        ElTightSF = lambda lepColl : self.lambda_ElectronTightSF(lepColl[0])
        MuTightSF = lambda lepColl : self.lambda_MuonTightSF(lepColl[0])
    else:
        ElLooseSF = lambda lepColl : []
        MuLooseSF = lambda lepColl : []
        ElTightSF = lambda lepColl : []
        MuTightSF = lambda lepColl : []

    
    #--- Preselection ---#
    selectionDict = {}

    if "Preselected" in sel_level:
        print('... PreSelected Lepton Selection ===>')
        ElPreSelObject = SelectionObject(sel          = baseSel,
                                         selName      = "HasElPreselected",
                                         yieldTitle   = "Preselected lepton (channel $e^+/e^-$)")
        MuPreSelObject = SelectionObject(sel          = baseSel,
                                         selName      = "HasMuPreselected",
                                         yieldTitle   = "Preselected lepton (channel $\mu^+/\mu^-$)")

        # PreSelectionObject # [PreSelIncl ????????]
        ElPreSelObject.refine(cut = [op.rng_len(self.electronsPreSel) >= 1])
        MuPreSelObject.refine(cut = [op.rng_len(self.muonsPreSel) >= 1])
        # Yield #
        if plot_yield:
            ElPreSelObject.makeYield(self.yieldPlots)
            MuPreSelObject.makeYield(self.yieldPlots)
        # Record if actual argument #
        if "Preselected" in save_level:
            selectionDict["Preselected"] = [ElPreSelObject,MuPreSelObject]

        #--- Fakeable ---#
        if "Fakeable" in sel_level:
            print('... Fakeable Lepton Selection ===>')
            ElFakeSelObject = SelectionObject(sel         = ElPreSelObject.sel,
                                              selName     = "HasElFakeable",
                                              yieldTitle  = "Fakeable lepton (channel $e^+/e^-$)")
            MuFakeSelObject = SelectionObject(sel         = MuPreSelObject.sel,
                                              selName     = "HasMuFakeable",
                                              yieldTitle  = "Fakeable lepton (channel $\mu^+/\mu^-$)")

            # Selection : at least one fakeable dilepton #
            ElFakeSelObject.refine(cut    = [op.rng_len(self.leadElectronFakeSel) == 1, lambdaHighPtCut(self.leadElectronFakeSel)],
                                   weight = ElLooseSF(self.leadElectronFakeSel))
            MuFakeSelObject.refine(cut    = [op.rng_len(self.leadMuonFakeSel) == 1, lambdaLowPtCut(self.leadMuonFakeSel)],
                                   weight = MuLooseSF(self.leadMuonFakeSel))
            
            # Yield #
            if plot_yield:
                ElFakeSelObject.makeYield(self.yieldPlots)
                MuFakeSelObject.makeYield(self.yieldPlots)


            # Selection : Mll cut + Z peak exclusion in **preselected** leptons # 
            ElFakeSelObject.selName += "PreMllCut"
            MuFakeSelObject.selName += "PreMllCut"
           
            ElFakeSelObject.yieldTitle += " + $M_{ll}>12$"
            MuFakeSelObject.yieldTitle += " + $M_{ll}>12$"
            
            mllCut = [lambda_lowMllCut(self.ElElDileptonPreSel),lambda_lowMllCut(self.MuMuDileptonPreSel),lambda_lowMllCut(self.ElMuDileptonPreSel)]

            print ('...... LowMass resonance veto')
            ElEl_MllCut = copy(mllCut)
            MuMu_MllCut = copy(mllCut)
            ElMu_MllCut = copy(mllCut)

            if not self.args.NoZVeto: # No Mz cut in OS leptons
                # All the preselected leptons outside Z peak region within 10 GeV
                print ('...... Z veto')
                ElFakeSelObject.selName += "OutZ"
                MuFakeSelObject.selName += "OutZ"

                ElFakeSelObject.yieldTitle += " + Z Veto $|M_{ll}^{Preselected SF}-M_Z|>10$"
                MuFakeSelObject.yieldTitle += " + Z Veto $|M_{ll}^{Preselected SF}-M_Z|>10$"

                outZCut = [lambda_outZ(self.OSElElDileptonPreSel),lambda_outZ(self.OSMuMuDileptonPreSel)]

                ElEl_MllCut.extend(outZCut)
                MuMu_MllCut.extend(outZCut)
                ElMu_MllCut.extend(outZCut)

            ElFakeSelObject.refine(cut = ElEl_MllCut)
            MuFakeSelObject.refine(cut = MuMu_MllCut)

            # Yield #
            if plot_yield:
                ElFakeSelObject.makeYield(self.yieldPlots)
                MuFakeSelObject.makeYield(self.yieldPlots)

                
            if not self.args.NoTauVeto: # nTau (isolated from fakeable leptons) = 0
                # All the preselected leptons outside Z peak region within 10 GeV
                print ('...... Tau veto')
                ElFakeSelObject.selName += "noTau"
                MuFakeSelObject.selName += "noTau"

                ElFakeSelObject.yieldTitle += " + Tau Veto"
                MuFakeSelObject.yieldTitle += " + Tau Veto"

                ElFakeSelObject.refine(cut = [op.rng_len(self.tauCleanSel) == 0])
                MuFakeSelObject.refine(cut = [op.rng_len(self.tauCleanSel) == 0])
                
                # Yield #
                if plot_yield:
                    ElFakeSelObject.makeYield(self.yieldPlots)
                    MuFakeSelObject.makeYield(self.yieldPlots)

            # Record if actual argument #
            if "Fakeable" in save_level:
                selectionDict["Fakeable"] = [ElFakeSelObject,MuFakeSelObject]
                
            #--- Tight ---#
            if "Tight" in sel_level:
                print ('... Tight Lepton Selection===>')
                ElTightSelObject = SelectionObject(sel          = ElFakeSelObject.sel,
                                                   selName      = "HasElTight",
                                                   yieldTitle   = "Tight lepton (channel $e^+/e^-$)")
                MuTightSelObject = SelectionObject(sel          = MuFakeSelObject.sel,
                                                   selName      = "HasMuTight",
                                                   yieldTitle   = "Tight lepton (channel $\mu^+/\mu^-$)")
                
                # Selection : at least one tight dilepton #
                if use_dd:
                    print ('...... Estimating fake contribution (DataDriven)')
                    ElTightSelObject.create(ddSuffix  = "FakeExtrapolation",
                                            cut       = [op.rng_len(self.leadElectronTightSel) == 1],
                                            weight    = ElTightSF(self.leadElectronTightSel),
                                            ddCut     = [op.rng_len(self.leadElectronFakeExtrapolationSel) == 1],
                                            ddWeight  = [self.ElFakeFactor(self.leadElectronFakeExtrapolationSel[0])]+ElTightSF(self.leadElectronFakeExtrapolationSel[0]),
                                            enable    = "FakeExtrapolation" in self.datadrivenContributions and 
                                            self.datadrivenContributions["FakeExtrapolation"].usesSample(self.sample, self.sampleCfg))
                    MuTightSelObject.create(ddSuffix  = "FakeExtrapolation",
                                            cut       = [op.rng_len(self.leadMuonTightSel) == 1],
                                            weight    = MuTightSF(self.leadMuonTightSel),
                                            ddCut     = [op.rng_len(self.leadMuonFakeExtrapolationSel) == 1],
                                            ddWeight  = [self.MuFakeFactor(self.leadMuonFakeExtrapolationSel[0])]+MuTightSF(self.leadMuonFakeExtrapolationSel[0]),
                                            enable    = "FakeExtrapolation" in self.datadrivenContributions 
                                            and self.datadrivenContributions["FakeExtrapolation"].usesSample(self.sample, self.sampleCfg))

                else:
                    ElTightSelObject.refine(cut    = [op.rng_len(self.leadElectronTightSel) == 1],
                                            weight = ElTightSF(self.leadElectronTightSel))
                    MuTightSelObject.refine(cut    = [op.rng_len(self.leadMuonTightSel) == 1],
                                            weight = MuTightSF(self.leadMuonTightSel))

                # Yield #
                if plot_yield:
                    ElTightSelObject.makeYield(self.yieldPlots)
                    MuTightSelObject.makeYield(self.yieldPlots)
                
                # Record if actual argument #
                if "Tight" in save_level:
                    selectionDict["Tight"] = [ElTightSelObject,MuTightSelObject]


            #--- Fake extrapolated ---#
            if "FakeExtrapolation" in sel_level:
                print ('... FakeExtrapolated Lepton Selection===>')
                ElFakeExtrapolationSelObject = SelectionObject(sel          = ElFakeSelObject.sel,
                                                               selName      = "HasElFakeExtrapolation",
                                                               yieldTitle   = "Fake extrapolated lepton (channel $e^+/e^-$)")
                MuFakeExtrapolationSelObject = SelectionObject(sel          = MuFakeSelObject.sel,
                                                               selName      = "HasMuFakeExtrapolation",
                                                               yieldTitle   = "Fake extrapolated lepton (channel $\mu^+/\mu^-$)")
                
                # Selection : at least one tight dilepton #
                ElFakeExtrapolationSelObject.refine(cut    = [op.rng_len(self.leadElectronFakeExtrapolationSel) == 1],
                                                    weight = ElTightSF(self.leadElectronFakeExtrapolationSel[0]))
                MuFakeExtrapolationSelObject.refine(cut    = [op.rng_len(self.leadMuonFakeExtrapolationSel) == 1],
                                                    weight = MuTightSF(self.leadMuonFakeExtrapolationSel[0]))
                
                # Yield #
                if plot_yield:
                    ElFakeExtrapolationSelObject.makeYield(self.yieldPlots)
                    MuFakeExtrapolationSelObject.makeYield(self.yieldPlots)
                    
                # Record if actual argument #
                if "FakeExtrapolation" in save_level:
                    selectionDict["FakeExtrapolation"] = [ElFakeExtrapolationSelObject,MuFakeExtrapolationSelObject]
    # Return # 
    return selectionDict

def makeCoarseResolvedSelection(self,selObject,nJet,copy_sel=False,plot_yield=False):
    """
    Produces the Ak4 jet selection
    inputs :
        - selObject : SelectionObject class objects (contains bamboo selection, its name and the yield title)
        - copy_sel : wether to copy and return a new selection object built on the one provided, if not will modify it
    Careful : if copy_sel is False, the selObject will be modified
    Selection :
         No ak8-bTagged Jets : To make exclusive from Boosted category
         TightResolvedSelection (nAk4 Jets ≥ 4) : Full Reco
         LooseResolvedSelection (nAk4 Jets < 4) : Missed
    """

    if copy_sel:
        selObject = copy(selObject)
    if nJet == 3: 
        selObject.selName += "Ak4JetsLoose"
        selObject.yieldTitle += " + Ak4 Jets Loose = 3"
        selObject.refine(cut  = [op.rng_len(self.ak4Jets) == nJet, 
                                 op.rng_len(self.ak8BJets) == 0])
    elif nJet == 4:
        selObject.selName += "Ak4JetsTight"
        selObject.yieldTitle += " + Ak4 Jets Tight $\geq 4$"
        selObject.refine(cut  = [op.rng_len(self.ak4Jets) >= nJet,
                                 op.rng_len(self.ak8BJets) == 0])
    if plot_yield:
        selObject.makeYield(self.yieldPlots)
    if copy_sel :
        return selObject 


def makeCoarseBoostedSelection(self,selObject,copy_sel=False,plot_yield=False):
    """
    Produces the Ak8b jet selection
    inputs :
        - selObject : SelectionObject class objects (contains bamboo selection, its name and the yield title)
        - copy_sel : wether to copy and return a new selection object built on the one provided, if not will modify it
    Careful : if copy_sel is False, the selObject will be modified
    Selection : at least one Ak8b jet
    """
    if copy_sel:
        selObject = copy(selObject)
    selObject.selName += "Ak8BJets"
    selObject.yieldTitle += " + Ak8 bJets $\geq 1$"
    selObject.refine(cut=[op.rng_len(self.ak8BJets) >= 1])
    if plot_yield:
        selObject.makeYield(self.yieldPlots)
    if copy_sel:
        return selObject


def makeExclusiveLooseResolvedJetComboSelection(self,selObject,nbJet,copy_sel=False,plot_yield=False):
    if copy_sel:
        selObject = copy(selObject)

    if nbJet == 0:
        AppliedSF = None
        selObject.selName += "ExclusiveResolved0b3j"
        selObject.yieldTitle += " + Exclusive Resolved (nbJet = 0, nAk4LightJet = 3)"
        selObject.refine(cut    = [op.rng_len(self.ak4BJets) == 0],
                         weight = AppliedSF)
        
    elif nbJet == 1:
        AppliedSF = None
        selObject.selName += "ExclusiveResolved1b2j"
        selObject.yieldTitle += " + Exclusive Resolved (nbJet = 1, nAk4LightJet = 2)"
        selObject.refine(cut    = [op.rng_len(self.ak4BJets) == 1],
                         weight = AppliedSF)

    elif nbJet == 2:
        AppliedSF = None
        selObject.selName += "ExclusiveResolved2b1j"
        selObject.yieldTitle += " + Exclusive Resolved (nbJet = 2, nAk4LightJet = 1)"
        selObject.refine(cut    = [op.rng_len(self.ak4BJets) == 2],
                         weight = None)
    else: raise RuntimeError ("Error in Jet Selection!!!")
    
    if plot_yield:
        selObject.makeYield(self.yieldPlots)
    if copy_sel:
        return selObject

def makeExclusiveTightResolvedJetComboSelection(self,selObject,nbJet,copy_sel=False,plot_yield=False):
    if copy_sel:
        selObject = copy(selObject)

    if nbJet == 0:
        AppliedSF = None
        selObject.selName += "ExclusiveResolved0b4j"
        selObject.yieldTitle += " + Exclusive Resolved (nbJet = 0, nAk4LightJets $\geq 4$)"
        selObject.refine(cut    = [op.rng_len(self.ak4BJets) == 0],
                         weight = AppliedSF)

    elif nbJet == 1:
        AppliedSF = None
        selObject.selName += "ExclusiveResolved1b3j"
        selObject.yieldTitle += " + Exclusive Resolved (nbjet = 1, nAk4LightJets $\geq 3$)"
        selObject.refine(cut    = [op.rng_len(self.ak4BJets) == 1],
                         weight = AppliedSF)

    elif nbJet == 2:
        AppliedSF = None
        selObject.selName += "ExclusiveResolved2b2j"
        selObject.yieldTitle += " + Exclusive Resolved (nbjet $\geq 2$, nAk4LightJets $\geq 2$)"
        selObject.refine(cut   = [op.rng_len(self.ak4BJets) >= 2,
                                  op.rng_len(self.ak4LightJetsByBtagScore) >= 2],
                        weight = AppliedSF)

    else: raise RuntimeError ("Error in Jet Selection!!!")

    if plot_yield:
        selObject.makeYield(self.yieldPlots)

    if copy_sel:
        return selObject


def makeSemiBoostedHbbSelection(self,selObject,nNonb,copy_sel=False,plot_yield=False):
    """
    Produces the inclusive boosted selection
    inputs :
        - selObject : SelectionObject class objects (contains bamboo selection, its name and the yield title)
        - copy_sel : wether to copy and return a new selection object built on the one provided, if not will modify it
    Careful : if copy_sel is False, the selObject will be modified
    Selection : At least one btagged Ak8 jet
    """
    if copy_sel:
        selObject = copy(selObject)
    if nNonb == 1:
        AppliedSF = None # TODO: correct at v7
        selObject.selName += "HbbBoostedWtoJ"
        selObject.yieldTitle += " + HbbBoostedWtoJ"
        selObject.refine(cut    = [op.rng_len(self.ak4JetsCleanedFromAk8b) == nNonb],
                         weight = AppliedSF)

    elif nNonb == 2:
        AppliedSF = None # TODO: correct at v7
        selObject.selName += "HbbBoostedWtoJJ"
        selObject.yieldTitle += " + HbbBoostedWtoJJ"
        selObject.refine(cut    = [op.rng_len(self.ak4JetsCleanedFromAk8b) >= nNonb],
                         weight = AppliedSF)

    if plot_yield:
        selObject.makeYield(self.yieldPlots)
    if copy_sel:
        return selObject



def makeDoubleLeptonSelection(self,baseSel,plot_yield=False,use_dd=True): 
    """
    Produces the requested lepton selection (encapsulated in SelectionObject class objects)
    Will produce a dict :
        - key = level required (Preselected, Fakeable, Tight and/or FakeExtrapolation)
        - value = list of SelectionObject class objects per channel [ElEl,MuMu,ElMu] 
    We start by leptons so no need to pass selObject
    Code organized such that selections are not repeated and hopefully optimized the RooDataFrame
    """
    # Select level #
    sel_level = []   # What lepton type selection to perform (e.g. Tight needs Fakeable and Preselected)
    lepton_level = ["Preselected","Fakeable","Tight","FakeExtrapolation"]
    # What lepton type selection to save for plotting
    save_level = [arg for (arg,boolean) in self.args.__dict__.items() if arg in lepton_level and boolean]  
        # Make list of requested lepton selection in the arguments
    if len(save_level) == 0: save_level = lepton_level # If nothing asked, all will be done

    if "Tight" in save_level:
        sel_level.extend(["Preselected","Fakeable","Tight"]) # Tight needs Fakeable and Preselected
    if "FakeExtrapolation" in save_level:
        sel_level.extend(["Preselected","Fakeable","FakeExtrapolation"]) # FakeExtrapolation needs Fakeable and Preselected
    if "Fakeable" in save_level:
        sel_level.extend(["Preselected","Fakeable"]) # Fakeable needs Preselected
    if "Preselected" in save_level:
        sel_level.extend(["Preselected"]) 
            # TODO : find cleaner-> if args.Tight And args.FakeExtrapolation : will be redundancies (not a problem though)
    #--- Lambdas ---#
    # Fakeable lambdas #
    lambdaOSDilepton = lambda dilep : dilep[0].charge != dilep[1].charge
    if self.args.POGID:
        lambdaLowPtCut = lambda dilep: op.AND(dilep[0].pt > 15 , dilep[1].pt > 15) # subleading above 15 GeV
        lambdaLeadingPtCut = lambda dilep: op.OR(dilep[0].pt > 25 , dilep[1].pt > 25) # leading above 25 GeV
    if self.args.TTHIDLoose or self.args.TTHIDTight:
        lambdaLowPtCutElEl = lambda dilep: op.AND( self.electron_conept[dilep[0].idx] > 15 , self.electron_conept[dilep[1].idx] > 15) # subleading above 15 GeV
        lambdaLowPtCutMuMu = lambda dilep: op.AND( self.muon_conept[dilep[0].idx] > 15 , self.muon_conept[dilep[1].idx] > 15) # subleading above 15 GeV
        lambdaLowPtCutElMu = lambda dilep: op.AND( self.electron_conept[dilep[0].idx] > 15 , self.muon_conept[dilep[1].idx] > 15) # subleading above 15 GeV
        lambdaLeadingPtCutElEl = lambda dilep: op.OR( self.electron_conept[dilep[0].idx] > 25 , self.electron_conept[dilep[1].idx] > 25) # leading above 25 GeV
        lambdaLeadingPtCutMuMu = lambda dilep: op.OR( self.muon_conept[dilep[0].idx] > 25 , self.muon_conept[dilep[1].idx] > 25) # leading above 25 GeV
        lambdaLeadingPtCutElMu = lambda dilep: op.OR( self.electron_conept[dilep[0].idx] > 25 , self.muon_conept[dilep[1].idx] > 25) # leading above 25 GeV
        # Need to split by flavour because conept definition needs it

    # Mll cut lambdas #
    ZMass = 91.1876
    lambda_lowMllCut    = lambda dileptons: op.NOT(op.rng_any(dileptons, lambda dilep : op.invariant_mass(dilep[0].p4, dilep[1].p4)<12.))
        # If any dilepton below 12GeV, returns False
    lambda_outZ         = lambda dileptons: op.NOT(op.rng_any(dileptons, lambda dilep : op.in_range(ZMass-10.,op.invariant_mass(dilep[0].p4, dilep[1].p4),ZMass+10.)))
        # If any dilepton within Z peak, returns False
    lambda_inZ          = lambda dilep : op.in_range(ZMass-10.,op.invariant_mass(dilep[0].p4, dilep[1].p4),ZMass+10.)
    lambda_highMllCut   = lambda dilep : op.invariant_mass(dilep[0].p4, dilep[1].p4) < 76.
    
    #--- Preselection ---#
    selectionDict = {}


    if self.args.POGID:
       # SF #
        if self.is_MC:
            ElElTightSF = lambda dilep : [op.defineOnFirstUse(self.lambda_ttH_doubleElectron_trigSF(dilep)),self.elPOGTight(dilep[0]),self.elPOGTight(dilep[1])]
            MuMuTightSF = lambda dilep : [op.defineOnFirstUse(self.lambda_ttH_doubleMuon_trigSF(dilep)),self.muPOGTightID(dilep[0]),self.muPOGTightISO(dilep[0]),self.muPOGTightID(dilep[1]),self.muPOGTightISO(dilep[1])]
            ElMuTightSF = lambda dilep : [op.defineOnFirstUse(self.lambda_ttH_electronMuon_trigSF(dilep)),self.elPOGTight(dilep[0]),self.muPOGTightID(dilep[1]),self.muPOGTightISO(dilep[1])]
        else:
            ElElTightSF = lambda dilep : [] 
            MuMuTightSF = lambda dilep : []
            ElMuTightSF = lambda dilep : []

        # Sel #
        ElElTightSelObject = SelectionObject(sel          = baseSel,
                                             selName      = "HasElElTight",
                                             yieldTitle   = "OS tight lepton pair (channel $e^+e^-$)")
        MuMuTightSelObject = SelectionObject(sel          = baseSel,
                                             selName      = "HasMuMuTight",
                                             yieldTitle   = "OS tight lepton pair (channel $\mu^+\mu^-$)")
        ElMuTightSelObject = SelectionObject(sel          = baseSel,
                                             selName      = "HasElMuTight",
                                             yieldTitle   = "OS tight lepton pair (channel $e^{\pm}\mu^{\mp}$)")
        # Tight selection #
        ElElTightSelObject.refine(cut    = [op.rng_len(self.ElElDileptonTightSel) >= 1,
                                            lambdaOSDilepton(self.ElElDileptonTightSel[0]),
                                            lambdaLowPtCut(self.ElElDileptonTightSel[0]),
                                            lambdaLeadingPtCut(self.ElElDileptonTightSel[0])],
                                  weight =  ElElTightSF(self.ElElDileptonTightSel[0]))
        MuMuTightSelObject.refine(cut    = [op.rng_len(self.MuMuDileptonTightSel) >= 1,
                                            lambdaOSDilepton(self.MuMuDileptonTightSel[0]),
                                            lambdaLowPtCut(self.MuMuDileptonTightSel[0]),
                                            lambdaLeadingPtCut(self.MuMuDileptonTightSel[0])],
                                  weight =  MuMuTightSF(self.MuMuDileptonTightSel[0]))
        ElMuTightSelObject.refine(cut    = [op.rng_len(self.ElMuDileptonTightSel) >= 1,
                                            lambdaOSDilepton(self.ElMuDileptonTightSel[0]),
                                            lambdaLowPtCut(self.ElMuDileptonTightSel[0]),
                                            lambdaLeadingPtCut(self.ElMuDileptonTightSel[0])],
                                  weight =  ElMuTightSF(self.ElMuDileptonTightSel[0]))

        # OS and PT cuts #

        ElElTightSelObject.selName += "OSPtCuts"
        MuMuTightSelObject.selName += "OSPtCuts"
        ElMuTightSelObject.selName += "OSPtCuts"
        ElElTightSelObject.yieldTitle += "OS + $P_T$ cuts"
        ElElTightSelObject.yieldTitle += "OS + $P_T$ cuts"
        ElMuTightSelObject.yieldTitle += "OS + $P_T$ cuts"

        # Yield #
        if plot_yield:
            ElElTightSelObject.makeYield(self.yieldPlots)
            MuMuTightSelObject.makeYield(self.yieldPlots)
            ElMuTightSelObject.makeYield(self.yieldPlots)

        # Mll cut #
        mllCut = [lambda_lowMllCut(self.ElElTightDileptonPairs),lambda_lowMllCut(self.MuMuTightDileptonPairs),lambda_lowMllCut(self.ElMuTightDileptonPairs)]
        ElElMllCut = copy(mllCut)
        MuMuMllCut = copy(mllCut)
        ElMuMllCut = copy(mllCut)

        ElElTightSelObject.selName += "LowMllCut"
        MuMuTightSelObject.selName += "LowMllCut"
        ElMuTightSelObject.selName += "LowMllCut"
        ElElTightSelObject.yieldTitle += " + $M_{ll} > 12 GeV$" 
        MuMuTightSelObject.yieldTitle += " + $M_{ll} > 12 GeV$"
        ElMuTightSelObject.yieldTitle += " + $M_{ll} > 12 GeV$"

        if not self.args.NoZVeto: # No Mz cut in OF leptons
            outZCut = [lambda_outZ(self.ElElTightDileptonOSPairs),lambda_outZ(self.MuMuTightDileptonOSPairs)]
            ElElMllCut.extend(outZCut)
            MuMuMllCut.extend(outZCut)
            ElMuMllCut.extend(outZCut)
            ElElTightSelObject.selName += "OutZ"
            MuMuTightSelObject.selName += "OutZ"
            ElMuTightSelObject.selName += "OutZ"
            ElElTightSelObject.yieldTitle += " + Z veto $|M_{ll}-M_Z|>10 GeV$" 
            MuMuTightSelObject.yieldTitle += " + Z veto $|M_{ll}-M_Z|>10 GeV$"
            ElMuTightSelObject.yieldTitle += " + Z veto $|M_{ll}-M_Z|>10 GeV$"

        if self.args.ZPeak:
            ElElTightSelObject.selName += "ZPeak"
            MuMuTightSelObject.selName += "ZPeak"
            ElElTightSelObject.yieldTitle+= " + Z peak $|M_{ll}-M_Z| < 10 GeV$" 
            MuMuTightSelObject.yieldTitle+= " + Z peak $|M_{ll}-M_Z| < 10 GeV$"
            ElElMllCut.append(lambda_inZ(self.ElElDileptonTightSel[0]))
            MuMuMllCut.append(lambda_inZ(self.MuMuDileptonTightSel[0]))
            
        ElElTightSelObject.refine(cut = ElElMllCut)
        MuMuTightSelObject.refine(cut = MuMuMllCut)
        ElMuTightSelObject.refine(cut = ElMuMllCut)

        # Yield #
        if plot_yield:
            ElElTightSelObject.makeYield(self.yieldPlots)
            MuMuTightSelObject.makeYield(self.yieldPlots)
            ElMuTightSelObject.makeYield(self.yieldPlots)
            
        # Record if actual argument #
        if "Tight" in save_level:
            selectionDict["Tight"] = [ElElTightSelObject,MuMuTightSelObject,ElMuTightSelObject]



    if self.args.TTHIDLoose or self.args.TTHIDTight:
        # SF #
        if self.is_MC:
            ElElLooseSF = lambda dilep : [op.defineOnFirstUse(self.lambda_ttH_doubleElectron_trigSF(dilep))] + self.lambda_ElectronLooseSF(dilep[0]) + self.lambda_ElectronLooseSF(dilep[1])
            MuMuLooseSF = lambda dilep : [op.defineOnFirstUse(self.lambda_ttH_doubleMuon_trigSF(dilep))] + self.lambda_MuonLooseSF(dilep[0]) + self.lambda_MuonLooseSF(dilep[1])
            ElMuLooseSF = lambda dilep : [op.defineOnFirstUse(self.lambda_ttH_electronMuon_trigSF(dilep))] + self.lambda_ElectronLooseSF(dilep[0]) + self.lambda_MuonLooseSF(dilep[1])
                # lepton SF are defined with get_scalefactors so by default op.defineOnFirstUse is used
                # DL trigger SF are using op.systematics so op.defineOnFirstUse must be expicitly used

            ElElTightSF = lambda dilep : self.lambda_ElectronTightSF(dilep[0]) + self.lambda_ElectronTightSF(dilep[1]) 
            MuMuTightSF = lambda dilep : self.lambda_MuonTightSF(dilep[0]) + self.lambda_MuonTightSF(dilep[1])
            ElMuTightSF = lambda dilep : self.lambda_ElectronTightSF(dilep[0]) + self.lambda_MuonTightSF(dilep[1])
        else:
            ElElLooseSF = lambda dilep : []
            MuMuLooseSF = lambda dilep : []
            ElMuLooseSF = lambda dilep : []

            ElElTightSF = lambda dilep : [] 
            MuMuTightSF = lambda dilep : []
            ElMuTightSF = lambda dilep : []

        # Sel #
        if "Preselected" in sel_level:
            ElElPreSelObject = SelectionObject(sel          = baseSel,
                                               selName      = "HasElElPreselected",
                                               yieldTitle   = "Preselected lepton pair (channel $e^+e^-$)")
            MuMuPreSelObject = SelectionObject(sel          = baseSel,
                                               selName      = "HasMuMuPreselected",
                                               yieldTitle   = "Preselected lepton pair (channel $\mu^+\mu^-$)")
            ElMuPreSelObject = SelectionObject(sel          = baseSel,
                                               selName      = "HasElMuPreselected",
                                               yieldTitle   = "Preselected lepton pair (channel $e^{\pm}\mu^{\mp}$)")

            # Selection #
            ElElPreSelObject.refine(cut     = [op.rng_len(self.ElElDileptonPreSel)>=1])
            MuMuPreSelObject.refine(cut     = [op.rng_len(self.MuMuDileptonPreSel)>=1])
            ElMuPreSelObject.refine(cut     = [op.rng_len(self.ElMuDileptonPreSel)>=1])

            # Yield #
            if plot_yield:
                ElElPreSelObject.makeYield(self.yieldPlots)
                MuMuPreSelObject.makeYield(self.yieldPlots)
                ElMuPreSelObject.makeYield(self.yieldPlots)

            # Record if actual argument #
            if "Preselected" in save_level:
                selectionDict["Preselected"] = [ElElPreSelObject,MuMuPreSelObject,ElMuPreSelObject]

            #--- Fakeable ---#
            if "Fakeable" in sel_level:
                ElElFakeSelObject = SelectionObject(sel         = ElElPreSelObject.sel,
                                                   selName      = "HasElElFakeable",
                                                   yieldTitle   = "OS fakeable lepton pair (channel $e^+e^-$)")
                MuMuFakeSelObject = SelectionObject(sel         = MuMuPreSelObject.sel,
                                                   selName      = "HasMuMuFakeable",
                                                   yieldTitle   = "OS fakeable lepton pair (channel $\mu^+\mu^-$)")
                ElMuFakeSelObject = SelectionObject(sel         = ElMuPreSelObject.sel,
                                                   selName      = "HasElMuFakeable",
                                                   yieldTitle   = "OS fakeable lepton pair (channel $e^{\pm}\mu^{\mp}$)")

                # Selection : at least one fakeable dilepton #
                # rng_len needs to be done on unsliced container to make sure it is not empty for sliced container afterwards
                ElElFakeSelObject.refine(cut = [op.rng_len(self.ElElDileptonFakeSel)>=1,
                                                lambdaOSDilepton(self.ElElDileptonFakeSel[0]),
                                                lambdaLowPtCutElEl(self.ElElDileptonFakeSel[0]),
                                                lambdaLeadingPtCutElEl(self.ElElDileptonFakeSel[0])],
                                         weight = ElElLooseSF(self.ElElDileptonFakeSel[0])) 
                MuMuFakeSelObject.refine(cut = [op.rng_len(self.MuMuDileptonFakeSel)>=1,
                                                lambdaOSDilepton(self.MuMuDileptonFakeSel[0]),
                                                lambdaLowPtCutMuMu(self.MuMuDileptonFakeSel[0]),
                                                lambdaLeadingPtCutMuMu(self.MuMuDileptonFakeSel[0])],
                                         weight = MuMuLooseSF(self.MuMuDileptonFakeSel[0]))
                ElMuFakeSelObject.refine(cut = [op.rng_len(self.ElMuDileptonFakeSel)>=1,
                                                lambdaOSDilepton(self.ElMuDileptonFakeSel[0]),
                                                lambdaLowPtCutElMu(self.ElMuDileptonFakeSel[0]),
                                                lambdaLeadingPtCutElMu(self.ElMuDileptonFakeSel[0])],
                                         weight = ElMuLooseSF(self.ElMuDileptonFakeSel[0]))

                # Yield #
                if plot_yield:
                    ElElFakeSelObject.makeYield(self.yieldPlots)
                    MuMuFakeSelObject.makeYield(self.yieldPlots)
                    ElMuFakeSelObject.makeYield(self.yieldPlots)

                # Selection : Mll cut + Z peak exclusion (charge already done in previous selection) in **preselected** leptons # 
                ElElFakeSelObject.selName += "PreMllCut"
                MuMuFakeSelObject.selName += "PreMllCut"
                ElMuFakeSelObject.selName += "PreMllCut"

                ElElFakeSelObject.yieldTitle += " + $M_{ll}>12$"
                MuMuFakeSelObject.yieldTitle += " + $M_{ll}>12$"
                ElMuFakeSelObject.yieldTitle += " + $M_{ll}>12$"

                mllCut = [lambda_lowMllCut(self.ElElDileptonPreSel),lambda_lowMllCut(self.MuMuDileptonPreSel),lambda_lowMllCut(self.ElMuDileptonPreSel)]

                ElEl_MllCut = copy(mllCut)
                MuMu_MllCut = copy(mllCut)
                ElMu_MllCut = copy(mllCut)

                if not self.args.NoZVeto: # No Mz cut in OF leptons
                    # All the preselected leptons outside Z peak region within 10 GeV
                    ElElFakeSelObject.selName += "OutZ"
                    MuMuFakeSelObject.selName += "OutZ"
                    ElMuFakeSelObject.selName += "OutZ"
                    ElElFakeSelObject.yieldTitle += " + Z Veto $|M_{ll}^{Preselected SF}-M_Z|>10$"
                    MuMuFakeSelObject.yieldTitle += " + Z Veto $|M_{ll}^{Preselected SF}-M_Z|>10$"
                    ElMuFakeSelObject.yieldTitle += " + Z Veto $|M_{ll}^{Preselected SF}-M_Z|>10$"
                    outZCut = [lambda_outZ(self.OSElElDileptonPreSel),lambda_outZ(self.OSMuMuDileptonPreSel)]
                    ElEl_MllCut.extend(outZCut)
                    MuMu_MllCut.extend(outZCut)
                    ElMu_MllCut.extend(outZCut)

                ElElFakeSelObject.refine(cut = ElEl_MllCut)
                MuMuFakeSelObject.refine(cut = MuMu_MllCut)
                ElMuFakeSelObject.refine(cut = ElMu_MllCut)

                # Yield #
                if plot_yield:
                    ElElFakeSelObject.makeYield(self.yieldPlots)
                    MuMuFakeSelObject.makeYield(self.yieldPlots)
                    ElMuFakeSelObject.makeYield(self.yieldPlots)
                                

                # Selection : tight cut veto, no more than two tight leptons in the event #
                ElElFakeSelObject.selName += "TightVeto"
                MuMuFakeSelObject.selName += "TightVeto"
                ElMuFakeSelObject.selName += "TightVeto"

                ElElFakeSelObject.yieldTitle += " + Two tights cut"
                MuMuFakeSelObject.yieldTitle += " + Two tights cut"
                ElMuFakeSelObject.yieldTitle += " + Two tights cut"

                ElElFakeSelObject.refine(cut = [op.rng_len(self.muonsTightSel)+op.rng_len(self.electronsTightSel) <= 2])
                MuMuFakeSelObject.refine(cut = [op.rng_len(self.muonsTightSel)+op.rng_len(self.electronsTightSel) <= 2])
                ElMuFakeSelObject.refine(cut = [op.rng_len(self.muonsTightSel)+op.rng_len(self.electronsTightSel) <= 2])

                # Yield #
                if plot_yield:
                    ElElFakeSelObject.makeYield(self.yieldPlots)
                    MuMuFakeSelObject.makeYield(self.yieldPlots)
                    ElMuFakeSelObject.makeYield(self.yieldPlots)

                # Record if actual argument #
                if "Fakeable" in save_level:
                    selectionDict["Fakeable"] = [ElElFakeSelObject,MuMuFakeSelObject,ElMuFakeSelObject]
                    
                #--- Tight ---#
                if "Tight" in sel_level:
                    ElElTightSelObject = SelectionObject(sel          = ElElFakeSelObject.sel,
                                                        selName      = "HasElElTight",
                                                        yieldTitle   = "OS tight lepton pair (channel $e^+e^-$)")
                    MuMuTightSelObject = SelectionObject(sel          = MuMuFakeSelObject.sel,
                                                        selName      = "HasMuMuTight",
                                                        yieldTitle   = "OS tight lepton pair (channel $\mu^+\mu^-$)")
                    ElMuTightSelObject = SelectionObject(sel          = ElMuFakeSelObject.sel,
                                                        selName      = "HasElMuTight",
                                                        yieldTitle   = "OS tight lepton pair (channel $e^{\pm}\mu^{\mp}$)")

                    # Selection : at least one tight dilepton #
                    if use_dd:
                        enable = "FakeExtrapolation" in self.datadrivenContributions and self.datadrivenContributions["FakeExtrapolation"].usesSample(self.sample, self.sampleCfg) 
                        ElElTightSelObject.create(ddSuffix  = "FakeExtrapolation",
                                                  cut       = [op.rng_len(self.ElElDileptonTightSel) >= 1],
                                                  weight    = ElElTightSF(self.ElElDileptonTightSel[0]),
                                                  ddCut     = [op.rng_len(self.ElElDileptonFakeExtrapolationSel) >= 1],
                                                  ddWeight  = [self.ElElFakeFactor(self.ElElDileptonFakeExtrapolationSel[0])]+ElElTightSF(self.ElElDileptonFakeExtrapolationSel[0]),
                                                  enable    = enable)
                        MuMuTightSelObject.create(ddSuffix  = "FakeExtrapolation",
                                                  cut       = [op.rng_len(self.MuMuDileptonTightSel) >= 1],
                                                  weight    = MuMuTightSF(self.MuMuDileptonTightSel[0]),
                                                  ddCut     = [op.rng_len(self.MuMuDileptonFakeExtrapolationSel) >= 1],
                                                  ddWeight  = [self.MuMuFakeFactor(self.MuMuDileptonFakeExtrapolationSel[0])]+MuMuTightSF(self.MuMuDileptonFakeExtrapolationSel[0]),
                                                  enable    = enable)
                        ElMuTightSelObject.create(ddSuffix  = "FakeExtrapolation",
                                                  cut       = [op.rng_len(self.ElMuDileptonTightSel) >= 1],
                                                  weight    = ElMuTightSF(self.ElMuDileptonTightSel[0]),
                                                  ddCut     = [op.rng_len(self.ElMuDileptonFakeExtrapolationSel) >= 1],
                                                  ddWeight  = [self.ElMuFakeFactor(self.ElMuDileptonFakeExtrapolationSel[0])]+ElMuTightSF(self.ElMuDileptonFakeExtrapolationSel[0]),
                                                  enable    = enable)
                    else:
                        ElElTightSelObject.refine(cut    = [op.rng_len(self.ElElDileptonTightSel) >= 1],
                                                  weight = ElElTightSF(self.ElElDileptonTightSel[0]))
                        MuMuTightSelObject.refine(cut    = [op.rng_len(self.MuMuDileptonTightSel) >= 1],
                                                  weight = MuMuTightSF(self.MuMuDileptonTightSel[0]))
                        ElMuTightSelObject.refine(cut    = [op.rng_len(self.ElMuDileptonTightSel) >= 1],
                                                  weight = ElMuTightSF(self.ElMuDileptonTightSel[0]))

                    # Yield #
                    if plot_yield:
                        ElElTightSelObject.makeYield(self.yieldPlots)
                        MuMuTightSelObject.makeYield(self.yieldPlots)
                        ElMuTightSelObject.makeYield(self.yieldPlots)

                    # Selection : Z peak for control region #
                    if self.args.ZPeak:
                        ElElTightSelObject.selName += "ZPeak"
                        MuMuTightSelObject.selName += "ZPeak"
                        ElElTightSelObject.yieldTitle+= " + Z peak $|M_{ll}-M_Z| < 10 GeV$" 
                        MuMuTightSelObject.yieldTitle+= " + Z peak $|M_{ll}-M_Z| < 10 GeV$"

                        ElElTightSelObject.refine(cut = [op.switch(op.rng_len(self.ElElDileptonTightSel)>=1,
                                                                   lambda_inZ(self.ElElDileptonTightSel[0]),
                                                                   lambda_inZ(self.ElElDileptonFakeExtrapolationSel[0]))])
                        MuMuTightSelObject.refine(cut = [op.switch(op.rng_len(self.MuMuDileptonTightSel)>=1,
                                                                   lambda_inZ(self.MuMuDileptonTightSel[0]),
                                                                   lambda_inZ(self.MuMuDileptonFakeExtrapolationSel[0]))])

                        if plot_yield:
                            ElElTightSelObject.makeYield(self.yieldPlots)
                            MuMuTightSelObject.makeYield(self.yieldPlots)
                    
                    # Record if actual argument #
                    if "Tight" in save_level:
                        selectionDict["Tight"] = [ElElTightSelObject,MuMuTightSelObject,ElMuTightSelObject]

                #--- Fake extrapolated ---#
                if "FakeExtrapolation" in sel_level:
                    ElElFakeExtrapolationSelObject = SelectionObject(sel          = ElElFakeSelObject.sel,
                                                                     selName      = "HasElElFakeExtrapolation",
                                                                     yieldTitle   = "OS fake extrapolated lepton pair (channel $e^+e^-$)")
                    MuMuFakeExtrapolationSelObject = SelectionObject(sel          = MuMuFakeSelObject.sel,
                                                                     selName      = "HasMuMuFakeExtrapolation",
                                                                     yieldTitle   = "OS fake extrapolated lepton pair (channel $\mu^+\mu^-$)")
                    ElMuFakeExtrapolationSelObject = SelectionObject(sel          = ElMuFakeSelObject.sel,
                                                                     selName      = "HasElMuFakeExtrapolation",
                                                                     yieldTitle   = "OS fake extrapolated lepton pair (channel $e^{\pm}\mu^{\mp}$)")
                    # Selection : at least one tight dilepton #
                    ElElFakeExtrapolationSelObject.refine(cut    = [op.rng_len(self.ElElDileptonFakeExtrapolationSel) >= 1],
                                                          weight = ElElTightSF(self.ElElDileptonFakeExtrapolationSel[0]))
                    MuMuFakeExtrapolationSelObject.refine(cut    = [op.rng_len(self.MuMuDileptonFakeExtrapolationSel) >= 1],
                                                          weight = MuMuTightSF(self.MuMuDileptonFakeExtrapolationSel[0]))
                    ElMuFakeExtrapolationSelObject.refine(cut    = [op.rng_len(self.ElMuDileptonFakeExtrapolationSel) >= 1],
                                                          weight = ElMuTightSF(self.ElMuDileptonFakeExtrapolationSel[0]))

                    # Yield #
                    if plot_yield:
                        ElElFakeExtrapolationSelObject.makeYield(self.yieldPlots)
                        MuMuFakeExtrapolationSelObject.makeYield(self.yieldPlots)
                        ElMuFakeExtrapolationSelObject.makeYield(self.yieldPlots)

                    # Record if actual argument #
                    if "FakeExtrapolation" in save_level:
                        selectionDict["FakeExtrapolation"] = [ElElFakeExtrapolationSelObject,MuMuFakeExtrapolationSelObject,ElMuFakeExtrapolationSelObject]
    # Return # 
    return selectionDict

def makeAtLeastTwoAk4JetSelection(self,selObject,copy_sel=False,plot_yield=False,use_dd=True):
    """
    Produces the Ak4 jet selection
    inputs :
        - selObject : SelectionObject class objects (contains bamboo selection, its name and the yield title)
        - copy_sel : wether to copy and return a new selection object built on the one provided, if not will modify it
    Careful : if copy_sel is False, the selObject will be modified
    Selection : at least two Ak4 jets
    """
    if copy_sel:
        selObject = copy(selObject)
    selObject.selName += "TwoAk4Jets"
    selObject.yieldTitle += " + Ak4 Jets $\geq 2$"
    selObject.refine(cut=[op.rng_len(self.ak4Jets)>=2])
    if plot_yield:
        selObject.makeYield(self.yieldPlots)
    if copy_sel :
        return selObject 

def makeAtLeastOneAk8JetSelection(self,selObject,copy_sel=False,plot_yield=False,use_dd=True):
    """
    Produces the Ak8 jet selection
    inputs :
        - selObject : SelectionObject class objects (contains bamboo selection, its name and the yield title)
        - copy_sel : wether to copy and return a new selection object built on the one provided, if not will modify it
    Careful : if copy_sel is False, the selObject will be modified
    Selection : at least one Ak8 jet
    """
    if copy_sel:
        selObject = copy(selObject)
    selObject.selName += "OneAk8Jet"
    selObject.yieldTitle += " + Ak8 Jets $\geq 1$"
    selObject.refine(cut=[op.rng_len(self.ak8Jets)>=1])
    if plot_yield:
        selObject.makeYield(self.yieldPlots)
    if copy_sel:
        return selObject
    
def makeExclusiveResolvedNoBtagSelection(self,selObject,copy_sel=False,plot_yield=False,use_dd=True):
    """
    Produces the exclusive resolved selection with no btagging
    inputs :
        - selObject : SelectionObject class objects (contains bamboo selection, its name and the yield title)
        - copy_sel : wether to copy and return a new selection object built on the one provided, if not will modify it
    Careful : if copy_sel is False, the selObject will be modified
    Selection : no btagged Ak8 jet (aka boosted), all Ak4 jets are non btagged
    """
    AppliedSF = None
    if copy_sel:
        selObject = copy(selObject)
    selObject.selName += "ExclusiveResolvedNoBtag"
    selObject.yieldTitle += " + Exclusive Resolved (0 bjet)"
    selObject.refine(cut=[op.rng_len(self.ak4BJets)==0,op.rng_len(self.ak8BJets)==0],weight=AppliedSF)
    if plot_yield:
        selObject.makeYield(self.yieldPlots)
    if copy_sel:
        return selObject

def makeExclusiveResolvedOneBtagSelection(self,selObject,copy_sel=False,plot_yield=False,use_dd=True):
    """
    Produces the exclusive resolved selection with only one btagging
    inputs :
        - selObject : SelectionObject class objects (contains bamboo selection, its name and the yield title)
        - copy_sel : wether to copy and return a new selection object built on the one provided, if not will modify it
    Careful : if copy_sel is False, the selObject will be modified
    Selection : no btagged Ak8 jet (aka boosted), only one btagged Ak4 jet
    """
    if copy_sel:
        selObject = copy(selObject)
    AppliedSF = None

    selObject.selName += "ExclusiveResolvedOneBtag"
    selObject.yieldTitle += " + Exclusive Resolved (1 bjet)"

    enable = use_dd and "DYEstimation" in self.datadrivenContributions and self.datadrivenContributions["DYEstimation"].usesSample(self.sample, self.sampleCfg) 
    if "ElEl" in selObject.selName:
        selObject.create(ddSuffix  = "DYEstimation",
                         cut       = [op.rng_len(self.ak4BJets)==1,op.rng_len(self.ak8BJets)==0],
                         weight    = AppliedSF,
                         ddCut     = [op.rng_len(self.ak4BJets)==0,op.rng_len(self.ak8BJets)==0],
                         ddWeight  = self.ResolvedDYReweighting1bElEl(self.ak4LightJetsByBtagScore[0]), 
                         enable    = enable)
    elif "MuMu" in selObject.selName:
        selObject.create(ddSuffix  = "DYEstimation",
                         cut       = [op.rng_len(self.ak4BJets)==1,op.rng_len(self.ak8BJets)==0],
                         weight    = AppliedSF,
                         ddCut     = [op.rng_len(self.ak4BJets)==0,op.rng_len(self.ak8BJets)==0],
                         ddWeight  = self.ResolvedDYReweighting1bMuMu(self.ak4LightJetsByBtagScore[0]), 
                         enable    = enable)
    elif "ElMu" in selObject.selName:
        selObject.refine(cut    = [op.rng_len(self.ak4BJets)==1,op.rng_len(self.ak8BJets)==0],
                         weight = AppliedSF)
    else:
        raise RuntimeError("Could not find the channel in selection name")

    if plot_yield:
        selObject.makeYield(self.yieldPlots)
    if copy_sel:
        return selObject

def makeExclusiveResolvedTwoBtagsSelection(self,selObject,copy_sel=False,plot_yield=False,use_dd=True):
    """
    Produces the exclusive resolved selection with no btagging
    inputs :
        - selObject : SelectionObject class objects (contains bamboo selection, its name and the yield title)
        - copy_sel : wether to copy and return a new selection object built on the one provided, if not will modify it
    Careful : if copy_sel is False, the selObject will be modified
    Selection : no btagged Ak8 jet (aka boosted), two Ak4 btagged jets 
    """
    AppliedSF = None
    if copy_sel:
        selObject = copy(selObject)

    selObject.selName += "ExclusiveResolvedTwoBtags"
    selObject.yieldTitle += " + Exclusive Resolved (2 bjets)"

    enable = use_dd and "DYEstimation" in self.datadrivenContributions and self.datadrivenContributions["DYEstimation"].usesSample(self.sample, self.sampleCfg) 
    if "ElEl" in selObject.selName:
        selObject.create(ddSuffix  = "DYEstimation",
                         cut       = [op.rng_len(self.ak4BJets)>=2,op.rng_len(self.ak8BJets)==0],
                         weight    = AppliedSF,
                         ddCut     = [op.rng_len(self.ak4BJets)==0,op.rng_len(self.ak8BJets)==0],
                         ddWeight  = self.ResolvedDYReweighting2bElEl(self.ak4LightJetsByBtagScore[0]), 
                         enable    = enable)
    elif "MuMu" in selObject.selName:
        selObject.create(ddSuffix  = "DYEstimation",
                         cut       = [op.rng_len(self.ak4BJets)>=2,op.rng_len(self.ak8BJets)==0],
                         weight    = AppliedSF,
                         ddCut     = [op.rng_len(self.ak4BJets)==0,op.rng_len(self.ak8BJets)==0],
                         ddWeight  = self.ResolvedDYReweighting2bMuMu(self.ak4LightJetsByBtagScore[0]), 
                         enable    = enable)
    elif "ElMu" in selObject.selName:
        selObject.refine(cut    = [op.rng_len(self.ak4BJets)>=2,op.rng_len(self.ak8BJets)==0],
                         weight = AppliedSF)
    else:
        raise RuntimeError("Could not find the channel in selection name")

    if plot_yield:
        selObject.makeYield(self.yieldPlots)
    if self.args.TTBarCR:
        selObject.selName += "MbbCut150"
        selObject.yieldTitle += " + $M_{bb}>150$"
        selObject.refine(cut = [op.invariant_mass(self.ak4BJets[0].p4,self.ak4BJets[1].p4)>150])
        if plot_yield:
            selObject.makeYield(self.yieldPlots)
    if copy_sel:
        return selObject

def makeInclusiveBoostedNoBtagSelection(self,selObject,copy_sel=False,plot_yield=False):
    """
    Produces the inclusive boosted selection(0 btag)
    inputs :
        - selObject : SelectionObject class objects (contains bamboo selection, its name and the yield title)
        - copy_sel : wether to copy and return a new selection object built on the one provided, if not will modify it
    Careful : if copy_sel is False, the selObject will be modified
    Selection : No btagged Ak8 jet
    """
    if copy_sel:
        selObject = copy(selObject)
    selObject.selName += "InclusiveBoostedNoBtag"
    selObject.yieldTitle += " + Inclusive Boosted 0 Btag"
    selObject.refine(cut = [op.rng_len(self.ak8BJets) == 0])

    if plot_yield:
        selObject.makeYield(self.yieldPlots)
    if copy_sel:
        return selObject

def makeInclusiveBoostedOneBtagSelection(self,selObject,copy_sel=False,plot_yield=False,use_dd=True):
    """
    Produces the inclusive boosted selection (1 btag)
    inputs :
        - selObject : SelectionObject class objects (contains bamboo selection, its name and the yield title)
        - copy_sel : wether to copy and return a new selection object built on the one provided, if not will modify it
    Careful : if copy_sel is False, the selObject will be modified
    Selection : At least one btagged Ak8 jet
    """
    AppliedSF = None
    if copy_sel:
        selObject = copy(selObject)
    selObject.selName += "InclusiveBoostedOneBtag"
    selObject.yieldTitle += " + Inclusive Boosted 1 Btag"

    enable = use_dd and "DYEstimation" in self.datadrivenContributions and self.datadrivenContributions["DYEstimation"].usesSample(self.sample, self.sampleCfg) 
    if "ElEl" in selObject.selName:
        selObject.create(ddSuffix  = "DYEstimation",
                         cut       = [op.rng_len(self.ak8BJets) >= 1],
                         weight    = AppliedSF,
                         ddCut     = [op.rng_len(self.ak8BJets) == 0],
                         ddWeight  = self.BoostedDYReweighting1bElEl(self.ak8Jets[0]), 
                         enable    = enable)
    elif "MuMu" in selObject.selName:
        selObject.create(ddSuffix  = "DYEstimation",
                         cut       = [op.rng_len(self.ak8BJets) >= 1],
                         weight    = AppliedSF,
                         ddCut     = [op.rng_len(self.ak8BJets) == 0],
                         ddWeight  = self.BoostedDYReweighting1bMuMu(self.ak8Jets[0]), 
                         enable    = enable)
    elif "ElMu" in selObject.selName:
        selObject.refine(cut    = [op.rng_len(self.ak8BJets)>=1],
                         weight = AppliedSF)
    else:
        raise RuntimeError("Could not find the channel in selection name")


    if plot_yield:
        selObject.makeYield(self.yieldPlots)
    if copy_sel:
        return selObject



def makeDNNOutputNodesSelections(self,selObject,output,plot_yield=False,suffix=""):
    selBaseObject = copy(selObject)
    selObjDict = {}
    maxIdx = op.rng_max_element_index(output)
    for i,node in enumerate(self.nodes):
        newSelObj = copy(selBaseObject)
        newSelObj.selName += suffix+'%snode'%node
        newSelObj.yieldTitle += " "+suffix+" in %s node"%node 
        newSelObj.refine(cut = [maxIdx == op.c_int(i)]) 
        selObjDict[node] = newSelObj
    return selObjDict

