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
    #--- Lambdas ---#
    # Fakeable lambdas #
    lambdaElPtCut = lambda lepColl: lepColl[0].pt > 32 
    lambdaMuPtCut = lambda lepColl: lepColl[0].pt > 25 

    # Tight #
    lambda_tight_ele = lambda ele : op.AND(self.lambda_is_matched(ele) , self.lambda_electronTightSel(ele))
    lambda_tight_mu  = lambda mu  : op.AND(self.lambda_is_matched(mu)  , self.lambda_muonTightSel(mu))
    
    lambda_fake_ele = lambda ele : op.AND(self.lambda_is_matched(ele) , op.NOT(self.lambda_electronTightSel(ele)))
    lambda_fake_mu  = lambda mu  : op.AND(self.lambda_is_matched(mu)  , op.NOT(self.lambda_muonTightSel(mu)))

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

    #---- Select fakeables -----#
    ElSelObj = SelectionObject(sel          = baseSel,
                               selName      = "Has1Fakeable",
                               yieldTitle   = "Fakeable lepton (channel $e^{\pm}$)")
    MuSelObj = SelectionObject(sel          = baseSel,
                               selName      = "Has1Fakeable",
                               yieldTitle   = "Fakeable lepton (channel $\mu^{\pm}$)")
    ElSelObj.refine(cut = [op.rng_len(self.electronsFakeSel) >= 1,
                           op.OR(op.rng_len(self.muonsFakeSel) == 0,
                                 self.electron_conept[self.electronsFakeSel[0].idx] > self.muon_conept[self.muonsFakeSel[0].idx])])
                        # With this selection the lead fakeable lepton is self.electronsFakeSel[0]
    MuSelObj.refine(cut = [op.rng_len(self.muonsFakeSel) >= 1,
                           op.OR(op.rng_len(self.electronsFakeSel) == 0,
                                 self.muon_conept[self.muonsFakeSel[0].idx] > self.electron_conept[self.electronsFakeSel[0].idx])])
                        # With this selection the lead fakeable lepton is self.muonsFakeSel[0]

    #---- Triggers ----#
    ElSelObj.selName += "WithTriggers"
    MuSelObj.selName += "WithTriggers"
    ElSelObj.yieldTitle += " + SE triggers"
    MuSelObj.yieldTitle += " + SM triggers"
    ElSelObj.refine(cut    = self.returnTriggers(["SingleElectron"]),
                    weight = ElLooseSF(self.electronsFakeSel))    
    MuSelObj.refine(cut    = self.returnTriggers(["SingleMuon"]),
                    weight = ElLooseSF(self.muonsFakeSel))    

    #---- PT cuts ----#
    ElSelObj.selName += "PtCuts"
    MuSelObj.selName += "PtCuts"
    ElSelObj.yieldTitle += " + $P_T > 32$"
    MuSelObj.yieldTitle += " + $P_T > 25$"
    ElSelObj.refine(cut = lambdaElPtCut(self.electronsFakeSel))
    MuSelObj.refine(cut = lambdaMuPtCut(self.muonsFakeSel))

    #---- Mll cuts ----#
    ElSelObject.selName += "PreMllCut"
    MuSelObject.selName += "PreMllCut"
    ElSelObject.yieldTitle += " + $M_{ll}>12$"
    MuSelObject.yieldTitle += " + $M_{ll}>12$"
            
    mllCut = [lambda_lowMllCut(self.ElElDileptonPreSel),lambda_lowMllCut(self.MuMuDileptonPreSel),lambda_lowMllCut(self.ElMuDileptonPreSel)]
    ElSelObject.refine(cut = mllCut)
    MuSelObject.refine(cut = mllCut)


    if not self.args.NoZVeto: 
        ElSelObject.selName += "OutZ"
        MuSelObject.selName += "OutZ"
        ElFakeSelObject.yieldTitle += " + Z Veto $|M_{ll}^{Preselected OS SF}-M_Z|>10$"
        MuFakeSelObject.yieldTitle += " + Z Veto $|M_{ll}^{Preselected OS SF}-M_Z|>10$"
        outZCut = [lambda_outZ(self.OSElElDileptonPreSel),lambda_outZ(self.OSMuMuDileptonPreSel)]
        ElSelObject.refine(cut = outZCut)
        MuSelObject.refine(cut = outZCut)
  
    #---- Tau veto ----#
    if not self.args.NoTauVeto: # nTau (isolated from fakeable leptons) = 0
        # All the preselected leptons outside Z peak region within 10 GeV
        ElSelObject.selName += "noTau"
        MuSelObject.selName += "noTau"
        ElSelObject.yieldTitle += " + Tau Veto"
        MuSelObject.yieldTitle += " + Tau Veto"

        ElSelObject.refine(cut = [op.rng_len(self.tauCleanSel) == 0])
        MuSelObject.refine(cut = [op.rng_len(self.tauCleanSel) == 0])


    #---- Tight selection ----#
    ElSelObject.selName += "TightSelected"
    MuSelObject.selName += "TightSelected"
    ElSelObject.yieldTitle += " + Tight selection"
    MuSelObject.yieldTitle += " + Tight selection"

    if use_dd:
        enable = "FakeExtrapolation" in self.datadrivenContributions and self.datadrivenContributions["FakeExtrapolation"].usesSample(self.sample, self.sampleCfg)
        ElSelObject.create(cut      = [lambda_tight_ele(self.electronsFakeSel[0]),
                                       op.rng_len(self.electronsTightSel) == 1,
                                       op.rng_len(self.muonsTightSel) == 0,
                                       self.electronsTightSel[0].idx == self.electronsFakeSel[0].idx],
                           weight   = ElTightSF(self.electronsTightSel),
                           ddSuffix = "FakeExtrapolation",
                           ddCut    = [lambda_fake_ele(self.electronsFakeSel[0])],
                           ddWeight = [self.ElFakeFactor(self.electronsFakeSel[0])]+ElTightSF(self.electronsFakeSel[0]),
                           enable   = enable)
        MuSelObject.create(cut      = [lambda_tight_mu(self.muonsFakeSel[0]),
                                       op.rng_len(self.muonsTightSel) == 1,
                                       op.rng_len(self.electronsTightSel) == 0,
                                       self.muonsTightSel[0].idx == self.muonsFakeSel[0].idx],
                           weight   = MuTightSF(self.muonsTightSel),
                           ddSuffix = "FakeExtrapolation",
                           ddCut    = [lambda_fake_mu(self.muonsFakeSel[0])],
                           ddWeight = [self.MuFakeFactor(self.muonsFakeSel[0])]+MuTightSF(self.muonsFakeSel[0]),
                           enable   = enable)
    else:
        ElSelObject.refine(cut    = [lambda_tight_ele(self.electronsFakeSel[0]),
                                     op.rng_len(self.electronsTightSel) == 1,
                                     op.rng_len(self.muonsTightSel) == 0,
                                     self.electronsTightSel[0].idx == self.electronsFakeSel[0].idx],
                           weight = ElTightSF(self.electronsTightSel))
        MuSelObject.refine(cut    = [lambda_tight_mu(self.muonsFakeSel[0]),
                                     op.rng_len(self.muonsTightSel) == 1,
                                     op.rng_len(self.electronsTightSel) == 0,
                                     self.muonsTightSel[0].idx == self.muonsFakeSel[0].idx],
                           weight = ElTightSF(self.electronsTightSel))
        # -> in SR : lead lepton is self.electronsTightSel[0] or self.muonsTightSel[0]
        # -> in Fake CR : lead lepton is self.electronsFakeSel[0] or self.muonsFakeSel[0]

    # Return # 
    return [ElSelObject,MuSelObject]

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
    #---- Common lambdas ----#
    lambdaOSDilepton = lambda dilep : dilep[0].charge != dilep[1].charge
    # Mll cut lambdas #
    ZMass = 91.1876
    lambda_lowMllCut    = lambda dileptons: op.NOT(op.rng_any(dileptons, lambda dilep : op.invariant_mass(dilep[0].p4, dilep[1].p4)<12.))
        # If any dilepton below 12GeV, returns False
    lambda_outZ         = lambda dileptons: op.NOT(op.rng_any(dileptons, lambda dilep : op.in_range(ZMass-10.,op.invariant_mass(dilep[0].p4, dilep[1].p4),ZMass+10.)))
        # If any dilepton within Z peak, returns False
    lambda_inZ          = lambda dilep : op.in_range(ZMass-10.,op.invariant_mass(dilep[0].p4, dilep[1].p4),ZMass+10.)

    ######### POG ID selection #########
    if self.args.POGID:
        #---- Lambdas ----#
        lambdaLowPtCut = lambda dilep: op.AND(dilep[0].pt > 15 , dilep[1].pt > 15) # subleading above 15 GeV
        lambdaLeadingPtCut = lambda dilep: op.OR(dilep[0].pt > 25 , dilep[1].pt > 25) # leading above 25 GeV
        #---- SF ----#
        if self.is_MC:
            ElElTightSF = lambda dilep : [op.defineOnFirstUse(self.lambda_ttH_doubleElectron_trigSF(dilep)),self.elPOGTight(dilep[0]),self.elPOGTight(dilep[1])]
            MuMuTightSF = lambda dilep : [op.defineOnFirstUse(self.lambda_ttH_doubleMuon_trigSF(dilep)),self.muPOGTightID(dilep[0]),self.muPOGTightISO(dilep[0]),self.muPOGTightID(dilep[1]),self.muPOGTightISO(dilep[1])]
            ElMuTightSF = lambda dilep : [op.defineOnFirstUse(self.lambda_ttH_electronMuon_trigSF(dilep)),self.elPOGTight(dilep[0]),self.muPOGTightID(dilep[1]),self.muPOGTightISO(dilep[1])]
        else:
            ElElTightSF = lambda dilep : [] 
            MuMuTightSF = lambda dilep : []
            ElMuTightSF = lambda dilep : []

        #---- Select tight leptons ----#
        ElElSelObj = SelectionObject(sel          = baseSel,
                                     selName      = "Has2Tight",
                                     yieldTitle   = "Tightdilepton (channel $e^+e^-$)")
        MuMuSelObj = SelectionObject(sel          = baseSel,
                                     selName      = "Has2Tight",
                                     yieldTitle   = "Tight dilepton (channel $\mu^+\mu^-$)")
        ElMuSelObj = SelectionObject(sel          = baseSel,
                                     selName      = "Has2Tight",
                                     yieldTitle   = "Tight dilepton (channel $e^{\pm}\mu^{\mp}$)")

        ElElSelObj.refine(cut    = [op.rng_len(self.ElElTightSel) == 1,
                                    op.rng_len(self.muonsTightSel) == 0],
                          weight = ElElTightSF(self.ElElTightSel[0]))  
        MuMuSelObj.refine(cut    = [op.rng_len(self.MuMuTightSel) == 1,
                                    op.rng_len(self.electronsTightSel) == 0],
                          weight = MuMuTightSF(self.MuMuTightSel[0]))  
        ElMuSelObj.refine(cut    = [op.rng_len(self.ElMuTightSel) == 1,
                                    op.rng_len(self.muonsTightSel) == 1,
                                    op.rng_len(self.electronsTightSel) == 1],
                          weight = ElMuTightSF(self.ElMuTightSel[0]))  


        #---- Opposite sign ----#
        ElElSelObj.selName += "OS"
        MuMuSelObj.selName += "OS"
        ElMuSelObj.selName += "OS"
        ElElSelObj.yieldTitle += " + OS"
        MuMuSelObj.yieldTitle += " + OS"
        ElMuSelObj.yieldTitle += " + OS"

        ElElSelObj.refine(cut = [lambdaOSDilepton(self.ElElTightSel[0])])
        MuMuSelObj.refine(cut = [lambdaOSDilepton(self.MuMuTightSel[0])])
        ElMuSelObj.refine(cut = [lambdaOSDilepton(self.ElMuTightSel[0])])

        #---- Triggers ----#
        ElElSelObj.selName += "WithTriggers"
        MuMuSelObj.selName += "WithTriggers"
        ElMuSelObj.selName += "WithTriggers"
        ElElSelObj.yieldTitle += " + SE/DE triggers"
        MuMuSelObj.yieldTitle += " + SM/DM triggers"
        ElMuSelObj.yieldTitle += " + SE/SM/EM triggers"

        ElElSelObj.refine(cut    = self.returnTriggers(["SingleElectron","DoubleEGamma"]),
                          weight = ElElLooseSF(self.ElElTightSel[0]))
        MuMuSelObj.refine(cut    = self.returnTriggers(["SingleMuon","DoubleMuon"]),
                          weight = MuMuLooseSF(self.MuMuTightSel[0]))
        ElMuSelObj.refine(cut    = self.returnTriggers(["SingleElectron","SingleMuon","MuonEG"]),
                          weight = ElMuLooseSF(self.ElMuTightSel[0]))

        #---- PT cuts ----#
        ElElSelObj.selName += "PtCuts"
        MuMuSelObj.selName += "PtCuts"
        ElMuSelObj.selName += "PtCuts"
        ElElSelObj.yieldTitle += " + $P_T^{leading} > 25$ + $P_T^{subleading} > 15$"
        MuMuSelObj.yieldTitle += " + $P_T^{leading} > 25$ + $P_T^{subleading} > 15$"
        ElMuSelObj.yieldTitle += " + $P_T^{leading} > 25$ + $P_T^{subleading} > 15$"

        ElElSelObj.refine(cut = [lambdaLowPtCutElEl(self.ElElTightSel[0]),lambdaHighPtCutElEl(self.ElElTightSel[0])])
        MuMuSelObj.refine(cut = [lambdaLowPtCutMuMu(self.MuMuTightSel[0]),lambdaHighPtCutMuMu(self.MuMuTightSel[0])])
        ElMuSelObj.refine(cut = [lambdaLowPtCutElMu(self.ElMuTightSel[0]),lambdaHighPtCutElMu(self.ElMuTightSel[0])])

        #---- Mll cuts ----#
        ElElSelObject.selName += "PreMllCut"
        ElElSelObject.yieldTitle += " + $M_{ll}>12$"
        MuMuSelObject.selName += "PreMllCut"
        MuMuSelObject.yieldTitle += " + $M_{ll}>12$"
        ElMuSelObject.selName += "PreMllCut"
        ElMuSelObject.yieldTitle += " + $M_{ll}>12$"

        mllCut = [lambda_lowMllCut(self.ElElDileptonPreSel),lambda_lowMllCut(self.MuMuDileptonPreSel)]
        ElElSelObject.refine(cut = mllCut)
        MuMuSelObject.refine(cut = mllCut)
        ElMuSelObject.refine(cut = mllCut)


        if not self.args.NoZVeto:
            ElElSelObject.selName += "OutZ"
            ElElSelObject.yieldTitle += " + Z Veto $|M_{ll}^{Preselected OS SF}-M_Z|>10$"
            MuMuSelObject.selName += "OutZ"
            MuMuSelObject.yieldTitle += " + Z Veto $|M_{ll}^{Preselected OS SF}-M_Z|>10$"
            ElMuSelObject.selName += "OutZ"
            ElMuSelObject.yieldTitle += " + Z Veto $|M_{ll}^{Preselected OS SF}-M_Z|>10$"

            outZCut = [lambda_outZ(self.OSElElDileptonPreSel),lambda_outZ(self.OSMuMuDileptonPreSel)]
            ElElSelObject.refine(cut = outZCut)
            MuMuSelObject.refine(cut = outZCut)
            ElMuSelObject.refine(cut = outZCut)

        if self.args.ZPeak:
            ElElSelObject.selName += "ZPeak"
            MuMuSelObject.selName += "ZPeak"
            ElElSelObject.yieldTitle+= " + Z peak $|M_{ll}-M_Z| < 10 GeV$" 
            MuMuSelObject.yieldTitle+= " + Z peak $|M_{ll}-M_Z| < 10 GeV$"

            ElElSelObject.refine(cut = [lambda_inZ(self.ElElTightSel[0])])
            MuMuSelObject.refine(cut = [lambda_inZ(self.MuMuTightSel[0])])



    ######### TTH ID selection #########
    if self.args.TTHIDLoose or self.args.TTHIDTight:
        #--- Lambdas ---#
        lambdaLowPtCutElEl = lambda dilep: op.AND( self.electron_conept[dilep[0].idx] > 15 , self.electron_conept[dilep[1].idx] > 15) # subleading above 15 GeV
        lambdaLowPtCutMuMu = lambda dilep: op.AND( self.muon_conept[dilep[0].idx] > 15 , self.muon_conept[dilep[1].idx] > 15) # subleading above 15 GeV
        lambdaLowPtCutElMu = lambda dilep: op.AND( self.electron_conept[dilep[0].idx] > 15 , self.muon_conept[dilep[1].idx] > 15) # subleading above 15 GeV
        lambdaLeadingPtCutElEl = lambda dilep: op.OR( self.electron_conept[dilep[0].idx] > 25 , self.electron_conept[dilep[1].idx] > 25) # leading above 25 GeV
        lambdaLeadingPtCutMuMu = lambda dilep: op.OR( self.muon_conept[dilep[0].idx] > 25 , self.muon_conept[dilep[1].idx] > 25) # leading above 25 GeV
        lambdaLeadingPtCutElMu = lambda dilep: op.OR( self.electron_conept[dilep[0].idx] > 25 , self.muon_conept[dilep[1].idx] > 25) # leading above 25 GeV
        # Need to split by flavour because conept definition needs it

        # Tight Dilepton : must also be Gen matched if MC #
        lambda_dilepton_matched = lambda dilep : op.AND(self.lambda_is_matched(dilep[0]),self.lambda_is_matched(dilep[1]))
        lambda_tightpair_ElEl = lambda dilep : op.AND(lambda_dilepton_matched(dilep),
                                                      self.lambda_electronTightSel(dilep[0]),
                                                      self.lambda_electronTightSel(dilep[1]))
        lambda_tightpair_MuMu = lambda dilep : op.AND(lambda_dilepton_matched(dilep),
                                                      self.lambda_muonTightSel(dilep[0]),
                                                      self.lambda_muonTightSel(dilep[1]))
        lambda_tightpair_ElMu = lambda dilep : op.AND(lambda_dilepton_matched(dilep),
                                                      self.lambda_electronTightSel(dilep[0]),
                                                      self.lambda_muonTightSel(dilep[1]))
        
        # Fake Extrapolation dilepton #
        lambda_fakepair_ElEl = lambda dilep : op.AND(lambda_dilepton_matched(dilep),
                                                     op.NOT(op.AND(self.lambda_electronTightSel(dilep[0]),
                                                                   self.lambda_electronTightSel(dilep[1]))))
        lambda_fakepair_MuMu = lambda dilep : op.AND(lambda_dilepton_matched(dilep),
                                                     op.NOT(op.AND(self.lambda_muonTightSel(dilep[0]),
                                                                   self.lambda_muonTightSel(dilep[1]))))
        lambda_fakepair_ElMu = lambda dilep : op.AND(lambda_dilepton_matched(dilep),
                                                     op.NOT(op.AND(self.lambda_electronTightSel(dilep[0]),
                                                                   self.lambda_muonTightSel(dilep[1]))))
        

        #---- SF ----#
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


        #---- Select fakeables -----# 
        ElElSelObj = SelectionObject(sel          = baseSel,
                                     selName      = "Has2Fakeable",
                                     yieldTitle   = "Fakeable dilepton (channel $e^+e^-$)")
        MuMuSelObj = SelectionObject(sel          = baseSel,
                                     selName      = "Has2Fakeable",
                                     yieldTitle   = "Fakeable dilepton (channel $\mu^+\mu^-$)")
        ElMuSelObj = SelectionObject(sel          = baseSel,
                                     selName      = "Has2Fakeable",
                                     yieldTitle   = "Fakeable dilepton (channel $e^{\pm}\mu^{\mp}$)")
        ElElSelObj.refine(cut = [op.rng_len(self.ElElFakeSel) >= 1,
                                 op.OR(op.rng_len(self.muonsFakeSel) == 0,
                                       op.AND(op.rng_len(self.muonsFakeSel) == 1,
                                              self.electron_conept[self.ElElFakeSel[0][0].idx] > self.muon_conept[self.muonsFakeSel[0].idx],
                                              self.electron_conept[self.ElElFakeSel[0][1].idx] > self.muon_conept[self.muonsFakeSel[0].idx]),
                                       op.AND(op.rng_len(self.muonsFakeSel) >= 2,
                                              self.electron_conept[self.ElElFakeSel[0][0].idx] > self.muon_conept[self.muonsFakeSel[0].idx],
                                              self.electron_conept[self.ElElFakeSel[0][1].idx] > self.muon_conept[self.muonsFakeSel[0].idx],
                                              self.electron_conept[self.ElElFakeSel[0][0].idx] > self.muon_conept[self.muonsFakeSel[1].idx],
                                              self.electron_conept[self.ElElFakeSel[0][1].idx] > self.muon_conept[self.muonsFakeSel[1].idx]))])
        MuMuSelObj.refine(cut = [op.rng_len(self.MuMuFakeSel) >= 1,
                                 op.OR(op.rng_len(self.electronsFakeSel) == 0,
                                       op.AND(op.rng_len(self.electronsFakeSel) == 1,
                                              self.muon_conept[self.MuMuFakeSel[0][0].idx] > self.electron_conept[self.electronsFakeSel[0].idx],
                                              self.muon_conept[self.MuMuFakeSel[0][1].idx] > self.electron_conept[self.electronsFakeSel[0].idx]),
                                       op.AND(op.rng_len(self.electronsFakeSel) >= 2,
                                              self.muon_conept[self.MuMuFakeSel[0][0].idx] > self.electron_conept[self.electronsFakeSel[0].idx],
                                              self.muon_conept[self.MuMuFakeSel[0][1].idx] > self.electron_conept[self.electronsFakeSel[0].idx],
                                              self.muon_conept[self.MuMuFakeSel[0][0].idx] > self.electron_conept[self.electronsFakeSel[1].idx],
                                              self.muon_conept[self.MuMuFakeSel[0][1].idx] > self.electron_conept[self.electronsFakeSel[1].idx]))])
        ElMuSelObj.refine(cut = [op.rng_len(self.ElMuFakeSel) >= 1,
                                 op.OR(op.rng_len(self.electronsFakeSel) + op.rng_len(self.muonsFakeSel) == 2,
                                       op.AND(op.rng_len(self.electronsFakeSel) >= 2,
                                              op.rng_len(self.muonsFakeSel) == 1,
                                              self.muon_conept[self.ElMuFakeSel[0][1].idx] > self.electron_conept[self.electronsFakeSel[1].idx]),
                                       op.AND(op.rng_len(self.muonsFakeSel) >= 2,
                                              op.rng_len(self.electronsFakeSel) == 1,
                                              self.electron_conept[self.ElMuFakeSel[0][0].idx] > self.muon_conept[self.muonsFakeSel[1].idx]),
                                       op.AND(op.rng_len(self.electronsFakeSel) >= 2,
                                              op.rng_len(self.muonsFakeSel) >= 2,
                                              self.muon_conept[self.ElMuFakeSel[0][1].idx] > self.electron_conept[self.electronsFakeSel[1].idx],
                                              self.electron_conept[self.ElMuFakeSel[0][0].idx] > self.muon_conept[self.muonsFakeSel[1].idx]))])

        #---- Opposite sign ----#
        ElElSelObj.selName += "OS"
        MuMuSelObj.selName += "OS"
        ElMuSelObj.selName += "OS"
        ElElSelObj.yieldTitle += " + OS"
        MuMuSelObj.yieldTitle += " + OS"
        ElMuSelObj.yieldTitle += " + OS"

        ElElSelObj.refine(cut = [lambdaOSDilepton(self.ElElFakeSel[0])])
        MuMuSelObj.refine(cut = [lambdaOSDilepton(self.MuMuFakeSel[0])])
        ElMuSelObj.refine(cut = [lambdaOSDilepton(self.ElMuFakeSel[0])])

        #---- Triggers ----#
        ElElSelObj.selName += "WithTriggers"
        MuMuSelObj.selName += "WithTriggers"
        ElMuSelObj.selName += "WithTriggers"
        ElElSelObj.yieldTitle += " + SE/DE triggers"
        MuMuSelObj.yieldTitle += " + SM/DM triggers"
        ElMuSelObj.yieldTitle += " + SE/SM/EM triggers"

        ElElSelObj.refine(cut    = self.returnTriggers(["SingleElectron","DoubleEGamma"]),
                          weight = ElElLooseSF(self.ElElFakeSel[0]))
        MuMuSelObj.refine(cut    = self.returnTriggers(["SingleMuon","DoubleMuon"]),
                          weight = MuMuLooseSF(self.MuMuFakeSel[0]))
        ElMuSelObj.refine(cut    = self.returnTriggers(["SingleElectron","SingleMuon","MuonEG"]),
                          weight = ElMuLooseSF(self.ElMuFakeSel[0]))

        #---- PT cuts ----#
        ElElSelObj.selName += "PtCuts"
        MuMuSelObj.selName += "PtCuts"
        ElMuSelObj.selName += "PtCuts"
        ElElSelObj.yieldTitle += " + $P_T^{leading} > 25$ + $P_T^{subleading} > 15$"
        MuMuSelObj.yieldTitle += " + $P_T^{leading} > 25$ + $P_T^{subleading} > 15$"
        ElMuSelObj.yieldTitle += " + $P_T^{leading} > 25$ + $P_T^{subleading} > 15$"

        ElElSelObj.refine(cut = [lambdaLowPtCutElEl(self.ElElFakeSel[0]),lambdaHighPtCutElEl(self.ElElFakeSel[0])])
        MuMuSelObj.refine(cut = [lambdaLowPtCutMuMu(self.MuMuFakeSel[0]),lambdaHighPtCutMuMu(self.MuMuFakeSel[0])])
        ElMuSelObj.refine(cut = [lambdaLowPtCutElMu(self.ElMuFakeSel[0]),lambdaHighPtCutElMu(self.ElMuFakeSel[0])])

        #---- Mll cuts ----#
        ElElSelObject.selName += "PreMllCut"
        ElElSelObject.yieldTitle += " + $M_{ll}>12$"
        MuMuSelObject.selName += "PreMllCut"
        MuMuSelObject.yieldTitle += " + $M_{ll}>12$"
        ElMuSelObject.selName += "PreMllCut"
        ElMuSelObject.yieldTitle += " + $M_{ll}>12$"

        mllCut = [lambda_lowMllCut(self.ElElDileptonPreSel),lambda_lowMllCut(self.MuMuDileptonPreSel),lambda_lowMllCut(self.ElMuDileptonPreSel)]
        ElElSelObject.refine(cut = mllCut)
        MuMuSelObject.refine(cut = mllCut)
        ElMuSelObject.refine(cut = mllCut)


        if not self.args.NoZVeto:
            ElElSelObject.selName += "OutZ"
            ElElSelObject.yieldTitle += " + Z Veto $|M_{ll}^{Preselected OS SF}-M_Z|>10$"
            MuMuSelObject.selName += "OutZ"
            MuMuSelObject.yieldTitle += " + Z Veto $|M_{ll}^{Preselected OS SF}-M_Z|>10$"
            ElMuSelObject.selName += "OutZ"
            ElMuSelObject.yieldTitle += " + Z Veto $|M_{ll}^{Preselected OS SF}-M_Z|>10$"

            outZCut = [lambda_outZ(self.OSElElDileptonPreSel),lambda_outZ(self.OSMuMuDileptonPreSel)]
            ElElSelObject.refine(cut = outZCut)
            MuMuSelObject.refine(cut = outZCut)
            ElMuSelObject.refine(cut = outZCut)

        if self.args.ZPeak:
            ElElSelObject.selName += "ZPeak"
            MuMuSelObject.selName += "ZPeak"
            ElElSelObject.yieldTitle+= " + Z peak $|M_{ll}-M_Z| < 10 GeV$" 
            MuMuSelObject.yieldTitle+= " + Z peak $|M_{ll}-M_Z| < 10 GeV$"

            ElElSelObject.refine(cut = [lambda_inZ(self.ElElFakeSel[0])])
            MuMuSelObject.refine(cut = [lambda_inZ(self.MuMuFakeSel[0])])

        
        #---- Tight selection ----#
        ElElSelObject.selName += "TightSelected"
        ElElSelObject.yieldTitle += " + Tight selection"
        MuMuSelObject.selName += "TightSelected"
        MuMuSelObject.yieldTitle += " + Tight selection"
        ElMuSelObject.selName += "TightSelected"
        ElMuSelObject.yieldTitle += " + Tight selection"
        
        if use_dd:
            enable = "FakeExtrapolation" in self.datadrivenContributions and self.datadrivenContributions["FakeExtrapolation"].usesSample(self.sample, self.sampleCfg)
            ElElSelObject.refine(cut      = [lambda_tightpair_ElEl(self.ElElFakeSel[0]),
                                             op.rng_len(self.electronsTightSel) == 2,
                                             op.rng_len(self.muonsTightSel) == 0,
                                             self.ElElTightSel[0][0].idx == self.ElElFakeSel[0][0].idx,
                                             self.ElElTightSel[0][1].idx == self.ElElFakeSel[0][1].idx],
                                 weight   = [ElElTightSF(self.ElElFakeSel[0])],
                                 ddSuffix = "FakeExtrapolation",
                                 ddCut    = [lambda_fakepair_ElEl(self.ElElFakeSel[0])],
                                 ddWeight = [self.ElElFakeFactor(ElElFakeSel[0])]+ElElTightSF(ElElFakeSel[0]),
                                 enable   = enable)
            MuMuSelObject.refine(cut      = [lambda_tightpair_MuMu(self.MuMuFakeSel[0]),
                                             op.rng_len(self.electronsTightSel) == 0,
                                             op.rng_len(self.muonsTightSel) == 2,
                                             self.MuMuTightSel[0][0].idx == self.MuMuFakeSel[0][0].idx,
                                             self.MuMuTightSel[0][1].idx == self.MuMuFakeSel[0][1].idx],
                                 weight   = [MuMuTightSF(self.MuMuFakeSel[0])],
                                 ddSuffix = "FakeExtrapolation",
                                 ddCut    = [lambda_fakepair_MuMu(self.MuMuFakeSel[0])],
                                 ddWeight = [self.MuMuFakeFactor(MuMuFakeSel[0])]+MuMuTightSF(MuMuFakeSel[0]),
                                 enable   = enable)
            ElMuSelObject.refine(cut      = [lambda_tightpair_ElMu(self.ElMuFakeSel[0]),
                                             op.rng_len(self.electronsTightSel) == 1,
                                             op.rng_len(self.muonsTightSel) == 1,
                                             self.ElMuTightSel[0][0].idx == self.ElMuFakeSel[0][0].idx,
                                             self.ElMuTightSel[0][1].idx == self.ElMuFakeSel[0][1].idx],
                                 weight   = [ElMuTightSF(self.ElMuFakeSel[0])],
                                 ddSuffix = "FakeExtrapolation",
                                 ddCut    = [lambda_fakepair_ElMu(self.ElMuFakeSel[0])],
                                 ddWeight = [self.ElMuFakeFactor(ElMuFakeSel[0])]+ElMuTightSF(ElMuFakeSel[0]),
                                 enable   = enable)
        else:
            ElElSelObject.refine(cut    = [lambda_tightpair_ElEl(self.ElElFakeSel[0]),
                                           op.rng_len(self.electronsTightSel) == 2,
                                           op.rng_len(self.muonsTightSel) == 0,
                                           self.ElElTightSel[0][0].idx == self.ElElFakeSel[0][0].idx,
                                           self.ElElTightSel[0][1].idx == self.ElElFakeSel[0][1].idx],
                                 weight = [ElElTightSF(self.ElElFakeSel[0])])
            MuMuSelObject.refine(cut    = [lambda_tightpair_MuMu(self.MuMuFakeSel[0]),
                                           op.rng_len(self.electronsTightSel) == 0,
                                           op.rng_len(self.muonsTightSel) == 2,
                                           self.MuMuTightSel[0][0].idx == self.MuMuFakeSel[0][0].idx,
                                           self.MuMuTightSel[0][1].idx == self.MuMuFakeSel[0][1].idx],
                                 weight = [MuMuTightSF(self.MuMuFakeSel[0])])
            ElMuSelObject.refine(cut    = [lambda_tightpair_ElMu(self.ElMuFakeSel[0]),
                                           op.rng_len(self.electronsTightSel) == 1,
                                           op.rng_len(self.muonsTightSel) == 1,
                                           self.ElMuTightSel[0][0].idx == self.ElMuFakeSel[0][0].idx,
                                           self.ElMuTightSel[0][1].idx == self.ElMuFakeSel[0][1].idx],
                                 weight = [ElMuTightSF(self.ElMuFakeSel[0])])

    return ElElSelObject,MuMuSelObject,ElMuSelObject
            


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

    if use_dd:
        enable = "DYEstimation" in self.datadrivenContributions and self.datadrivenContributions["DYEstimation"].usesSample(self.sample, self.sampleCfg) 
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
    else:
        selObject.refine(cut    = [op.rng_len(self.ak4BJets)==1,op.rng_len(self.ak8BJets)==0], 
                         weight = AppliedSF)

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

    if use_dd:
        enable = "DYEstimation" in self.datadrivenContributions and self.datadrivenContributions["DYEstimation"].usesSample(self.sample, self.sampleCfg) 
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
    else:
        selObject.refine(cut    = [op.rng_len(self.ak4BJets)>=2,op.rng_len(self.ak8BJets)==0], 
                         weight = AppliedSF)

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

    if use_dd:
        enable = "DYEstimation" in self.datadrivenContributions and self.datadrivenContributions["DYEstimation"].usesSample(self.sample, self.sampleCfg) 
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
    else:
        selObject.refine(cut    = [op.rng_len(self.ak8BJets)>=1],
                         weight = AppliedSF)
        


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

