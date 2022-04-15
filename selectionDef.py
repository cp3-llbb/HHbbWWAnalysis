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
    def __init__(self,sel,selName,yieldTitle,yieldObj=None,record_yields=False):
        self.sel           = sel
        self.selName       = selName
        self.yieldTitle    = yieldTitle
        self.yieldObj      = yieldObj
        self.record_yields = record_yields
    
    def refine(self,cut=None,weight=None,autoSyst=True):
        """
        Overload the refine function to avoid repeating the selName arg 
        """
        self.sel = self.sel.refine(name    = self.selName,
                                   cut     = cut,
                                   weight  = weight,
                                   autoSyst= autoSyst)
        if self.record_yields:
            self._makeYield()

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
        if self.record_yields:
            self._makeYield()

    def _makeYield(self):
        """
        Record an entry in the yield table 
        """
        self.yieldObj.add(self.sel,title=self.yieldTitle)

#===============================================================================================#
#                                        Selections                                             #
#===============================================================================================#
def makeSingleLeptonSelection(self,baseSel,plot_yield=False,use_dd=True,fake_selection=False): 
    """
    Produces the requested lepton selection (encapsulated in SelectionObject class objects)
    Will produce a dict :
        - key = level required (Preselected, Fakeable, Tight and/or FakeExtrapolation)
        - value = list of SelectionObject class objects per channel [El,Mu] 
    We start by leptons so no need to pass selObject
    Code organized such that selections are not repeated and hopefully optimzed the RooDataFrame
    """
    #--- Lambdas ---#
    # Fakeable lambdas #
    #lambdaElPtCut = lambda lepColl: lepColl[0].pt > 32 
    #lambdaMuPtCut = lambda lepColl: lepColl[0].pt > 25
    lambdaElPtCut = lambda lepColl: self.electron_conept[lepColl[0].idx] > 32.0
    lambdaMuPtCut = lambda lepColl: self.muon_conept[lepColl[0].idx] > 25.0

    # Tight #
    lambda_tight_ele = lambda ele : op.AND(self.lambda_is_matched(ele) , self.lambda_electronTightSel(ele))
    lambda_tight_mu  = lambda mu  : op.AND(self.lambda_is_matched(mu)  , self.lambda_muonTightSel(mu))
    
    if self.args.FakeRateNonClosureMCFakes:
        # MC fakes are events in the SR (tight) but non-prompt (the actual MC fakes, very unreliable, hence our datadriven estimate) #
        lambda_fake_ele = lambda ele : op.AND(op.NOT(self.lambda_is_matched(ele)) , self.lambda_electronTightSel(ele))
        lambda_fake_mu  = lambda mu  : op.AND(op.NOT(self.lambda_is_matched(mu))  , self.lambda_muonTightSel(mu))
        # We use unit-weight Fake Factors, because the events are actually taken from SR #
        lambda_FF_ele   = lambda ele : [op.c_float(1.)]
        lambda_FF_mu    = lambda mu  : [op.c_float(1.)]
        # Override the selections #
        fake_selection = True
        use_dd = False
    elif self.args.FakeRateNonClosureMCClosure:
        # MC closure are events in the CR (non-tight) but non-prompt -> used with MC QCD fake rates to evaluate closure #
        lambda_fake_ele = lambda ele : op.AND(op.NOT(self.lambda_is_matched(ele)) , op.NOT(self.lambda_electronTightSel(ele)))
        lambda_fake_mu  = lambda mu  : op.AND(op.NOT(self.lambda_is_matched(mu))  , op.NOT(self.lambda_muonTightSel(mu)))
        # We use the mc-driven Fake Factors #
        lambda_FF_ele   = lambda ele : [self.ElFakeFactorNonClosure(ele)]
        lambda_FF_mu    = lambda mu  : [self.MuFakeFactorNonClosure(mu)]
        # Override the selections #
        fake_selection = True
        use_dd = False
    else:
        # SR events used for Fake estimation are in the Fake CR (non-tight) and prompt for MC (for data, both prompt and non prompt, indistinguishable) #
        lambda_fake_ele = lambda ele : op.AND(self.lambda_is_matched(ele) , op.NOT(self.lambda_electronTightSel(ele)))
        lambda_fake_mu  = lambda mu  : op.AND(self.lambda_is_matched(mu)  , op.NOT(self.lambda_muonTightSel(mu)))
        # We use the datadriven Fake Factors  #
        lambda_FF_ele   = lambda ele : [self.ElFakeFactor(ele)]
        lambda_FF_mu    = lambda mu  : [self.MuFakeFactor(mu)]


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
        ElLooseSF = lambda lep : [op.defineOnFirstUse(self.ttH_singleElectron_trigSF(lep))] + self.lambda_ElectronLooseSF(lep)
        MuLooseSF = lambda lep : [op.defineOnFirstUse(self.ttH_singleMuon_trigSF(lep))] + self.lambda_MuonLooseSF(lep)
        # Tight SF #
        ElTightSF = lambda lep : self.lambda_ElectronTightSF(lep)
        MuTightSF = lambda lep : self.lambda_MuonTightSF(lep)
    else:
        ElLooseSF = lambda lep : []
        MuLooseSF = lambda lep : []
        ElTightSF = lambda lep : []
        MuTightSF = lambda lep : []


    #---- Select fakeables -----#
    ElSelObject = SelectionObject(sel          = baseSel,
                                  selName      = "Has1FakeableEl",
                                  yieldTitle   = "Fakeable lepton (channel $e^{\pm}$)",
                                  yieldObj     = self.yields,
                                  record_yields = self.args.PrintYield or self.args.OnlyYield)
    MuSelObject = SelectionObject(sel          = baseSel,
                                  selName      = "Has1FakeableMu",
                                  yieldTitle   = "Fakeable lepton (channel $\mu^{\pm}$)",
                                  yieldObj     = self.yields,
                                  record_yields = self.args.PrintYield or self.args.OnlyYield)
    ElSelObject.refine(cut = [op.rng_len(self.electronsFakeSel) >= 1,
                              op.OR(op.rng_len(self.muonsFakeSel) == 0,
                                    self.electron_conept[self.electronsFakeSel[0].idx] > self.muon_conept[self.muonsFakeSel[0].idx])])
    # With this selection the lead fakeable lepton is self.electronsFakeSel[0]
    MuSelObject.refine(cut = [op.rng_len(self.muonsFakeSel) >= 1,
                              op.OR(op.rng_len(self.electronsFakeSel) == 0,
                                    self.muon_conept[self.muonsFakeSel[0].idx] > self.electron_conept[self.electronsFakeSel[0].idx])])
    # With this selection the lead fakeable lepton is self.muonsFakeSel[0]

    #---- Triggers ----#
    ElSelObject.selName += "WithTriggers"
    MuSelObject.selName += "WithTriggers"
    ElSelObject.yieldTitle += " + SE triggers"
    MuSelObject.yieldTitle += " + SM triggers"
    ElSelObject.refine(cut    = self.returnTriggers(["SingleElectron"]),
                       weight = ElLooseSF(self.electronsFakeSel[0]))    
    MuSelObject.refine(cut    = self.returnTriggers(["SingleMuon"]),
                       weight = MuLooseSF(self.muonsFakeSel[0]))    

    #---- PT cuts ----#
    ElSelObject.selName += "PtCuts"
    MuSelObject.selName += "PtCuts"
    ElSelObject.yieldTitle += " + $P_T > 32$"
    MuSelObject.yieldTitle += " + $P_T > 25$"
    ElSelObject.refine(cut = lambdaElPtCut(self.electronsFakeSel))
    MuSelObject.refine(cut = lambdaMuPtCut(self.muonsFakeSel))

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
        ElSelObject.yieldTitle += " + Z Veto $|M_{ll}^{Preselected OS SF}-M_Z|>10$"
        MuSelObject.yieldTitle += " + Z Veto $|M_{ll}^{Preselected OS SF}-M_Z|>10$"
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
        # Use the datadriven in parallel of the SR selection #
        enable = "FakeExtrapolation" in self.datadrivenContributions and self.datadrivenContributions["FakeExtrapolation"].usesSample(self.sample, self.sampleCfg)
        ElSelObject.create(cut      = [lambda_tight_ele(self.electronsFakeSel[0]),
                                       op.rng_len(self.electronsTightSel) == 1,
                                       op.rng_len(self.muonsTightSel) == 0,
                                       self.electronsTightSel[0].idx == self.electronsFakeSel[0].idx],
                           weight   = ElTightSF(self.electronsTightSel[0]),
                           ddSuffix = "FakeExtrapolation",
                           ddCut    = [lambda_fake_ele(self.electronsFakeSel[0]), op.rng_len(self.electronsTightSel)+op.rng_len(self.muonsTightSel)<=1],
                           ddWeight = lambda_FF_ele(self.electronsFakeSel[0])+ElTightSF(self.electronsFakeSel[0]),
                           enable   = enable)
        MuSelObject.create(cut      = [lambda_tight_mu(self.muonsFakeSel[0]),
                                       op.rng_len(self.muonsTightSel) == 1,
                                       op.rng_len(self.electronsTightSel) == 0,
                                       self.muonsTightSel[0].idx == self.muonsFakeSel[0].idx],
                           weight   = MuTightSF(self.muonsTightSel[0]),
                           ddSuffix = "FakeExtrapolation",
                           ddCut    = [lambda_fake_mu(self.muonsFakeSel[0]), op.rng_len(self.electronsTightSel)+op.rng_len(self.muonsTightSel)<=1],
                           ddWeight = lambda_FF_mu(self.muonsFakeSel[0])+MuTightSF(self.muonsFakeSel[0]),
                           enable   = enable)
    elif fake_selection:
        # Only return the datadriven CR selection (eg, skimmer or non closure) #
        ElSelObject.refine(cut    = [lambda_fake_ele(self.electronsFakeSel[0]), op.rng_len(self.electronsTightSel)+op.rng_len(self.muonsTightSel)<=1],
                           weight = lambda_FF_ele(self.electronsFakeSel[0]))
        MuSelObject.refine(cut    = [lambda_fake_mu(self.muonsFakeSel[0]), op.rng_len(self.electronsTightSel)+op.rng_len(self.muonsTightSel)<=1],
                           weight = lambda_FF_mu(self.muonsFakeSel[0]))
        # -> in SR : lead lepton is self.electronsTightSel[0] or self.muonsTightSel[0]
        # -> in Fake CR : lead lepton is self.electronsFakeSel[0] or self.muonsFakeSel[0]
    else:
        # Only return the SR selection #
        ElSelObject.refine(cut    = [lambda_tight_ele(self.electronsFakeSel[0]),
                                     op.rng_len(self.electronsTightSel) == 1,
                                     op.rng_len(self.muonsTightSel) == 0,
                                     self.electronsTightSel[0].idx == self.electronsFakeSel[0].idx],
                           weight = ElTightSF(self.electronsTightSel[0]))
        MuSelObject.refine(cut    = [lambda_tight_mu(self.muonsFakeSel[0]),
                                     op.rng_len(self.muonsTightSel) == 1,
                                     op.rng_len(self.electronsTightSel) == 0,
                                     self.muonsTightSel[0].idx == self.muonsFakeSel[0].idx],
                           weight = MuTightSF(self.muonsTightSel[0]))

    # Return # 
    return [ElSelObject,MuSelObject]

def makeResolvedSelection(self,selObject,copy_sel=False):
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
    selObject.selName += "Resolved"
    selObject.yieldTitle += " + Ak4 Jets $\geq 3$"
    selObject.refine(cut  = [op.rng_len(self.ak4Jets)  >= 3,
                             op.rng_len(self.ak4BJets) >= 1,
                             op.rng_len(self.ak8BJets) == 0])
    #if plot_yield:
    #    selObject.makeYield(self.yieldPlots)
    if copy_sel :
        return selObject 

def makeBoostedSelection(self,selObject,copy_sel=False):
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
    selObject.selName += "Boosted"
    selObject.yieldTitle += " + Ak8 jets $\geq 1$ Ak4jets $\geq 1$"
    selObject.refine(cut=[op.rng_len(self.ak8BJets) >= 1, op.rng_len(self.ak4JetsCleanedFromAk8b) >= 1])
    #selObject.refine(cut=[op.rng_len(self.ak8BJets) >= 1])
    #if plot_yield:
    #    selObject.makeYield(self.yieldPlots)
    if copy_sel:
        return selObject

def makeCoarseResolvedSelection(self,selObject,nJet,copy_sel=False):
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
    if copy_sel :
        return selObject 


def makeCoarseBoostedSelection(self,selObject,copy_sel=False):
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
    if copy_sel:
        return selObject

def makeExclusiveLooseResolvedJetComboSelection(self,selObject,nbJet,copy_sel=False):
    if copy_sel:
        selObject = copy(selObject)

    if nbJet == 0:
        AppliedSF = None
        selObject.selName += "has0bJets"
        selObject.yieldTitle += " + Exclusive Resolved (nbJet = 0)"
        selObject.refine(cut    = [op.rng_len(self.ak4BJets) == 0],
                         weight = AppliedSF)
        
    elif nbJet == 1:
        AppliedSF = None
        #selObject.selName += "ExclusiveResolved1b2j"
        #selObject.yieldTitle += " + Exclusive Resolved (nbJet = 1, nAk4LightJet = 2)"
        selObject.selName += "has1bJets"
        selObject.yieldTitle += " + Exclusive Resolved (nbJet = 1)"
        selObject.refine(cut    = [op.rng_len(self.ak4BJets) == 1],
                         weight = AppliedSF)

    elif nbJet == 2:
        AppliedSF = None
        selObject.selName += "has2bJets"
        selObject.yieldTitle += " + Exclusive Resolved (nbJet = 2)"
        selObject.refine(cut    = [op.rng_len(self.ak4BJets) >= 2],
                         weight = None)
    else: raise RuntimeError ("Error in Jet Selection!!!")
    
    if copy_sel:
        return selObject


def makeExclusiveTightResolvedJetComboSelection(self,selObject,nbJet,copy_sel=False):
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

    if copy_sel:
        return selObject


def makeExclusiveResolvedJetComboSelection(self,selObject,nbJet,nJet,copy_sel=False):
    if copy_sel:
        selObject = copy(selObject)

    if nJet == 4 : 
        if nbJet == 0:
            AppliedSF = None
            selObject.selName += "ExclusiveResolved0b4j"
            selObject.yieldTitle += " + Exclusive Resolved (nbJet = 0, nAk4LightJets $\geq 4$)"
            selObject.refine(cut    = [op.rng_len(self.ak4BJets) == 0, op.rng_len(self.ak4LightJetsByPt) >= 4],
                             weight = AppliedSF)
            
        elif nbJet == 1:
            AppliedSF = None
            selObject.selName += "ExclusiveResolved1b3j"
            selObject.yieldTitle += " + Exclusive Resolved (nbjet = 1, nAk4LightJets $\geq 3$)"
            selObject.refine(cut    = [op.rng_len(self.ak4BJets) == 1, op.rng_len(self.ak4LightJetsByPt) >= 3],
                             weight = AppliedSF)
            
        elif nbJet == 2:
            AppliedSF = None
            selObject.selName += "ExclusiveResolved2b2j"
            selObject.yieldTitle += " + Exclusive Resolved (nbjet $\geq 2$, nAk4LightJets $\geq 2$)"
            selObject.refine(cut   = [op.rng_len(self.ak4BJets) >= 2,op.rng_len(self.ak4LightJetsByPt) >= 2],
                             weight = AppliedSF)

        else: raise RuntimeError ("Error in Jet Selection!!!")

    if nJet == 3 : 
        if nbJet == 0:
            AppliedSF = None
            selObject.selName += "ExclusiveResolved0b3j"
            selObject.yieldTitle += " + Exclusive Resolved (nbJet = 0, nAk4LightJet = 3)"
            selObject.refine(cut    = [op.rng_len(self.ak4BJets) == 0, op.rng_len(self.ak4LightJetsByPt) == 3],
                             weight = AppliedSF)
        
        elif nbJet == 1:
            AppliedSF = None
            selObject.selName += "ExclusiveResolved1b2j"
            selObject.yieldTitle += " + Exclusive Resolved (nbJet = 1, nAk4LightJet = 2)"
            selObject.refine(cut    = [op.rng_len(self.ak4BJets) == 1, op.rng_len(self.ak4LightJetsByPt) == 2],
                             weight = AppliedSF)
            
        elif nbJet == 2:
            AppliedSF = None
            selObject.selName += "ExclusiveResolved2b1j"
            selObject.yieldTitle += " + Exclusive Resolved (nbJet = 2, nAk4LightJet = 1)"
            selObject.refine(cut    = [op.rng_len(self.ak4BJets) == 2, op.rng_len(self.ak4LightJetsByPt) == 1],
                             weight = None)

        else: raise RuntimeError ("Error in Jet Selection!!!")

def makeExclusiveResolvedSelection(self,selObject,nbJet,copy_sel=False):
    if copy_sel:
        selObject = copy(selObject)
        
    if nbJet == 1:    
        AppliedSF = None
        selObject.selName += "ExclusiveResolved1b"
        selObject.yieldTitle += " + Exclusive Resolved (nbjet = 1)"
        selObject.refine(cut    = [op.rng_len(self.ak4BJets) == 1],
                         weight = AppliedSF)
        
    elif nbJet == 2:
        AppliedSF = None
        selObject.selName += "ExclusiveResolved2b"
        selObject.yieldTitle += " + Exclusive Resolved (nbjet $\geq 2$)"
        selObject.refine(cut   = [op.rng_len(self.ak4BJets) >= 2],
                         weight = AppliedSF)
        
    else: raise RuntimeError ("Error in Jet Selection!!!")
    
    if copy_sel:
        return selObject

def makeSemiBoostedHbbSelection(self,selObject,nNonb,copy_sel=False):
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

    if copy_sel:
        return selObject

def makeVBFSelection(self,selObject,vbfJetPairs,copy_sel=False):
    if copy_sel:
        selObject = copy(selObject)
    selObject.selName += "HasVBFjets"
    selObject.yieldTitle += " + VBF Jet Pairs $\geq 1$"
    selObject.refine(cut=[op.rng_len(vbfJetPairs) >= 1])
    if copy_sel:
        return selObject

def makeNoVBFSelection(self,selObject,vbfJetPairs,copy_sel=False):
    if copy_sel:
        selObject = copy(selObject)
    selObject.selName += "HasNoVBFjets"
    selObject.yieldTitle += " + VBF Jet Pairs = 0"
    selObject.refine(cut=[op.rng_len(vbfJetPairs) == 0])
    if copy_sel:
        return selObject


def makeDoubleLeptonSelection(self,baseSel,use_dd=True,fake_selection=False): 
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

    # Pt cuts #
    lambdaLowPtCutElEl = lambda dilep: op.AND( self.electron_conept[dilep[0].idx] > 15 , self.electron_conept[dilep[1].idx] > 15) # subleading above 15 GeV
    lambdaLowPtCutMuMu = lambda dilep: op.AND( self.muon_conept[dilep[0].idx] > 15 , self.muon_conept[dilep[1].idx] > 15) # subleading above 15 GeV
    lambdaLowPtCutElMu = lambda dilep: op.AND( self.electron_conept[dilep[0].idx] > 15 , self.muon_conept[dilep[1].idx] > 15) # subleading above 15 GeV
    lambdaLeadingPtCutElEl = lambda dilep: op.OR( self.electron_conept[dilep[0].idx] > 25 , self.electron_conept[dilep[1].idx] > 25) # leading above 25 GeV
    lambdaLeadingPtCutMuMu = lambda dilep: op.OR( self.muon_conept[dilep[0].idx] > 25 , self.muon_conept[dilep[1].idx] > 25) # leading above 25 GeV
    lambdaLeadingPtCutElMu = lambda dilep: op.OR( self.electron_conept[dilep[0].idx] > 25 , self.muon_conept[dilep[1].idx] > 25) # leading above 25 GeV
    # Need to split by flavour because conept definition needs it

    ######### TTH ID selection #########
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
    ElElSelObj = SelectionObject(sel           = baseSel,
                                 selName       = "Has2FakeableElEl",
                                 yieldTitle    = "Fakeable dilepton (channel $e^+e^-$)          ",
                                 yieldObj      = self.yields,
                                 record_yields = self.args.PrintYield or self.args.OnlyYield)
    MuMuSelObj = SelectionObject(sel           = baseSel,
                                 selName       = "Has2FakeableMuMu",
                                 yieldTitle    = "Fakeable dilepton (channel $\mu^+\mu^-$)      ",
                                 yieldObj      = self.yields,
                                 record_yields = self.args.PrintYield or self.args.OnlyYield)
    ElMuSelObj = SelectionObject(sel           = baseSel,
                                 selName       = "Has2FakeableElMu",
                                 yieldTitle    = "Fakeable dilepton (channel $e^{\pm}\mu^{\mp}$)",
                                 yieldObj      = self.yields,
                                 record_yields = self.args.PrintYield or self.args.OnlyYield)
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
                             op.OR(op.AND(op.rng_len(self.electronsFakeSel) == 1,
                                          op.rng_len(self.muonsFakeSel) == 1),
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

    ElElSelObj.refine(cut = [lambdaLowPtCutElEl(self.ElElFakeSel[0]),lambdaLeadingPtCutElEl(self.ElElFakeSel[0])])
    MuMuSelObj.refine(cut = [lambdaLowPtCutMuMu(self.MuMuFakeSel[0]),lambdaLeadingPtCutMuMu(self.MuMuFakeSel[0])])
    ElMuSelObj.refine(cut = [lambdaLowPtCutElMu(self.ElMuFakeSel[0]),lambdaLeadingPtCutElMu(self.ElMuFakeSel[0])])

    #---- Mll cuts ----#
    ElElSelObj.selName += "PreMllCut"
    ElElSelObj.yieldTitle += " + $M_{ll}>12$"
    MuMuSelObj.selName += "PreMllCut"
    MuMuSelObj.yieldTitle += " + $M_{ll}>12$"
    ElMuSelObj.selName += "PreMllCut"
    ElMuSelObj.yieldTitle += " + $M_{ll}>12$"

    mllCut = [lambda_lowMllCut(self.ElElDileptonPreSel),lambda_lowMllCut(self.MuMuDileptonPreSel),lambda_lowMllCut(self.ElMuDileptonPreSel)]
    ElElSelObj.refine(cut = mllCut)
    MuMuSelObj.refine(cut = mllCut)
    ElMuSelObj.refine(cut = mllCut)


    if not self.args.NoZVeto:
        ElElSelObj.selName += "OutZ"
        ElElSelObj.yieldTitle += " + Z Veto $|M_{ll}^{Preselected OS SF}-M_Z|>10$"
        MuMuSelObj.selName += "OutZ"
        MuMuSelObj.yieldTitle += " + Z Veto $|M_{ll}^{Preselected OS SF}-M_Z|>10$"
        ElMuSelObj.selName += "OutZ"
        ElMuSelObj.yieldTitle += " + Z Veto $|M_{ll}^{Preselected OS SF}-M_Z|>10$"

        outZCut = [lambda_outZ(self.OSElElDileptonPreSel),lambda_outZ(self.OSMuMuDileptonPreSel)]
        ElElSelObj.refine(cut = outZCut)
        MuMuSelObj.refine(cut = outZCut)
        ElMuSelObj.refine(cut = outZCut)

    if self.args.ZPeak:
        ElElSelObj.selName += "ZPeak"
        MuMuSelObj.selName += "ZPeak"
        ElElSelObj.yieldTitle+= " + Z peak $|M_{ll}-M_Z| < 10 GeV$" 
        MuMuSelObj.yieldTitle+= " + Z peak $|M_{ll}-M_Z| < 10 GeV$"

        ElElSelObj.refine(cut = [lambda_inZ(self.ElElFakeSel[0])])
        MuMuSelObj.refine(cut = [lambda_inZ(self.MuMuFakeSel[0])])
    
    #---- Tight selection ----#
    ElElSelObj.selName += "TightSelected"
    ElElSelObj.yieldTitle += " + Tight selection"
    MuMuSelObj.selName += "TightSelected"
    MuMuSelObj.yieldTitle += " + Tight selection"
    ElMuSelObj.selName += "TightSelected"
    ElMuSelObj.yieldTitle += " + Tight selection"
    
    if use_dd:
        # Use datadriven in parrallel of actual selection #
        enable = "FakeExtrapolation" in self.datadrivenContributions and self.datadrivenContributions["FakeExtrapolation"].usesSample(self.sample, self.sampleCfg)
        ElElSelObj.create(cut      = [self.lambda_tightpair_ElEl(self.ElElFakeSel[0]),
                                      op.rng_len(self.electronsTightSel) == 2,
                                      op.rng_len(self.muonsTightSel) == 0,
                                      self.ElElTightSel[0][0].idx == self.ElElFakeSel[0][0].idx,
                                      self.ElElTightSel[0][1].idx == self.ElElFakeSel[0][1].idx],
                          weight   = ElElTightSF(self.ElElFakeSel[0]),
                          ddSuffix = "FakeExtrapolation",
                          ddCut    = [self.lambda_fakepair_ElEl(self.ElElFakeSel[0]),
                                      op.rng_len(self.electronsTightSel) + op.rng_len(self.muonsTightSel)<=2],
                          ddWeight = [self.ElElFakeFactor(self.ElElFakeSel[0])]+ElElTightSF(self.ElElFakeSel[0]),
                          enable   = enable)
        MuMuSelObj.create(cut      = [self.lambda_tightpair_MuMu(self.MuMuFakeSel[0]),
                                      op.rng_len(self.electronsTightSel) == 0,
                                      op.rng_len(self.muonsTightSel) == 2,
                                      self.MuMuTightSel[0][0].idx == self.MuMuFakeSel[0][0].idx,
                                      self.MuMuTightSel[0][1].idx == self.MuMuFakeSel[0][1].idx],
                          weight   = MuMuTightSF(self.MuMuFakeSel[0]),
                          ddSuffix = "FakeExtrapolation",
                          ddCut    = [self.lambda_fakepair_MuMu(self.MuMuFakeSel[0]),
                                      op.rng_len(self.electronsTightSel) + op.rng_len(self.muonsTightSel)<=2],
                          ddWeight = [self.MuMuFakeFactor(self.MuMuFakeSel[0])]+MuMuTightSF(self.MuMuFakeSel[0]),
                          enable   = enable)
        ElMuSelObj.create(cut      = [self.lambda_tightpair_ElMu(self.ElMuFakeSel[0]),
                                      op.rng_len(self.electronsTightSel) == 1,
                                      op.rng_len(self.muonsTightSel) == 1,
                                      self.ElMuTightSel[0][0].idx == self.ElMuFakeSel[0][0].idx,
                                      self.ElMuTightSel[0][1].idx == self.ElMuFakeSel[0][1].idx],
                          weight   = ElMuTightSF(self.ElMuFakeSel[0]),
                          ddSuffix = "FakeExtrapolation",
                          ddCut    = [self.lambda_fakepair_ElMu(self.ElMuFakeSel[0]),
                                      op.rng_len(self.electronsTightSel) + op.rng_len(self.muonsTightSel)<=2],
                          ddWeight = [self.ElMuFakeFactor(self.ElMuFakeSel[0])]+ElMuTightSF(self.ElMuFakeSel[0]),
                          enable   = enable)
    elif fake_selection:
        # Return only the selection fot the fake datadriven (eg for skimmer) #
        ElElSelObj.refine(cut    = [self.lambda_fakepair_ElEl(self.ElElFakeSel[0]),
                                    op.rng_len(self.electronsTightSel) + op.rng_len(self.muonsTightSel)<=2],
                          weight = ElElTightSF(self.ElElFakeSel[0]))
        MuMuSelObj.refine(cut    = [self.lambda_fakepair_MuMu(self.MuMuFakeSel[0]),
                                    op.rng_len(self.electronsTightSel) + op.rng_len(self.muonsTightSel)<=2],
                          weight = MuMuTightSF(self.MuMuFakeSel[0]))
        ElMuSelObj.refine(cut    = [self.lambda_fakepair_ElMu(self.ElMuFakeSel[0]),
                                    op.rng_len(self.electronsTightSel) + op.rng_len(self.muonsTightSel)<=2],
                          weight = ElMuTightSF(self.ElMuFakeSel[0]))
    else:
        # Return only the selection for the SR #
        ElElSelObj.refine(cut    = [self.lambda_tightpair_ElEl(self.ElElFakeSel[0]),
                                       op.rng_len(self.electronsTightSel) == 2,
                                       op.rng_len(self.muonsTightSel) == 0,
                                       self.ElElTightSel[0][0].idx == self.ElElFakeSel[0][0].idx,
                                       self.ElElTightSel[0][1].idx == self.ElElFakeSel[0][1].idx],
                          weight = ElElTightSF(self.ElElFakeSel[0]))
        MuMuSelObj.refine(cut    = [self.lambda_tightpair_MuMu(self.MuMuFakeSel[0]),
                                       op.rng_len(self.electronsTightSel) == 0,
                                       op.rng_len(self.muonsTightSel) == 2,
                                       self.MuMuTightSel[0][0].idx == self.MuMuFakeSel[0][0].idx,
                                       self.MuMuTightSel[0][1].idx == self.MuMuFakeSel[0][1].idx],
                          weight = MuMuTightSF(self.MuMuFakeSel[0]))
        ElMuSelObj.refine(cut    = [self.lambda_tightpair_ElMu(self.ElMuFakeSel[0]),
                                       op.rng_len(self.electronsTightSel) == 1,
                                       op.rng_len(self.muonsTightSel) == 1,
                                       self.ElMuTightSel[0][0].idx == self.ElMuFakeSel[0][0].idx,
                                       self.ElMuTightSel[0][1].idx == self.ElMuFakeSel[0][1].idx],
                          weight = ElMuTightSF(self.ElMuFakeSel[0]))



    return ElElSelObj,MuMuSelObj,ElMuSelObj
            


def makeAtLeastTwoAk4JetSelection(self,selObject,copy_sel=False,use_dd=True):
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
    if copy_sel :
        return selObject 

def makeAtLeastOneAk8JetSelection(self,selObject,copy_sel=False,use_dd=True):
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
    if copy_sel:
        return selObject
    
def makeExclusiveResolvedNoBtagSelection(self,selObject,copy_sel=False,use_dd=True):
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
    selObject.refine(cut=[op.rng_len(self.ak4BJets)==0,op.rng_len(self.ak8BJets)==0,op.rng_len(self.ak8Jets)==0],weight=AppliedSF)
    if copy_sel:
        return selObject

def makeExclusiveResolvedOneBtagSelection(self,selObject,copy_sel=False,use_dd=True,dy_selection=False):
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
        # Use the datadriven selection in parallel of the SR one #
        enable = "DYEstimation" in self.datadrivenContributions and self.datadrivenContributions["DYEstimation"].usesSample(self.sample, self.sampleCfg) 
        selObject.create(ddSuffix  = "DYEstimation",
                         cut       = [op.rng_len(self.ak4BJets)==1,op.rng_len(self.ak8BJets)==0], 
                         weight    = AppliedSF,
                         ddCut     = [op.rng_len(self.ak4BJets)==0,
                                      op.rng_len(self.ak8BJets)==0,
                                      op.rng_len(self.ak8Jets)==0,
                                      (self.tree.event//5)%2==0],
                         ddWeight  = 2*self.ResolvedDYReweighting1b(self.ak4Jets),
                         enable    = enable)
    elif dy_selection:
        # Only return the datadriven selection (eg for Skimmer) #
        selObject.refine(cut    = [op.rng_len(self.ak4BJets)==0,op.rng_len(self.ak8BJets)==0,op.rng_len(self.ak8Jets)==0,(self.tree.event//5)%2==0],
                         weight = 2*self.ResolvedDYReweighting1b(self.ak4Jets))
    else:
        # Only return the SR selection #
        selObject.refine(cut    = [op.rng_len(self.ak4BJets)==1,op.rng_len(self.ak8BJets)==0], 
                         weight = AppliedSF)

    if copy_sel:
        return selObject

def makeExclusiveResolvedTwoBtagsSelection(self,selObject,copy_sel=False,use_dd=True,dy_selection=False):
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
        # Use the datadriven selection in parallel of the SR one #
        enable = "DYEstimation" in self.datadrivenContributions and self.datadrivenContributions["DYEstimation"].usesSample(self.sample, self.sampleCfg) 
        selObject.create(ddSuffix  = "DYEstimation",
                         cut       = [op.rng_len(self.ak4BJets)>=2,op.rng_len(self.ak8BJets)==0],
                         weight    = AppliedSF,
                         ddCut     = [op.rng_len(self.ak4BJets)==0,
                                      op.rng_len(self.ak8BJets)==0,
                                      op.rng_len(self.ak8Jets)==0,
                                      (self.tree.event//5)%2==1],
                         ddWeight  = 2*self.ResolvedDYReweighting2b(self.ak4Jets),
                         enable    = enable)

    elif dy_selection:
        # Only return the datadriven selection (eg for Skimmer) #
        selObject.refine(cut    = [op.rng_len(self.ak4BJets)==0,op.rng_len(self.ak8BJets)==0,op.rng_len(self.ak8Jets)==0,(self.tree.event//5)%2==1],
                         weight = 2*self.ResolvedDYReweighting2b(self.ak4Jets))
    else:
        # Only return the SR selection #
        selObject.refine(cut    = [op.rng_len(self.ak4BJets)>=2,op.rng_len(self.ak8BJets)==0], 
                         weight = AppliedSF)

    if self.args.TTBarCR:
        selObject.selName += "MbbCut150"
        selObject.yieldTitle += " + $M_{bb}>150$"
        selObject.refine(cut = [op.invariant_mass(self.ak4BJets[0].p4,self.ak4BJets[1].p4)>150])
    if copy_sel:
        return selObject

def makeInclusiveBoostedNoBtagSelection(self,selObject,copy_sel=False):
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
    selObject.refine(cut = [op.rng_len(self.ak8BJets) == 0,op.rng_len(self.ak4BJets) == 0])
    if copy_sel:
        return selObject

def makeInclusiveBoostedOneBtagSelection(self,selObject,copy_sel=False,use_dd=True,dy_selection=False):
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
        # Use the datadriven selection in parallel of the SR one #
        enable = "DYEstimation" in self.datadrivenContributions and self.datadrivenContributions["DYEstimation"].usesSample(self.sample, self.sampleCfg) 
        selObject.create(ddSuffix  = "DYEstimation",
                         cut       = [op.rng_len(self.ak8BJets) >= 1],
                         weight    = AppliedSF,
                         ddCut     = [op.rng_len(self.ak8BJets) == 0,op.rng_len(self.ak4BJets) == 0],
                         ddWeight  = self.BoostedDYReweighting1b(self.ak8Jets[0]), 
                         enable    = enable)
    elif dy_selection:
        # Only return the datadriven selection (eg for Skimmer) #
        selObject.refine(cut    = [op.rng_len(self.ak8BJets) == 0,op.rng_len(self.ak4BJets) == 0],
                         weight = self.BoostedDYReweighting1b(self.ak8Jets[0]))
    else:
        # Only return the SR selection #
        selObject.refine(cut    = [op.rng_len(self.ak8BJets)>=1],
                         weight = AppliedSF)
        


    if copy_sel:
        return selObject



def makeDNNOutputNodesSelections(self,selObject,output,suffix=""):
    selBaseObject = copy(selObject)
    selObjDict = {}
    maxIdx = op.rng_max_element_index(output)
    for i,node in enumerate(self.nodes):
        newSelObj = copy(selBaseObject)
        newSelObj.selName += f"{suffix}{node}"
        newSelObj.yieldTitle += f"{suffix} in {node}"
        newSelObj.refine(cut = [maxIdx == op.c_int(i)]) 
        selObjDict[node] = newSelObj
    return selObjDict

