import os
import sys
from copy import copy
from bamboo import treefunctions as op

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

def makeLeptonSelection(self,baseSel,plot_yield=False): 
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
    lambdaLowPtCut = lambda dilep: op.AND( dilep[0].pt > 15 , dilep[1].pt > 15) # subleading above 15 GeV
    lambdaLeadingPtCut = lambda dilep: op.OR( dilep[0].pt > 25 , dilep[1].pt > 25) # leading above 25 GeV

    # Mll cut lambdas #
    lambda_lowMllCut    = lambda dileptons: op.NOT(op.rng_any(dileptons, lambda dilep : op.invariant_mass(dilep[0].p4, dilep[1].p4)<12.))
        # If any dilepton preselected below 12GeV, returns False
    lambda_outZ         = lambda dileptons: op.NOT(op.rng_any(dileptons, lambda dilep : op.in_range(80.,op.invariant_mass(dilep[0].p4, dilep[1].p4),100.)))
        # If any dilepton preselected within Z peak, returns False
    lambda_inZ          = lambda dilep : op.in_range(80.,op.invariant_mass(dilep[0].p4, dilep[1].p4),100.)
    

    # Loose SF lambdas #
    MuonLooseSF = lambda mu : [self.muLooseId(mu)]
    ElectronLooseSF = lambda el : [self.elLooseId(el) , self.elLooseEff(el), op.switch(el.pt>20 , self.elLooseRecoPtGt20(el) , self.elLooseRecoPtLt20(el))] 

    ElElLooseSF = lambda dilep : [self.ttH_doubleElectron_trigSF]+ElectronLooseSF(dilep[0])+ElectronLooseSF(dilep[1]) if self.is_MC else None
    MuMuLooseSF = lambda dilep : [self.ttH_doubleMuon_trigSF]+MuonLooseSF(dilep[0])+MuonLooseSF(dilep[1]) if self.is_MC else None
    ElMuLooseSF = lambda dilep : [self.ttH_electronMuon_trigSF]+ElectronLooseSF(dilep[0])+MuonLooseSF(dilep[1]) if self.is_MC else None
    
    # Tight SF lambdas #
    ElElTightSF = lambda dilep : [self.elTightMVA(dilep[0]),self.elTightMVA(dilep[1])] if self.is_MC else None
    MuMuTightSF = lambda dilep : [self.muTightMVA(dilep[0]),self.muTightMVA(dilep[1])] if self.is_MC else None
    ElMuTightSF = lambda dilep : [self.elTightMVA(dilep[0]),self.muTightMVA(dilep[1])] if self.is_MC else None

    #--- Preselection ---#
    selectionDict = {}

    if "Preselected" in sel_level:
        ElElPreSelObject = SelectionObject(sel          = baseSel,
                                           selName      = "HasElElPreselected",
                                           yieldTitle   = "OS preselected lepton pair (channel $e^+e^-$)")
        MuMuPreSelObject = SelectionObject(sel          = baseSel,
                                           selName      = "HasMuMuPreselected",
                                           yieldTitle   = "OS preselected lepton pair (channel $\mu^+\mu^-$)")
        ElMuPreSelObject = SelectionObject(sel          = baseSel,
                                           selName      = "HasElMuPreselected",
                                           yieldTitle   = "OS preselected lepton pair (channel $e^{\pm}\mu^{\mp}$)")

        # Selection #
        ElElPreSelObject.refine(cut     = [op.rng_len(self.OSElElDileptonPreSel)>=1],
                                weight  = ElElLooseSF(self.OSElElDileptonPreSel[0]))
        MuMuPreSelObject.refine(cut     = [op.rng_len(self.OSMuMuDileptonPreSel)>=1],
                                weight  = MuMuLooseSF(self.OSMuMuDileptonPreSel[0]))
        ElMuPreSelObject.refine(cut     = [op.rng_len(self.OSElMuDileptonPreSel)>=1],
                                weight  = ElMuLooseSF(self.OSElMuDileptonPreSel[0]))

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
            ElElFakeSelObject.refine(cut = [op.rng_len(self.OSElElDileptonFakeSel)>=1,lambdaLowPtCut(self.OSElElDileptonFakeSel[0]),lambdaLeadingPtCut(self.OSElElDileptonFakeSel[0])])
            MuMuFakeSelObject.refine(cut = [op.rng_len(self.OSMuMuDileptonFakeSel)>=1,lambdaLowPtCut(self.OSMuMuDileptonFakeSel[0]),lambdaLeadingPtCut(self.OSMuMuDileptonFakeSel[0])])
            ElMuFakeSelObject.refine(cut = [op.rng_len(self.OSElMuDileptonFakeSel)>=1,lambdaLowPtCut(self.OSElMuDileptonFakeSel[0]),lambdaLeadingPtCut(self.OSElMuDileptonFakeSel[0])])

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

            ElEl_MllCut = [lambda_lowMllCut(self.OSElElDileptonPreSel)]
            MuMu_MllCut = [lambda_lowMllCut(self.OSMuMuDileptonPreSel)]
            ElMu_MllCut = [lambda_lowMllCut(self.OSElMuDileptonPreSel)]

            if not self.args.NoZVeto: # No Mz cut in OF leptons
                ElElFakeSelObject.selName += "OutZ"
                MuMuFakeSelObject.selName += "OutZ"
                ElElFakeSelObject.yieldTitle += " + Z Veto $|M_{ll}-M_Z|>10$"
                MuMuFakeSelObject.yieldTitle += " + Z Veto $|M_{ll}-M_Z|>10$"
                ElEl_MllCut.append(lambda_outZ(self.OSElElDileptonPreSel))
                MuMu_MllCut.append(lambda_outZ(self.OSMuMuDileptonPreSel))

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

            ElElFakeSelObject.yieldTitle += " + Two tights veto"
            MuMuFakeSelObject.yieldTitle += " + Two tights veto"
            ElMuFakeSelObject.yieldTitle += " + Two tights veto"

                    # TODO : that or N(tight dileptons) <= 1 ?
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
                ElElTightSelObject.refine(cut    = [op.rng_len(self.OSElElDileptonTightSel) >= 1],
                                         weight = ElElTightSF(self.OSElElDileptonTightSel[0]))
                MuMuTightSelObject.refine(cut    = [op.rng_len(self.OSMuMuDileptonTightSel) >= 1],
                                         weight = MuMuTightSF(self.OSMuMuDileptonTightSel[0]))
                ElMuTightSelObject.refine(cut    = [op.rng_len(self.OSElMuDileptonTightSel) >= 1],
                                         weight = ElMuTightSF(self.OSElMuDileptonTightSel[0]))
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

                    ElElTightSelObject.refine(cut = [lambda_inZ(self.OSElElDileptonTightSel[0])])
                    MuMuTightSelObject.refine(cut = [lambda_inZ(self.OSMuMuDileptonTightSel[0])])

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
                ElElFakeExtrapolationSelObject.refine(cut    = [op.rng_len(self.OSElElDileptonTightSel) == 0,
                                                                op.rng_len(self.OSElElDileptonFakeExtrapolationSel) >= 1],
                                                      weight = ElElTightSF(self.OSElElDileptonFakeExtrapolationSel[0])) # TODO : should I ?
                MuMuFakeExtrapolationSelObject.refine(cut    = [op.rng_len(self.OSMuMuDileptonTightSel) == 0,
                                                                op.rng_len(self.OSMuMuDileptonFakeExtrapolationSel) >= 1],
                                                      weight = MuMuTightSF(self.OSMuMuDileptonFakeExtrapolationSel[0])) # TODO : should I ?
                ElMuFakeExtrapolationSelObject.refine(cut    = [op.rng_len(self.OSElMuDileptonTightSel) == 0,
                                                                op.rng_len(self.OSElMuDileptonFakeExtrapolationSel) >= 1],
                                                      weight = ElMuTightSF(self.OSElMuDileptonFakeExtrapolationSel[0])) # TODO : should I ?

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

def makeAk4JetSelection(self,selObject,copy_sel=False,plot_yield=False):
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
    selObject.sel = selObject.sel.refine(selObject.selName,cut=[op.rng_len(self.ak4Jets)>=2])
    if plot_yield:
        selObject.makeYield(self.yieldPlots)
    if copy_sel :
        return selObject 

def makeAk8JetSelection(self,selObject,copy_sel=False,plot_yield=False):
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
    
def makeExclusiveResolvedNoBtagSelection(self,selObject,copy_sel=False,plot_yield=False):
    """
    Produces the exclusive resolved selection with no btagging
    inputs :
        - selObject : SelectionObject class objects (contains bamboo selection, its name and the yield title)
        - copy_sel : wether to copy and return a new selection object built on the one provided, if not will modify it
    Careful : if copy_sel is False, the selObject will be modified
    Selection : no btagged Ak8 jet (aka boosted), all Ak4 jets are non btagged
    """
    if copy_sel:
        selObject = copy(selObject)
    selObject.selName += "ExclusiveResolvedNoBtag"
    selObject.yieldTitle += " + Exclusive Resolved (0 bjet)"
    selObject.refine(cut=[op.rng_len(self.ak4BJets)==0,op.rng_len(self.ak8BJets)==0])
    if plot_yield:
        selObject.makeYield(self.yieldPlots)
    if copy_sel:
        return selObject

def makeExclusiveResolvedOneBtagSelection(self,selObject,copy_sel=False,plot_yield=False):
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
    AppliedSF = [self.DeepJetMediumSF(self.ak4BJets[0])] if self.is_MC else None
    selObject.selName += "ExclusiveResolvedOneBtag"
    selObject.yieldTitle += " + Exclusive Resolved (1 bjet)"
    selObject.refine(cut    = [op.rng_len(self.ak4BJets)==1,op.rng_len(self.ak8BJets)==0],
                     weight = AppliedSF)
    if plot_yield:
        selObject.makeYield(self.yieldPlots)
    if copy_sel:
        return selObject

def makeExclusiveResolvedTwoBtagsSelection(self,selObject,copy_sel=False,plot_yield=False):
    """
    Produces the exclusive resolved selection with no btagging
    inputs :
        - selObject : SelectionObject class objects (contains bamboo selection, its name and the yield title)
        - copy_sel : wether to copy and return a new selection object built on the one provided, if not will modify it
    Careful : if copy_sel is False, the selObject will be modified
    Selection : no btagged Ak8 jet (aka boosted), two Ak4 btagged jets 
    """
    if copy_sel:
        selObject = copy(selObject)
    AppliedSF = [self.DeepJetMediumSF(self.ak4BJets[0]),self.DeepJetMediumSF(self.ak4BJets[1])] if self.is_MC else None
    selObject.selName += "ExclusiveResolvedTwoBtags"
    selObject.yieldTitle += " + Exclusive Resolved (2 bjets)"
    selObject.refine(cut    = [op.rng_len(self.ak4BJets)==2,op.rng_len(self.ak8BJets)==0],
                     weight = AppliedSF)
    if plot_yield:
        selObject.makeYield(self.yieldPlots)
    if copy_sel:
        return selObject

def makeInclusiveBoostedSelection(self,selObject,copy_sel=False,plot_yield=False):
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
    AppliedSF = None # TODO: correct at v7
    selObject.selName += "InclusiveBoosted"
    selObject.yieldTitle += " + Inclusive Boosted"
    selObject.refine(cut    = [op.rng_len(self.ak8BJets)==1],
                     weight = AppliedSF)
    if plot_yield:
        selObject.makeYield(self.yieldPlots)
    if copy_sel:
        return selObject
    

