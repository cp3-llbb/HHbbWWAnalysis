import os
import sys
from copy import copy, deepcopy
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

def makeDoubleLeptonSelection(self,baseSel,plot_yield=False): 
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
        # If any dilepton preselected below 12GeV, returns False
    lambda_outZ         = lambda dileptons: op.NOT(op.rng_any(dileptons, lambda dilep : op.in_range(ZMass-10.,op.invariant_mass(dilep[0].p4, dilep[1].p4),ZMass+10.)))
        # If any dilepton preselected within Z peak, returns False
    lambda_inZ          = lambda dilep : op.in_range(ZMass-10.,op.invariant_mass(dilep[0].p4, dilep[1].p4),ZMass+10.)
    lambda_highMllCut   = lambda dilep : op.invariant_mass(dilep[0].p4, dilep[1].p4) < 76.
    
    # SF #
    if self.is_MC:
        ElElLooseSF = lambda dilep : [self.lambda_ttH_doubleElectron_trigSF(dilep)] + self.lambda_ElectronLooseSF(dilep[0]) + self.lambda_ElectronLooseSF(dilep[1])
        MuMuLooseSF = lambda dilep : [self.lambda_ttH_doubleMuon_trigSF(dilep)] + self.lambda_MuonLooseSF(dilep[0]) + self.lambda_MuonLooseSF(dilep[1])
        ElMuLooseSF = lambda dilep : [self.lambda_ttH_electronMuon_trigSF(dilep)] + self.lambda_ElectronLooseSF(dilep[0]) + self.lambda_MuonLooseSF(dilep[1])

        ElElTightSF = lambda dilep : self.lambda_ElectronTightSF(dilep[0]) + self.lambda_ElectronTightSF(dilep[1]) 
        MuMuTightSF = lambda dilep : self.lambda_MuonTightSF(dilep[0]) + self.lambda_MuonTightSF(dilep[1])
        ElMuTightSF = lambda dilep : self.lambda_ElectronTightSF(dilep[0]) + self.lambda_MuonTightSF(dilep[1])
    else:
        ElElLooseSF = lambda dilep : None
        MuMuLooseSF = lambda dilep : None
        ElMuLooseSF = lambda dilep : None

        ElElTightSF = lambda dilep : None
        MuMuTightSF = lambda dilep : None
        ElMuTightSF = lambda dilep : None


        
    #--- Preselection ---#
    selectionDict = {}

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
                ## M_ll(fakeable leptons) < 76
                #ElEl_MllCut.append(lambda_highMllCut(self.OSElElDileptonFakeSel[0]))
                #MuMu_MllCut.append(lambda_highMllCut(self.OSMuMuDileptonFakeSel[0]))
                #ElMu_MllCut.append(lambda_highMllCut(self.OSElMuDileptonFakeSel[0]))
                #ElElFakeSelObject.yieldTitle += " + M_{ll}^{Fakeable} < 76"
                #MuMuFakeSelObject.yieldTitle += " + M_{ll}^{Fakeable} < 76"
                #ElMuFakeSelObject.yieldTitle += " + M_{ll}^{Fakeable} < 76"

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

                    ElElTightSelObject.refine(cut = [lambda_inZ(self.ElElDileptonTightSel[0])])
                    MuMuTightSelObject.refine(cut = [lambda_inZ(self.MuMuDileptonTightSel[0])])

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
                ElElFakeExtrapolationSelObject.refine(cut    = [op.rng_len(self.ElElDileptonFakeExtrapolationSel) >= 1])
                MuMuFakeExtrapolationSelObject.refine(cut    = [op.rng_len(self.MuMuDileptonFakeExtrapolationSel) >= 1])
                ElMuFakeExtrapolationSelObject.refine(cut    = [op.rng_len(self.ElMuDileptonFakeExtrapolationSel) >= 1])

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

def makeAtLeastTwoAk4JetSelection(self,selObject,copy_sel=False,plot_yield=False):
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

def makeAtLeastOneAk8JetSelection(self,selObject,copy_sel=False,plot_yield=False):
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
    AppliedSF = None
    #----- DY estimation from data -----#
    if self.args.DYDataEstimation1Btag or self.args.DYDataEstimation2Btag:
        if self.is_MC:
            print ("Warning : the DY reweighting is not applied on MC, will ignore that step")
        else:
            if "ElEl" in selObject.selName:
                if self.DYReweightingElEl is None:
                    raise RuntimeError('DY reweighting for ElEl is not initialized')
                AppliedDYReweighting = self.DYReweightingElEl(self.OSElElDileptonTightSel[0][0]) # Only on highest PT electron
            elif "MuMu" in selObject.selName:
                if self.DYReweightingMuMu is None:
                    raise RuntimeError('DY reweighting for MuMu is not initialized')
                AppliedDYReweighting = self.DYReweightingMuMu(self.OSMuMuDileptonTightSel[0][0]) # Only on highest PT muon
            else: # No DY reweighting in ElMu (marginal contribution)
                AppliedDYReweighting = None
            AppliedSF = [AppliedDYReweighting] if AppliedDYReweighting is not None else None
                # weight = None works, but not weight = [None]
    if copy_sel:
        selObject = copy(selObject)
    selObject.selName += "ExclusiveResolvedNoBtag"
    selObject.yieldTitle += " + Exclusive Resolved (0 bjet)"
    selObject.refine(cut=[op.rng_len(self.ak4BJets)==0,op.rng_len(self.ak8BJets)==0],weight=AppliedSF)
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
    AppliedSF = None
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
    AppliedSF = None
    if copy_sel:
        selObject = copy(selObject)
    selObject.selName += "ExclusiveResolvedTwoBtags"
    selObject.yieldTitle += " + Exclusive Resolved (2 bjets)"
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
    selObject.refine(cut    = [op.rng_len(self.ak8BJets)>=1],
                     weight = AppliedSF)
    if plot_yield:
        selObject.makeYield(self.yieldPlots)
    if copy_sel:
        return selObject
    

