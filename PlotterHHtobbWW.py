import os
import sys
from copy import copy

from bamboo.analysismodules import HistogramsModule
from bamboo import treefunctions as op
from bamboo.plots import Plot, EquidistantBinning, SummedPlot

sys.path.append('/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/') # Add scripts in this directory -- TODO : make cleaner
from BaseHHtobbWW import BaseNanoHHtobbWW
from plotDef import *
from scalefactorsbbWW import ScaleFactorsbbWW


#===============================================================================================#
#                                       SelectionObject                                         #
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
#                                       PlotterHHtobbWW                                         #
#===============================================================================================#
class PlotterNanoHHtobbWW(BaseNanoHHtobbWW,HistogramsModule):
    """ Plotter module: HH->bbW(->e/µ nu)W(->e/µ nu) histograms from NanoAOD """
    def __init__(self, args):
        super(PlotterNanoHHtobbWW, self).__init__(args)

    #-------------------------------------------------------------------------------------------#
    #                                       Selections                                          #
    #-------------------------------------------------------------------------------------------#
    def makeLeptonSelection(self,baseSel): # We start by leptons so no need to pass selObject
        # Select level #
        level = []
        if self.args.Preselected:
            level.extend(["Preselected"])
        if self.args.Fakeable:
            level.extend(["Preselected","Fakeable"]) # Fakeable needs the preselected 
        if self.args.Tight:
            level.extend(["Preselected","Fakeable","Tight"]) # Tight needs fakeable and preselected
        if self.args.FakeExtrapolation:
            level.extend(["Preselected","Fakeable","FakeExtrapolation"]) # FakeExtrapolation needs fakeable and preselected
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

        if "Preselected" in level:
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
            ElElPreSelObject.makeYield(self.yieldPlots)
            MuMuPreSelObject.makeYield(self.yieldPlots)
            ElMuPreSelObject.makeYield(self.yieldPlots)

            # Record if actual argument #
            if self.args.Preselected:
                selectionDict["Preselected"] = [ElElPreSelObject,MuMuPreSelObject,ElMuPreSelObject]

            #--- Fakeable ---#
            if "Fakeable" in level:
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
                ElElFakeSelObject.makeYield(self.yieldPlots)
                MuMuFakeSelObject.makeYield(self.yieldPlots)
                ElMuFakeSelObject.makeYield(self.yieldPlots)

                # Selection : Mll cut + Z peak exclusion (charge already done in previous selection) in **preselected** leptons # 
                ElElFakeSelObject.selName += "PreMllCutOutZ"
                MuMuFakeSelObject.selName += "PreMllCutOutZ"
                ElMuFakeSelObject.selName += "PreMllCut"

                ElElFakeSelObject.yieldTitle += " + $M_{ll}$ cut"
                MuMuFakeSelObject.yieldTitle += " + $M_{ll}$ cut"
                ElMuFakeSelObject.yieldTitle += " + $M_{ll}$ cut"

                ElElFakeSelObject.refine(cut = [lambda_lowMllCut(self.OSElElDileptonPreSel),lambda_outZ(self.OSElElDileptonPreSel)])
                MuMuFakeSelObject.refine(cut = [lambda_lowMllCut(self.OSMuMuDileptonPreSel),lambda_outZ(self.OSMuMuDileptonPreSel)])
                ElMuFakeSelObject.refine(cut = [lambda_lowMllCut(self.OSElMuDileptonPreSel)]) # No Mz cut in OF leptons

                # Yield #
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
                ElElFakeSelObject.makeYield(self.yieldPlots)
                MuMuFakeSelObject.makeYield(self.yieldPlots)
                ElMuFakeSelObject.makeYield(self.yieldPlots)

                # Record if actual argument #
                if self.args.Fakeable:
                    selectionDict["Fakeable"] = [ElElFakeSelObject,MuMuFakeSelObject,ElMuFakeSelObject]
                    
                #--- Tight ---#
                if "Tight" in level:
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
                    ElElTightSelObject.makeYield(self.yieldPlots)
                    MuMuTightSelObject.makeYield(self.yieldPlots)
                    ElMuTightSelObject.makeYield(self.yieldPlots)
                    
                    # Record if actual argument #
                    if self.args.Tight:
                        selectionDict["Tight"] = [ElElTightSelObject,MuMuTightSelObject,ElMuTightSelObject]

                #--- Fake extrapolated ---#
                if "FakeExtrapolation" in level:
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
                    ElElFakeExtrapolationSelObject.makeYield(self.yieldPlots)
                    MuMuFakeExtrapolationSelObject.makeYield(self.yieldPlots)
                    ElMuFakeExtrapolationSelObject.makeYield(self.yieldPlots)

                    # Record if actual argument #
                    if self.args.FakeExtrapolation:
                        selectionDict["FakeExtrapolation"] = [ElElFakeExtrapolationSelObject,MuMuFakeExtrapolationSelObject,ElMuFakeExtrapolationSelObject]
        # Return # 
        return selectionDict

    def makeAk4JetSelection(self,selObject):
        selObject = copy(selObject)
        selObject.selName += "TwoAk4Jets"
        selObject.yieldTitle += " + Ak4 Jets $\geq 2$"
        selObject.sel = selObject.sel.refine(selObject.selName,cut=[op.rng_len(self.ak4Jets)>=2])
        selObject.makeYield(self.yieldPlots)
        return selObject 

    def makeAk8JetSelection(self,selObject):
        selObject = copy(selObject)
        selObject.selName += "OneAk8Jet"
        selObject.yieldTitle += " + Ak8 Jets $\geq 1$"
        selObject.refine(cut=[op.rng_len(self.ak8Jets)>=1])
        selObject.makeYield(self.yieldPlots)
        return selObject
        
    def makeExclusiveResolvedNoBtagSelection(self,selObject):
        selObject = copy(selObject)
        selObject.selName += "ExclusiveResolvedNoBtag"
        selObject.yieldTitle += " + Exclusive Resolved (0 bjet)"
        selObject.refine(cut=[op.rng_len(self.ak4BJets)==0,op.rng_len(self.ak8BJets)==0])
        selObject.makeYield(self.yieldPlots)
        return selObject

    def makeExclusiveResolvedOneBtagSelection(self,selObject):
        selObject = copy(selObject)
        AppliedSF = [self.DeepJetMediumSF(self.ak4BJets[0])] if self.is_MC else None
        selObject.selName += "ExclusiveResolvedOneBtag"
        selObject.yieldTitle += " + Exclusive Resolved (1 bjet)"
        selObject.refine(cut    = [op.rng_len(self.ak4BJets)==1,op.rng_len(self.ak8BJets)==0],
                         weight = AppliedSF)
        selObject.makeYield(self.yieldPlots)
        return selObject

    def makeExclusiveResolvedTwoBtagsSelection(self,selObject):
        selObject = copy(selObject)
        AppliedSF = [self.DeepJetMediumSF(self.ak4BJets[0]),self.DeepJetMediumSF(self.ak4BJets[1])] if self.is_MC else None
        selObject.selName += "ExclusiveResolvedTwoBtags"
        selObject.yieldTitle += " + Exclusive Resolved (2 bjets)"
        selObject.refine(cut    = [op.rng_len(self.ak4BJets)==2,op.rng_len(self.ak8BJets)==0],
                         weight = AppliedSF)
        selObject.makeYield(self.yieldPlots)
        return selObject

    def makeInclusiveBoostedSelection(self,selObject):
        selObject = copy(selObject)
        AppliedSF = None # TODO: correct at v7
        selObject.selName += "InclusiveBoosted"
        selObject.yieldTitle += " + Inclusive Boosted"
        selObject.refine(cut    = [op.rng_len(self.ak8BJets)>=1],
                         weight = AppliedSF)
        selObject.makeYield(self.yieldPlots)
        return selObject
        
    #-------------------------------------------------------------------------------------------#
    #                                       Plots                                               #
    #-------------------------------------------------------------------------------------------#

    def definePlots(self, t, noSel, sample=None, sampleCfg=None): 
        noSel = super(PlotterNanoHHtobbWW,self).prepareObjects(t, noSel, sample, sampleCfg)

        plots = []

        era = sampleCfg['era']

        self.is_MC = self.isMC(sample)

        self.yieldPlots = makeYieldPlots()

        #############################################################################
        #                             Scalefactors                                  #
        #############################################################################
        if self.is_MC:
            # Initialize scalefactors class #
            SF = ScaleFactorsbbWW()    
            ####  Muons ####
            self.muLooseId = SF.get_scalefactor("lepton", 'muon_loose_{}'.format(era), combine="weight", systName="mu_loose")
            self.muTightMVA = SF.get_scalefactor("lepton", 'muon_tightMVA_{}'.format(era), combine="weight", systName="mu_tightmva")

            ####  Electrons ####
            self.elLooseRecoPtLt20 = SF.get_scalefactor("lepton", ('electron_loosereco_{}'.format(era) , 'electron_loosereco_ptgt20'), combine="weight", systName="el_looserecoptlt20")
            self.elLooseRecoPtGt20 = SF.get_scalefactor("lepton", ('electron_loosereco_{}'.format(era) , 'electron_loosereco_ptlt20'), combine="weight", systName="el_looserecoptgt20")
            self.elLooseId = SF.get_scalefactor("lepton", 'electron_looseid_{}'.format(era) , combine="weight", systName="el_looseid")
            self.elLooseEff = SF.get_scalefactor("lepton", 'electron_looseeff_{}'.format(era) , combine="weight", systName="el_looseeff")
            self.elTightMVA = SF.get_scalefactor("lepton", 'electron_tightMVA_{}'.format(era) , combine="weight", systName="el_tightmva")

            #### Triggers (from TTH) ####
            self.ttH_doubleMuon_trigSF = op.systematic(op.c_float(1.010), name="ttH_doubleMuon_trigSF", up=op.c_float(1.020), down=op.c_float(1.000))
            self.ttH_doubleElectron_trigSF = op.systematic(op.c_float(1.020), name="ttH_doubleElectron_trigSF", up=op.c_float(1.040), down=op.c_float(1.000))
            self.ttH_electronMuon_trigSF = op.systematic(op.c_float(1.020), name="ttH_electronMuon_trigSF", up=op.c_float(1.030), down=op.c_float(1.010))
            #### Ak4 Btagging ####
            DeepJetTag_discriVar = {"BTagDiscri": lambda j : j.btagDeepFlavB}
            self.DeepJetMediumSF = SF.get_scalefactor("jet", ("btag_"+era+"_"+self.sfTag, "DeepJet_medium"), additionalVariables=DeepJetTag_discriVar, systName="deepjet") # For RESOLVED
            #### Ak8 Btagging ####
#            DeepCSVTag_discriVar = {"BTagDiscri": lambda j : j.btagDeepB}
#            self.DeepCSVMediumSF = SF.get_scalefactor("jet", ("subjet_btag_"+era+"_"+self.sfTag, "DeepCSV_medium"), additionalVariables=DeepCSVTag_discriVar, systName="deepcsv") # For BOOSTED (btag on subjet)
                # TODO : add Boosted SF (need NanoAODv7)

        #############################################################################
        #                           Selection and plot                              #
        #############################################################################
        #----- Dileptons -----#
        selObjectDict = self.makeLeptonSelection(noSel)
        # selObjectDict : keys -> level (str)
        #                 values -> [ElEl,MuMu,ElMu] x Selection object
        for selectionType, selectionList in selObjectDict.items():
            print ("... Processing %s lepton type"%selectionType)
            #----- Select correct dilepton -----#
            if selectionType == "Preselected":  
                OSElElDilepton = self.OSElElDileptonPreSel
                OSMuMuDilepton = self.OSMuMuDileptonPreSel
                OSElMuDilepton = self.OSElMuDileptonPreSel
            elif selectionType == "Fakeable":
                OSElElDilepton = self.OSElElDileptonFakeSel
                OSMuMuDilepton = self.OSMuMuDileptonFakeSel
                OSElMuDilepton = self.OSElMuDileptonFakeSel
            elif selectionType == "Tight":
                OSElElDilepton = self.OSElElDileptonTightSel
                OSMuMuDilepton = self.OSMuMuDileptonTightSel
                OSElMuDilepton = self.OSElMuDileptonTightSel
            elif selectionType == "FakeExtrapolation":
                OSElElDilepton = self.OSElElDileptonFakeExtrapolationSel
                OSMuMuDilepton = self.OSMuMuDileptonFakeExtrapolationSel
                OSElMuDilepton = self.OSElMuDileptonFakeExtrapolationSel

            #----- Separate selections ------#
            ElElSelObjectDilepton = selectionList[0]
            MuMuSelObjectDilepton = selectionList[1]
            ElMuSelObjectDilepton = selectionList[2]

            if not self.args.OnlyYield:
                #----- Trigger plots -----#
                plots.extend(triggerPlots(sel         = ElElSelObjectDilepton.sel,
                                          triggerDict = self.triggersPerPrimaryDataset,
                                          suffix      = ElElSelObjectDilepton.selName,
                                          channel     = "ElEl"))
                plots.extend(triggerPlots(sel         = MuMuSelObjectDilepton.sel,
                                          triggerDict = self.triggersPerPrimaryDataset,
                                          suffix      = MuMuSelObjectDilepton.selName,
                                          channel     = "MuMu"))
                plots.extend(triggerPlots(sel         = ElMuSelObjectDilepton.sel,
                                          triggerDict = self.triggersPerPrimaryDataset,
                                          suffix      = ElMuSelObjectDilepton.selName,
                                          channel     = "ElMu"))

                #----- Lepton plots -----#
                # Dilepton channel plots #
                plots.extend(channelPlot(sel        = ElElSelObjectDilepton.sel,
                                         DilepElEl  = OSElElDilepton,
                                         DilepMuMu  = OSMuMuDilepton,
                                         DilepElMu  = OSElMuDilepton,
                                         suffix     = ElElSelObjectDilepton.selName,
                                         channel    = "ElEl"))
                plots.extend(channelPlot(sel        = MuMuSelObjectDilepton.sel,
                                         DilepElEl  = OSElElDilepton,
                                         DilepMuMu  = OSMuMuDilepton,
                                         DilepElMu  = OSElMuDilepton,
                                         suffix     = MuMuSelObjectDilepton.selName,
                                         channel    = "MuMu"))
                plots.extend(channelPlot(sel        = ElMuSelObjectDilepton.sel,
                                         DilepElEl  = OSElElDilepton,
                                         DilepMuMu  = OSMuMuDilepton,
                                         DilepElMu  = OSElMuDilepton,
                                         suffix     = ElMuSelObjectDilepton.selName,
                                         channel    = "ElMu"))
        
            #----- Ak4 jets selection -----#
            ElElSelObjectDileptonAk4Jets = self.makeAk4JetSelection(ElElSelObjectDilepton)
            MuMuSelObjectDileptonAk4Jets = self.makeAk4JetSelection(MuMuSelObjectDilepton)
            ElMuSelObjectDileptonAk4Jets = self.makeAk4JetSelection(ElMuSelObjectDilepton)

            # Jet and lepton plots #
            ChannelDictList = []
            if self.args.Ak4 and not self.args.OnlyYield:
                ChannelDictList.append({'channel':'ElEl','sel':ElElSelObjectDileptonAk4Jets.sel,'dilepton':OSElElDilepton[0],'leadjet':self.ak4Jets[0],'subleadjet':self.ak4Jets[1],'lead_is_b':False,'sublead_is_b':False,'suffix':ElElSelObjectDileptonAk4Jets.selName})
                ChannelDictList.append({'channel':'MuMu','sel':MuMuSelObjectDileptonAk4Jets.sel,'dilepton':OSMuMuDilepton[0],'leadjet':self.ak4Jets[0],'subleadjet':self.ak4Jets[1],'lead_is_b':False,'sublead_is_b':False,'suffix':MuMuSelObjectDileptonAk4Jets.selName})
                ChannelDictList.append({'channel':'ElMu','sel':ElMuSelObjectDileptonAk4Jets.sel,'dilepton':OSElMuDilepton[0],'leadjet':self.ak4Jets[0],'subleadjet':self.ak4Jets[1],'lead_is_b':False,'sublead_is_b':False,'suffix':ElMuSelObjectDileptonAk4Jets.selName})

            LeptonKeys = ['channel','sel','dilepton','suffix']
            JetKeys    = ['channel','sel','leadjet','subleadjet','lead_is_b','sublead_is_b','suffix']
            commonItems = ['channel','sel','suffix']
            JetsN = {'objName':'Ak4Jets','objCont':self.ak4Jets,'Nmax':10,'xTitle':'N(Ak4 jets)'}
                                                        
            for channelDict in ChannelDictList:
                # Dilepton #
                plots.extend(makeDileptonPlots(**{k:channelDict[k] for k in LeptonKeys},isMC=self.is_MC))
                # Number of jets #
                plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**JetsN))
                # Ak4 Jets #
                plots.extend(makeAk4JetsPlots(**{k:channelDict[k] for k in JetKeys}))
                # MET #
                plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))

            ##### Ak8 jets selection #####
            ElElSelObjectDileptonAk8Jets = self.makeAk8JetSelection(ElElSelObjectDilepton)
            MuMuSelObjectDileptonAk8Jets = self.makeAk8JetSelection(MuMuSelObjectDilepton)
            ElMuSelObjectDileptonAk8Jets = self.makeAk8JetSelection(ElMuSelObjectDilepton)

            # Fatjets plots #
            ChannelDictList = []
            if self.args.Ak8 and not self.args.OnlyYield:
                ChannelDictList.append({'channel':'ElEl','sel':ElElSelObjectDileptonAk8Jets.sel,'dilepton':OSElElDilepton[0],'fatjet':self.ak8Jets[0],'suffix':ElElSelObjectDileptonAk8Jets.selName})
                ChannelDictList.append({'channel':'MuMu','sel':MuMuSelObjectDileptonAk8Jets.sel,'dilepton':OSMuMuDilepton[0],'fatjet':self.ak8Jets[0],'suffix':MuMuSelObjectDileptonAk8Jets.selName})
                ChannelDictList.append({'channel':'ElMu','sel':ElMuSelObjectDileptonAk8Jets.sel,'dilepton':OSElMuDilepton[0],'fatjet':self.ak8Jets[0],'suffix':ElMuSelObjectDileptonAk8Jets.selName})

            FatJetKeys = ['channel','sel','fatjet','suffix']
            FatJetsN = {'objName':'Ak8Jets','objCont':self.ak8Jets,'Nmax':5,'xTitle':'N(Ak8 jets)'}
            
            for channelDict in ChannelDictList:
                # Dilepton #
                plots.extend(makeDileptonPlots(**{k:channelDict[k] for k in LeptonKeys},isMC=self.is_MC))
                # Number of jets #
                plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**FatJetsN))
                # Ak8 Jets #
                plots.extend(makeAk8JetsPlots(**{k:channelDict[k] for k in FatJetKeys}))
                # MET #
                plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))

                         
            #----- Resolved selection : 0 Btag -----#
            ElElSelObjectDileptonAk4JetsExclusiveResolvedNoBtag = self.makeExclusiveResolvedNoBtagSelection(ElElSelObjectDileptonAk4Jets)
            MuMuSelObjectDileptonAk4JetsExclusiveResolvedNoBtag = self.makeExclusiveResolvedNoBtagSelection(MuMuSelObjectDileptonAk4Jets)
            ElMuSelObjectDileptonAk4JetsExclusiveResolvedNoBtag = self.makeExclusiveResolvedNoBtagSelection(ElMuSelObjectDileptonAk4Jets)

            #----  Resolved selection : 1 Btag  -----#
            ElElSelObjectDileptonAk4JetsExclusiveResolvedOneBtag = self.makeExclusiveResolvedOneBtagSelection(ElElSelObjectDileptonAk4Jets)
            MuMuSelObjectDileptonAk4JetsExclusiveResolvedOneBtag = self.makeExclusiveResolvedOneBtagSelection(MuMuSelObjectDileptonAk4Jets)
            ElMuSelObjectDileptonAk4JetsExclusiveResolvedOneBtag = self.makeExclusiveResolvedOneBtagSelection(ElMuSelObjectDileptonAk4Jets)

            #----- Resolved selection : 2 Btags -----#
            ElElSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags = self.makeExclusiveResolvedTwoBtagsSelection(ElElSelObjectDileptonAk4Jets)
            MuMuSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags = self.makeExclusiveResolvedTwoBtagsSelection(MuMuSelObjectDileptonAk4Jets)
            ElMuSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags = self.makeExclusiveResolvedTwoBtagsSelection(ElMuSelObjectDileptonAk4Jets)

            # Lepton + jet Plots #
            ChannelDictList = []
            if self.args.Resolved0Btag and not self.args.OnlyYield:
               ChannelDictList.append({'channel':'ElEl','sel':ElElSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.sel,'dilepton':OSElElDilepton[0],'leadjet':self.ak4LightJets[0],'subleadjet':self.ak4LightJets[1],'lead_is_b':False,'sublead_is_b':False,'suffix':ElElSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.selName})
               ChannelDictList.append({'channel':'MuMu','sel':MuMuSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.sel,'dilepton':OSMuMuDilepton[0],'leadjet':self.ak4LightJets[0],'subleadjet':self.ak4LightJets[1],'lead_is_b':False,'sublead_is_b':False,'suffix':MuMuSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.selName})
               ChannelDictList.append({'channel':'ElMu','sel':ElMuSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.sel,'dilepton':OSElMuDilepton[0],'leadjet':self.ak4LightJets[0],'subleadjet':self.ak4LightJets[1],'lead_is_b':False,'sublead_is_b':False,'suffix':ElMuSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.selName})
            if self.args.Resolved1Btag and not self.args.OnlyYield:
               ChannelDictList.append({'channel':'ElEl','sel':ElElSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.sel,'dilepton':OSElElDilepton[0],'leadjet':self.ak4BJets[0],'subleadjet':self.ak4LightJets[0],'lead_is_b':True,'sublead_is_b':False,'suffix':ElElSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.selName})
               ChannelDictList.append({'channel':'MuMu','sel':MuMuSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.sel,'dilepton':OSMuMuDilepton[0],'leadjet':self.ak4BJets[0],'subleadjet':self.ak4LightJets[0],'lead_is_b':True,'sublead_is_b':False,'suffix':MuMuSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.selName})
               ChannelDictList.append({'channel':'ElMu','sel':ElMuSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.sel,'dilepton':OSElMuDilepton[0],'leadjet':self.ak4BJets[0],'subleadjet':self.ak4LightJets[0],'lead_is_b':True,'sublead_is_b':False,'suffix':ElMuSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.selName})
            if self.args.Resolved2Btag and not self.args.OnlyYield:
               ChannelDictList.append({'channel':'ElEl','sel':ElElSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.sel,'dilepton':OSElElDilepton[0],'leadjet':self.ak4BJets[0],'subleadjet':self.ak4BJets[1],'lead_is_b':True,'sublead_is_b':True,'suffix':ElElSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.selName})
               ChannelDictList.append({'channel':'MuMu','sel':MuMuSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.sel,'dilepton':OSMuMuDilepton[0],'leadjet':self.ak4BJets[0],'subleadjet':self.ak4BJets[1],'lead_is_b':True,'sublead_is_b':True,'suffix':MuMuSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.selName})
               ChannelDictList.append({'channel':'ElMu','sel':ElMuSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.sel,'dilepton':OSElMuDilepton[0],'leadjet':self.ak4BJets[0],'subleadjet':self.ak4BJets[1],'lead_is_b':True,'sublead_is_b':True,'suffix':ElMuSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.selName})

            ResolvedJetsN = {'objName':'Ak4BJets','objCont':self.ak4BJets,'Nmax':5,'xTitle':'N(Ak4 Bjets)'}
    
            for channelDict in ChannelDictList:
                # Dilepton #
                plots.extend(makeDileptonPlots(**{k:channelDict[k] for k in LeptonKeys},isMC=self.is_MC))
                # Number of jets #
                plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**ResolvedJetsN))
                # Ak4 Jets #
                plots.extend(makeAk4JetsPlots(**{k:channelDict[k] for k in JetKeys}))
                # MET #
                plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))
        
            #----- Boosted selection -----#
            ElElSelObjectDileptonAk8JetsInclusiveBoosted = self.makeInclusiveBoostedSelection(ElElSelObjectDileptonAk8Jets)
            MuMuSelObjectDileptonAk8JetsInclusiveBoosted = self.makeInclusiveBoostedSelection(MuMuSelObjectDileptonAk8Jets)
            ElMuSelObjectDileptonAk8JetsInclusiveBoosted = self.makeInclusiveBoostedSelection(ElMuSelObjectDileptonAk8Jets)

            # Lepton + jet Plots #
            ChannelDictList = []
            if self.args.Boosted and not self.args.OnlyYield:
                ChannelDictList.append({'channel':'ElEl','sel':ElElSelObjectDileptonAk8JetsInclusiveBoosted.sel,'dilepton':OSElElDilepton[0],'fatjet':self.ak8BJets[0],'suffix':ElElSelObjectDileptonAk8JetsInclusiveBoosted.selName})
                ChannelDictList.append({'channel':'MuMu','sel':MuMuSelObjectDileptonAk8JetsInclusiveBoosted.sel,'dilepton':OSMuMuDilepton[0],'fatjet':self.ak8BJets[0],'suffix':MuMuSelObjectDileptonAk8JetsInclusiveBoosted.selName})
                ChannelDictList.append({'channel':'ElMu','sel':ElMuSelObjectDileptonAk8JetsInclusiveBoosted.sel,'dilepton':OSElMuDilepton[0],'fatjet':self.ak8BJets[0],'suffix':ElMuSelObjectDileptonAk8JetsInclusiveBoosted.selName})

            BoostedJetsN = {'objName':'Ak8BJets','objCont':self.ak8BJets,'Nmax':5,'xTitle':'N(Ak8 Bjets)'}
            
            for channelDict in ChannelDictList:
                # Dilepton #
                plots.extend(makeDileptonPlots(**{k:channelDict[k] for k in LeptonKeys},isMC=self.is_MC))
                # Number of jets #
                plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**BoostedJetsN))
                # Ak8 Jets #
                plots.extend(makeAk8JetsPlots(**{k:channelDict[k] for k in FatJetKeys}))
                # MET #
                plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))
            
            #----- High-level combinations -----#
            ChannelDictList = []
            # Resolved No Btag #
            if self.args.Resolved0Btag and not self.args.OnlyYield:
                ChannelDictList.append({'channel': 'ElEl','met': self.corrMET,'l1':OSElElDilepton[0][0],'l2':OSElElDilepton[0][1],'j1':self.ak4LightJets[0],'j2':self.ak4LightJets[1],'sel':ElElSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.sel,'suffix':ElElSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.selName})
                ChannelDictList.append({'channel': 'MuMu','met': self.corrMET,'l1':OSMuMuDilepton[0][0],'l2':OSMuMuDilepton[0][1],'j1':self.ak4LightJets[0],'j2':self.ak4LightJets[1],'sel':MuMuSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.sel,'suffix':MuMuSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.selName})
                ChannelDictList.append({'channel': 'ElMu','met': self.corrMET,'l1':OSElMuDilepton[0][0],'l2':OSElMuDilepton[0][1],'j1':self.ak4LightJets[0],'j2':self.ak4LightJets[1],'sel':ElMuSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.sel,'suffix':ElMuSelObjectDileptonAk4JetsExclusiveResolvedNoBtag.selName})
            # Resolved One Btag #
            if self.args.Resolved1Btag and not self.args.OnlyYield:
                ChannelDictList.append({'channel': 'ElEl','met': self.corrMET,'l1':OSElElDilepton[0][0],'l2':OSElElDilepton[0][1],'j1':self.ak4BJets[0],'j2':self.ak4LightJets[0],'sel':ElElSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.sel,'suffix':ElElSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.selName})
                ChannelDictList.append({'channel': 'MuMu','met': self.corrMET,'l1':OSMuMuDilepton[0][0],'l2':OSMuMuDilepton[0][1],'j1':self.ak4BJets[0],'j2':self.ak4LightJets[0],'sel':MuMuSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.sel,'suffix':MuMuSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.selName})
                ChannelDictList.append({'channel': 'ElMu','met': self.corrMET,'l1':OSElMuDilepton[0][0],'l2':OSElMuDilepton[0][1],'j1':self.ak4BJets[0],'j2':self.ak4LightJets[0],'sel':ElMuSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.sel,'suffix':ElMuSelObjectDileptonAk4JetsExclusiveResolvedOneBtag.selName})
            # Resolved Two Btags  #
            if self.args.Resolved2Btag and not self.args.OnlyYield:
                ChannelDictList.append({'channel': 'ElEl','met': self.corrMET,'l1':OSElElDilepton[0][0],'l2':OSElElDilepton[0][1],'j1':self.ak4BJets[0],'j2':self.ak4BJets[1],'sel':ElElSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.sel,'suffix':ElElSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.selName})
                ChannelDictList.append({'channel': 'MuMu','met': self.corrMET,'l1':OSMuMuDilepton[0][0],'l2':OSMuMuDilepton[0][1],'j1':self.ak4BJets[0],'j2':self.ak4BJets[1],'sel':MuMuSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.sel,'suffix':MuMuSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.selName})
                ChannelDictList.append({'channel': 'ElMu','met': self.corrMET,'l1':OSElMuDilepton[0][0],'l2':OSElMuDilepton[0][1],'j1':self.ak4BJets[0],'j2':self.ak4BJets[1],'sel':ElMuSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.sel,'suffix':ElMuSelObjectDileptonAk4JetsExclusiveResolvedTwoBtags.selName})
            # Boosted #
            if self.args.Boosted and not self.args.OnlyYield:
                ChannelDictList.append({'channel': 'ElEl','met': self.corrMET,'l1':OSElElDilepton[0][0],'l2':OSElElDilepton[0][1],'j1':self.ak8BJets[0].subJet1,'j2':self.ak8BJets[0].subJet2,'sel':ElElSelObjectDileptonAk8JetsInclusiveBoosted.sel,'suffix':ElElSelObjectDileptonAk8JetsInclusiveBoosted.selName})
                ChannelDictList.append({'channel': 'MuMu','met': self.corrMET,'l1':OSMuMuDilepton[0][0],'l2':OSMuMuDilepton[0][1],'j1':self.ak8BJets[0].subJet1,'j2':self.ak8BJets[0].subJet2,'sel':MuMuSelObjectDileptonAk8JetsInclusiveBoosted.sel,'suffix':MuMuSelObjectDileptonAk8JetsInclusiveBoosted.selName})
                ChannelDictList.append({'channel': 'ElMu','met': self.corrMET,'l1':OSElMuDilepton[0][0],'l2':OSElMuDilepton[0][1],'j1':self.ak8BJets[0].subJet1,'j2':self.ak8BJets[0].subJet2,'sel':ElMuSelObjectDileptonAk8JetsInclusiveBoosted.sel,'suffix':ElMuSelObjectDileptonAk8JetsInclusiveBoosted.selName})

            for channelDict in ChannelDictList:
                plots.extend(makeHighLevelQuantities(**channelDict))

        #----- Add the Yield plots -----#
        plots.extend(self.yieldPlots.returnPlots())

        return plots

