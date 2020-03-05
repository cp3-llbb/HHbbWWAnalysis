import os
import sys

from bamboo.analysismodules import HistogramsModule
from bamboo import treefunctions as op
from bamboo.plots import Plot, EquidistantBinning, SummedPlot

sys.path.append('/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/') # Add scripts in this directory -- TODO : make cleaner
from BaseHHtobbWW import BaseNanoHHtobbWW
from plotDef import *
from scalefactorsbbWW import ScaleFactorsbbWW

#===============================================================================================#
#                                 PlotterHHtobbWW                                               #
#===============================================================================================#
class PlotterNanoHHtobbWW(BaseNanoHHtobbWW,HistogramsModule):
    """ Plotter module: HH->bbW(->e/µ nu)W(->e/µ nu) histograms from NanoAOD """
    def __init__(self, args):
        super(PlotterNanoHHtobbWW, self).__init__(args)

    def definePlots(self, t, noSel, sample=None, sampleCfg=None): 
        noSel = super(PlotterNanoHHtobbWW,self).prepareObjects(t, noSel, sample, sampleCfg)

        plots = []

        era = sampleCfg['era']

        isMC = self.isMC(sample)

        yieldPlots = makeYieldPlots()

        #############################################################################
        #                             Scalefactors                                  #
        #############################################################################
        if isMC:
            # Initialize scalefactors class #
            SF = ScaleFactorsbbWW()    
            ####  Muons ####
            muTightIDSF = SF.get_scalefactor("lepton", ("muon_{0}_{1}".format(era, self.sfTag), "id_tight"), combine="weight", systName="muid")
            muTightISOSF = SF.get_scalefactor("lepton", ("muon_{0}_{1}".format(era, self.sfTag), "iso_tight_id_medium"), combine="weight", systName="muiso")
            ####  Electrons ####
            elTightIDSF = SF.get_scalefactor("lepton", ("electron_{0}_{1}".format(era,self.sfTag), "id_tight"), systName="elid")
            #### Triggers (from TTH) ####
            ttH_doubleMuon_trigSF = op.systematic(op.c_float(1.010), name="ttH_doubleMuon_trigSF", up=op.c_float(1.020), down=op.c_float(1.000))
            ttH_doubleElectron_trigSF = op.systematic(op.c_float(1.020), name="ttH_doubleElectron_trigSF", up=op.c_float(1.040), down=op.c_float(1.000))
            ttH_electronMuon_trigSF = op.systematic(op.c_float(1.020), name="ttH_electronMuon_trigSF", up=op.c_float(1.030), down=op.c_float(1.010))
            #### Ak4 Btagging ####
            DeepJetTag_discriVar = {"BTagDiscri": lambda j : j.btagDeepFlavB}
            DeepJetMediumSF = SF.get_scalefactor("jet", ("btag_"+era+"_"+self.sfTag, "DeepJet_medium"), additionalVariables=DeepJetTag_discriVar, systName="deepjet") # For RESOLVED
            #### Ak8 Btagging ####
            DeepCSVTag_discriVar = {"BTagDiscri": lambda j : j.btagDeepB}
            DeepCSVMediumSF = SF.get_scalefactor("jet", ("subjet_btag_"+era+"_"+self.sfTag, "DeepCSV_medium"), additionalVariables=DeepCSVTag_discriVar, systName="deepcsv") # For BOOSTED (btag on subjet)

        #############################################################################
        #                               Dilepton                                    #
        #############################################################################
        # Make dilepton fakeable pairs wrt to channel #
        lambdaOS  = lambda l1,l2 : l1.charge != l2.charge           # Opposite sign charges
        lambdaPTCut = lambda l1,l2: op.OR( l1.pt > 25 , l2.pt > 25) # PT cut 
        lambdaFakeableDilepton = lambda l1, l2 : op.AND( lambdaOS(l1,l2) , lambdaPTCut(l1,l2) )
        OsElElFakeable = op.combine(self.electronsFakeSel, N=2, pred=lambdaFakeableDilepton)
        OsMuMuFakeable = op.combine(self.muonsFakeSel, N=2, pred=lambdaFakeableDilepton)
        OsElMuFakeable = op.combine((self.electronsFakeSel, self.muonsFakeSel), pred=lambdaFakeableDilepton)

        # Applied SF #
        ElElTriggerSFApplied = [ttH_doubleElectron_trigSF] if isMC else None
        MuMuTriggerSFApplied = [ttH_doubleMuon_trigSF] if isMC else None
        ElMuTriggerSFApplied = [ttH_electronMuon_trigSF] if isMC else None

        # At least 2 fakeables leptons #
        ElElHasFakeableDilepton = noSel.refine("ElElHasFakeableDilepton",
                                        cut = [op.rng_len(OsElElFakeable)>0],
                                        weight = ElElTriggerSFApplied)
        MuMuHasFakeableDilepton = noSel.refine("MuMuHasFakeableDilepton",
                                        cut = [op.rng_len(OsMuMuFakeable)>0],
                                        weight = MuMuTriggerSFApplied)
        ElMuHasFakeableDilepton = noSel.refine("ElMuHasFakeableDilepton",
                                        cut = [op.rng_len(OsElMuFakeable)>0],
                                        weight = ElMuTriggerSFApplied)

        # Yield #
        plots.append(yieldPlots.addYield(ElElHasFakeableDilepton,"ElElHasFakeableDilepton","OS fakeable leptons (channel $e^+e^-$)"))
        plots.append(yieldPlots.addYield(MuMuHasFakeableDilepton,"MuMuHasFakeableDilepton","OS fakeable leptons (channel $\mu^+\mu^-$)"))
        plots.append(yieldPlots.addYield(ElMuHasFakeableDilepton,"ElMuHasFakeableDilepton","OS fakeable leptons (channel $e^{\pm}\mu^{\mp}$)"))
    
        # Mll cut + Z peak exclusion (charge already done in previous selection) in **preselected** leptons # 
        OsElElPresel = op.combine(self.electronsPreSel, N=2, pred=lambdaOS)
        OsMuMuPresel = op.combine(self.muonsPreSel, N=2, pred=lambdaOS)
        OsElMuPresel = op.combine((self.electronsPreSel, self.muonsPreSel), pred=lambdaOS)
        lambda_lowMllCut    = lambda dileptons: op.NOT(op.rng_any(dileptons, lambda dilep : op.invariant_mass(dilep[0].p4, dilep[1].p4)<12.))
            # If any dilepton preselected below 12GeV, returns False
        lambda_outZ         = lambda dileptons: op.NOT(op.rng_any(dileptons, lambda dilep : op.in_range(80.,op.invariant_mass(dilep[0].p4, dilep[1].p4),100.)))
            # If any dilepton preselected within Z peak, returns False

        ElElHasFakeableDileptonPreMllCutOutZ = ElElHasFakeableDilepton.refine("ElElHasFakeableDileptonPreMllCutOutZ", 
                                cut=[lambda_lowMllCut(OsElElPresel), lambda_outZ(OsElElPresel)])
        MuMuHasFakeableDileptonPreMllCutOutZ = MuMuHasFakeableDilepton.refine("MuMuHasFakeableDileptonPreMllCutOutZ", 
                                cut=[lambda_lowMllCut(OsMuMuPresel),lambda_outZ(OsMuMuPresel)])
        ElMuHasFakeableDileptonPreMllCut     = ElMuHasFakeableDilepton.refine("ElMuHasFakeableDileptonPreMllCut", 
                                cut=[lambda_lowMllCut(OsElMuPresel)])

        # Yield #
        plots.append(yieldPlots.addYield(ElElHasFakeableDileptonPreMllCutOutZ,"ElElHasFakeableDileptonPreMllCutOutZ","OS fakeable leptons + $M_{ll}$ cut (channel $e^+e^-$)"))
        plots.append(yieldPlots.addYield(MuMuHasFakeableDileptonPreMllCutOutZ,"MuMuHasFakeableDileptonPreMllCutOutZ","OS fakeable leptons + $M_{ll}$ cut (channel $\mu^+\mu^-$)"))
        plots.append(yieldPlots.addYield(ElMuHasFakeableDileptonPreMllCut,"ElMuHasFakeableDileptonPreMllCut","OS fakeable leptons + $M_{ll}$ cut (channel $e^{\pm}\mu^{\mp}$)"))

        # Tight cut veto : no more than two tight leptons in the event #
        ElElHasFakeableDileptonPreMllCutOutZTightVeto = ElElHasFakeableDileptonPreMllCutOutZ.refine("ElElHasFakeableDileptonPreMllCutOutZTightVeto", 
                                cut=[op.rng_len(self.muonsTightSel)+op.rng_len(self.electronsTightSel) <= 2])
        MuMuHasFakeableDileptonPreMllCutOutZTightVeto = MuMuHasFakeableDileptonPreMllCutOutZ.refine("MuMuHasFakeableDileptonPreMllCutOutZTightVeto", 
                                cut=[op.rng_len(self.muonsTightSel)+op.rng_len(self.electronsTightSel) <= 2])
        ElMuHasFakeableDileptonPreMllCutTightVeto     = ElMuHasFakeableDileptonPreMllCut.refine("ElMuHasFakeableDileptonPreMllCutTightVeto", 
                                cut=[op.rng_len(self.muonsTightSel)+op.rng_len(self.electronsTightSel) <= 2])

        # Yield #
        plots.append(yieldPlots.addYield(ElElHasFakeableDileptonPreMllCutOutZTightVeto,"ElElHasFakeableDileptonPreMllCutOutZTightVeto","OS fakeable leptons + $M_{ll}$ cut + Two tights veto (channel $e^+e^-$)"))
        plots.append(yieldPlots.addYield(MuMuHasFakeableDileptonPreMllCutOutZTightVeto,"MuMuHasFakeableDileptonPreMllCutOutZTightVeto","OS fakeable leptons + $M_{ll}$ cut + Two tights veto (channel $\mu^+\mu^-$)"))
        plots.append(yieldPlots.addYield(ElMuHasFakeableDileptonPreMllCutTightVeto,"ElMuHasFakeableDileptonPreMllCutTightVeto","OS fakeable leptons + $M_{ll}$ cut + Two tights veto (channel $e^{\pm}\mu^{\mp}$)"))
        # Trigger Plots #
        plots.extend(triggerPlots(sel         = ElElHasFakeableDileptonPreMllCutOutZTightVeto,
                                  triggerDict = self.triggersPerPrimaryDataset,
                                  suffix      = "ElElHasFakeableDileptonPreMllCutOutZTightVeto",
                                  channel     = "ElEl"))
        plots.extend(triggerPlots(sel         = MuMuHasFakeableDileptonPreMllCutOutZTightVeto,
                                  triggerDict = self.triggersPerPrimaryDataset,
                                  suffix      = "MuMuHasFakeableDileptonPreMllCutOutZTightVeto",
                                  channel     = "MuMu"))
        plots.extend(triggerPlots(sel         = ElMuHasFakeableDileptonPreMllCutTightVeto,
                                  triggerDict = self.triggersPerPrimaryDataset,
                                  suffix      = "ElMuHasFakeableDileptonPreMllCutTightVeto",
                                  channel     = "ElMu"))

        # Dilepton channel plots #
        plots.extend(channelPlot(sel        = ElElHasFakeableDileptonPreMllCutOutZTightVeto,
                                 DilepElEl  = OsElElFakeable,
                                 DilepMuMu  = OsMuMuFakeable,
                                 DilepElMu  = OsElMuFakeable,
                                 suffix     = "HasFakeableDileptonPreMllCutOutZTightVeto",
                                 channel    = "ElEl"))
        plots.extend(channelPlot(sel        = MuMuHasFakeableDileptonPreMllCutOutZTightVeto,
                                 DilepElEl  = OsElElFakeable,
                                 DilepMuMu  = OsMuMuFakeable,
                                 DilepElMu  = OsElMuFakeable,
                                 suffix     = "HasFakeableDileptonPreMllCutOutZTightVeto",
                                 channel    = "MuMu"))
        plots.extend(channelPlot(sel        = ElMuHasFakeableDileptonPreMllCutTightVeto,
                                 DilepElEl  = OsElElFakeable,
                                 DilepMuMu  = OsMuMuFakeable,
                                 DilepElMu  = OsElMuFakeable,
                                 suffix     = "HasFakeableDileptonPreMllCutTightVeto",
                                 channel    = "ElMu"))


#    If both fakeable leptons are electrons, the event needs to pass either the single electron or the double electron trigger; if both fakeable leptons are muons, the event needs to pass either the single muon or the double muon trigger; if one fakeable lepton is an electron and the other fakeable lepton is a muon, the event needs to pass either the single electron or the single muon or the muon+electron trigger
        #############################################################################
        #                              AK4 Jets                                     #
        #############################################################################
        # Jet Selection #
        ElElHasFakeableDileptonPreMllCutOutZTightVetoTwoAk4Jets = ElElHasFakeableDileptonPreMllCutOutZTightVeto.refine("ElElHasFakeableDileptonPreMllCutOutZTightVetoTwoAk4Jets", 
                                cut=[op.rng_len(self.ak4Jets)>=2])
        MuMuHasFakeableDileptonPreMllCutOutZTightVetoTwoAk4Jets = MuMuHasFakeableDileptonPreMllCutOutZTightVeto.refine("MuMuHasFakeableDileptonPreMllCutOutZTightVetoTwoAk4Jets", 
                                cut=[op.rng_len(self.ak4Jets)>=2])
        ElMuHasFakeableDileptonPreMllCutTightVetoTwoAk4Jets     = ElMuHasFakeableDileptonPreMllCutTightVeto.refine("ElMuHasFakeableDileptonPreMllCutTightVetoTwoAk4Jets", 
                                cut=[op.rng_len(self.ak4Jets)>=2])

        # Yield 
        plots.append(yieldPlots.addYield(ElElHasFakeableDileptonPreMllCutOutZTightVetoTwoAk4Jets,
                                         "ElElHasFakeableDileptonPreMllCutOutZTightVetoTwoAk4Jets",
                                         "OS leptons + Ak4 Jets $\geq 2$ (channel $e^+e^-$)"))
        plots.append(yieldPlots.addYield(MuMuHasFakeableDileptonPreMllCutOutZTightVetoTwoAk4Jets,
                                         "MuMuHasFakeableDileptonPreMllCutOutZTightVetoTwoAk4Jets",
                                         "OS leptons + Ak4 Jets $\geq 2$ (channel $\mu^+\mu^-$)"))
        plots.append(yieldPlots.addYield(ElMuHasFakeableDileptonPreMllCutTightVetoTwoAk4Jets,
                                         "ElMuHasFakeableDileptonPreMllCutTightVetoTwoAk4Jets",
                                         "OS leptons + Ak4 Jets $\geq 2$ (channel $e^{\pm}\mu^{\mp}$)"))

        # Plots #
        HasFakeableDileptonTwoAk4JetsChannelList = [
                             {'channel':'ElEl','sel':ElElHasFakeableDileptonPreMllCutOutZTightVetoTwoAk4Jets,'dilepton':OsElElFakeable[0],'suffix':'HasFakeableDileptonPreMllCutOutZTightVetoTwoAk4Jets'},
                             {'channel':'MuMu','sel':MuMuHasFakeableDileptonPreMllCutOutZTightVetoTwoAk4Jets,'dilepton':OsMuMuFakeable[0],'suffix':'HasFakeableDileptonPreMllCutOutZTightVetoTwoAk4Jets'},
                             {'channel':'ElMu','sel':ElMuHasFakeableDileptonPreMllCutTightVetoTwoAk4Jets,    'dilepton':OsElMuFakeable[0],'suffix':'HasFakeableDileptonPreMllCutTightVetoTwoAk4Jets'},
                                                  ]
        commonItems = ['channel','sel','suffix']
        HasFakeableDileptonTwoAk4JetsBaseDictJetsPlots = {'leadjet':self.ak4Jets[0],'subleadjet':self.ak4Jets[1],'lead_is_b':False,'sublead_is_b':False}
        HasFakeableDileptonTwoAk4JetsBaseDictJetsN = {'objName':'Ak4Jets','objCont':self.ak4Jets,'Nmax':10,'xTitle':'N(Ak4 jets)'}
                                                    
        for channelDict in HasFakeableDileptonTwoAk4JetsChannelList:
            # Dilepton #
            plots.extend(makeDileptonPlots(**channelDict,isMC=isMC))
            # Ak4 Jets #
            plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**HasFakeableDileptonTwoAk4JetsBaseDictJetsN))
            plots.extend(makeAk4JetsPlots(**{k:channelDict[k] for k in commonItems},**HasFakeableDileptonTwoAk4JetsBaseDictJetsPlots))
            # MET #
            plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))


        #############################################################################
        #                              Ak4 Btagging                                 #
        #############################################################################
        # Resolved selection #
        # Inclusive resolved #
#        ElElHasFakeableDileptonPreMllCutOutZTightVetoInclusiveResolved = ElElHasFakeableDileptonPreMllCutOutZTightVetoTwoAk4Jets.refine(
#                                "ElElHasFakeableDileptonPreMllCutOutZTightVetoInclusiveResolved", 
#                                cut=[op.rng_len(self.ak4BJets) >= 1])
#        MuMuHasFakeableDileptonPreMllCutOutZTightVetoInclusiveResolved = MuMuHasFakeableDileptonPreMllCutOutZTightVetoTwoAk4Jets.refine(
#                                "MuMuHasFakeableDileptonPreMllCutOutZTightVetoInclusiveResolved", 
#                                cut=[op.rng_len(self.ak4BJets) >= 1])
#        ElMuHasFakeableDileptonPreMllCutTightVetoInclusiveResolved = ElMuHasFakeableDileptonPreMllCutTightVetoTwoAk4Jets.refine("
#                                ElMuHasFakeableDileptonPreMllCutTightVetoInclusiveResolved", 
#                                cut=[op.rng_len(self.ak4BJets) >= 1])
        # Applied SF #
        Ak4BtaggingOneBtag  = [DeepJetMediumSF(self.ak4BJets[0])] if isMC else None
        Ak4BtaggingTwoBtags = [DeepJetMediumSF(self.ak4BJets[0]),DeepJetMediumSF(self.ak4BJets[1])] if isMC else None

        # Exclusive Resolved No Btag #
        ElElHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedNoBtag = ElElHasFakeableDileptonPreMllCutOutZTightVetoTwoAk4Jets.refine(
                                "ElElHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedNoBtag", 
                                cut=[op.rng_len(self.ak4BJets)==0,op.rng_len(self.ak8BJets)==0])
        MuMuHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedNoBtag = MuMuHasFakeableDileptonPreMllCutOutZTightVetoTwoAk4Jets.refine(
                                "MuMuHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedNoBtag", 
                                cut=[op.rng_len(self.ak4BJets)==0,op.rng_len(self.ak8BJets)==0])
        ElMuHasFakeableDileptonPreMllCutTightVetoExclusiveResolvedNoBtag = ElMuHasFakeableDileptonPreMllCutTightVetoTwoAk4Jets.refine(
                                "ElMuHasFakeableDileptonPreMllCutTightVetoExclusiveResolvedNoBtag", 
                                cut=[op.rng_len(self.ak4BJets)==0,op.rng_len(self.ak8BJets)==0])
        # Exclusive resolved 1 Btag #
        ElElHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedOneBtag = ElElHasFakeableDileptonPreMllCutOutZTightVetoTwoAk4Jets.refine(
                                "ElElHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedOneBtag", 
                                cut=[op.rng_len(self.ak4BJets)==1,op.rng_len(self.ak8BJets)==0],
                                weight=Ak4BtaggingOneBtag)
        MuMuHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedOneBtag = MuMuHasFakeableDileptonPreMllCutOutZTightVetoTwoAk4Jets.refine(
                                "MuMuHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedOneBtag", 
                                cut=[op.rng_len(self.ak4BJets)==1,op.rng_len(self.ak8BJets)==0],
                                weight=Ak4BtaggingOneBtag)
        ElMuHasFakeableDileptonPreMllCutTightVetoExclusiveResolvedOneBtag = ElMuHasFakeableDileptonPreMllCutTightVetoTwoAk4Jets.refine(
                                "ElMuHasFakeableDileptonPreMllCutTightVetoExclusiveResolvedOneBtag", 
                                cut=[op.rng_len(self.ak4BJets)==1,op.rng_len(self.ak8BJets)==0],
                                weight=Ak4BtaggingOneBtag)
        # Exclusive resolved 2 Btags #
        ElElHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedTwoBtags = ElElHasFakeableDileptonPreMllCutOutZTightVetoTwoAk4Jets.refine(
                                "ElElHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedTwoBtags", 
                                cut=[op.rng_len(self.ak4BJets)==2,op.rng_len(self.ak8BJets)==0],
                                weight=Ak4BtaggingTwoBtags)
        MuMuHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedTwoBtags = MuMuHasFakeableDileptonPreMllCutOutZTightVetoTwoAk4Jets.refine(
                                "MuMuHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedTwoBtags", 
                                cut=[op.rng_len(self.ak4BJets)==2,op.rng_len(self.ak8BJets)==0],
                                weight=Ak4BtaggingTwoBtags)
        ElMuHasFakeableDileptonPreMllCutTightVetoExclusiveResolvedTwoBtags = ElMuHasFakeableDileptonPreMllCutTightVetoTwoAk4Jets.refine(
                                "ElMuHasFakeableDileptonPreMllCutTightVetoExclusiveResolvedTwoBtags", 
                                cut=[op.rng_len(self.ak4BJets)==2,op.rng_len(self.ak8BJets)==0],
                                weight=Ak4BtaggingTwoBtags)

                                
        # Yield 
#        plots.append(yieldPlots.addYield(ElElHasFakeableDileptonPreMllCutOutZTightVetoInclusiveResolved,
#                                         "ElElHasFakeableDileptonPreMllCutOutZTightVetoInclusiveResolved",
#                                         "OS leptons + Inclusive Resolved Jets (channel $e^+e^-$)"))
#        plots.append(yieldPlots.addYield(MuMuHasFakeableDileptonPreMllCutOutZTightVetoInclusiveResolved,
#                                         "MuMuHasFakeableDileptonPreMllCutOutZTightVetoInclusiveResolved",
#                                         "OS leptons + Inclusive Resolved Jets (channel $\mu^+\mu^-$)"))
#        plots.append(yieldPlots.addYield(ElMuHasFakeableDileptonPreMllCutTightVetoInclusiveResolved,
#                                         "ElMuHasFakeableDileptonPreMllCutTightVetoInclusiveResolved",
#                                         "OS leptons + Inclusive Resolved Jets (channel $e^{\pm}\mu^{\mp}$)"))

        plots.append(yieldPlots.addYield(ElElHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedNoBtag,
                                         "ElElHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedNoBtag",
                                         "OS leptons + Exclusive Resolved Jets (0 bjet) (channel $e^+e^-$)"))
        plots.append(yieldPlots.addYield(MuMuHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedNoBtag,
                                         "MuMuHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedNoBtag",
                                         "OS leptons + Exclusive Resolved Jets (0 bjet) (channel $\mu^+\mu^-$)"))
        plots.append(yieldPlots.addYield(ElMuHasFakeableDileptonPreMllCutTightVetoExclusiveResolvedNoBtag,
                                         "ElMuHasFakeableDileptonPreMllCutTightVetoExclusiveResolvedNoBtag",
                                         "OS leptons + Exclusive Resolved Jets (0 bjet) (channel $e^{\pm}\mu^{\mp}$)"))

        plots.append(yieldPlots.addYield(ElElHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedOneBtag,
                                         "ElElHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedOneBtag",
                                         "OS leptons + Exclusive Resolved Jets (1 bjet) (channel $e^+e^-$)"))
        plots.append(yieldPlots.addYield(MuMuHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedOneBtag,
                                         "MuMuHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedOneBtag",
                                         "OS leptons + Exclusive Resolved Jets (1 bjet) (channel $\mu^+\mu^-$)"))
        plots.append(yieldPlots.addYield(ElMuHasFakeableDileptonPreMllCutTightVetoExclusiveResolvedOneBtag,
                                         "ElMuHasFakeableDileptonPreMllCutTightVetoExclusiveResolvedOneBtag",
                                         "OS leptons + Exclusive Resolved Jets (1 bjet) (channel $e^{\pm}\mu^{\mp}$)"))

        plots.append(yieldPlots.addYield(ElElHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedTwoBtags,
                                         "ElElHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedTwoBtags",
                                         "OS leptons + Exclusive Resolved Jets (2 bjets) (channel $e^+e^-$)"))
        plots.append(yieldPlots.addYield(MuMuHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedTwoBtags,
                                         "MuMuHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedTwoBtags",
                                         "OS leptons + Exclusive Resolved Jets (2 bjets) (channel $\mu^+\mu^-$)"))
        plots.append(yieldPlots.addYield(ElMuHasFakeableDileptonPreMllCutTightVetoExclusiveResolvedTwoBtags,
                                         "ElMuHasFakeableDileptonPreMllCutTightVetoExclusiveResolvedTwoBtags",
                                         "OS leptons + Exclusive Resolved Jets (2 bjets) (channel $e^{\pm}\mu^{\mp}$)"))

        # Dilepton plots #
        HasFakeableDileptonResolvedJetsChannelListLeptons = [
#             {'channel':'ElEl','sel':ElElHasFakeableDileptonPreMllCutOutZTightVetoInclusiveResolved,'dilepton':OsElElFakeable[0],'suffix':'HasFakeableDileptonPreMllCutOutZTightVetoInclusiveResolved'},
#             {'channel':'MuMu','sel':MuMuHasFakeableDileptonPreMllCutOutZTightVetoInclusiveResolved,'dilepton':OsMuMuFakeable[0],'suffix':'HasFakeableDileptonPreMllCutOutZTightVetoInclusiveResolved'},
#             {'channel':'ElMu','sel':ElMuHasFakeableDileptonPreMllCutTightVetoInclusiveResolved,    'dilepton':OsElMuFakeable[0],'suffix':'HasFakeableDileptonPreMllCutTightVetoInclusiveResolved'},
             {'channel':'ElEl','sel':ElElHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedNoBtag,'dilepton':OsElElFakeable[0],'suffix':'HasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedNoBtag'},
             {'channel':'MuMu','sel':MuMuHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedNoBtag,'dilepton':OsMuMuFakeable[0],'suffix':'HasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedNoBtag'},
             {'channel':'ElMu','sel':ElMuHasFakeableDileptonPreMllCutTightVetoExclusiveResolvedNoBtag,    'dilepton':OsElMuFakeable[0],'suffix':'HasFakeableDileptonPreMllCutTightVetoExclusiveResolvedNoBtag'},
             {'channel':'ElEl','sel':ElElHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedOneBtag,'dilepton':OsElElFakeable[0],'suffix':'HasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedOneBtag'},
             {'channel':'MuMu','sel':MuMuHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedOneBtag,'dilepton':OsMuMuFakeable[0],'suffix':'HasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedOneBtag'},
             {'channel':'ElMu','sel':ElMuHasFakeableDileptonPreMllCutTightVetoExclusiveResolvedOneBtag,    'dilepton':OsElMuFakeable[0],'suffix':'HasFakeableDileptonPreMllCutTightVetoExclusiveResolvedOneBtag'},
             {'channel':'ElEl','sel':ElElHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedTwoBtags,'dilepton':OsElElFakeable[0],'suffix':'HasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedTwoBtags'},
             {'channel':'MuMu','sel':MuMuHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedTwoBtags,'dilepton':OsMuMuFakeable[0],'suffix':'HasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedTwoBtags'},
             {'channel':'ElMu','sel':ElMuHasFakeableDileptonPreMllCutTightVetoExclusiveResolvedTwoBtags,    'dilepton':OsElMuFakeable[0],'suffix':'HasFakeableDileptonPreMllCutTightVetoExclusiveResolvedTwoBtags'},
                                                            ]

        commonItems = ['channel','sel','suffix']
        HasFakeableDileptonResolvedJetsBaseDictJetsN = {'objName':'Ak4BJets','objCont':self.ak4BJets,'Nmax':10,'xTitle':'N(Ak4 Bjets)'}

        for channelDict in HasFakeableDileptonResolvedJetsChannelListLeptons:
            # Dilepton #
            plots.extend(makeDileptonPlots(**channelDict,isMC=isMC))
            # Number of jets #
            plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**HasFakeableDileptonResolvedJetsBaseDictJetsN))
            # MET #
            plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))
 

        # Dijet plots #
        HasFakeableDileptonResolvedJetsChannelListJets = [
            {'channel':'ElEl','sel':ElElHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedNoBtag,'suffix':'ElElHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedNoBtag','leadjet':self.ak4LightJets[0],'subleadjet':self.ak4LightJets[1],'lead_is_b':False,'sublead_is_b':False},
            {'channel':'MuMu','sel':MuMuHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedNoBtag,'suffix':'MuMuHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedNoBtag','leadjet':self.ak4LightJets[0],'subleadjet':self.ak4LightJets[1],'lead_is_b':False,'sublead_is_b':False},
            {'channel':'ElMu','sel':ElMuHasFakeableDileptonPreMllCutTightVetoExclusiveResolvedNoBtag,'suffix':'ElMuHasFakeableDileptonPreMllCutTightVetoExclusiveResolvedNoBtag','leadjet':self.ak4LightJets[0],'subleadjet':self.ak4LightJets[1],'lead_is_b':False,'sublead_is_b':False},

            {'channel':'ElEl','sel':ElElHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedOneBtag,'suffix':'ElElHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedOneBtag','leadjet':self.ak4BJets[0],'subleadjet':self.ak4LightJets[0],'lead_is_b':True,'sublead_is_b':False},
            {'channel':'MuMu','sel':MuMuHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedOneBtag,'suffix':'MuMuHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedOneBtag','leadjet':self.ak4BJets[0],'subleadjet':self.ak4LightJets[0],'lead_is_b':True,'sublead_is_b':False},
            {'channel':'ElMu','sel':ElMuHasFakeableDileptonPreMllCutTightVetoExclusiveResolvedOneBtag,'suffix':'ElMuHasFakeableDileptonPreMllCutTightVetoExclusiveResolvedOneBtag','leadjet':self.ak4BJets[0],'subleadjet':self.ak4LightJets[0],'lead_is_b':True,'sublead_is_b':False},

            {'channel':'ElEl','sel':ElElHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedTwoBtags,'suffix':'ElElHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedTwoBtags','leadjet':self.ak4BJets[0],'subleadjet':self.ak4BJets[1],'lead_is_b':True,'sublead_is_b':True},
            {'channel':'MuMu','sel':MuMuHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedTwoBtags,'suffix':'MuMuHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedTwoBtags','leadjet':self.ak4BJets[0],'subleadjet':self.ak4BJets[1],'lead_is_b':True,'sublead_is_b':True},
            {'channel':'ElMu','sel':ElMuHasFakeableDileptonPreMllCutTightVetoExclusiveResolvedTwoBtags,'suffix':'ElMuHasFakeableDileptonPreMllCutTightVetoExclusiveResolvedTwoBtags','leadjet':self.ak4BJets[0],'subleadjet':self.ak4BJets[1],'lead_is_b':True,'sublead_is_b':True},
                                                          ]
        for channelDict in HasFakeableDileptonResolvedJetsChannelListJets:
            plots.extend(makeAk4JetsPlots(**channelDict))

        #############################################################################
        #                              AK8 Jets                                     #
        #############################################################################
        # Jet Selection #
        ElElHasFakeableDileptonPreMllCutOutZTightVetoOneAk8Jet = ElElHasFakeableDileptonPreMllCutOutZTightVeto.refine("ElElHasFakeableDileptonPreMllCutOutZTightVetoOneAk8Jet", 
                                cut=[op.rng_len(self.ak8Jets)>=1])
        MuMuHasFakeableDileptonPreMllCutOutZTightVetoOneAk8Jet = MuMuHasFakeableDileptonPreMllCutOutZTightVeto.refine("MuMuHasFakeableDileptonPreMllCutOutZTightVetoOneAk8Jet", 
                                cut=[op.rng_len(self.ak8Jets)>=1])
        ElMuHasFakeableDileptonPreMllCutTightVetoOneAk8Jet     = ElMuHasFakeableDileptonPreMllCutTightVeto.refine("ElMuHasFakeableDileptonPreMllCutTightVetoOneAk8Jet", 
                                cut=[op.rng_len(self.ak8Jets)>=1])

        # Yield 
        plots.append(yieldPlots.addYield(ElElHasFakeableDileptonPreMllCutOutZTightVetoOneAk8Jet,
                                         "ElElHasFakeableDileptonPreMllCutOutZTightVetoOneAk8Jet",
                                         "OS leptons + Ak8 Jets $\geq 1$ (channel $e^+e^-$)"))
        plots.append(yieldPlots.addYield(MuMuHasFakeableDileptonPreMllCutOutZTightVetoOneAk8Jet,
                                         "MuMuHasFakeableDileptonPreMllCutOutZTightVetoOneAk8Jet",
                                         "OS leptons + Ak8 Jets $\geq 1$ (channel $\mu^+\mu^-$)"))
        plots.append(yieldPlots.addYield(ElMuHasFakeableDileptonPreMllCutTightVetoOneAk8Jet,
                                         "ElMuHasFakeableDileptonPreMllCutTightVetoOneAk8Jet",
                                         "OS leptons + Ak8 Jets $\geq 1$ (channel $e^{\pm}\mu^{\mp}$)"))

        # Plots #
        HasFakeableDileptonOneAk8JetJetsChannelList = [
                             {'channel':'ElEl','sel':ElElHasFakeableDileptonPreMllCutOutZTightVetoOneAk8Jet,'dilepton':OsElElFakeable[0],'suffix':'HasFakeableDileptonPreMllCutOutZTightVetoOneAk8Jet'},
                             {'channel':'MuMu','sel':MuMuHasFakeableDileptonPreMllCutOutZTightVetoOneAk8Jet,'dilepton':OsMuMuFakeable[0],'suffix':'HasFakeableDileptonPreMllCutOutZTightVetoOneAk8Jet'},
                             {'channel':'ElMu','sel':ElMuHasFakeableDileptonPreMllCutTightVetoOneAk8Jet,    'dilepton':OsElMuFakeable[0],'suffix':'HasFakeableDileptonPreMllCutTightVetoOneAk8Jet'},
                                                  ]
        commonItems = ['channel','sel','suffix']
        HasFakeableDileptonOneAk8JetsBaseDictFatjetsN = {'objName':'Ak8Jets','objCont':self.ak8Jets,'Nmax':10,'xTitle':'N(Ak8 jets)'}

                                                    
        for channelDict in HasFakeableDileptonOneAk8JetJetsChannelList:
            # Dilepton #
            plots.extend(makeDileptonPlots(**channelDict,isMC=isMC))
            # Ak8 Jets #
            plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**HasFakeableDileptonOneAk8JetsBaseDictFatjetsN))
            plots.extend(makeAk8Plots(fatjet=self.ak8Jets[0],**{k:channelDict[k] for k in commonItems}))
            # MET #
            plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))


        #############################################################################
        #                             Boosted jets                                  #
        #############################################################################
        # Boosted selection #
        ElElHasFakeableDileptonPreMllCutOutZTightVetoInclusiveBoosted = ElElHasFakeableDileptonPreMllCutOutZTightVetoOneAk8Jet.refine("ElElHasFakeableDileptonPreMllCutOutZTightVetoInclusiveBoosted", 
                                cut=[op.rng_len(self.ak8BJets)>=1])
        MuMuHasFakeableDileptonPreMllCutOutZTightVetoInclusiveBoosted = MuMuHasFakeableDileptonPreMllCutOutZTightVetoOneAk8Jet.refine("MuMuHasFakeableDileptonPreMllCutOutZTightVetoInclusiveBoosted", 
                                cut=[op.rng_len(self.ak8BJets)>=1])
        ElMuHasFakeableDileptonPreMllCutTightVetoInclusiveBoosted     = ElMuHasFakeableDileptonPreMllCutTightVetoOneAk8Jet.refine("ElMuHasFakeableDileptonPreMllCutTightVetoInclusiveBoosted", 
                                cut=[op.rng_len(self.ak8BJets)>=1])

#        ElElHasFakeableDileptonPreMllCutOutZTightVetoExclusiveBoosted = ElElHasFakeableDileptonPreMllCutOutZTightVetoOneAk8Jet.refine("ElElHasFakeableDileptonPreMllCutOutZTightVetoExclusiveBoosted", 
#                                cut=[op.rng_len(self.ak8BJets)>=1, op.OR(op.rng_len(self.ak4Jets)<=1 , op.rng_len(self.ak4BJets)==0)])
#        MuMuHasFakeableDileptonPreMllCutOutZTightVetoExclusiveBoosted = MuMuHasFakeableDileptonPreMllCutOutZTightVetoOneAk8Jet.refine("MuMuHasFakeableDileptonPreMllCutOutZTightVetoExclusiveBoosted", 
#                                cut=[op.rng_len(self.ak8BJets)>=1, op.OR(op.rng_len(self.ak4Jets)<=1 , op.rng_len(self.ak4BJets)==0)])
#        ElMuHasFakeableDileptonPreMllCutTightVetoExclusiveBoosted     = ElMuHasFakeableDileptonPreMllCutTightVetoOneAk8Jet.refine("ElMuHasFakeableDileptonPreMllCutTightVetoExclusiveBoosted", 
#                                cut=[op.rng_len(self.ak8BJets)>=1, op.OR(op.rng_len(self.ak4Jets)<=1 , op.rng_len(self.ak4BJets)==0)])

        # Yield 
        plots.append(yieldPlots.addYield(ElElHasFakeableDileptonPreMllCutOutZTightVetoInclusiveBoosted,
                                         "ElElHasFakeableDileptonPreMllCutOutZTightVetoInclusiveBoosted",
                                         "OS leptons + Inclusive Boosted Jets (channel $e^+e^-$)"))
        plots.append(yieldPlots.addYield(MuMuHasFakeableDileptonPreMllCutOutZTightVetoInclusiveBoosted,
                                         "MuMuHasFakeableDileptonPreMllCutOutZTightVetoInclusiveBoosted",
                                         "OS leptons + Inclusive Boosted Jets (channel $\mu^+\mu^-$)"))
        plots.append(yieldPlots.addYield(ElMuHasFakeableDileptonPreMllCutTightVetoInclusiveBoosted,
                                         "ElMuHasFakeableDileptonPreMllCutTightVetoInclusiveBoosted",
                                         "OS leptons + Inclusive Boosted Jets (channel $e^{\pm}\mu^{\mp}$)"))
#        plots.append(yieldPlots.addYield(ElElHasFakeableDileptonPreMllCutOutZTightVetoExclusiveBoosted,
#                                         "ElElHasFakeableDileptonPreMllCutOutZTightVetoExclusiveBoosted",
#                                         "OS leptons + Exclusive Boosted Jets (channel $e^+e^-$)"))
#        plots.append(yieldPlots.addYield(MuMuHasFakeableDileptonPreMllCutOutZTightVetoExclusiveBoosted,
#                                         "MuMuHasFakeableDileptonPreMllCutOutZTightVetoExclusiveBoosted",
#                                         "OS leptons + Exclusive Boosted Jets (channel $\mu^+\mu^-$)"))
#        plots.append(yieldPlots.addYield(ElMuHasFakeableDileptonPreMllCutTightVetoExclusiveBoosted,
#                                         "ElMuHasFakeableDileptonPreMllCutTightVetoExclusiveBoosted",
#                                         "OS leptons + Exclusive Boosted Jets (channel $e^{\pm}\mu^{\mp}$)"))

        # Plots #
        HasFakeableDileptonBoostedJetsChannelList = [
             {'channel':'ElEl','sel':ElElHasFakeableDileptonPreMllCutOutZTightVetoInclusiveBoosted,'dilepton':OsElElFakeable[0],'suffix':'HasFakeableDileptonPreMllCutOutZTightVetoInclusiveBoosted'},
             {'channel':'MuMu','sel':MuMuHasFakeableDileptonPreMllCutOutZTightVetoInclusiveBoosted,'dilepton':OsMuMuFakeable[0],'suffix':'HasFakeableDileptonPreMllCutOutZTightVetoInclusiveBoosted'},
             {'channel':'ElMu','sel':ElMuHasFakeableDileptonPreMllCutTightVetoInclusiveBoosted,    'dilepton':OsElMuFakeable[0],'suffix':'HasFakeableDileptonPreMllCutTightVetoInclusiveBoosted'},
#             {'channel':'ElEl','sel':ElElHasFakeableDileptonPreMllCutOutZTightVetoExclusiveBoosted,'dilepton':OsElElFakeable[0],'suffix':'HasFakeableDileptonPreMllCutOutZTightVetoExclusiveBoosted'},
#             {'channel':'MuMu','sel':MuMuHasFakeableDileptonPreMllCutOutZTightVetoExclusiveBoosted,'dilepton':OsMuMuFakeable[0],'suffix':'HasFakeableDileptonPreMllCutOutZTightVetoExclusiveBoosted'},
#             {'channel':'ElMu','sel':ElMuHasFakeableDileptonPreMllCutTightVetoExclusiveBoosted,    'dilepton':OsElMuFakeable[0],'suffix':'HasFakeableDileptonPreMllCutTightVetoExclusiveBoosted'},
                                                  ]
        commonItems = ['channel','sel','suffix']
        HasFakeableDileptonBoostedJetsBaseDictJetsN = {'objName':'Ak8BJets','objCont':self.ak8BJets,'Nmax':10,'xTitle':'N(Ak8 Bjets)'}

        for channelDict in HasFakeableDileptonBoostedJetsChannelList:
            # Dilepton #
            plots.extend(makeDileptonPlots(**channelDict,isMC=isMC))
            # Ak8 Jets #
            plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**HasFakeableDileptonBoostedJetsBaseDictJetsN))
            plots.extend(makeAk8Plots(fatjet=self.ak8BJets[0],**{k:channelDict[k] for k in commonItems}))
            # MET #
            plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))



        #############################################################################
        #                     High-level combinations                               #
        #############################################################################
        HasOsPreMllCutHighLevelVariablesChannelList = [
            # Resolved No Btag #
            {'channel'      : 'ElEl',
             'met'          : self.corrMET,
             'l1'           : OsElElFakeable[0][0],
             'l2'           : OsElElFakeable[0][1],
             'j1'           : self.ak4Jets[0],
             'j2'           : self.ak4Jets[1],
             'sel'          : ElElHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedNoBtag ,
             'suffix'       : 'ElElHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedNoBtag'},
            {'channel'      : 'MuMu',
             'met'          : self.corrMET,
             'l1'           : OsMuMuFakeable[0][0],
             'l2'           : OsMuMuFakeable[0][1],
             'j1'           : self.ak4Jets[0],
             'j2'           : self.ak4Jets[1],
             'sel'          : MuMuHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedNoBtag,
             'suffix'       : 'MuMuHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedNoBtag'},
            {'channel'      : 'ElMu',
             'met'          : self.corrMET,
             'l1'           : OsElMuFakeable[0][0],
             'l2'           : OsElMuFakeable[0][1],
             'j1'           : self.ak4Jets[0],
             'j2'           : self.ak4Jets[1],
             'sel'          : ElMuHasFakeableDileptonPreMllCutTightVetoExclusiveResolvedNoBtag,
             'suffix'       : 'ElMuHasFakeableDileptonPreMllCutTightVetoExclusiveResolvedNoBtag'},
            # Resolved One Btag #
            {'channel'      : 'ElEl',
             'met'          : self.corrMET,
             'l1'           : OsElElFakeable[0][0],
             'l2'           : OsElElFakeable[0][1],
             'j1'           : self.ak4BJets[0],
             'j2'           : self.ak4Jets[0],
             'sel'          : ElElHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedOneBtag,
             'suffix'       : 'ElElHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedOneBtag'},
            {'channel'      : 'MuMu',
             'met'          : self.corrMET,
             'l1'           : OsMuMuFakeable[0][0],
             'l2'           : OsMuMuFakeable[0][1],
             'j1'           : self.ak4BJets[0],
             'j2'           : self.ak4Jets[0],
             'sel'          : MuMuHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedOneBtag, 
             'suffix'       : 'MuMuHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedOneBtag'},
            {'channel'      : 'ElMu',
             'met'          : self.corrMET,
             'l1'           : OsElMuFakeable[0][0],
             'l2'           : OsElMuFakeable[0][1],
             'j1'           : self.ak4BJets[0],
             'j2'           : self.ak4Jets[0],
             'sel'          : ElMuHasFakeableDileptonPreMllCutTightVetoExclusiveResolvedOneBtag, 
             'suffix'       : 'ElMuHasFakeableDileptonPreMllCutTightVetoExclusiveResolvedOneBtag'},
            # Resolved Two Btags  #
            {'channel'      : 'ElEl',
             'met'          : self.corrMET,
             'l1'           : OsElElFakeable[0][0],
             'l2'           : OsElElFakeable[0][1],
             'j1'           : self.ak4BJets[0],
             'j2'           : self.ak4BJets[1],
             'sel'          : ElElHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedTwoBtags, 
             'suffix'       : 'ElElHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedTwoBtags'},
            {'channel'      : 'MuMu',
             'met'          : self.corrMET,
             'l1'           : OsMuMuFakeable[0][0],
             'l2'           : OsMuMuFakeable[0][1],
             'j1'           : self.ak4BJets[0],
             'j2'           : self.ak4BJets[1],
             'sel'          : MuMuHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedTwoBtags, 
             'suffix'       : 'MuMuHasFakeableDileptonPreMllCutOutZTightVetoExclusiveResolvedTwoBtags'},
            {'channel'      : 'ElMu',
             'met'          : self.corrMET,
             'l1'           : OsElMuFakeable[0][0],
             'l2'           : OsElMuFakeable[0][1],
             'j1'           : self.ak4BJets[0],
             'j2'           : self.ak4BJets[1],
             'sel'          : ElMuHasFakeableDileptonPreMllCutTightVetoExclusiveResolvedTwoBtags,
             'suffix'       : 'ElMuHasFakeableDileptonPreMllCutTightVetoExclusiveResolvedTwoBtags'},
            # Boosted #
            {'channel'      : 'ElEl',
             'met'          : self.corrMET,
             'l1'           : OsElElFakeable[0][0],
             'l2'           : OsElElFakeable[0][1],
             'j1'           : self.ak8BJets[0].subJet1,
             'j2'           : self.ak8BJets[0].subJet2,
             'sel'          : ElElHasFakeableDileptonPreMllCutOutZTightVetoInclusiveBoosted, 
             'suffix'       : 'ElElHasFakeableDileptonPreMllCutOutZTightVetoInclusiveBoosted'},
            {'channel'      : 'MuMu',
             'met'          : self.corrMET,
             'l1'           : OsMuMuFakeable[0][0],
             'l2'           : OsMuMuFakeable[0][1],
             'j1'           : self.ak8BJets[0].subJet1,
             'j2'           : self.ak8BJets[0].subJet2,
             'sel'          : MuMuHasFakeableDileptonPreMllCutOutZTightVetoInclusiveBoosted, 
             'suffix'       : 'MuMuHasFakeableDileptonPreMllCutOutZTightVetoInclusiveBoosted'},
            {'channel'      : 'ElMu',
             'met'          : self.corrMET,
             'l1'           : OsElMuFakeable[0][0],
             'l2'           : OsElMuFakeable[0][1],
             'j1'           : self.ak8BJets[0].subJet1,
             'j2'           : self.ak8BJets[0].subJet2,
             'sel'          : ElMuHasFakeableDileptonPreMllCutTightVetoInclusiveBoosted, 
             'suffix'       : 'ElMuHasFakeableDileptonPreMllCutTightVetoInclusiveBoosted'},
                                      ]

        for channelDict in HasOsPreMllCutHighLevelVariablesChannelList:
            plots.extend(makeHighLevelQuantities(**channelDict))

        return plots

