import os
import sys

from bamboo.analysismodules import HistogramsModule
from bamboo import treefunctions as op
from bamboo.plots import Plot, EquidistantBinning, SummedPlot

sys.path.append('/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/') # Add scripts in this directory -- TODO : make cleaner
from BaseHHtobbWW import BaseNanoHHtobbWW
from plotDef import *

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
        #                               Dilepton                                    #
        #############################################################################
        # Make dilepton fakeable pairs wrt to channel #
        lambdaOS  = lambda l1,l2 : l1.charge != l2.charge           # Opposite sign charges
        lambdaPTCut = lambda l1,l2: op.OR( l1.pt > 25 , l2.pt > 25) # PT cut 
        lambdaFakeableDilepton = lambda l1, l2 : op.AND( lambdaOS(l1,l2) , lambdaPTCut(l1,l2) )
        OsElElFakeable = op.combine(self.electronsFakeSel, N=2, pred=lambdaFakeableDilepton)
        OsMuMuFakeable = op.combine(self.muonsFakeSel, N=2, pred=lambdaFakeableDilepton)
        OsElMuFakeable = op.combine((self.electronsFakeSel, self.muonsFakeSel), pred=lambdaFakeableDilepton)

        # At least 2 fakeables leptons #
        ElElHasFakeableDilepton = noSel.refine("ElElHasFakeableDilepton",
                                        cut = [op.rng_len(OsElElFakeable)>0])
        MuMuHasFakeableDilepton = noSel.refine("MuMuHasFakeableDilepton",
                                        cut = [op.rng_len(OsMuMuFakeable)>0])
        ElMuHasFakeableDilepton = noSel.refine("ElMuHasFakeableDilepton",
                                        cut = [op.rng_len(OsElMuFakeable)>0])

#        plots.append(Plot.make1D("test1",
#                                 op.rng_len(self.electronsFakeSel),
#                                 ElElHasFakeableDilepton,
#                                 EquidistantBinning(1, 0., 1.),
#                                 xTitle = 'test'))
#        plots.append(Plot.make1D("test2",
#                                 self.muonsFakeSel[0].pt,
#                                 noSel.refine('test',cut=[op.rng_len(self.muonsFakeSel)>0]),
#                                 EquidistantBinning(1, 0., 1.),
#                                 xTitle = 'test'))
#        plots.append(Plot.make1D("test3",
#                                 op.rng_len(OsMuMuFakeable),
#                                 noSel,
#                                 EquidistantBinning(1, 0., 1.),
#                                 xTitle = 'test'))
        # Yield #
        plots.append(yieldPlots.addYield(ElElHasFakeableDilepton,"ElElHasFakeableDilepton","OS fakeable leptons (channel $e^+e^-$)"))
        plots.append(yieldPlots.addYield(MuMuHasFakeableDilepton,"MuMuHasFakeableDilepton","OS fakeable leptons (channel $\mu^+\mu^-$)"))
        plots.append(yieldPlots.addYield(ElMuHasFakeableDilepton,"ElMuHasFakeableDilepton","OS fakeable leptons (channel $e^{\pm}\mu^{\mp}$)"))
    
        # Dilepton Z peak exclusion (charge already done in previous selection) # 
        lambda_lowMllCut    = lambda dilep: op.invariant_mass(dilep[0].p4, dilep[1].p4)>12.
        lambda_outZ         = lambda dilep: op.NOT(op.in_range(80.,op.invariant_mass(dilep[0].p4, dilep[1].p4),100.))

        ElElHasFakeableDileptonMllCutOutZ = ElElHasFakeableDilepton.refine("ElElHasFakeableDileptonMllCutOutZ", 
                                cut=[lambda_lowMllCut(OsElElFakeable[0]),lambda_outZ(OsElElFakeable[0])])
        MuMuHasFakeableDileptonMllCutOutZ = MuMuHasFakeableDilepton.refine("MuMuHasFakeableDileptonMllCutOutZ", 
                                cut=[lambda_lowMllCut(OsMuMuFakeable[0]),lambda_outZ(OsMuMuFakeable[0])])
        ElMuHasFakeableDileptonMllCut     = ElMuHasFakeableDilepton.refine("ElMuHasFakeableDileptonMllCutOutZ", 
                                cut=[lambda_lowMllCut(OsElMuFakeable[0])])

        # Yield #
        plots.append(yieldPlots.addYield(ElElHasFakeableDileptonMllCutOutZ,"ElElHasFakeableDileptonMllCutOutZ","OS fakeable leptons + M_{ll} cut (channel $e^+e^-$)"))
        plots.append(yieldPlots.addYield(MuMuHasFakeableDileptonMllCutOutZ,"MuMuHasFakeableDileptonMllCutOutZ","OS fakeable leptons + M_{ll} cut (channel $\mu^+\mu^-$)"))
        plots.append(yieldPlots.addYield(ElMuHasFakeableDileptonMllCut,"ElMuHasFakeableDileptonMllCut","OS fakeable leptons + M_{ll} cut (channel $e^{\pm}\mu^{\mp}$)"))

        # Tight cut veto : no more than two tight leptons in the event #
        ElElHasFakeableDileptonMllCutOutZTightVeto = ElElHasFakeableDileptonMllCutOutZ.refine("ElElHasFakeableDileptonMllCutOutZTightVeto", 
                                cut=[op.rng_len(self.muonsTightSel)+op.rng_len(self.electronsTightSel) <= 2])
        MuMuHasFakeableDileptonMllCutOutZTightVeto = MuMuHasFakeableDileptonMllCutOutZ.refine("MuMuHasFakeableDileptonMllCutOutZTightVeto", 
                                cut=[op.rng_len(self.muonsTightSel)+op.rng_len(self.electronsTightSel) <= 2])
        ElMuHasFakeableDileptonMllCutTightVeto     = ElMuHasFakeableDileptonMllCut.refine("ElMuHasFakeableDileptonMllCutTightVeto", 
                                cut=[op.rng_len(self.muonsTightSel)+op.rng_len(self.electronsTightSel) <= 2])

        # Yield #
        plots.append(yieldPlots.addYield(ElElHasFakeableDileptonMllCutOutZTightVeto,"ElElHasFakeableDileptonMllCutOutZTightVeto","OS fakeable leptons + M_{ll} cut + Two tights veto (channel $e^+e^-$)"))
        plots.append(yieldPlots.addYield(MuMuHasFakeableDileptonMllCutOutZTightVeto,"MuMuHasFakeableDileptonMllCutOutZTightVeto","OS fakeable leptons + M_{ll} cut + Two tights veto (channel $\mu^+\mu^-$)"))
        plots.append(yieldPlots.addYield(ElMuHasFakeableDileptonMllCutTightVeto,"ElMuHasFakeableDileptonMllCutTightVeto","OS fakeable leptons + M_{ll} cut + Two tights veto (channel $e^{\pm}\mu^{\mp}$)"))


        # Dilepton channel plots #
        plots.extend(channelPlot(sel        = ElElHasFakeableDileptonMllCutOutZTightVeto,
                                 DilepElEl  = OsElElFakeable,
                                 DilepMuMu  = OsMuMuFakeable,
                                 DilepElMu  = OsElMuFakeable,
                                 suffix     = "HasFakeableDileptonMllCutOutZTightVeto",
                                 channel    = "ElEl"))
        plots.extend(channelPlot(sel        = MuMuHasFakeableDileptonMllCutOutZTightVeto,
                                 DilepElEl  = OsElElFakeable,
                                 DilepMuMu  = OsMuMuFakeable,
                                 DilepElMu  = OsElMuFakeable,
                                 suffix     = "HasFakeableDileptonMllCutOutZTightVeto",
                                 channel    = "MuMu"))
        plots.extend(channelPlot(sel        = ElMuHasFakeableDileptonMllCutTightVeto,
                                 DilepElEl  = OsElElFakeable,
                                 DilepMuMu  = OsMuMuFakeable,
                                 DilepElMu  = OsElMuFakeable,
                                 suffix     = "HasFakeableDileptonMllCutTightVeto",
                                 channel    = "ElMu"))

#    If both fakeable leptons are electrons, the event needs to pass either the single electron or the double electron trigger; if both fakeable leptons are muons, the event needs to pass either the single muon or the double muon trigger; if one fakeable lepton is an electron and the other fakeable lepton is a muon, the event needs to pass either the single electron or the single muon or the muon+electron trigger
#
#        # Scalefactors #
#        if self.isMC(sample):
#            doubleEleTrigSF = SF.get_scalefactor("dilepton", ("doubleEleLeg_HHMoriond17_2016"), systName="eleltrig")     
#            doubleMuTrigSF  = SF.get_scalefactor("dilepton", ("doubleMuLeg_HHMoriond17_2016"), systName="mumutrig")    
#            elemuTrigSF     = SF.get_scalefactor("dilepton", ("elemuLeg_HHMoriond17_2016"), systName="elmutrig")
#            mueleTrigSF     = SF.get_scalefactor("dilepton", ("mueleLeg_HHMoriond17_2016"), systName="mueltrig")
#
#            # From https://gitlab.cern.ch/ttH_leptons/doc/blob/master/Legacy/data_to_mc_corrections.md#trigger-efficiency-scale-factors
#        ttH_doubleMuon_trigSF = op.systematic(op.c_float(1.010), name="ttH_doubleMuon_trigSF", up=op.c_float(1.020), down=op.c_float(1.000))
#        ttH_doubleElectron_trigSF = op.systematic(op.c_float(1.020), name="ttH_doubleElectron_trigSF", up=op.c_float(1.040), down=op.c_float(1.000))
#        ttH_electronMuon_trigSF = op.systematic(op.c_float(1.020), name="ttH_electronMuon_trigSF", up=op.c_float(1.030), down=op.c_float(1.010))
#
#        llSF =  {
#            "ElEl" : (lambda ll : [ elTightIDSF(ll[0][0]),                                                                     # First lepton SF
#                                    elTightIDSF(ll[0][1]),                                                                     # Second lepton SF
#sFakeableDileptonMllCutOutZTightVeto                                    #doubleEleTrigSF(ll[0]),                                                                  # Dilepton SF
#                                    ttH_doubleElectron_trigSF,
#                                  ]),
#            "MuMu" : (lambda ll : [ muTightIDSF(ll[0][0]), muTightISOSF(ll[0][0]),                                             # First lepton SF
#                                    muTightIDSF(ll[0][1]), muTightISOSF(ll[0][1]),                                             # Second lepton SF
#                                    #doubleMuTrigSF(ll[0]),                                                                   # Dilepton SF
#                                    ttH_doubleMuon_trigSF,
#                                  ]),
#            "ElMu" : (lambda ll : [ elTightIDSF(ll[0][0]),                                                                     # First lepton SF
#                                    muTightIDSF(ll[0][1]), muTightISOSF(ll[0][1]),                                             # Second lepton SF
#                                    #elemuTrigSF(ll[0]),                                                                      # Dilepton SF
#                                    ttH_electronMuon_trigSF,
#                                  ]),
#                # ll is a proxy list of dileptons 
#                # ll[0] is the first dilepton 
#                # ll[0][0] is the first lepton and ll[0][1] the second in the dilepton
#                }
#
#
#        self.llSFApplied = {
#            "ElEl": llSF["ElEl"](self.OsElEl) if isMC else None,
#            "MuMu": llSF["MuMu"](self.OsMuMu) if isMC else None,
#            "ElMu": llSF["ElMu"](self.OsElMu) if isMC else None,
#                      }

        #############################################################################
        #                              AK4 Jets                                     #
        #############################################################################
        # Jet Selection #
        ElElHasFakeableDileptonMllCutOutZTightVetoTwoAk4Jets = ElElHasFakeableDileptonMllCutOutZTightVeto.refine("ElElHasFakeableDileptonMllCutOutZTightVetoTwoAk4Jets", 
                                cut=[op.rng_len(self.ak4Jets)>=2])
        MuMuHasFakeableDileptonMllCutOutZTightVetoTwoAk4Jets = MuMuHasFakeableDileptonMllCutOutZTightVeto.refine("MuMuHasFakeableDileptonMllCutOutZTightVetoTwoAk4Jets", 
                                cut=[op.rng_len(self.ak4Jets)>=2])
        ElMuHasFakeableDileptonMllCutTightVetoTwoAk4Jets     = ElMuHasFakeableDileptonMllCutTightVeto.refine("ElMuHasFakeableDileptonMllCutTightVetoTwoAk4Jets", 
                                cut=[op.rng_len(self.ak4Jets)>=2])

        # Yield 
        plots.append(yieldPlots.addYield(ElElHasFakeableDileptonMllCutOutZTightVetoTwoAk4Jets,
                                         "ElElHasFakeableDileptonMllCutOutZTightVetoTwoAk4Jets",
                                         "OS leptons + Ak4 Jets (channel $e^+e^-$)"))
        plots.append(yieldPlots.addYield(MuMuHasFakeableDileptonMllCutOutZTightVetoTwoAk4Jets,
                                         "MuMuHasFakeableDileptonMllCutOutZTightVetoTwoAk4Jets",
                                         "OS leptons + Ak4 Jets (channel $\mu^+\mu^-$)"))
        plots.append(yieldPlots.addYield(ElMuHasFakeableDileptonMllCutTightVetoTwoAk4Jets,
                                         "ElMuHasFakeableDileptonMllCutTightVetoTwoAk4Jets",
                                         "OS leptons + Ak4 Jets (channel $e^{\pm}\mu^{\mp}$)"))

        # Plots #
        HasFakeableDileptonTwoAk4JetsChannelList = [
                             {'channel':'ElEl','sel':ElElHasFakeableDileptonMllCutOutZTightVetoTwoAk4Jets,'dilepton':OsElElFakeable[0],'suffix':'HasFakeableDileptonMllCutOutZTightVetoTwoAk4Jets'},
                             {'channel':'MuMu','sel':MuMuHasFakeableDileptonMllCutOutZTightVetoTwoAk4Jets,'dilepton':OsMuMuFakeable[0],'suffix':'HasFakeableDileptonMllCutOutZTightVetoTwoAk4Jets'},
                             {'channel':'ElMu','sel':ElMuHasFakeableDileptonMllCutTightVetoTwoAk4Jets,    'dilepton':OsElMuFakeable[0],'suffix':'HasFakeableDileptonMllCutTightVetoTwoAk4Jets'},
                                                  ]
        commonItems = ['channel','sel','suffix']
        HasFakeableDileptonTwoAk4JetsBaseDictJetsPlots = {'leadjet':self.ak4Jets[0],'subleadjet':self.ak4Jets[1],'plot_type':'NoBtag'}
        HasFakeableDileptonTwoAk4JetsBaseDictJetsN = {'objName':'Ak4Jets','objCont':self.ak4Jets,'Nmax':5,'xTitle':'N(Ak4 jets)'}
                                                    
        for channelDict in HasFakeableDileptonTwoAk4JetsChannelList:
            # Dilepton #
            plots.extend(makeDileptonPlots(**channelDict))
            # Ak4 Jets #
            plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**HasFakeableDileptonTwoAk4JetsBaseDictJetsN))
            plots.extend(makeSeparateJetsPlots(**{k:channelDict[k] for k in commonItems},**HasFakeableDileptonTwoAk4JetsBaseDictJetsPlots))
            # MET #
            plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))


        #############################################################################
        #                            Resolved jets                                  #
        #############################################################################
        # Resolved selection #
        ElElHasFakeableDileptonMllCutOutZTightVetoInclusiveResolved = ElElHasFakeableDileptonMllCutOutZTightVetoTwoAk4Jets.refine("ElElHasFakeableDileptonMllCutOutZTightVetoInclusiveResolved", 
                                cut=[op.rng_len(self.ak4BJets) >= 1])
        MuMuHasFakeableDileptonMllCutOutZTightVetoInclusiveResolved = MuMuHasFakeableDileptonMllCutOutZTightVetoTwoAk4Jets.refine("MuMuHasFakeableDileptonMllCutOutZTightVetoInclusiveResolved", 
                                cut=[op.rng_len(self.ak4BJets) >= 1])
        ElMuHasFakeableDileptonMllCutTightVetoInclusiveResolved = ElMuHasFakeableDileptonMllCutTightVetoTwoAk4Jets.refine("ElMuHasFakeableDileptonMllCutTightVetoInclusiveResolved", 
                                cut=[op.rng_len(self.ak4BJets) >= 1])

        ElElHasFakeableDileptonMllCutOutZTightVetoExclusiveResolved = ElElHasFakeableDileptonMllCutOutZTightVetoTwoAk4Jets.refine("ElElHasFakeableDileptonMllCutOutZTightVetoExclusiveResolved", 
                                cut=[op.rng_len(self.ak4BJets)>=1,op.rng_len(self.ak8BJets)==0])
        MuMuHasFakeableDileptonMllCutOutZTightVetoExclusiveResolved = MuMuHasFakeableDileptonMllCutOutZTightVetoTwoAk4Jets.refine("MuMuHasFakeableDileptonMllCutOutZTightVetoExclusiveResolved", 
                                cut=[op.rng_len(self.ak4BJets)>=1,op.rng_len(self.ak8BJets)==0])
        ElMuHasFakeableDileptonMllCutTightVetoExclusiveResolved = ElMuHasFakeableDileptonMllCutTightVetoTwoAk4Jets.refine("ElMuHasFakeableDileptonMllCutTightVetoExclusiveResolved", 
                                cut=[op.rng_len(self.ak4BJets)>=1,op.rng_len(self.ak8BJets)==0])
                                
        # Yield 
        plots.append(yieldPlots.addYield(ElElHasFakeableDileptonMllCutOutZTightVetoInclusiveResolved,
                                         "ElElHasFakeableDileptonMllCutOutZTightVetoInclusiveResolved",
                                         "OS leptons + Inclusive Resolved Jets (channel $e^+e^-$)"))
        plots.append(yieldPlots.addYield(MuMuHasFakeableDileptonMllCutOutZTightVetoInclusiveResolved,
                                         "MuMuHasFakeableDileptonMllCutOutZTightVetoInclusiveResolved",
                                         "OS leptons + Inclusive Resolved Jets (channel $\mu^+\mu^-$)"))
        plots.append(yieldPlots.addYield(ElMuHasFakeableDileptonMllCutTightVetoInclusiveResolved,
                                         "ElMuHasFakeableDileptonMllCutTightVetoInclusiveResolved",
                                         "OS leptons + Inclusive Resolved Jets (channel $e^{\pm}\mu^{\mp}$)"))
        plots.append(yieldPlots.addYield(ElElHasFakeableDileptonMllCutOutZTightVetoExclusiveResolved,
                                         "ElElHasFakeableDileptonMllCutOutZTightVetoExclusiveResolved",
                                         "OS leptons + Exclusive Resolved Jets (channel $e^+e^-$)"))
        plots.append(yieldPlots.addYield(MuMuHasFakeableDileptonMllCutOutZTightVetoExclusiveResolved,
                                         "MuMuHasFakeableDileptonMllCutOutZTightVetoExclusiveResolved",
                                         "OS leptons + Exclusive Resolved Jets (channel $\mu^+\mu^-$)"))
        plots.append(yieldPlots.addYield(ElMuHasFakeableDileptonMllCutTightVetoExclusiveResolved,
                                         "ElMuHasFakeableDileptonMllCutTightVetoExclusiveResolved",
                                         "OS leptons + Exclusive Resolved Jets (channel $e^{\pm}\mu^{\mp}$)"))

        # Plots #
        HasFakeableDileptonResolvedJetsChannelList = [
             {'channel':'ElEl','sel':ElElHasFakeableDileptonMllCutOutZTightVetoInclusiveResolved,'dilepton':OsElElFakeable[0],'suffix':'HasFakeableDileptonMllCutOutZTightVetoInclusiveResolved'},
             {'channel':'MuMu','sel':MuMuHasFakeableDileptonMllCutOutZTightVetoInclusiveResolved,'dilepton':OsMuMuFakeable[0],'suffix':'HasFakeableDileptonMllCutOutZTightVetoInclusiveResolved'},
             {'channel':'ElMu','sel':ElMuHasFakeableDileptonMllCutTightVetoInclusiveResolved,    'dilepton':OsElMuFakeable[0],'suffix':'HasFakeableDileptonMllCutTightVetoInclusiveResolved'},
             {'channel':'ElEl','sel':ElElHasFakeableDileptonMllCutOutZTightVetoExclusiveResolved,'dilepton':OsElElFakeable[0],'suffix':'HasFakeableDileptonMllCutOutZTightVetoExclusiveResolved'},
             {'channel':'MuMu','sel':MuMuHasFakeableDileptonMllCutOutZTightVetoExclusiveResolved,'dilepton':OsMuMuFakeable[0],'suffix':'HasFakeableDileptonMllCutOutZTightVetoExclusiveResolved'},
             {'channel':'ElMu','sel':ElMuHasFakeableDileptonMllCutTightVetoExclusiveResolved,    'dilepton':OsElMuFakeable[0],'suffix':'HasFakeableDileptonMllCutTightVetoExclusiveResolved'},
                                                  ]
        commonItems = ['channel','sel','suffix']
        HasFakeableDileptonResolvedJetsBaseDictJetsPlots = {'bjets':self.ak4BJets, 'lightjets':self.ak4LightJets, 'alljets': self.ak4Jets}
        HasFakeableDileptonResolvedJetsBaseDictJetsN = {'objName':'Ak4BJets','objCont':self.ak4BJets,'Nmax':5,'xTitle':'N(Ak4 Bjets)'}

        for channelDict in HasFakeableDileptonResolvedJetsChannelList:
            # Dilepton #
            plots.extend(makeDileptonPlots(**channelDict))
            # Ak4 Jets #
            plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**HasFakeableDileptonResolvedJetsBaseDictJetsN))
            plots.extend(makeAk4BJetsPlots(**{k:channelDict[k] for k in commonItems},**HasFakeableDileptonResolvedJetsBaseDictJetsPlots))
            # MET #
            plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))
 
        #############################################################################
        #                              AK8 Jets                                     #
        #############################################################################
        # Jet Selection #
        ElElHasFakeableDileptonMllCutOutZTightVetoOneAk8Jet = ElElHasFakeableDileptonMllCutOutZTightVeto.refine("ElElHasFakeableDileptonMllCutOutZTightVetoOneAk8Jet", 
                                cut=[op.rng_len(self.ak8Jets)>=1])
        MuMuHasFakeableDileptonMllCutOutZTightVetoOneAk8Jet = MuMuHasFakeableDileptonMllCutOutZTightVeto.refine("MuMuHasFakeableDileptonMllCutOutZTightVetoOneAk8Jet", 
                                cut=[op.rng_len(self.ak8Jets)>=1])
        ElMuHasFakeableDileptonMllCutTightVetoOneAk8Jet     = ElMuHasFakeableDileptonMllCutTightVeto.refine("ElMuHasFakeableDileptonMllCutTightVetoOneAk8Jet", 
                                cut=[op.rng_len(self.ak8Jets)>=1])

        # Yield 
        plots.append(yieldPlots.addYield(ElElHasFakeableDileptonMllCutOutZTightVetoOneAk8Jet,
                                         "ElElHasFakeableDileptonMllCutOutZTightVetoOneAk8Jet",
                                         "OS leptons + Ak4 Jets (channel $e^+e^-$)"))
        plots.append(yieldPlots.addYield(MuMuHasFakeableDileptonMllCutOutZTightVetoOneAk8Jet,
                                         "MuMuHasFakeableDileptonMllCutOutZTightVetoOneAk8Jet",
                                         "OS leptons + Ak4 Jets (channel $\mu^+\mu^-$)"))
        plots.append(yieldPlots.addYield(ElMuHasFakeableDileptonMllCutTightVetoOneAk8Jet,
                                         "ElMuHasFakeableDileptonMllCutTightVetoOneAk8Jet",
                                         "OS leptons + Ak4 Jets (channel $e^{\pm}\mu^{\mp}$)"))

        # Plots #
        HasFakeableDileptonOneAk8JetJetsChannelList = [
                             {'channel':'ElEl','sel':ElElHasFakeableDileptonMllCutOutZTightVetoOneAk8Jet,'dilepton':OsElElFakeable[0],'suffix':'HasFakeableDileptonMllCutOutZTightVetoOneAk8Jet'},
                             {'channel':'MuMu','sel':MuMuHasFakeableDileptonMllCutOutZTightVetoOneAk8Jet,'dilepton':OsMuMuFakeable[0],'suffix':'HasFakeableDileptonMllCutOutZTightVetoOneAk8Jet'},
                             {'channel':'ElMu','sel':ElMuHasFakeableDileptonMllCutTightVetoOneAk8Jet,    'dilepton':OsElMuFakeable[0],'suffix':'HasFakeableDileptonMllCutTightVetoOneAk8Jet'},
                                                  ]
        commonItems = ['channel','sel','suffix']
        HasFakeableDileptonOneAk8JetsBaseDictFatjetsN = {'objName':'Ak8Jets','objCont':self.ak8Jets,'Nmax':5,'xTitle':'N(Ak8 jets)'}

                                                    
        for channelDict in HasFakeableDileptonOneAk8JetJetsChannelList:
            # Dilepton #
            plots.extend(makeDileptonPlots(**channelDict))
            # Ak8 Jets #
            plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**HasFakeableDileptonOneAk8JetsBaseDictFatjetsN))
            plots.extend(makeAk8Plots(fatjet=self.ak8Jets[0],**{k:channelDict[k] for k in commonItems}))
            # MET #
            plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))


        #############################################################################
        #                             Boosted jets                                  #
        #############################################################################
        # Boosted selection #
        ElElHasFakeableDileptonMllCutOutZTightVetoInclusiveBoosted = ElElHasFakeableDileptonMllCutOutZTightVetoOneAk8Jet.refine("ElElHasFakeableDileptonMllCutOutZTightVetoInclusiveBoosted", 
                                cut=[op.rng_len(self.ak8BJets)>=1])
        MuMuHasFakeableDileptonMllCutOutZTightVetoInclusiveBoosted = MuMuHasFakeableDileptonMllCutOutZTightVetoOneAk8Jet.refine("MuMuHasFakeableDileptonMllCutOutZTightVetoInclusiveBoosted", 
                                cut=[op.rng_len(self.ak8BJets)>=1])
        ElMuHasFakeableDileptonMllCutTightVetoInclusiveBoosted     = ElMuHasFakeableDileptonMllCutTightVetoOneAk8Jet.refine("ElMuHasFakeableDileptonMllCutTightVetoInclusiveBoosted", 
                                cut=[op.rng_len(self.ak8BJets)>=1])

        ElElHasFakeableDileptonMllCutOutZTightVetoExclusiveBoosted = ElElHasFakeableDileptonMllCutOutZTightVetoOneAk8Jet.refine("ElElHasFakeableDileptonMllCutOutZTightVetoExclusiveBoosted", 
                                cut=[op.rng_len(self.ak8BJets)>=1, op.OR(op.rng_len(self.ak4Jets)<=1 , op.rng_len(self.ak4BJets)==0)])
        MuMuHasFakeableDileptonMllCutOutZTightVetoExclusiveBoosted = MuMuHasFakeableDileptonMllCutOutZTightVetoOneAk8Jet.refine("MuMuHasFakeableDileptonMllCutOutZTightVetoExclusiveBoosted", 
                                cut=[op.rng_len(self.ak8BJets)>=1, op.OR(op.rng_len(self.ak4Jets)<=1 , op.rng_len(self.ak4BJets)==0)])
        ElMuHasFakeableDileptonMllCutTightVetoExclusiveBoosted     = ElMuHasFakeableDileptonMllCutTightVetoOneAk8Jet.refine("ElMuHasFakeableDileptonMllCutTightVetoExclusiveBoosted", 
                                cut=[op.rng_len(self.ak8BJets)>=1, op.OR(op.rng_len(self.ak4Jets)<=1 , op.rng_len(self.ak4BJets)==0)])

        # Yield 
        plots.append(yieldPlots.addYield(ElElHasFakeableDileptonMllCutOutZTightVetoInclusiveBoosted,
                                         "ElElHasFakeableDileptonMllCutOutZTightVetoInclusiveBoosted",
                                         "OS leptons + Inclusive Boosted Jets (channel $e^+e^-$)"))
        plots.append(yieldPlots.addYield(MuMuHasFakeableDileptonMllCutOutZTightVetoInclusiveBoosted,
                                         "MuMuHasFakeableDileptonMllCutOutZTightVetoInclusiveBoosted",
                                         "OS leptons + Inclusive Boosted Jets (channel $\mu^+\mu^-$)"))
        plots.append(yieldPlots.addYield(ElMuHasFakeableDileptonMllCutTightVetoInclusiveBoosted,
                                         "ElMuHasFakeableDileptonMllCutTightVetoInclusiveBoosted",
                                         "OS leptons + Inclusive Boosted Jets (channel $e^{\pm}\mu^{\mp}$)"))
        plots.append(yieldPlots.addYield(ElElHasFakeableDileptonMllCutOutZTightVetoExclusiveBoosted,
                                         "ElElHasFakeableDileptonMllCutOutZTightVetoExclusiveBoosted",
                                         "OS leptons + Exclusive Boosted Jets (channel $e^+e^-$)"))
        plots.append(yieldPlots.addYield(MuMuHasFakeableDileptonMllCutOutZTightVetoExclusiveBoosted,
                                         "MuMuHasFakeableDileptonMllCutOutZTightVetoExclusiveBoosted",
                                         "OS leptons + Exclusive Boosted Jets (channel $\mu^+\mu^-$)"))
        plots.append(yieldPlots.addYield(ElMuHasFakeableDileptonMllCutTightVetoExclusiveBoosted,
                                         "ElMuHasFakeableDileptonMllCutTightVetoExclusiveBoosted",
                                         "OS leptons + Exclusive Boosted Jets (channel $e^{\pm}\mu^{\mp}$)"))

        # Plots #
        HasFakeableDileptonBoostedJetsChannelList = [
             {'channel':'ElEl','sel':ElElHasFakeableDileptonMllCutOutZTightVetoInclusiveBoosted,'dilepton':OsElElFakeable[0],'suffix':'HasFakeableDileptonMllCutOutZTightVetoInclusiveBoosted'},
             {'channel':'MuMu','sel':MuMuHasFakeableDileptonMllCutOutZTightVetoInclusiveBoosted,'dilepton':OsMuMuFakeable[0],'suffix':'HasFakeableDileptonMllCutOutZTightVetoInclusiveBoosted'},
             {'channel':'ElMu','sel':ElMuHasFakeableDileptonMllCutTightVetoInclusiveBoosted,    'dilepton':OsElMuFakeable[0],'suffix':'HasFakeableDileptonMllCutTightVetoInclusiveBoosted'},
             {'channel':'ElEl','sel':ElElHasFakeableDileptonMllCutOutZTightVetoExclusiveBoosted,'dilepton':OsElElFakeable[0],'suffix':'HasFakeableDileptonMllCutOutZTightVetoExclusiveBoosted'},
             {'channel':'MuMu','sel':MuMuHasFakeableDileptonMllCutOutZTightVetoExclusiveBoosted,'dilepton':OsMuMuFakeable[0],'suffix':'HasFakeableDileptonMllCutOutZTightVetoExclusiveBoosted'},
             {'channel':'ElMu','sel':ElMuHasFakeableDileptonMllCutTightVetoExclusiveBoosted,    'dilepton':OsElMuFakeable[0],'suffix':'HasFakeableDileptonMllCutTightVetoExclusiveBoosted'},
                                                  ]
        commonItems = ['channel','sel','suffix']
        HasFakeableDileptonBoostedJetsBaseDictJetsN = {'objName':'Ak8BJets','objCont':self.ak8BJets,'Nmax':5,'xTitle':'N(Ak8 Bjets)'}

        for channelDict in HasFakeableDileptonBoostedJetsChannelList:
            # Dilepton #
            plots.extend(makeDileptonPlots(**channelDict))
            # Ak8 Jets #
            plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**HasFakeableDileptonBoostedJetsBaseDictJetsN))
            plots.extend(makeAk8Plots(fatjet=self.ak8BJets[0],**{k:channelDict[k] for k in commonItems}))
            # MET #
            plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))



        #############################################################################
        #                     High-level combinations                               #
        #############################################################################
        highLevelBaseQuantities = { # No need to repeat for each channelDict
             'met'          : self.corrMET,
             'jets'         : self.ak4Jets, 
             'resolvedjets' : self.ak4BJets,
             'lightjets'    : self.ak4LightJets,
             'boostedjets'  : self.ak8BJets,
                                  }
 
        HasOsMllCutHighLevelVariablesChannelList = [
            # Resolved #
            {'channel'      :'ElEl',
             'dilepton'     : OsElElFakeable[0],
             'sel'          : ElElHasFakeableDileptonMllCutOutZTightVetoExclusiveResolved,
             'suffix'       : 'HasFakeableDileptonMllCutOutZTightVetoExclusiveResolved'},
            {'channel'      :'MuMu',
             'dilepton'     : OsMuMuFakeable[0],
             'sel'          : MuMuHasFakeableDileptonMllCutOutZTightVetoExclusiveResolved,
             'suffix'       : 'MuMuHasFakeableDileptonMllCutOutZTightVetoExclusiveResolved'},
            {'channel'      :'ElMu',
             'dilepton'     : OsElMuFakeable[0],
             'sel'          : ElMuHasFakeableDileptonMllCutTightVetoExclusiveResolved,
             'suffix'       : 'ElMuHasFakeableDileptonMllCutTightVetoExclusiveResolved'},
            # Boosted #
            {'channel'      :'ElEl',
             'dilepton'     : OsElElFakeable[0],
             'sel'          : ElElHasFakeableDileptonMllCutOutZTightVetoInclusiveBoosted,
             'suffix'       : 'HasFakeableDileptonMllCutOutZTightVetoInclusiveBoosted'},
            {'channel'      :'MuMu',
             'dilepton'     : OsMuMuFakeable[0],
             'sel'          : MuMuHasFakeableDileptonMllCutOutZTightVetoInclusiveBoosted,
             'suffix'       : 'MuMuHasFakeableDileptonMllCutOutZTightVetoInclusiveBoosted'},
            {'channel'      :'ElMu',
             'dilepton'     : OsElMuFakeable[0],
             'sel'          : ElMuHasFakeableDileptonMllCutTightVetoInclusiveBoosted,
             'suffix'       : 'ElMuHasFakeableDileptonMllCutTightVetoInclusiveBoosted'}
                                                    ]

        for channelDict in HasOsMllCutHighLevelVariablesChannelList:
            plots.extend(makeHighLevelQuantities(**highLevelBaseQuantities,**channelDict))

        return plots

