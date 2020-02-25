import os
import sys

from bamboo import treefunctions as op
from bamboo.plots import Plot, EquidistantBinning, SummedPlot

sys.path.append('/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/') # Add scripts in this directory -- TODO : make cleaner
from BaseHHtobbWW import BaseNanoHHtobbWW
from plotDef import makeDileptonPlots, makeJetsPlots, makeFatJetPlots, makeMETPlots, makeDeltaRPlots, makeYieldPlot, makeHighLevelQuantities

#===============================================================================================#
#                                 PlotterHHtobbWW                                               #
#===============================================================================================#
class PlotterNanoHHtobbWW(BaseNanoHHtobbWW):
    """ Plotter module: HH->bbW(->e/µ nu)W(->e/µ nu) histograms from NanoAOD """
    def __init__(self, args):
        super(PlotterNanoHHtobbWW, self).__init__(args)

    def prepareTree(self, tree, sample=None, sampleCfg=None, enableSystematics=None):
        return super(PlotterNanoHHtobbWW,self).prepareTree(tree, sample, sampleCfg, enableSystematics)

    def definePlots(self, t, noSel, sample=None, sampleCfg=None): 
        noSel = super(PlotterNanoHHtobbWW,self).preparePlots(t, noSel, sample, sampleCfg)

        plots = []

        era = sampleCfg['era']

        isMC = self.isMC(sample)
        ############################################################
        #                    Dilepton plots                        #
        ############################################################

        # Plots Numbers of dilepton in each channel #
        plots.append(Plot.make1D("ElEl_channel",op.rng_len(self.OsElEl),noSel,EquidistantBinning(5,0,5.),title='Number of dilepton events in ElEl channel',xTitle='N_{dilepton} (ElEl channel)'))
        plots.append(Plot.make1D("MuMu_channel",op.rng_len(self.OsMuMu),noSel,EquidistantBinning(5,0,5.),title='Number of dilepton events in MuMu channel',xTitle='N_{dilepton} (MuMu channel)'))
        plots.append(Plot.make1D("ElMu_channel",op.rng_len(self.OsElMu),noSel,EquidistantBinning(5,0,5.),title='Number of dilepton events in ElMu channel',xTitle='N_{dilepton} (ElMu channel)'))

        # Selection #
        hasOsElEl = noSel.refine("hasOsElEl",
                                 cut = [op.rng_len(self.OsElEl) >= 1,                        # Require at least one dilepton ElEl
                                        (op.rng_len(self.electrons)+op.rng_len(self.muons))<=2],  # Not more than two tight leptons
                                 weight = self.llSFApplied["ElEl"])
        hasOsMuMu = noSel.refine("hasOsMuMu",
                                 cut = [op.rng_len(self.OsMuMu) >= 1,                        # Require at least one dilepton MuMu
                                        (op.rng_len(self.electrons)+op.rng_len(self.muons))<=2],  # Not more than two tight leptons
                                 weight = self.llSFApplied["MuMu"])
        hasOsElMu = noSel.refine("hasOsElMu",
                                 cut = [op.rng_len(self.OsElMu) >= 1,                        # Require at least one dilepton ElMu
                                        (op.rng_len(self.electrons)+op.rng_len(self.muons))<=2],  # Not more than two tight leptons
                                 weight = self.llSFApplied["ElMu"])

        # Yield plots #
        plots.append(makeYieldPlot(self,hasOsElEl,"OSElEl","OS leptons (channel $e^+e^-$)",0))
        plots.append(makeYieldPlot(self,hasOsMuMu,"OSMuMu","OS leptons (channel : $\mu^+\mu^-$)",1))
        plots.append(makeYieldPlot(self,hasOsElMu,"OSMuEl","OS leptons (channel $e^{\pm}\mu^{\mp}$)",2))

        # Dilepton channel plot #

        hasOsCutChannelList = [
                             {'channel':'ElEl','sel':hasOsElEl,'dilepton':self.OsElEl[0],'suffix':'hasOsElEl'},
                             {'channel':'MuMu','sel':hasOsMuMu,'dilepton':self.OsMuMu[0],'suffix':'hasOsMuMu'},
                             {'channel':'ElMu','sel':hasOsElMu,'dilepton':self.OsElMu[0],'suffix':'hasOsElMu'},
                                       ]
        for channelDict in hasOsCutChannelList:
            # Dilepton plots #
            plots.extend(makeDileptonPlots(self, **channelDict))
            # MET plots #
            plots.extend(makeMETPlots(self      = self,
                                      sel       = channelDict['sel'],
                                      met       = self.corrMET,
                                      suffix    = channelDict['suffix'],
                                      channel   = channelDict['channel']))

        # Dilepton Z peak exclusion (charge already done in previous selection) # 
        lambda_lowMllCut    = lambda dilep: op.invariant_mass(dilep[0].p4, dilep[1].p4)>12.
        lambda_outZ         = lambda dilep: op.NOT(op.in_range(80.,op.invariant_mass(dilep[0].p4, dilep[1].p4),100.))

        hasOsElElLowMllCutOutZ = hasOsElEl.refine("hasOsElElLowMllCutOutZ",cut=[lambda_lowMllCut(self.OsElEl[0]),lambda_outZ(self.OsElEl[0])])
        hasOsMuMuLowMllCutOutZ = hasOsMuMu.refine("hasOsMuMuLowMllCutOutZ",cut=[lambda_lowMllCut(self.OsMuMu[0]),lambda_outZ(self.OsMuMu[0])])
        hasOsElMuLowMllCut     = hasOsElMu.refine("hasOsElMuLowMllCut",cut=[lambda_lowMllCut(self.OsElMu[0])]) # Z peak cut not needed because Opposite Flavour

        # Yield plots #
        plots.append(makeYieldPlot(self,hasOsElElLowMllCutOutZ,"OSElElMllCutOutZ","OS leptons + $M_{ll}$ (channel : $e^+e^-$)",3))
        plots.append(makeYieldPlot(self,hasOsMuMuLowMllCutOutZ,"OSMuMuMllCutOutZ","OS leptons + $M_{ll}$ (channel : $\mu^+\mu^-$)",4))
        plots.append(makeYieldPlot(self,hasOsElMuLowMllCut,"OSMuMuMllCut","OS leptons + $M_{ll}$ (channel : $e^{\pm}\mu^{\mp}$)",5))

        # Dilepton plots #

        hasOsMllCutChannelList = [
                             {'channel':'ElEl','sel':hasOsElElLowMllCutOutZ,'dilepton':self.OsElEl[0],'suffix':'hasOsElElLowMllCutOutZ'},
                             {'channel':'MuMu','sel':hasOsMuMuLowMllCutOutZ,'dilepton':self.OsMuMu[0],'suffix':'hasOsMuMuLowMllCutOutZ'},
                             {'channel':'ElMu','sel':hasOsElMuLowMllCut,    'dilepton':self.OsElMu[0],'suffix':'hasOsElMuLowMllCut'},
                           ]

        for channelDict in hasOsMllCutChannelList:
            # Dilepton plots #
            plots.extend(makeDileptonPlots(self, **channelDict))
            # MET plots #
            plots.extend(makeMETPlots(self      = self,
                                      sel       = channelDict['sel'],
                                      met       = self.corrMET,
                                      suffix    = channelDict['suffix'],
                                      channel   = channelDict['channel']))

        ############################################################
        #                        Jets plots                        #
        ############################################################
        # Plot lepton/jet angle separartion
        DeltaRChannelList = [
                     {'isMC':isMC,'channel':'ElEl','sel':hasOsElElLowMllCutOutZ,'cont1':self.electrons,'cont2':self.jetsSel,'suffix':'hasOsElElLowMllCutOutZ_ElectronJet'},
                     {'isMC':isMC,'channel':'MuMu','sel':hasOsMuMuLowMllCutOutZ,'cont1':self.electrons,'cont2':self.jetsSel,'suffix':'hasOsMuMuLowMllCutOutZ_ElectronJet'},
                     {'isMC':isMC,'channel':'ElMu','sel':hasOsElMuLowMllCut,    'cont1':self.electrons,'cont2':self.jetsSel,'suffix':'hasOsElMuLowMllCut_ElectronJet'},
                     {'isMC':isMC,'channel':'ElEl','sel':hasOsElElLowMllCutOutZ,'cont1':self.electrons,'cont2':self.fatjetsSel,'suffix':'hasOsElElLowMllCutOutZ_ElectronFatjet'},
                     {'isMC':isMC,'channel':'MuMu','sel':hasOsMuMuLowMllCutOutZ,'cont1':self.electrons,'cont2':self.fatjetsSel,'suffix':'hasOsMuMuLowMllCutOutZ_ElectronFatjet'},
                     {'isMC':isMC,'channel':'ElMu','sel':hasOsElMuLowMllCut,    'cont1':self.electrons,'cont2':self.fatjetsSel,'suffix':'hasOsElMuLowMllCut_ElectronFatjet'},
                     {'isMC':isMC,'channel':'ElEl','sel':hasOsElElLowMllCutOutZ,'cont1':self.muons,'cont2':self.jetsSel,'suffix':'hasOsElElLowMllCutOutZ_MuonJet'},
                     {'isMC':isMC,'channel':'MuMu','sel':hasOsMuMuLowMllCutOutZ,'cont1':self.muons,'cont2':self.jetsSel,'suffix':'hasOsMuMuLowMllCutOutZ_MuonJet'},
                     {'isMC':isMC,'channel':'ElMu','sel':hasOsElMuLowMllCut,    'cont1':self.muons,'cont2':self.jetsSel,'suffix':'hasOsElMuLowMllCut_MuonJet'},
                     {'isMC':isMC,'channel':'ElEl','sel':hasOsElElLowMllCutOutZ,'cont1':self.muons,'cont2':self.fatjetsSel,'suffix':'hasOsElElLowMllCutOutZ_MuonFatjet'},
                     {'isMC':isMC,'channel':'MuMu','sel':hasOsMuMuLowMllCutOutZ,'cont1':self.muons,'cont2':self.fatjetsSel,'suffix':'hasOsMuMuLowMllCutOutZ_MuonFatjet'},
                     {'isMC':isMC,'channel':'ElMu','sel':hasOsElMuLowMllCut,    'cont1':self.muons,'cont2':self.fatjetsSel,'suffix':'hasOsElMuLowMllCut_MuonFatjet'},
                   ]
        for channelDict in DeltaRChannelList:
            plots.extend(makeDeltaRPlots(self, **channelDict))

        # Yield plots #
        hasOsElElLowMllCutOutZAtLeast2Jets = hasOsElElLowMllCutOutZ.refine("hasOsElElLowMllCutOutZAtLeast2Jets",cut=[op.rng_len(self.jets)>=2])
        hasOsMuMuLowMllCutOutZAtLeast2Jets = hasOsMuMuLowMllCutOutZ.refine("hasOsMuMuLowMllCutOutZAtLeast2Jets",cut=[op.rng_len(self.jets)>=2])
        hasOsElMuLowMllCutAtLeast2Jets     = hasOsElMuLowMllCut.refine("hasOsElMuLowMllCutOutZAtLeast2Jets",cut=[op.rng_len(self.jets)>=2])

        hasOsElElLowMllCutOutZAtLeast1Fatjet = hasOsElElLowMllCutOutZ.refine("hasOsElElLowMllCutOutZAtLeast1Fatjet",cut=[op.rng_len(self.fatjets)>=1])
        hasOsMuMuLowMllCutOutZAtLeast1Fatjet = hasOsMuMuLowMllCutOutZ.refine("hasOsMuMuLowMllCutOutZAtLeast1Fatjet",cut=[op.rng_len(self.fatjets)>=1])
        hasOsElMuLowMllCutAtLeast1Fatjet     = hasOsElMuLowMllCut.refine("hasOsElMuLowMllCutOutZAtLeast1Fatjet",cut=[op.rng_len(self.fatjets)>=1])
        
        plots.append(makeYieldPlot(self,hasOsElElLowMllCutOutZAtLeast2Jets,"OSElElMllCutOutZAtLeast2Jets","OS leptons + $M_{ll}$ + Jets $\geq 2$ (channel : $e^+e^-$)",6))
        plots.append(makeYieldPlot(self,hasOsMuMuLowMllCutOutZAtLeast2Jets,"OSMuMuMllCutOutZAtLeast2Jets","OS leptons + $M_{ll}$ + Jets $\geq 2$ (channel : $\mu^+\mu^-$)",7))
        plots.append(makeYieldPlot(self,hasOsElMuLowMllCutAtLeast2Jets,"OSElMuMllCutAtLeast2Jets","OS leptons + $M_{ll}$ + Jets $\geq 2$ (channel : $e^{\pm}\mu^{\mp}$",8))

        plots.append(makeYieldPlot(self,hasOsElElLowMllCutOutZAtLeast1Fatjet,"OSElElMllCutOutZAtLeast1Fatjet","OS leptons + $M_{ll}$ + Fatjets $\geq 1$ (channel : $e^+e^-$)",9))
        plots.append(makeYieldPlot(self,hasOsMuMuLowMllCutOutZAtLeast1Fatjet,"OSMuMuMllCutOutZAtLeast1Fatjet","OS leptons + $M_{ll}$ + Fatjets $\geq 1$ (channel : $\mu^+\mu^-$)",10))
        plots.append(makeYieldPlot(self,hasOsElMuLowMllCutAtLeast1Fatjet,"OSElMuMllCutAtLeast1Fatjet","OS leptons + $M_{ll}$ + Fatjets $\geq 1$ (channel : $e^{\pm}\mu^{\mp}$)",11))

   
        # Define the boosted and Resolved (+exclusive) selections #
        hasBoostedJets = noSel.refine("hasBoostedJets",
                               cut=[op.rng_len(self.bjetsBoosted)>=1],
                               weight = self.DeepCSVMediumSFApplied)
        hasResolvedJets = noSel.refine("hasResolvedJets",
                               cut=[op.rng_len(self.jets)>=2,op.rng_len(self.bjetsResolved)>=1],
                               weight = self.DeepJetMediumSFApplied)
        hasExclusiveResolvedJets = noSel.refine("hasExclusiveResolved",
                               cut=[op.rng_len(self.jets)>=2, op.rng_len(self.bjetsResolved)>=1,op.rng_len(self.bjetsBoosted)==0],
                               weight = self.DeepJetMediumSFApplied)
        hasExclusiveBoostedJets = noSel.refine("hasExclusiveBoostedJets",
                               cut=[op.rng_len(self.bjetsBoosted)>=1,op.OR(op.rng_len(self.jets)<=1,op.rng_len(self.bjetsResolved)==0)],
                               weight = self.DeepCSVMediumSFApplied)

        hasNotBoostedJets = noSel.refine("hasNotBoostedJets",
                               cut=[op.rng_len(self.bjetsBoosted)==0])
        hasNotResolvedJets = noSel.refine("hasNotResolvedJets",
                               cut=[op.OR(op.rng_len(self.jets)<=1,op.rng_len(self.bjetsResolved)==0)])
        hasBoostedAndResolvedJets = noSel.refine("hasBoostedAndResolvedJets", 
                               cut=[op.rng_len(self.bjetsBoosted)>=1,op.rng_len(self.jets)>=2,op.rng_len(self.bjetsResolved)>=1])
        hasNotBoostedAndResolvedJets = noSel.refine("hasNotBoostedAndResolvedJets", 
                               cut=[op.OR(op.rng_len(self.bjetsBoosted)==0,op.rng_len(self.jets)<=1,op.rng_len(self.bjetsResolved)==0)])
        hasNotExclusiveResolvedJets = noSel.refine("hasNotExclusiveResolved",
                               cut=[op.OR(op.OR(op.rng_len(self.jets)<=1,op.rng_len(self.bjetsResolved)==0),op.AND(op.rng_len(self.bjetsBoosted)>=1,op.rng_len(self.jets)>=2,op.rng_len(self.bjetsResolved)>=1))])
        hasNotExclusiveBoostedJets = noSel.refine("hasNotExclusiveBoostedJets",
                               cut=[op.OR(op.rng_len(self.bjetsBoosted)==0,op.AND(op.rng_len(self.jets)>=2,op.rng_len(self.bjetsResolved)>=1))])
        # Note that these selection should be done the same (and SF) when combining the OS lepton selection #
        # Because so far there is no way to "concatenate" two refine's #


        # Scalefactors plot #
#        plots.append(Plot.make1D("DeepJetMediumSF_ResolvedJets",
#                                self.DeepJetMediumSF(self.bjetsResolved[0]),
#                                hasResolvedJets,
#                                EquidistantBinning(100,0.,2.),
#                                title='DeepJetMediumSF',
#                                xTitle='DeepJetMediumSF'))
#        plots.append(Plot.make1D("DeepJetMediumSF_ExclusiveResolvedJets",
#                                self.DeepJetMediumSF(self.bjetsResolved[0]),
#                                hasExclusiveResolvedJets,
#                                EquidistantBinning(100,0.,2.),
#                                title='DeepJetMediumSF',
#                                xTitle='DeepJetMediumSF'))

        # Counting events from different selections for debugging #
        # Passing Boosted selection #
        PassedBoosted = Plot.make1D("PassedBoosted",
                                    op.c_int(1),
                                    hasBoostedJets,
                                    EquidistantBinning(2,0.,2.),
                                    title='Passed Boosted',
                                    xTitle='Passed Boosted')
        FailedBoosted = Plot.make1D("FailedBoosted",
                                    op.c_int(0),
                                    hasNotBoostedJets,
                                    EquidistantBinning(2,0.,2.),
                                    title='Failed Boosted',
                                    xTitle='Failed Boosted')
        plots.append(SummedPlot("BoostedCase",
                                [FailedBoosted,PassedBoosted],
                                xTitle="Boosted selection"))

        # Passing Resolved selection #
        PassedResolved = Plot.make1D("PassedResolved",
                                     op.c_int(1),
                                     hasResolvedJets,
                                     EquidistantBinning(2,0.,2.),
                                     title='Passed Resolved',
                                     xTitle='Passed Resolved')
        FailedResolved = Plot.make1D("FailedResolved",
                                     op.c_int(0),
                                     hasNotResolvedJets,
                                     EquidistantBinning(2,0.,2.),
                                     title='Failed Resolved',
                                     xTitle='Failed Resolved')
        plots.append(SummedPlot("ResolvedCase",
                                [FailedResolved,PassedResolved],
                                xTitle="Resolved selection"))

        # Passing Exclusive Resolved (Resolved AND NOT Boosted) #
        PassedExclusiveResolved = Plot.make1D("PassedExclusiveResolved",
                                              op.c_int(1),
                                              hasExclusiveResolvedJets,
                                              EquidistantBinning(2,0.,2.),
                                              title='Passed Exclusive Resolved',
                                              xTitle='Passed Exclusive Resolved')
        FailedExclusiveResolved = Plot.make1D("FailedExclusiveResolved",
                                              op.c_int(0),
                                              hasNotExclusiveResolvedJets,
                                              EquidistantBinning(2,0.,2.),
                                              title='Failed Exclusive Resolved',
                                              xTitle='Failed Exclusive Resolved')
        plots.append(SummedPlot("ExclusiveResolvedCase",
                                [FailedExclusiveResolved,PassedExclusiveResolved],
                                xTitle="Exclusive Resolved selection"))

        # Passing Exclusive Boosted (Boosted AND NOT Resolved) #
        PassedExclusiveBoosted = Plot.make1D("PassedExclusiveBoosted",
                                             op.c_int(1),
                                             hasExclusiveBoostedJets,
                                             EquidistantBinning(2,0.,2.),
                                             title='Passed Exclusive Boosted',
                                             xTitle='Passed Exclusive Boosted')
        FailedExclusiveBoosted = Plot.make1D("FailedExclusiveBoosted",
                                             op.c_int(0),
                                             hasNotExclusiveBoostedJets,
                                             EquidistantBinning(2,0.,2.),
                                             title='Failed Exclusive Boosted',
                                             xTitle='Failed Exclusive Boosted')
        plots.append(SummedPlot("ExclusiveBoostedCase",
                                [FailedExclusiveBoosted,PassedExclusiveBoosted],
                                xTitle="Exclusive Boosted selection"))

        # Passing Boosted AND Resolved #
        PassedBoth = Plot.make1D("PassedBoth",
                                 op.c_int(1),
                                 hasBoostedAndResolvedJets,
                                 EquidistantBinning(2,0.,2.),
                                 title='Passed Both Boosted and Resolved',
                                 xTitle='Passed Boosted and Resolved')
        FailedBoth = Plot.make1D("FailedBoth", # Means failed the (Boosted AND Resolved) = either one or the other 
                                 op.c_int(0),
                                 hasNotBoostedAndResolvedJets,
                                 EquidistantBinning(2,0.,2.),
                                 title='Failed combination Boosted and Resolved',
                                 xTitle='Failed combination')
        plots.append(SummedPlot("BoostedAndResolvedCase",[FailedBoth,PassedBoth],xTitle="Boosted and Resolved selection"))

        # Plot number of subjets in the boosted fatjets #
        lambda_noSubjet  = lambda fatjet : op.AND(fatjet.subJet1._idx.result == -1, op.AND(fatjet.subJet2._idx.result == -1 )) 
        lambda_oneSubjet = lambda fatjet : op.AND(fatjet.subJet1._idx.result != -1, op.AND(fatjet.subJet2._idx.result == -1 ))
        lambda_twoSubjet = lambda fatjet : op.AND(fatjet.subJet1._idx.result != -1, op.AND(fatjet.subJet2._idx.result != -1 ))

        hasNoSubjet   = hasBoostedJets.refine("hasNoSubjet",
                       cut=[lambda_noSubjet(self.bjetsBoosted[0])])
        hasOneSubjet  = hasBoostedJets.refine("hasOneSubjet",
                       cut=[lambda_oneSubjet(self.bjetsBoosted[0])])
        hasTwoSubjet  = hasBoostedJets.refine("hasTwoSubjet",
                       cut=[lambda_twoSubjet(self.bjetsBoosted[0])])

        plot_hasNoSubjet = Plot.make1D("plot_hasNoSubjet", # Fill bin 0
                                       op.c_int(0),
                                       hasNoSubjet,
                                       EquidistantBinning(3,0.,3.),
                                       title='Boosted jet without subjet')
        plot_hasOneSubjet = Plot.make1D("plot_hasOneSubjet", # Fill bin 1
                                        op.c_int(1),
                                        hasOneSubjet,
                                        EquidistantBinning(3,0.,3.),
                                        title='Boosted jet with one subjet')
        plot_hasTwoSubjet = Plot.make1D("plot_hasTwoSubjet", # Fill bin 2
                                        op.c_int(2),
                                        hasTwoSubjet,
                                        EquidistantBinning(3,0.,3.),title='Boosted jet with two subjets')
        plots.append(SummedPlot("NumberOfSubjets",
                                [plot_hasNoSubjet,plot_hasOneSubjet,plot_hasTwoSubjet],
                                xTitle="Number of subjets in boosted jet"))


        #############################################################################
        #                     Jets + Dilepton combination                           #
        #############################################################################

        ##### BOOSTED #####
        # Combine dilepton and Boosted selections #
        hasOsElElLowMllCutOutZBoostedJets = hasOsElElLowMllCutOutZ.refine("hasOsElElLowMllCutOutZBoostedJets", 
                                      cut=[op.rng_len(self.bjetsBoosted)>=1],
                                      weight = self.DeepCSVMediumSFApplied)
        hasOsMuMuLowMllCutOutZBoostedJets = hasOsMuMuLowMllCutOutZ.refine("hasOsMuMuLowMllCutOutZBoostedJets", 
                                      cut=[op.rng_len(self.bjetsBoosted)>=1],
                                      weight = self.DeepCSVMediumSFApplied)
        hasOsElMuLowMllCutBoostedJets = hasOsElMuLowMllCut.refine("hasOsElMuLowMllCutBoostedJets", 
                                      cut=[op.rng_len(self.bjetsBoosted)>=1],
                                      weight = self.DeepCSVMediumSFApplied)

        hasOsElElLowMllCutOutZExclusiveBoostedJets = hasOsElElLowMllCutOutZ.refine("hasOsElElLowMllCutOutZExclusiveBoostedJets", 
                                      cut=[op.rng_len(self.bjetsBoosted)>=1,op.OR(op.rng_len(self.jets)<=1,op.rng_len(self.bjetsResolved)==0)],
                                      weight = self.DeepCSVMediumSFApplied)
        hasOsMuMuLowMllCutOutZExclusiveBoostedJets = hasOsMuMuLowMllCutOutZ.refine("hasOsMuMuLowMllCutOutZExclusiveBoostedJets", 
                                      cut=[op.rng_len(self.bjetsBoosted)>=1,op.OR(op.rng_len(self.jets)<=1,op.rng_len(self.bjetsResolved)==0)],
                                      weight = self.DeepCSVMediumSFApplied)
        hasOsElMuLowMllCutExclusiveBoostedJets = hasOsElMuLowMllCut.refine("hasOsElMuLowMllCutExclusiveBoostedJets", 
                                      cut=[op.rng_len(self.bjetsBoosted)>=1,op.OR(op.rng_len(self.jets)<=1,op.rng_len(self.bjetsResolved)==0)],
                                      weight = self.DeepCSVMediumSFApplied)

        # Yield plots #
        plots.append(makeYieldPlot(self,hasOsElElLowMllCutOutZBoostedJets,"OSElElMllCutOutZBoosted","OS leptons + $M_{ll}$ + Boosted (channel : $e^+e^-$)",12))
        plots.append(makeYieldPlot(self,hasOsMuMuLowMllCutOutZBoostedJets,"OSMuMuMllCutOutZBoosted","OS leptons + $M_{ll}$ + Boosted (channel : $\mu^+\mu^-$)",13))
        plots.append(makeYieldPlot(self,hasOsElMuLowMllCutBoostedJets,"OSMuMuMllCutBoosted","OS leptons + $M_{ll}$ + Boosted (channel : $e^{\pm}\mu^{\mp}$)",14))

        plots.append(makeYieldPlot(self,hasOsElElLowMllCutOutZExclusiveBoostedJets,"OSElElMllCutOutZExclusiveBoosted","OS leptons + $M_{ll}$ + Exclusive Boosted (channel : $e^+e^-$)",15))
        plots.append(makeYieldPlot(self,hasOsMuMuLowMllCutOutZExclusiveBoostedJets,"OSMuMuMllCutOutZExclusiveBoosted","OS leptons + $M_{ll}$ + Exclusive Boosted (channel : $\mu^+\mu^-$)",16))
        plots.append(makeYieldPlot(self,hasOsElMuLowMllCutExclusiveBoostedJets,"OSMuMuMllCutExclusiveBoosted","OS leptons + $M_{ll}$ + Exclusive Boosted (channel : $e^{\pm}\mu^{\mp}$)",17))

        # Boosted + OS dilepton plots #
        hasOsMllCutBoostedChannelList = [
                             {'channel':'ElEl','sel':hasOsElElLowMllCutOutZBoostedJets,'dilepton':self.OsElEl[0],'suffix':'hasOsElElLowMllCutOutZBoostedJets'},
                             {'channel':'MuMu','sel':hasOsMuMuLowMllCutOutZBoostedJets,'dilepton':self.OsMuMu[0],'suffix':'hasOsMuMuLowMllCutOutZBoostedJets'},
                             {'channel':'ElMu','sel':hasOsElMuLowMllCutBoostedJets,    'dilepton':self.OsElMu[0],'suffix':'hasOsElMuLowMllCutBoostedJets'},

        #                     {'channel':'ElEl','sel':hasOsElElLowMllCutOutZExclusiveBoostedJets,'dilepton':self.OsElEl[0],'suffix':'hasOsElElLowMllCutOutZExclusiveBoostedJets'},
        #                     {'channel':'MuMu','sel':hasOsMuMuLowMllCutOutZExclusiveBoostedJets,'dilepton':self.OsMuMu[0],'suffix':'hasOsMuMuLowMllCutOutZExclusiveBoostedJets'},
        #                     {'channel':'ElMu','sel':hasOsElMuLowMllCutExclusiveBoostedJets,    'dilepton':self.OsElMu[0],'suffix':'hasOsElMuLowMllCutExclusiveBoostedJets'},
                           ]

        for channelDict in hasOsMllCutBoostedChannelList:
            # Dilepton plots #
            plots.extend(makeDileptonPlots(self, **channelDict))
            # MET plots #
            plots.extend(makeMETPlots(self      = self,
                                      sel       = channelDict['sel'],
                                      met       = self.corrMET,
                                      suffix    = channelDict['suffix'],
                                      channel   = channelDict['channel']))
            # Fatjet plots #
            plots.extend(makeFatJetPlots(self,
                                         sel       = channelDict['sel'],
                                         fatjets   = self.bjetsBoosted,
                                         suffix    = channelDict['suffix'],
                                         channel   = channelDict['channel']))
  
        ##### RESOLVED #####
        # Combine dilepton and Exclusive Resolved (Exclusive = NOT Boosted) selections #
        hasOsElElLowMllCutOutZResolvedJets = hasOsElElLowMllCutOutZ.refine("hasOsElElLowMllCutOutZResolvedJets", 
                                      cut=[op.rng_len(self.jets)>=2, op.rng_len(self.bjetsResolved)>=1],
                                      weight = self.DeepJetMediumSFApplied)
        hasOsMuMuLowMllCutOutZResolvedJets = hasOsMuMuLowMllCutOutZ.refine("hasOsMuMuLowMllCutOutZResolvedJets", 
                                      cut=[op.rng_len(self.jets)>=2, op.rng_len(self.bjetsResolved)>=1], 
                                      weight = self.DeepJetMediumSFApplied)
        hasOsElMuLowMllCutResolvedJets = hasOsElMuLowMllCut.refine("hasOsElMuLowMllCutResolvedJets", 
                                      cut=[op.rng_len(self.jets)>=2, op.rng_len(self.bjetsResolved)>=1],
                                      weight = self.DeepJetMediumSFApplied)

        hasOsElElLowMllCutOutZExclusiveResolvedJets = hasOsElElLowMllCutOutZ.refine("hasOsElElLowMllCutOutZExclusiveResolvedJets", 
                                      cut=[op.rng_len(self.jets)>=2, op.rng_len(self.bjetsResolved)>=1,op.rng_len(self.bjetsBoosted)==0],
                                      weight = self.DeepJetMediumSFApplied)
        hasOsMuMuLowMllCutOutZExclusiveResolvedJets = hasOsMuMuLowMllCutOutZ.refine("hasOsMuMuLowMllCutOutZExclusiveResolvedJets", 
                                      cut=[op.rng_len(self.jets)>=2, op.rng_len(self.bjetsResolved)>=1,op.rng_len(self.bjetsBoosted)==0], 
                                      weight = self.DeepJetMediumSFApplied)
        hasOsElMuLowMllCutExclusiveResolvedJets = hasOsElMuLowMllCut.refine("hasOsElMuLowMllCutExclusiveResolvedJets", 
                                      cut=[op.rng_len(self.jets)>=2, op.rng_len(self.bjetsResolved)>=1,op.rng_len(self.bjetsBoosted)==0],
                                      weight = self.DeepJetMediumSFApplied)

        # Yield plots #
        plots.append(makeYieldPlot(self,hasOsElElLowMllCutOutZResolvedJets,"OSElElMllCutOutZResolved","OS leptons + $M_{ll}$ + Resolved (channel : $e^+e^-$)",18))
        plots.append(makeYieldPlot(self,hasOsMuMuLowMllCutOutZResolvedJets,"OSMuMuMllCutOutZResolved","OS leptons + $M_{ll}$ + Resolved (channel : $\mu^+\mu^-$)",19))
        plots.append(makeYieldPlot(self,hasOsElMuLowMllCutResolvedJets,"OSMuMuMllCutResolved","OS leptons + $M_{ll}$ + Resolved (channel : $e^{\pm}\mu^{\mp}$)",20))

        plots.append(makeYieldPlot(self,hasOsElElLowMllCutOutZExclusiveResolvedJets,"OSElElMllCutOutZExclusiveResolved","OS leptons + $M_{ll}$ + Exclusive Resolved (channel : $e^+e^-$)",21))
        plots.append(makeYieldPlot(self,hasOsMuMuLowMllCutOutZExclusiveResolvedJets,"OSMuMuMllCutOutZExclusiveResolved","OS leptons + $M_{ll}$ + Exclusive Resolved (channel : $\mu^+\mu^-$)",22))
        plots.append(makeYieldPlot(self,hasOsElMuLowMllCutExclusiveResolvedJets,"OSMuMuMllCutExclusiveResolved","OS leptons + $M_{ll}$ + Exclusive Resolved (channel : $e^{\pm}\mu^{\mp}$)",23))

        # ExclusiveResolved + OS dilepton plots #
        hasOsMllCutExclusiveResolvedChannelList = [
        #                     {'channel':'ElEl','sel':hasOsElElLowMllCutOutZResolvedJets,'dilepton':self.OsElEl[0],'suffix':'hasOsElElLowMllCutOutZResolvedJets'},
        #                     {'channel':'MuMu','sel':hasOsMuMuLowMllCutOutZResolvedJets,'dilepton':self.OsMuMu[0],'suffix':'hasOsMuMuLowMllCutOutZResolvedJets'},
        #                     {'channel':'ElMu','sel':hasOsElMuLowMllCutResolvedJets,    'dilepton':self.OsElMu[0],'suffix':'hasOsElMuLowMllCutResolvedJets'},

                             {'channel':'ElEl','sel':hasOsElElLowMllCutOutZExclusiveResolvedJets,'dilepton':self.OsElEl[0],'suffix':'hasOsElElLowMllCutOutZExclusiveResolvedJets'},
                             {'channel':'MuMu','sel':hasOsMuMuLowMllCutOutZExclusiveResolvedJets,'dilepton':self.OsMuMu[0],'suffix':'hasOsMuMuLowMllCutOutZExclusiveResolvedJets'},
                             {'channel':'ElMu','sel':hasOsElMuLowMllCutExclusiveResolvedJets,    'dilepton':self.OsElMu[0],'suffix':'hasOsElMuLowMllCutExclusiveResolvedJets'},
                           ]

        for channelDict in hasOsMllCutExclusiveResolvedChannelList:
            # Dilepton plots #
            plots.extend(makeDileptonPlots(self, **channelDict))
            # MET plots #
            plots.extend(makeMETPlots(self      = self,
                                      sel       = channelDict['sel'],
                                      met       = self.corrMET,
                                      suffix    = channelDict['suffix'],
                                      channel   = channelDict['channel']))
            # Jets plots #
            plots.extend(makeJetsPlots(self,
                                       sel         = channelDict['sel'],
                                       bjets       = self.bjetsResolved,
                                       lightjets   = self.lightjetsResolved,
                                       alljets     = self.jets,
                                       suffix      = channelDict['suffix'],
                                       channel     = channelDict['channel']))

        #############################################################################
        #                     High-level combinations                               #
        #############################################################################
        highLevelBaseQuantities = { # No need to repeat for each channelDict
             'met'          : self.corrMET,
             'jets'         : self.jets, 
             'resolvedjets' : self.bjetsResolved,
             'lightjets'    : self.lightjetsResolved,
             'boostedjets'  : self.bjetsBoosted,
                                  }
 
        hasOsMllCutHighLevelVariablesChannelList = [
            # Boosted #
            {'channel'      :'ElEl',
             'dilepton'     : self.OsElEl[0],
             'sel'          : hasOsElElLowMllCutOutZBoostedJets,
             'suffix'       : 'hasOsElElLowMllCutOutZBoostedJets'},
            {'channel'      :'MuMu',
             'dilepton'     : self.OsMuMu[0],
             'sel'          : hasOsMuMuLowMllCutOutZBoostedJets,
             'suffix'       : 'hasOsMuMuLowMllCutOutZBoostedJets'},
            {'channel'      :'ElMu',
             'dilepton'     : self.OsElMu[0],
             'sel'          : hasOsElMuLowMllCutBoostedJets,
             'suffix'       : 'hasOsElMuLowMllCutBoostedJets'},
            # Resolved #
            {'channel'      :'ElEl',
             'dilepton'     : self.OsElEl[0],
             'sel'          : hasOsElElLowMllCutOutZExclusiveResolvedJets,
             'suffix'       : 'hasOsElElLowMllCutOutZExclusiveResolvedJets'},
            {'channel'      :'MuMu',
             'dilepton'     : self.OsMuMu[0],
             'sel'          : hasOsMuMuLowMllCutOutZExclusiveResolvedJets,
             'suffix'       : 'hasOsMuMuLowMllCutOutZExclusiveResolvedJets'},
            {'channel'      :'ElMu',
             'dilepton'     : self.OsElMu[0],
             'sel'          : hasOsElMuLowMllCutExclusiveResolvedJets,
             'suffix'       : 'hasOsElMuLowMllCutExclusiveResolvedJets'},
                                                   ]

        for channelDict in hasOsMllCutHighLevelVariablesChannelList:
            plots.extend(makeHighLevelQuantities(self,**highLevelBaseQuantities,**channelDict))


        ## helper selection (OR) to make sure jet calculations are only done once
        #hasOSLL = noSel.refine("hasOSLL", cut=op.OR(*( hasOSLL_cmbRng(rng) for rng in osLLRng.values())))
        #forceDefine(t._Jet.calcProd, hasOSLL)
        #for varNm in t._Jet.available:
        #    forceDefine(t._Jet[varNm], hasOSLL)
        return plots

