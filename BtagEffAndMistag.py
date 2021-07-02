import os
import sys
import json

from bamboo.analysismodules import HistogramsModule
from bamboo import treefunctions as op
from bamboo.plots import Plot, EquidistantBinning, VariableBinning, SummedPlot
from bamboo.analysisutils import forceDefine

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)))) # Add scripts in this directory
from BaseHHtobbWW import BaseNanoHHtobbWW

from scalefactorsbbWW import ScaleFactorsbbWW


#===============================================================================================#
#                                       PlotterHHtobbWW                                         #
#===============================================================================================#
class BtagEffAndMistagNano(BaseNanoHHtobbWW,HistogramsModule):
    """ Plotter module: HH->bbW(->e/µ nu)W(->e/µ nu) histograms from NanoAOD """
    def __init__(self, args):
        super(BtagEffAndMistagNano, self).__init__(args)

    def initialize(self):
        super(BtagEffAndMistagNano, self).initialize(forSkimmer=True)

    def definePlots(self, t, noSel, sample=None, sampleCfg=None): 
        noSel = super(BtagEffAndMistagNano,self).prepareObjects(t, noSel, sample, sampleCfg, channel='DL', forSkimmer=True)

        era = sampleCfg['era']
        plots = []

        forceDefine(self.tree._Jet.calcProd, noSel)  # Jets for configureJets
        forceDefine(self.tree._FatJet.calcProd, noSel)  # FatJets for configureJets
        forceDefine(getattr(self.tree, "_{0}".format("MET" if self.era != "2017" else "METFixEE2017")).calcProd,noSel) # MET for configureMET

        # protection against data #
        if not self.is_MC:
            return []

        #############################################################################
        #                                   AK4                                     #
        #############################################################################
        ak4_truth_lightjets = op.select(self.ak4Jets, lambda j : j.hadronFlavour == 0)
        ak4_truth_cjets     = op.select(self.ak4Jets, lambda j : j.hadronFlavour == 4)
        ak4_truth_bjets     = op.select(self.ak4Jets, lambda j : j.hadronFlavour == 5)

        N_ak4_truth_lightjets = Plot.make3D('N_ak4_truth_lightjets',
                                        [op.map(ak4_truth_lightjets, lambda j : op.abs(j.eta)),
                                         op.map(ak4_truth_lightjets, lambda j : j.pt),
                                         op.map(ak4_truth_lightjets, lambda j : j.btagDeepFlavB)],
                                        noSel,
                                        [EquidistantBinning(100,0.,2.5),
                                         EquidistantBinning(100,0.,1000),
                                         EquidistantBinning(100,0.,1.)],
                                        xTitle = 'lightjet #eta',
                                        yTitle = 'lightjet P_{T}',
                                        zTitle = 'lightjet Btagging score')

        N_ak4_truth_cjets = Plot.make3D('N_ak4_truth_cjets',
                                        [op.map(ak4_truth_cjets, lambda j : op.abs(j.eta)),
                                         op.map(ak4_truth_cjets, lambda j : j.pt),
                                         op.map(ak4_truth_cjets, lambda j : j.btagDeepFlavB)],
                                        noSel,
                                        [EquidistantBinning(100,0.,2.5),
                                         EquidistantBinning(100,0.,1000),
                                         EquidistantBinning(100,0.,1.)],
                                        xTitle = 'cjet #eta',
                                        yTitle = 'cjet P_{T}',
                                        zTitle = 'cjet Btagging score')

        N_ak4_truth_bjets = Plot.make3D('N_ak4_truth_bjets',
                                        [op.map(ak4_truth_bjets, lambda j : op.abs(j.eta)),
                                         op.map(ak4_truth_bjets, lambda j : j.pt),
                                         op.map(ak4_truth_bjets, lambda j : j.btagDeepFlavB)],
                                        noSel,
                                        [EquidistantBinning(100,0.,2.5),
                                         EquidistantBinning(100,0.,1000),
                                         EquidistantBinning(100,0.,1.)],
                                        xTitle = 'bjet #eta',
                                        yTitle = 'bjet P_{T}',
                                        zTitle = 'bjet Btagging score')

        plots.extend([N_ak4_truth_lightjets,N_ak4_truth_cjets,N_ak4_truth_bjets])

        ak4_btagged_lightjets = op.select(ak4_truth_lightjets, self.lambda_ak4Btag)
        ak4_btagged_cjets     = op.select(ak4_truth_cjets,     self.lambda_ak4Btag)
        ak4_btagged_bjets     = op.select(ak4_truth_bjets,     self.lambda_ak4Btag)

        N_ak4_btagged_lightjets = Plot.make3D('N_ak4_btagged_lightjets',
                                        [op.map(ak4_btagged_lightjets, lambda j : op.abs(j.eta)),
                                         op.map(ak4_btagged_lightjets, lambda j : j.pt),
                                         op.map(ak4_btagged_lightjets, lambda j : j.btagDeepFlavB)],
                                        noSel,
                                        [EquidistantBinning(100,0.,2.5),
                                         EquidistantBinning(100,0.,1000),
                                         EquidistantBinning(100,0.,1.)],
                                        xTitle = 'lightjet #eta',
                                        yTitle = 'lightjet P_{T}',
                                        zTitle = 'lightjet Btagging score')

        N_ak4_btagged_cjets = Plot.make3D('N_ak4_btagged_cjets',
                                        [op.map(ak4_btagged_cjets, lambda j : op.abs(j.eta)),
                                         op.map(ak4_btagged_cjets, lambda j : j.pt),
                                         op.map(ak4_btagged_cjets, lambda j : j.btagDeepFlavB)],
                                        noSel,
                                        [EquidistantBinning(100,0.,2.5),
                                         EquidistantBinning(100,0.,1000),
                                         EquidistantBinning(100,0.,1.)],
                                        xTitle = 'cjet #eta',
                                        yTitle = 'cjet P_{T}',
                                        zTitle = 'cjet Btagging score')

        N_ak4_btagged_bjets = Plot.make3D('N_ak4_btagged_bjets',
                                        [op.map(ak4_btagged_bjets, lambda j : op.abs(j.eta)),
                                         op.map(ak4_btagged_bjets, lambda j : j.pt),
                                         op.map(ak4_btagged_bjets, lambda j : j.btagDeepFlavB)],
                                        noSel,
                                        [EquidistantBinning(100,0.,2.5),
                                         EquidistantBinning(100,0.,1000),
                                         EquidistantBinning(100,0.,1.)],
                                        xTitle = 'bjet #eta',
                                        yTitle = 'bjet P_{T}',
                                        zTitle = 'bjet Btagging score')

        plots.extend([N_ak4_btagged_lightjets,N_ak4_btagged_cjets,N_ak4_btagged_bjets])


        #############################################################################
        #                                   AK8                                     #
        #############################################################################

        # Truth MC object hadron flavour #
        # subjet.nBHadrons>0 : bjet
        # subjet.nCHadrons>0 : cjet
        # else               : lightjet
        ak8subJet1_truth_lightjets = op.select(self.ak8Jets, lambda j : op.AND(j.subJet1.nBHadrons==0,j.subJet1.nCHadrons==0))
        ak8subJet2_truth_lightjets = op.select(self.ak8Jets, lambda j : op.AND(j.subJet2.nBHadrons==0,j.subJet2.nCHadrons==0))
        ak8subJet1_truth_cjets     = op.select(self.ak8Jets, lambda j : j.subJet1.nCHadrons>0)
        ak8subJet2_truth_cjets     = op.select(self.ak8Jets, lambda j : j.subJet2.nCHadrons>0)
        ak8subJet1_truth_bjets     = op.select(self.ak8Jets, lambda j : j.subJet1.nBHadrons>0)
        ak8subJet2_truth_bjets     = op.select(self.ak8Jets, lambda j : j.subJet2.nBHadrons>0)

        N_ak8_truth_subJet1_lightjets = Plot.make3D('N_ak8_truth_subJet1_lightjets',
                                                [op.map(ak8subJet1_truth_lightjets, lambda j : op.abs(j.subJet1.eta)),
                                                 op.map(ak8subJet1_truth_lightjets, lambda j : j.subJet1.pt),
                                                 op.map(ak8subJet1_truth_lightjets, lambda j : j.subJet1.btagDeepB)],
                                                noSel,
                                                [EquidistantBinning(100,0.,2.5),
                                                 EquidistantBinning(100,0.,1000),
                                                 EquidistantBinning(100,0.,1.)],
                                                xTitle = 'lightjet subJet1 #eta',
                                                yTitle = 'lightjet subJet1 P_{T}',
                                                zTitle = 'lightjet subJet1 Btagging score')
        N_ak8_truth_subJet2_lightjets = Plot.make3D('N_ak8_truth_subJet2_lightjets',
                                                [op.map(ak8subJet2_truth_lightjets, lambda j : op.abs(j.subJet2.eta)),
                                                 op.map(ak8subJet2_truth_lightjets, lambda j : j.subJet2.pt),
                                                 op.map(ak8subJet2_truth_lightjets, lambda j : j.subJet2.btagDeepB)],
                                                noSel,
                                                [EquidistantBinning(100,0.,2.5),
                                                 EquidistantBinning(100,0.,1000),
                                                 EquidistantBinning(100,0.,1.)],
                                                xTitle = 'lightjet subJet2 #eta',
                                                yTitle = 'lightjet subJet2 P_{T}',
                                                zTitle = 'lightjet subJet2 Btagging score')
        N_ak8_truth_subJet1_cjets = Plot.make3D('N_ak8_truth_subJet1_cjets',
                                                [op.map(ak8subJet1_truth_cjets, lambda j : op.abs(j.subJet1.eta)),
                                                 op.map(ak8subJet1_truth_cjets, lambda j : j.subJet1.pt),
                                                 op.map(ak8subJet1_truth_cjets, lambda j : j.subJet1.btagDeepB)],
                                                noSel,
                                                [EquidistantBinning(100,0.,2.5),
                                                 EquidistantBinning(100,0.,1000),
                                                 EquidistantBinning(100,0.,1.)],
                                                xTitle = 'cjet subJet1 #eta',
                                                yTitle = 'cjet subJet1 P_{T}',
                                                zTitle = 'cjet subJet1 Btagging score')
        N_ak8_truth_subJet2_cjets = Plot.make3D('N_ak8_truth_subJet2_cjets',
                                                [op.map(ak8subJet2_truth_cjets, lambda j : op.abs(j.subJet2.eta)),
                                                 op.map(ak8subJet2_truth_cjets, lambda j : j.subJet2.pt),
                                                 op.map(ak8subJet2_truth_cjets, lambda j : j.subJet2.btagDeepB)],
                                                noSel,
                                                [EquidistantBinning(100,0.,2.5),
                                                 EquidistantBinning(100,0.,1000),
                                                 EquidistantBinning(100,0.,1.)],
                                                xTitle = 'cjet subJet2 #eta',
                                                yTitle = 'cjet subJet2 P_{T}',
                                                zTitle = 'cjet subJet2 Btagging score')
        N_ak8_truth_subJet1_bjets = Plot.make3D('N_ak8_truth_subJet1_bjets',
                                                [op.map(ak8subJet1_truth_bjets, lambda j : op.abs(j.subJet1.eta)),
                                                 op.map(ak8subJet1_truth_bjets, lambda j : j.subJet1.pt),
                                                 op.map(ak8subJet1_truth_bjets, lambda j : j.subJet1.btagDeepB)],
                                                noSel,
                                                [EquidistantBinning(100,0.,2.5),
                                                 EquidistantBinning(100,0.,1000),
                                                 EquidistantBinning(100,0.,1.)],
                                                xTitle = 'bjet subJet1 #eta',
                                                yTitle = 'bjet subJet1 P_{T}',
                                                zTitle = 'bjet subJet1 Btagging score')
        N_ak8_truth_subJet2_bjets = Plot.make3D('N_ak8_truth_subJet2_bjets',
                                                [op.map(ak8subJet2_truth_bjets, lambda j : op.abs(j.subJet2.eta)),
                                                 op.map(ak8subJet2_truth_bjets, lambda j : j.subJet2.pt),
                                                 op.map(ak8subJet2_truth_bjets, lambda j : j.subJet2.btagDeepB)],
                                                noSel,
                                                [EquidistantBinning(100,0.,2.5),
                                                 EquidistantBinning(100,0.,1000),
                                                 EquidistantBinning(100,0.,1.)],
                                                xTitle = 'bjet subJet2 #eta',
                                                yTitle = 'bjet subJet2 P_{T}',
                                                zTitle = 'bjet subJet2 Btagging score')

        plots.extend([N_ak8_truth_subJet1_lightjets,N_ak8_truth_subJet2_lightjets,N_ak8_truth_subJet1_cjets,N_ak8_truth_subJet2_cjets,N_ak8_truth_subJet1_bjets,N_ak8_truth_subJet2_bjets])
        plots.append(SummedPlot('N_ak8_truth_lightjets',
                                [N_ak8_truth_subJet1_lightjets,N_ak8_truth_subJet2_lightjets],
                                xTitle = 'lightjet #eta',
                                yTitle = 'lightjet P_{T}',
                                zTitle = 'lightjet Btagging score'))
        plots.append(SummedPlot('N_ak8_truth_cjets',
                                [N_ak8_truth_subJet1_cjets,N_ak8_truth_subJet2_cjets],
                                xTitle = 'cjet #eta',
                                yTitle = 'cjet P_{T}',
                                zTitle = 'cjet Btagging score'))
        plots.append(SummedPlot('N_ak8_truth_bjets',
                                [N_ak8_truth_subJet1_bjets,N_ak8_truth_subJet2_bjets],
                                xTitle = 'bjet #eta',
                                yTitle = 'bjet P_{T}',
                                zTitle = 'bjet Btagging score'))

                                    
        # Btagged objects per flavour #
        ak8subJet1_btagged_lightjets = op.select(ak8subJet1_truth_lightjets, lambda j: self.lambda_subjetBtag(j.subJet1))
        ak8subJet2_btagged_lightjets = op.select(ak8subJet2_truth_lightjets, lambda j: self.lambda_subjetBtag(j.subJet2))
        ak8subJet1_btagged_cjets     = op.select(ak8subJet1_truth_cjets, lambda j: self.lambda_subjetBtag(j.subJet1))
        ak8subJet2_btagged_cjets     = op.select(ak8subJet2_truth_cjets, lambda j: self.lambda_subjetBtag(j.subJet2))
        ak8subJet1_btagged_bjets     = op.select(ak8subJet1_truth_bjets, lambda j: self.lambda_subjetBtag(j.subJet1))
        ak8subJet2_btagged_bjets     = op.select(ak8subJet2_truth_bjets, lambda j: self.lambda_subjetBtag(j.subJet2))

        N_ak8_btagged_subJet1_lightjets = Plot.make3D('N_ak8_btagged_subJet1_lightjets',
                                                [op.map(ak8subJet1_btagged_lightjets, lambda j : op.abs(j.subJet1.eta)),
                                                 op.map(ak8subJet1_btagged_lightjets, lambda j : j.subJet1.pt),
                                                 op.map(ak8subJet1_btagged_lightjets, lambda j : j.subJet1.btagDeepB)],
                                                noSel,
                                                [EquidistantBinning(100,0.,2.5),
                                                 EquidistantBinning(100,0.,1000),
                                                 EquidistantBinning(100,0.,1.)],
                                                xTitle = 'lightjet subJet1 #eta',
                                                yTitle = 'lightjet subJet1 P_{T}',
                                                zTitle = 'lightjet subJet1 Btagging score')
        N_ak8_btagged_subJet2_lightjets = Plot.make3D('N_ak8_btagged_subJet2_lightjets',
                                                [op.map(ak8subJet2_btagged_lightjets, lambda j : op.abs(j.subJet2.eta)),
                                                 op.map(ak8subJet2_btagged_lightjets, lambda j : j.subJet2.pt),
                                                 op.map(ak8subJet2_btagged_lightjets, lambda j : j.subJet2.btagDeepB)],
                                                noSel,
                                                [EquidistantBinning(100,0.,2.5),
                                                 EquidistantBinning(100,0.,1000),
                                                 EquidistantBinning(100,0.,1.)],
                                                xTitle = 'lightjet subJet2 #eta',
                                                yTitle = 'lightjet subJet2 P_{T}',
                                                zTitle = 'lightjet subJet2 Btagging score')
        N_ak8_btagged_subJet1_cjets = Plot.make3D('N_ak8_btagged_subJet1_cjets',
                                                [op.map(ak8subJet1_btagged_cjets, lambda j : op.abs(j.subJet1.eta)),
                                                 op.map(ak8subJet1_btagged_cjets, lambda j : j.subJet1.pt),
                                                 op.map(ak8subJet1_btagged_cjets, lambda j : j.subJet1.btagDeepB)],
                                                noSel,
                                                [EquidistantBinning(100,0.,2.5),
                                                 EquidistantBinning(100,0.,1000),
                                                 EquidistantBinning(100,0.,1.)],
                                                xTitle = 'cjet subJet1 #eta',
                                                yTitle = 'cjet subJet1 P_{T}',
                                                zTitle = 'cjet subJet1 Btagging score')
        N_ak8_btagged_subJet2_cjets = Plot.make3D('N_ak8_btagged_subJet2_cjets',
                                                [op.map(ak8subJet2_btagged_cjets, lambda j : op.abs(j.subJet2.eta)),
                                                 op.map(ak8subJet2_btagged_cjets, lambda j : j.subJet2.pt),
                                                 op.map(ak8subJet2_btagged_cjets, lambda j : j.subJet2.btagDeepB)],
                                                noSel,
                                                [EquidistantBinning(100,0.,2.5),
                                                 EquidistantBinning(100,0.,1000),
                                                 EquidistantBinning(100,0.,1.)],
                                                xTitle = 'cjet subJet2 #eta',
                                                yTitle = 'cjet subJet2 P_{T}',
                                                zTitle = 'cjet subJet2 Btagging score')
        N_ak8_btagged_subJet1_bjets = Plot.make3D('N_ak8_btagged_subJet1_bjets',
                                                [op.map(ak8subJet1_btagged_bjets, lambda j : op.abs(j.subJet1.eta)),
                                                 op.map(ak8subJet1_btagged_bjets, lambda j : j.subJet1.pt),
                                                 op.map(ak8subJet1_btagged_bjets, lambda j : j.subJet1.btagDeepB)],
                                                noSel,
                                                [EquidistantBinning(100,0.,2.5),
                                                 EquidistantBinning(100,0.,1000),
                                                 EquidistantBinning(100,0.,1.)],
                                                xTitle = 'bjet subJet1 #eta',
                                                yTitle = 'bjet subJet1 P_{T}',
                                                zTitle = 'bjet subJet1 Btagging score')
        N_ak8_btagged_subJet2_bjets = Plot.make3D('N_ak8_btagged_subJet2_bjets',
                                                [op.map(ak8subJet2_btagged_bjets, lambda j : op.abs(j.subJet2.eta)),
                                                 op.map(ak8subJet2_btagged_bjets, lambda j : j.subJet2.pt),
                                                 op.map(ak8subJet2_btagged_bjets, lambda j : j.subJet2.btagDeepB)],
                                                noSel,
                                                [EquidistantBinning(100,0.,2.5),
                                                 EquidistantBinning(100,0.,1000),
                                                 EquidistantBinning(100,0.,1.)],
                                                xTitle = 'bjet subJet2 #eta',
                                                yTitle = 'bjet subJet2 P_{T}',
                                                zTitle = 'bjet subJet2 Btagging score')

        plots.extend([N_ak8_btagged_subJet1_lightjets,N_ak8_btagged_subJet2_lightjets,N_ak8_btagged_subJet1_cjets,N_ak8_btagged_subJet2_cjets,N_ak8_btagged_subJet1_bjets,N_ak8_btagged_subJet2_bjets])
        plots.append(SummedPlot('N_ak8_btagged_lightjets',
                                [N_ak8_btagged_subJet1_lightjets,N_ak8_btagged_subJet2_lightjets],
                                xTitle = 'lightjet #eta',
                                yTitle = 'lightjet P_{T}',
                                zTitle = 'lightjet Btagging score'))
        plots.append(SummedPlot('N_ak8_btagged_cjets',
                                [N_ak8_btagged_subJet1_cjets,N_ak8_btagged_subJet2_cjets],
                                xTitle = 'cjet #eta',
                                yTitle = 'cjet P_{T}',
                                zTitle = 'cjet Btagging score'))
        plots.append(SummedPlot('N_ak8_btagged_bjets',
                                [N_ak8_btagged_subJet1_bjets,N_ak8_btagged_subJet2_bjets],
                                xTitle = 'bjet #eta',
                                yTitle = 'bjet P_{T}',
                                zTitle = 'bjet Btagging score'))


        return plots

