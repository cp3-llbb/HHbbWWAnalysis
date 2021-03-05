import os
import sys
import json

from bamboo.analysismodules import HistogramsModule
from bamboo import treefunctions as op
<<<<<<< HEAD
from bamboo.plots import Plot, EquidistantBinning, VariableBinning, SummedPlot
from bamboo.analysisutils import forceDefine
=======
from bamboo.plots import Plot, EquidistantBinning, VariableBinning
>>>>>>> eed57b2aa0f195c36370ef8af9adb2772972dcc4

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

<<<<<<< HEAD
    def initialize(self):
        super(BtagEffAndMistagNano, self).initialize(forSkimmer=True)

    def definePlots(self, t, noSel, sample=None, sampleCfg=None): 
        noSel = super(BtagEffAndMistagNano,self).prepareObjects(t, noSel, sample, sampleCfg, channel='DL', forSkimmer=True)
=======
    def definePlots(self, t, noSel, sample=None, sampleCfg=None): 
        noSel = super(BtagEffAndMistagNano,self).prepareObjects(t, noSel, sample, sampleCfg)
>>>>>>> eed57b2aa0f195c36370ef8af9adb2772972dcc4

        era = sampleCfg['era']
        plots = []

<<<<<<< HEAD
        forceDefine(self.tree._Jet.calcProd, noSel)  # Jets for configureJets
        forceDefine(self.tree._FatJet.calcProd, noSel)  # FatJets for configureJets
        forceDefine(getattr(self.tree, "_{0}".format("MET" if self.era != "2017" else "METFixEE2017")).calcProd,noSel) # MET for configureMET

=======
>>>>>>> eed57b2aa0f195c36370ef8af9adb2772972dcc4
        # protection against data #
        if not self.is_MC:
            return []

<<<<<<< HEAD
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

        N_truth_subJet1_lightjets = Plot.make3D('N_truth_subJet1_lightjets',
                                                [op.map(ak8subJet1_truth_lightjets, lambda j : j.subJet1.eta),
                                                 op.map(ak8subJet1_truth_lightjets, lambda j : j.subJet1.pt),
                                                 op.map(ak8subJet1_truth_lightjets, lambda j : j.subJet1.btagDeepB)],
                                                noSel,
                                                [EquidistantBinning(100,0.,2.5),
                                                 EquidistantBinning(100,0.,1000),
                                                 EquidistantBinning(100,0.,1.)],
                                                xTitle = 'lightjet subJet1 #eta',
                                                yTitle = 'lightjet subJet1 P_{T}',
                                                zTitle = 'lightjet subJet1 Btagging score')
        N_truth_subJet2_lightjets = Plot.make3D('N_truth_subJet2_lightjets',
                                                [op.map(ak8subJet2_truth_lightjets, lambda j : j.subJet2.eta),
                                                 op.map(ak8subJet2_truth_lightjets, lambda j : j.subJet2.pt),
                                                 op.map(ak8subJet2_truth_lightjets, lambda j : j.subJet2.btagDeepB)],
                                                noSel,
                                                [EquidistantBinning(100,0.,2.5),
                                                 EquidistantBinning(100,0.,1000),
                                                 EquidistantBinning(100,0.,1.)],
                                                xTitle = 'lightjet subJet2 #eta',
                                                yTitle = 'lightjet subJet2 P_{T}',
                                                zTitle = 'lightjet subJet2 Btagging score')
        N_truth_subJet1_cjets = Plot.make3D('N_truth_subJet1_cjets',
                                                [op.map(ak8subJet1_truth_cjets, lambda j : j.subJet1.eta),
                                                 op.map(ak8subJet1_truth_cjets, lambda j : j.subJet1.pt),
                                                 op.map(ak8subJet1_truth_cjets, lambda j : j.subJet1.btagDeepB)],
                                                noSel,
                                                [EquidistantBinning(100,0.,2.5),
                                                 EquidistantBinning(100,0.,1000),
                                                 EquidistantBinning(100,0.,1.)],
                                                xTitle = 'cjet subJet1 #eta',
                                                yTitle = 'cjet subJet1 P_{T}',
                                                zTitle = 'cjet subJet1 Btagging score')
        N_truth_subJet2_cjets = Plot.make3D('N_truth_subJet2_cjets',
                                                [op.map(ak8subJet2_truth_cjets, lambda j : j.subJet2.eta),
                                                 op.map(ak8subJet2_truth_cjets, lambda j : j.subJet2.pt),
                                                 op.map(ak8subJet2_truth_cjets, lambda j : j.subJet2.btagDeepB)],
                                                noSel,
                                                [EquidistantBinning(100,0.,2.5),
                                                 EquidistantBinning(100,0.,1000),
                                                 EquidistantBinning(100,0.,1.)],
                                                xTitle = 'cjet subJet2 #eta',
                                                yTitle = 'cjet subJet2 P_{T}',
                                                zTitle = 'cjet subJet2 Btagging score')
        N_truth_subJet1_bjets = Plot.make3D('N_truth_subJet1_bjets',
                                                [op.map(ak8subJet1_truth_bjets, lambda j : j.subJet1.eta),
                                                 op.map(ak8subJet1_truth_bjets, lambda j : j.subJet1.pt),
                                                 op.map(ak8subJet1_truth_bjets, lambda j : j.subJet1.btagDeepB)],
                                                noSel,
                                                [EquidistantBinning(100,0.,2.5),
                                                 EquidistantBinning(100,0.,1000),
                                                 EquidistantBinning(100,0.,1.)],
                                                xTitle = 'bjet subJet1 #eta',
                                                yTitle = 'bjet subJet1 P_{T}',
                                                zTitle = 'bjet subJet1 Btagging score')
        N_truth_subJet2_bjets = Plot.make3D('N_truth_subJet2_bjets',
                                                [op.map(ak8subJet2_truth_bjets, lambda j : j.subJet2.eta),
                                                 op.map(ak8subJet2_truth_bjets, lambda j : j.subJet2.pt),
                                                 op.map(ak8subJet2_truth_bjets, lambda j : j.subJet2.btagDeepB)],
                                                noSel,
                                                [EquidistantBinning(100,0.,2.5),
                                                 EquidistantBinning(100,0.,1000),
                                                 EquidistantBinning(100,0.,1.)],
                                                xTitle = 'bjet subJet2 #eta',
                                                yTitle = 'bjet subJet2 P_{T}',
                                                zTitle = 'bjet subJet2 Btagging score')

        plots.extend([N_truth_subJet1_lightjets,N_truth_subJet2_lightjets,N_truth_subJet1_cjets,N_truth_subJet2_cjets,N_truth_subJet1_bjets,N_truth_subJet2_bjets])
        plots.append(SummedPlot('N_truth_lightjets',
                                [N_truth_subJet1_lightjets,N_truth_subJet2_lightjets],
                                xTitle = 'lightjet #eta',
                                yTitle = 'lightjet P_{T}',
                                zTitle = 'lightjet Btagging score'))
        plots.append(SummedPlot('N_truth_cjets',
                                [N_truth_subJet1_cjets,N_truth_subJet2_cjets],
                                xTitle = 'cjet #eta',
                                yTitle = 'cjet P_{T}',
                                zTitle = 'cjet Btagging score'))
        plots.append(SummedPlot('N_truth_bjets',
                                [N_truth_subJet1_bjets,N_truth_subJet2_bjets],
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

        N_btagged_subJet1_lightjets = Plot.make3D('N_btagged_subJet1_lightjets',
                                                [op.map(ak8subJet1_btagged_lightjets, lambda j : j.subJet1.eta),
                                                 op.map(ak8subJet1_btagged_lightjets, lambda j : j.subJet1.pt),
                                                 op.map(ak8subJet1_btagged_lightjets, lambda j : j.subJet1.btagDeepB)],
                                                noSel,
                                                [EquidistantBinning(100,0.,2.5),
                                                 EquidistantBinning(100,0.,1000),
                                                 EquidistantBinning(100,0.,1.)],
                                                xTitle = 'lightjet subJet1 #eta',
                                                yTitle = 'lightjet subJet1 P_{T}',
                                                zTitle = 'lightjet subJet1 Btagging score')
        N_btagged_subJet2_lightjets = Plot.make3D('N_btagged_subJet2_lightjets',
                                                [op.map(ak8subJet2_btagged_lightjets, lambda j : j.subJet2.eta),
                                                 op.map(ak8subJet2_btagged_lightjets, lambda j : j.subJet2.pt),
                                                 op.map(ak8subJet2_btagged_lightjets, lambda j : j.subJet2.btagDeepB)],
                                                noSel,
                                                [EquidistantBinning(100,0.,2.5),
                                                 EquidistantBinning(100,0.,1000),
                                                 EquidistantBinning(100,0.,1.)],
                                                xTitle = 'lightjet subJet2 #eta',
                                                yTitle = 'lightjet subJet2 P_{T}',
                                                zTitle = 'lightjet subJet2 Btagging score')
        N_btagged_subJet1_cjets = Plot.make3D('N_btagged_subJet1_cjets',
                                                [op.map(ak8subJet1_btagged_cjets, lambda j : j.subJet1.eta),
                                                 op.map(ak8subJet1_btagged_cjets, lambda j : j.subJet1.pt),
                                                 op.map(ak8subJet1_btagged_cjets, lambda j : j.subJet1.btagDeepB)],
                                                noSel,
                                                [EquidistantBinning(100,0.,2.5),
                                                 EquidistantBinning(100,0.,1000),
                                                 EquidistantBinning(100,0.,1.)],
                                                xTitle = 'cjet subJet1 #eta',
                                                yTitle = 'cjet subJet1 P_{T}',
                                                zTitle = 'cjet subJet1 Btagging score')
        N_btagged_subJet2_cjets = Plot.make3D('N_btagged_subJet2_cjets',
                                                [op.map(ak8subJet2_btagged_cjets, lambda j : j.subJet2.eta),
                                                 op.map(ak8subJet2_btagged_cjets, lambda j : j.subJet2.pt),
                                                 op.map(ak8subJet2_btagged_cjets, lambda j : j.subJet2.btagDeepB)],
                                                noSel,
                                                [EquidistantBinning(100,0.,2.5),
                                                 EquidistantBinning(100,0.,1000),
                                                 EquidistantBinning(100,0.,1.)],
                                                xTitle = 'cjet subJet2 #eta',
                                                yTitle = 'cjet subJet2 P_{T}',
                                                zTitle = 'cjet subJet2 Btagging score')
        N_btagged_subJet1_bjets = Plot.make3D('N_btagged_subJet1_bjets',
                                                [op.map(ak8subJet1_btagged_bjets, lambda j : j.subJet1.eta),
                                                 op.map(ak8subJet1_btagged_bjets, lambda j : j.subJet1.pt),
                                                 op.map(ak8subJet1_btagged_bjets, lambda j : j.subJet1.btagDeepB)],
                                                noSel,
                                                [EquidistantBinning(100,0.,2.5),
                                                 EquidistantBinning(100,0.,1000),
                                                 EquidistantBinning(100,0.,1.)],
                                                xTitle = 'bjet subJet1 #eta',
                                                yTitle = 'bjet subJet1 P_{T}',
                                                zTitle = 'bjet subJet1 Btagging score')
        N_btagged_subJet2_bjets = Plot.make3D('N_btagged_subJet2_bjets',
                                                [op.map(ak8subJet2_btagged_bjets, lambda j : j.subJet2.eta),
                                                 op.map(ak8subJet2_btagged_bjets, lambda j : j.subJet2.pt),
                                                 op.map(ak8subJet2_btagged_bjets, lambda j : j.subJet2.btagDeepB)],
                                                noSel,
                                                [EquidistantBinning(100,0.,2.5),
                                                 EquidistantBinning(100,0.,1000),
                                                 EquidistantBinning(100,0.,1.)],
                                                xTitle = 'bjet subJet2 #eta',
                                                yTitle = 'bjet subJet2 P_{T}',
                                                zTitle = 'bjet subJet2 Btagging score')

        plots.extend([N_btagged_subJet1_lightjets,N_btagged_subJet2_lightjets,N_btagged_subJet1_cjets,N_btagged_subJet2_cjets,N_btagged_subJet1_bjets,N_btagged_subJet2_bjets])
        plots.append(SummedPlot('N_btagged_lightjets',
                                [N_btagged_subJet1_lightjets,N_btagged_subJet2_lightjets],
                                xTitle = 'lightjet #eta',
                                yTitle = 'lightjet P_{T}',
                                zTitle = 'lightjet Btagging score'))
        plots.append(SummedPlot('N_btagged_cjets',
                                [N_btagged_subJet1_cjets,N_btagged_subJet2_cjets],
                                xTitle = 'cjet #eta',
                                yTitle = 'cjet P_{T}',
                                zTitle = 'cjet Btagging score'))
        plots.append(SummedPlot('N_btagged_bjets',
                                [N_btagged_subJet1_bjets,N_btagged_subJet2_bjets],
                                xTitle = 'bjet #eta',
                                yTitle = 'bjet P_{T}',
                                zTitle = 'bjet Btagging score'))
=======
        # Get ScaleFactor binning #
        instance = ScaleFactorsbbWW()
        all_SF = instance.all_scalefactors

        tuple_medium_DeepJet_json = all_SF['btag_'+era]['DeepJet_medium'] 
        variables = []
        binning = []
        index = ['x','y','z']
        for json_file in tuple_medium_DeepJet_json:
            with open(json_file,'r') as handle:
                content = json.load(handle)
            var  = content['variables']
            bins = content['binning']
            variables.append(var)
            binning.append(bins)

            print ('Variables in file %s'%json_file)
            for var,bin_name in zip(var,index):
                print (('... Variable : %s'%var).ljust(30,' ')+'Binning (index %s) :'%bin_name,bins[bin_name])

        # Find best binning (finest)  #
        opt_binning = {} 
        for var, bins in zip(variables,binning):
            for v, bn in zip(var,index):
                if v == 'AbsEta':
                    continue
                b = bins[bn]
                if not v in opt_binning.keys():
                    opt_binning[v] = b
                else:
                    opt_binning[v].extend([a for a in b if a not in opt_binning[v]])

        # Sort #
        for var,bins in opt_binning.items():
            opt_binning[var] = sorted(opt_binning[var])

        print ("Optimal Binning chosen")
        for var,bins in opt_binning.items():
            print (('... Variable : %s'%var).ljust(30,' ')+'Binning :',bins)


        # Truth MC object hadron flavour #
        ak4_truth_lightjets = op.select(self.ak4Jets, lambda j : op.abs(j.hadronFlavour) == 0)
        ak4_truth_cjets = op.select(self.ak4Jets, lambda j : op.abs(j.hadronFlavour) == 4)
        ak4_truth_bjets = op.select(self.ak4Jets, lambda j : op.abs(j.hadronFlavour) == 5)

        plots.append(Plot.make3D('N_total_lightjets',
                                 [op.map(ak4_truth_lightjets, lambda j : j.eta),
                                  op.map(ak4_truth_lightjets, lambda j : j.pt),
                                  op.map(ak4_truth_lightjets, lambda j : j.btagDeepFlavB)],
                                 noSel,
                                 [VariableBinning(opt_binning['Eta']),
                                  VariableBinning(opt_binning['Pt']),
                                  VariableBinning(opt_binning['BTagDiscri'])],
                                 xTitle = 'lightjet #eta',
                                 yTitle = 'lightjet P_{T}',
                                 zTitle = 'lightjet Btagging score'))
        plots.append(Plot.make3D('N_total_cjets',
                                 [op.map(ak4_truth_cjets, lambda j : j.eta),
                                  op.map(ak4_truth_cjets, lambda j : j.pt),
                                  op.map(ak4_truth_cjets, lambda j : j.btagDeepFlavB)],
                                 noSel,
                                 [VariableBinning(opt_binning['Eta']),
                                  VariableBinning(opt_binning['Pt']),
                                  VariableBinning(opt_binning['BTagDiscri'])],
                                 xTitle = 'cjet #eta',
                                 yTitle = 'cjet P_{T}',
                                 zTitle = 'cjet Btagging score'))
        plots.append(Plot.make3D('N_total_bjets',
                                 [op.map(ak4_truth_bjets, lambda j : j.eta),
                                  op.map(ak4_truth_bjets, lambda j : j.pt),
                                  op.map(ak4_truth_bjets, lambda j : j.btagDeepFlavB)],
                                 noSel,
                                 [VariableBinning(opt_binning['Eta']),
                                  VariableBinning(opt_binning['Pt']),
                                  VariableBinning(opt_binning['BTagDiscri'])],
                                 xTitle = 'bjet #eta',
                                 yTitle = 'bjet P_{T}',
                                 zTitle = 'bjet Btagging score'))

                                    
        # Btagged objects per flavour #
        ak4_btagged_lightjets = op.select(ak4_truth_lightjets, self.lambda_ak4Btag)
        ak4_btagged_cjets = op.select(ak4_truth_cjets, self.lambda_ak4Btag) 
        ak4_btagged_bjets = op.select(ak4_truth_bjets, self.lambda_ak4Btag) 

        plots.append(Plot.make3D('N_btagged_lightjets',
                                 [op.map(ak4_btagged_lightjets, lambda j : j.eta),
                                  op.map(ak4_btagged_lightjets, lambda j : j.pt),
                                  op.map(ak4_btagged_lightjets, lambda j : j.btagDeepFlavB)],
                                 noSel,
                                 [VariableBinning(opt_binning['Eta']),
                                  VariableBinning(opt_binning['Pt']),
                                  VariableBinning(opt_binning['BTagDiscri'])],
                                 xTitle = 'lightjet #eta',
                                 yTitle = 'lightjet P_{T}',
                                 zTitle = 'lightjet Btagging score'))
        plots.append(Plot.make3D('N_btagged_cjets',
                                 [op.map(ak4_btagged_cjets, lambda j : j.eta),
                                  op.map(ak4_btagged_cjets, lambda j : j.pt),
                                  op.map(ak4_btagged_cjets, lambda j : j.btagDeepFlavB)],
                                 noSel,
                                 [VariableBinning(opt_binning['Eta']),
                                  VariableBinning(opt_binning['Pt']),
                                  VariableBinning(opt_binning['BTagDiscri'])],
                                 xTitle = 'cjet #eta',
                                 yTitle = 'cjet P_{T}',
                                 zTitle = 'cjet Btagging score'))
        plots.append(Plot.make3D('N_btagged_bjets',
                                 [op.map(ak4_btagged_bjets, lambda j : j.eta),
                                  op.map(ak4_btagged_bjets, lambda j : j.pt),
                                  op.map(ak4_btagged_bjets, lambda j : j.btagDeepFlavB)],
                                 noSel,
                                 [VariableBinning(opt_binning['Eta']),
                                  VariableBinning(opt_binning['Pt']),
                                  VariableBinning(opt_binning['BTagDiscri'])],
                                 xTitle = 'bjet #eta',
                                 yTitle = 'bjet P_{T}',
                                 zTitle = 'bjet Btagging score'))

>>>>>>> eed57b2aa0f195c36370ef8af9adb2772972dcc4


        return plots

