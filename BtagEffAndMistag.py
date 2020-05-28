import os
import sys
import json

from bamboo.analysismodules import HistogramsModule
from bamboo import treefunctions as op
from bamboo.plots import Plot, EquidistantBinning, VariableBinning

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

    def definePlots(self, t, noSel, sample=None, sampleCfg=None): 
        noSel = super(BtagEffAndMistagNano,self).prepareObjects(t, noSel, sample, sampleCfg)

        era = sampleCfg['era']
        plots = []

        # protection against data #
        if not self.is_MC:
            return []

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



        return plots

