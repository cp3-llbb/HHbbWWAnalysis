import os
import sys
import yaml
import copy
import glob
import json
import ROOT
import argparse
import numpy as np

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)


class ComputeBtagEff:
    def __init__(self,path,suffix,dict_names,era,rebinX=None,rebinY=None,rebinZ=None):
        self.path = os.path.abspath(path)
        self.suffix = suffix
        self.dict_names = dict_names
        self.era = era
        self.rebinX = rebinX
        self.rebinY = rebinY
        self.rebinZ = rebinZ

        self.recoverHistograms()
        self.rebinHistograms()
        self.makeRatio()

    def recoverHistograms(self):
        self.dict_hist = {cat:{typ:None for typ in self.dict_names[cat].keys()} for cat in self.dict_names.keys()} 
        files = sorted(glob.glob(os.path.join(self.path,'results','*root')))
        if len(files) == 0:
            raise RuntimeError('Could not find files in path')
        for f in files:
            if '__skeleton__' in f:
                continue
            r = ROOT.TFile(f)
            print (f'Looking at{f}')
            for cat in self.dict_names.keys():
                # Get hist #
                h_truth = r.Get(self.dict_names[cat]['truth'])
                h_btagged = r.Get(self.dict_names[cat]['btagged'])

                # Generate subdict #
                if cat not in self.dict_hist.keys():
                    self.dict_hist[cat] = {}

                # Total hist #
                if self.dict_hist[cat]['truth'] is None:
                    self.dict_hist[cat]['truth'] = copy.deepcopy(h_truth)
                else:
                    self.dict_hist[cat]['truth'].Add(h_truth)

                # Btagged hist #
                if self.dict_hist[cat]['btagged'] is None:
                    self.dict_hist[cat]['btagged'] = copy.deepcopy(h_btagged)
                else:
                    self.dict_hist[cat]['btagged'].Add(h_btagged)
            r.Close()

    def rebinHistograms(self):
        for cat in self.dict_names.keys():
            print (f"Rebinning category {cat}")
            assert  self.dict_hist[cat]['btagged'].GetNbinsX() == self.dict_hist[cat]['truth'].GetNbinsX()
            assert  self.dict_hist[cat]['btagged'].GetNbinsY() == self.dict_hist[cat]['truth'].GetNbinsY()
            assert  self.dict_hist[cat]['btagged'].GetNbinsZ() == self.dict_hist[cat]['truth'].GetNbinsZ()
            NbinsX = self.dict_hist[cat]['btagged'].GetNbinsX()
            NbinsY = self.dict_hist[cat]['btagged'].GetNbinsY()
            NbinsZ = self.dict_hist[cat]['btagged'].GetNbinsZ()

            if self.rebinX is not None:
                assert NbinsX % self.rebinX == 0
                groupX = NbinsX // self.rebinX
            else:
                groupX = 1

            if self.rebinY is not None:
                assert NbinsY % self.rebinY == 0
                groupY = NbinsY // self.rebinY
            else:
                groupY = 1

            if self.rebinZ is not None:
                assert NbinsZ % self.rebinZ == 0
                groupZ = NbinsZ // self.rebinZ
            else:
                groupZ = 1

            self.dict_hist[cat]['btagged'] = self.dict_hist[cat]['btagged'].Rebin3D(groupX,groupY,groupZ)
            self.dict_hist[cat]['truth'] = self.dict_hist[cat]['truth'].Rebin3D(groupX,groupY,groupZ)

    def makeRatio(self):
        # Produce ratio for each category #
        for cat in self.dict_hist.keys():
            ratio = copy.deepcopy(self.dict_hist[cat]['btagged'])
            ratio.Divide(self.dict_hist[cat]['truth'])
            self.dict_hist[cat]['ratio'] = ratio

        

    def produceROOT(self):
        root_path = f"Efficiency/Btag_{self.suffix}_{self.era}.root"
        f = ROOT.TFile(root_path,'recreate')
        for cat in self.dict_hist.keys():
            print (f'Producing root file for category {cat}')
            for name in ['truth','btagged','ratio']:
                self.dict_hist[cat][name].Write(f'{cat}_{name}')
        f.Write()
        f.Close()
        print (f"Produced file {root_path}")

    def producePlots(self):
        for cat in self.dict_hist.keys():
            print (f'Plotting category {cat}')
            truthx = self.dict_hist[cat]['truth'].Project3D("x")
            truthy = self.dict_hist[cat]['truth'].Project3D("y")
            truthz = self.dict_hist[cat]['truth'].Project3D("z")

            truthx.SetTitle(f'True {cat};#eta')
            truthy.SetTitle(f'True {cat};P_{{T}}')
            truthz.SetTitle(f'True {cat};Btag discriminant')

            btagx = self.dict_hist[cat]['btagged'].Project3D("x")
            btagy = self.dict_hist[cat]['btagged'].Project3D("y")
            btagz = self.dict_hist[cat]['btagged'].Project3D("z")

            btagx.SetTitle(f'Btagged {cat};#eta')
            btagy.SetTitle(f'Btagged {cat};P_{{T}}')
            btagz.SetTitle(f'Btagged {cat};Btag discriminant')

            # /!\ for efficiency profile : cannot project ratio !!! -> Need to do ratio of profiles
    
            ratiox = btagx.Clone("ratiox") #self.dict_hist[cat]['ratio'].Project3D("x")
            ratiox.Divide(truthx)
            ratioy = btagy.Clone("ratioy") #self.dict_hist[cat]['ratio'].Project3D("y")
            ratioy.Divide(truthy)
            ratioz = btagz.Clone("ratioz") #self.dict_hist[cat]['ratio'].Project3D("z")
            ratioz.Divide(truthz)

            ratiox.SetTitle(f'Efficiency {cat};#eta')
            ratioy.SetTitle(f'Efficiency {cat};P_{{T}}')
            ratioz.SetTitle(f'Efficiency {cat};Btag discriminant')

            ratiox.GetYaxis().SetRangeUser(0.,1.) 
            ratioy.GetYaxis().SetRangeUser(0.,1.)
            ratioz.GetYaxis().SetRangeUser(0.,1.)

            list_hist = [truthx,truthy,truthz,
                         btagx,btagy,btagz,
                         ratiox,ratioy,ratioz]

            C = ROOT.TCanvas('c','c',2100,2100)
            C.Divide(3,3)

            for idx, h in enumerate(list_hist,1):
                C.cd(idx)
                h.Draw("H")
                h.SetLineWidth(2)
                h.GetXaxis().SetTitleOffset(1.2)

            C.Print(f"Efficiency/BtagEff_{self.suffix}_{cat}_{self.era}.pdf")

    def produceJson(self):
        for cat in self.dict_hist.keys():
            print (f'Producing json file for category {cat}')
            json_dict = {'dimension': 3, 'variables': ['AbsEta','Pt','BTagDiscri'], 'error_type': 'absolute'}
            ratio = self.dict_hist[cat]['ratio']
            xAxis = ratio.GetXaxis()
            yAxis = ratio.GetYaxis()
            zAxis = ratio.GetZaxis()
            data = []
            etaBinning = []
            PtBinning = []
            DiscBinning = []
            # Loop over X = eta #
            for x in range(1,ratio.GetNbinsX()+1):
                etaLow = round(xAxis.GetBinLowEdge(x),8)
                etaHigh = round(xAxis.GetBinUpEdge(x),8)
                etaDict = {'bin':[etaLow,etaHigh],'values':[]}
                if len(etaBinning) == 0:
                    etaBinning.append(etaLow)
                etaBinning.append(etaHigh)
                # Loop over Y = PT #
                for y in range(1,ratio.GetNbinsY()+1):
                    PtLow = round(yAxis.GetBinLowEdge(y),8)
                    PtHigh = round(yAxis.GetBinUpEdge(y),8)
                    PtDict = {'bin':[PtLow,PtHigh],'values':[]}
                    if len(PtBinning) == 0:
                        PtBinning.append(PtLow)
                    if PtHigh not in PtBinning:
                        PtBinning.append(PtHigh)
                    # Loop over Z = btag discriminant #
                    for z in range(1,ratio.GetNbinsZ()+1):
                        DiscLow = round(zAxis.GetBinLowEdge(z),8)
                        DiscHigh = round(zAxis.GetBinUpEdge(z),8)
                        if len(DiscBinning) == 0:
                            DiscBinning.append(DiscLow)
                        if DiscHigh not in DiscBinning:
                            DiscBinning.append(DiscHigh)
                        DiscDict = {'bin' : [DiscLow,DiscHigh],
                                    'value' : ratio.GetBinContent(x,y,z),
                                    'error_low' : 0., 
                                    'error_high' : 0} 
                        PtDict['values'].append(DiscDict)
                    etaDict['values'].append(PtDict)
                data.append(etaDict)
            json_dict['binning'] = {'x':etaBinning,'y':PtBinning, 'z':DiscBinning}
            json_dict['data'] = data

            # Save to json #
            json_path = f"Efficiency/BtagEff_{self.suffix}_{cat}_{self.era}.json"
            with open(json_path,'w') as f:
                json.dump(json_dict,f,indent=2)
            print (f"Saved json to {json_path}")


parser = argparse.ArgumentParser(description='Produce plots and json files of tagging efficiencies')
parser.add_argument('--path', action='store', required=True, type=str, 
                    help='Path to bamboo output')
parser.add_argument('--era', action='store', required=True, type=str, 
                    help='Era : 2016 | 2017 | 2018 ')
parser.add_argument('--plots', action='store_true', required=False, default=False,
                    help='Produce plots')
parser.add_argument('--json', action='store_true', required=False, default=False,
                    help='Produce json')
parser.add_argument('--root', action='store_true', required=False, default=False,
                    help='Produce root')
args = parser.parse_args()


hist_names_ak4 = {'bjets'    : {'truth' : 'N_ak4_truth_bjets',      'btagged': 'N_ak4_btagged_bjets'},
                  'cjets'    : {'truth' : 'N_ak4_truth_cjets',      'btagged': 'N_ak4_btagged_cjets'},
                  'lightjets': {'truth' : 'N_ak4_truth_lightjets',  'btagged': 'N_ak4_btagged_lightjets'}}
hist_names_ak8 = {'bjets'    : {'truth' : 'N_ak8_truth_bjets',      'btagged': 'N_ak8_btagged_bjets'},
                  'cjets'    : {'truth' : 'N_ak8_truth_cjets',      'btagged': 'N_ak8_btagged_cjets'},
                  'lightjets': {'truth' : 'N_ak8_truth_lightjets',  'btagged': 'N_ak8_btagged_lightjets'}}

instances = [
#    ComputeBtagEff(path=args.path,dict_names=hist_names_ak4,suffix='ak4_raw',era=args.era),
#    ComputeBtagEff(path=args.path,dict_names=hist_names_ak8,suffix='ak8_raw',era=args.era),
    ComputeBtagEff(path=args.path,dict_names=hist_names_ak4,suffix='ak4',era=args.era,rebinX=20,rebinY=20,rebinZ=20),
    ComputeBtagEff(path=args.path,dict_names=hist_names_ak8,suffix='ak8',era=args.era,rebinX=20,rebinY=20,rebinZ=20),
]


for instance in instances:
    if args.plots:
        instance.producePlots()
    if args.json:
        instance.produceJson()
    if args.root:
        instance.produceROOT()

