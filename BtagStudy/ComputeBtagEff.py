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
    def __init__(self,path,dict_names,era,rebinX=None,rebinY=None,rebinZ=None):
        self.path = os.path.abspath(path)
        self.dict_names = dict_names
        self.era = era
        self.rebinX = rebinX
        self.rebinY = rebinY
        self.rebinZ = rebinZ

        self.recoverHistograms()
        self.rebinHistograms()
        self.makeRatio()

    def recoverHistograms(self):
        self.dict_hist = {}
        for f in glob.glob(os.path.join(self.path,'results','*root')):
            if '__skeleton__' in f:
                continue
            r = ROOT.TFile(f)
            print ('Looking at ',f)
            for cat in self.dict_names.keys():
                # Get hist #
                h_tot = r.Get(self.dict_names[cat]['truth'])
                h_btag = r.Get(self.dict_names[cat]['btagged'])

                # Generate subdict #
                if cat not in self.dict_hist.keys():
                    self.dict_hist[cat] = {}

                # Total hist #
                if 'truth' not in self.dict_hist[cat]:
                    self.dict_hist[cat]['truth'] = copy.deepcopy(h_tot.Clone(cat+'_truth'))
                else:
                    self.dict_hist[cat]['truth'].Add(h_tot)

                # Btagged hist #
                if 'btagged' not in self.dict_hist[cat]:
                    self.dict_hist[cat]['btagged'] = copy.deepcopy(h_btag.Clone(cat+'_btagged'))
                else:
                    self.dict_hist[cat]['btagged'].Add(h_btag)
            r.Close()
                    
    def rebinHistograms(self):
        for cat in self.dict_names.keys():
            print ("Rebinning category %s"%cat)
            assert  self.dict_hist[cat]['btagged'].GetNbinsX() == self.dict_hist[cat]['truth'].GetNbinsX()
            assert  self.dict_hist[cat]['btagged'].GetNbinsY() == self.dict_hist[cat]['truth'].GetNbinsY()
            assert  self.dict_hist[cat]['btagged'].GetNbinsZ() == self.dict_hist[cat]['truth'].GetNbinsZ()
            NbinsX = self.dict_hist[cat]['btagged'].GetNbinsX()
            NbinsY = self.dict_hist[cat]['btagged'].GetNbinsY()
            NbinsZ = self.dict_hist[cat]['btagged'].GetNbinsZ()

            assert NbinsX % self.rebinX == 0
            assert NbinsY % self.rebinY == 0
            assert NbinsZ % self.rebinZ == 0
            
            groupX = NbinsX // self.rebinX
            groupY = NbinsY // self.rebinY
            groupZ = NbinsZ // self.rebinZ

            self.dict_hist[cat]['btagged'] = self.dict_hist[cat]['btagged'].Rebin3D(groupX,groupY,groupZ)
            self.dict_hist[cat]['truth'] = self.dict_hist[cat]['truth'].Rebin3D(groupX,groupY,groupZ)

    def makeRatio(self):
#        # Add btag mistag rate #
#        self.dict_hist['mistag'] = {'truth' : copy.deepcopy(self.dict_hist['cjets']['truth']),
#                                    'btagged' : copy.deepcopy(self.dict_hist['cjets']['btagged'])}
#        self.dict_hist['mistag']['truth'].Add(self.dict_hist['lightjets']['truth'])
#        self.dict_hist['mistag']['btagged'].Add(self.dict_hist['lightjets']['btagged'])

        # Produce ratio for each category #
        for cat in self.dict_hist.keys():
            ratio = self.dict_hist[cat]['btagged'].Clone(cat+'_ratio')
            ratio.Divide(self.dict_hist[cat]['truth'])
            self.dict_hist[cat]['ratio'] = ratio

        

    def produceROOT(self):
        f = ROOT.TFile("Btag_%.root"%self.era,'recreate')
        for cat in self.dict_hist.keys():
            print ('Producing root file for category %s'%cat)
            for name in ['truth','btagged','ratio']:
                self.dict_hist[cat][name].Write(cat+name)
        f.Write()
        f.Close()
        print ("Produced file Btag_%s.root"%self.era)

    def producePlots(self):
        for cat in self.dict_hist.keys():
            print ('Plotting category %s'%cat)
            truthx = self.dict_hist[cat]['truth'].Project3D("x")
            truthy = self.dict_hist[cat]['truth'].Project3D("y")
            truthz = self.dict_hist[cat]['truth'].Project3D("z")

            truthx.SetTitle('True %s;#eta'%cat)
            truthy.SetTitle('True %s;P_{T}'%cat)
            truthz.SetTitle('True %s;Btag discriminant'%cat)

            btagx = self.dict_hist[cat]['btagged'].Project3D("x")
            btagy = self.dict_hist[cat]['btagged'].Project3D("y")
            btagz = self.dict_hist[cat]['btagged'].Project3D("z")

            btagx.SetTitle('Btagged %s;#eta'%cat)
            btagy.SetTitle('Btagged %s;P_{T}'%cat)
            btagz.SetTitle('Btagged %s;Btag discriminant'%cat)

            # /!\ for efficiency profile : cannot project ratio !!! -> Need to do ratio of profiles
    
            ratiox = btagx.Clone("ratiox") #self.dict_hist[cat]['ratio'].Project3D("x")
            ratiox.Divide(truthx)
            ratioy = btagy.Clone("ratioy") #self.dict_hist[cat]['ratio'].Project3D("y")
            ratioy.Divide(truthy)
            ratioz = btagz.Clone("ratioz") #self.dict_hist[cat]['ratio'].Project3D("z")
            ratioz.Divide(truthz)

            ratiox.SetTitle('Efficiency %s;#eta'%cat)
            ratioy.SetTitle('Efficiency %s;P_{T}'%cat)
            ratioz.SetTitle('Efficiency %s;Btag discriminant'%cat)

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

            C.Print("BtagEff_%s_%s.pdf"%(cat,self.era))

    def produceJson(self):
        for cat in self.dict_hist.keys():
            print ('Producing json file for category %s'%cat)
            json_dict = {'dimension': 3, 'variables': ['Eta','Pt','BTagDiscri'], 'error_type': 'absolute'}
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
                                    'error_low' : 0., #ratio.GetBinError(x,y,z),
                                    'error_high' : 0} #ratio.GetBinError(x,y,z)}
                        PtDict['values'].append(DiscDict)
                    etaDict['values'].append(PtDict)
                data.append(etaDict)
            json_dict['binning'] = {'x':etaBinning,'y':PtBinning, 'z':DiscBinning}
            json_dict['data'] = data

            # Save to json #
            with open("BtagEff_%s_%s.json"%(cat,self.era),'w') as f:
                json.dump(json_dict,f,indent=2)
            print ("Saved json to BtagEff_%s_%s.json"%(cat,self.era))


parser = argparse.ArgumentParser(description='Produce plots and json files of btagging efficiency and mistage rate')
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


hist_names = {'bjets'   : {'truth' : 'N_truth_bjets', 'btagged': 'N_btagged_bjets'},
              'cjets'    : {'truth' : 'N_truth_cjets', 'btagged': 'N_btagged_cjets'},
              'lightjets': {'truth' : 'N_truth_lightjets', 'btagged': 'N_btagged_lightjets'}}

instance = ComputeBtagEff(args.path,hist_names,era=args.era,rebinX=10,rebinY=20,rebinZ=1)

if args.plots:
    instance.producePlots()
if args.json:
    instance.produceJson()
if args.root:
    instance.produceROOT()

