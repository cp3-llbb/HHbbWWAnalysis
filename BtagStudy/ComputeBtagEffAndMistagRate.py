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


class ComputeBtagEffAndMistag:
    def __init__(self,path,dict_names):
        self.path = os.path.abspath(path)
        self.dict_names = dict_names

        self.recoverHistograms()
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
                h_tot = r.Get(self.dict_names[cat]['total'])
                h_btag = r.Get(self.dict_names[cat]['btagged'])

                # Underflow and overflow bins #
                h_tot.SetBinContent(1,h_tot.GetBinContent(0)+h_tot.GetBinContent(1))
                h_btag.SetBinContent(1,h_btag.GetBinContent(0)+h_btag.GetBinContent(1))
                h_tot.SetBinContent(h_tot.GetNbinsX(),h_tot.GetBinContent(h_tot.GetNbinsX())+h_tot.GetBinContent(h_tot.GetNbinsX()+1))
                h_btag.SetBinContent(h_btag.GetNbinsX(),h_btag.GetBinContent(h_btag.GetNbinsX())+h_btag.GetBinContent(h_btag.GetNbinsX()+1))

                # Generate subdict #
                if cat not in self.dict_hist.keys():
                    self.dict_hist[cat] = {}

                # Total hist #
                if 'total' not in self.dict_hist[cat]:
                    self.dict_hist[cat]['total'] = copy.deepcopy(h_tot.Clone(cat+'_total'))
                    self.dict_hist[cat]['total'].Sumw2()
                else:
                    self.dict_hist[cat]['total'].Add(h_tot)

                # Btagged hist #
                if 'btagged' not in self.dict_hist[cat]:
                    self.dict_hist[cat]['btagged'] = copy.deepcopy(h_btag.Clone(cat+'_btagged'))
                    self.dict_hist[cat]['btagged'].Sumw2()
                else:
                    self.dict_hist[cat]['btagged'].Add(h_btag)
            r.Close()
                    
    def makeRatio(self):
        # Add btag mistag rate #
        self.dict_hist['mistag'] = {'total' : copy.deepcopy(self.dict_hist['cjets']['total']),
                                    'btagged' : copy.deepcopy(self.dict_hist['cjets']['btagged'])}
        self.dict_hist['mistag']['total'].Add(self.dict_hist['lightjets']['total'])
        self.dict_hist['mistag']['btagged'].Add(self.dict_hist['lightjets']['btagged'])

        # Produce ratio for each category #
        for cat in self.dict_hist.keys():
            ratio = self.dict_hist[cat]['btagged'].Clone(cat+'_ratio')
            ratio.Divide(self.dict_hist[cat]['total'])
            self.dict_hist[cat]['ratio'] = ratio

        

    def produceROOT(self):
        f = ROOT.TFile("Btag.root",'recreate')
        for cat in self.dict_hist.keys():
            for name in ['total','btagged','ratio']:
                self.dict_hist[cat][name].Write(cat+name)
        f.Write()
        f.Close()

    def producePlots(self):
        for cat in self.dict_hist.keys():
            C = ROOT.TCanvas('c','c',2100,2100)
            C.Divide(3,3)

            totx = self.dict_hist[cat]['total'].Project3D("x")
            toty = self.dict_hist[cat]['total'].Project3D("y")
            totz = self.dict_hist[cat]['total'].Project3D("z")

            totx.SetTitle('Total %s;#eta'%cat)
            toty.SetTitle('Total %s;P_{T}'%cat)
            totz.SetTitle('Total %s;Btag discriminant'%cat)

            btagx = self.dict_hist[cat]['btagged'].Project3D("x")
            btagy = self.dict_hist[cat]['btagged'].Project3D("y")
            btagz = self.dict_hist[cat]['btagged'].Project3D("z")

            btagx.SetTitle('Btagged %s;#eta'%cat)
            btagy.SetTitle('Btagged %s;P_{T}'%cat)
            btagz.SetTitle('Btagged %s;Btag discriminant'%cat)

            # /!\ for efficiency profile : cannot project ratio !!! -> Need to do ratio of profiles
    
            ratiox = btagx.Clone("ratiox")
            ratiox.Divide(totx)
            ratioy = btagy.Clone("ratioy")
            ratioy.Divide(toty)
            ratioz = btagz.Clone("ratioz")
            ratioz.Divide(totz)

            ratiox.SetTitle('Efficiency %s;#eta'%cat)
            ratioy.SetTitle('Efficiency %s;P_{T}'%cat)
            ratioz.SetTitle('Efficiency %s;Btag discriminant'%cat)

            list_hist = [totx,toty,totz,
                         btagx,btagy,btagz,
                         ratiox,ratioy,ratioz]

            for idx, h in enumerate(list_hist,1):
                C.cd(idx)
                h.Draw("H")
                h.SetLineWidth(2)
                h.GetXaxis().SetTitleOffset(1.2)

            C.Print("BtagEff_%s.pdf"%cat)

    def produceJson(self):
        for cat in self.dict_hist.keys():
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
                etaLow = xAxis.GetBinLowEdge(x) 
                etaHigh = xAxis.GetBinUpEdge(x)
                etaDict = {'bin':[etaLow,etaHigh],'values':[]}
                if len(etaBinning) == 0:
                    etaBinning.append(etaLow)
                etaBinning.append(etaHigh)
                # Loop over Y = PT #
                for y in range(1,ratio.GetNbinsY()+1):
                    PtLow = yAxis.GetBinLowEdge(y) 
                    PtHigh = yAxis.GetBinUpEdge(y)
                    PtDict = {'bin':[PtLow,PtHigh],'values':[]}
                    if len(PtBinning) == 0:
                        PtBinning.append(PtLow)
                    PtBinning.append(PtHigh)
                    # Loop over Z = btag discriminant #
                    for z in range(1,ratio.GetNbinsZ()+1):
                        DiscLow = zAxis.GetBinLowEdge(z)
                        DiscHigh = zAxis.GetBinUpEdge(z)
                        if len(DiscBinning) == 0:
                            DiscBinning.append(DiscLow)
                        if DiscHigh not in DiscBinning:
                            DiscBinning.append(DiscHigh)
                        DiscDict = {'bin' : [DiscLow,DiscHigh],
                                    'value' : ratio.GetBinContent(x,y,z),
                                    'error_low' : ratio.GetBinError(x,y,z),
                                    'error_high' : ratio.GetBinError(x,y,z)}
                        PtDict['values'].append(DiscDict)
                    etaDict['values'].append(PtDict)
                data.append(etaDict)
            json_dict['binning'] = {'x':etaBinning,'y':PtBinning, 'z':DiscBinning}
            json_dict['data'] = data

            # Save to json #
            with open("BtagEff_%s.json"%cat,'w') as f:
                json.dump(json_dict,f,indent=2)
            print ("Saved json to BtagEff_%s.json"%cat)


parser = argparse.ArgumentParser(description='Produce plots and json files of btagging efficiency and mistage rate')
parser.add_argument('--path', action='store', required=True, type=str, 
                    help='Path to bamboo output')
parser.add_argument('--plots', action='store_true', required=False, default=False,
                    help='Produce plots')
parser.add_argument('--json', action='store_true', required=False, default=False,
                    help='Produce json')
parser.add_argument('--root', action='store_true', required=False, default=False,
                    help='Produce root')
args = parser.parse_args()


hist_names = {'bjets' : {'total' : 'N_total_bjets', 'btagged': 'N_btagged_bjets'},
             'cjets' : {'total' : 'N_total_cjets', 'btagged': 'N_btagged_cjets'},
             'lightjets' : {'total' : 'N_total_lightjets', 'btagged': 'N_btagged_lightjets'}}

instance = ComputeBtagEffAndMistag(args.path,hist_names)

if args.plots:
    instance.producePlots()
if args.json:
    instance.produceJson()
if args.root:
    instance.produceROOT()

