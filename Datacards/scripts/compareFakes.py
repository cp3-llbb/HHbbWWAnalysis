import os
import sys
import glob
import copy
import json 
import math
import argparse
from array import array
import numpy as np
import ROOT

from IPython import embed

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)

parser = argparse.ArgumentParser(description='Comparing DY')
parser.add_argument('--mc', action='store', required=True, type=str,
                    help='Path to Fake MC shapes')
parser.add_argument('--cl', action='store', required=True, type=str,
                    help='Path to Fake Closure shapes')
parser.add_argument('--name', action='store', required=True, type=str,
                    help='Name of the output files')
args = parser.parse_args()

#if '{era}' not in args.mc:
#    raise RuntimeError(f'Missing {{era}} in {args.mc}')

factors = {}

plotDict_MC = {}
plotDict_CL = {}

path_MC  = args.mc
path_CL  = args.cl

plotNames = [os.path.basename(name) for name in glob.glob(os.path.join(path_MC,'*root'))]


for plotName in sorted(plotNames):
    f_MC = ROOT.TFile(os.path.join(path_MC,plotName))
    f_CL = ROOT.TFile(os.path.join(path_CL,plotName))

    plotName = plotName.replace('.root','')

    h_MC = f_MC.Get("TT")
    h_CL = f_CL.Get("TT")

    #xn = [0.,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.9,0.95,1.]
    xn = np.linspace(0,1,50)

    h_MC = h_MC.Rebin(len(xn)-1,'rebin',array('d',xn))
    h_CL = h_CL.Rebin(len(xn)-1,'rebin',array('d',xn))

    for i in range(1,h_MC.GetNbinsX()+1):
        if h_MC.GetBinContent(i)<1e-3:
            h_MC.SetBinContent(i,0.)
    for i in range(1,h_CL.GetNbinsX()+1):
        if h_CL.GetBinContent(i)<1e-3:
            h_CL.SetBinContent(i,0.)
        if h_MC.GetBinContent(i) == 0.:
            h_CL.SetBinContent(i,0.)

    if plotName not in plotDict_MC.keys():
        plotDict_MC[plotName] = copy.deepcopy(h_MC)
    else:
        plotDict_MC[plotName].Add(h_MC)

    if plotName not in plotDict_CL.keys():
        plotDict_CL[plotName] = copy.deepcopy(h_CL)
    else:
        plotDict_CL[plotName].Add(h_CL)

    f_MC.Close()
    f_CL.Close()

def getY(h):
    return np.array([h.GetBinContent(i) for i in range(1,h.GetNbinsX()+1)])
def getX(h):
    return np.array([h.GetXaxis().GetBinCenter(i) for i in range(1,h.GetNbinsX()+1)])
def getEdges(h):
    return np.array([h.GetXaxis().GetBinLowEdge(i) for i in range(1,h.GetNbinsX()+2)])

categories = sorted(list(plotDict_MC.keys()))

C = ROOT.TCanvas("C","C",700,800)
pdfName = f"CompareFakes/ComparisonFakes_{args.name}.pdf"
C.Print(pdfName+"[")
data = {}
for cat in categories:
    C.Clear()

    h_MC = plotDict_MC[cat]
    h_CL = plotDict_CL[cat]


    N_MC = h_MC.Integral()
    N_CL = h_CL.Integral()

    h_MC.SetTitle(f"{cat};;Events / {h_MC.GetXaxis().GetBinWidth(1):0.2f}")
    hmax = max(h_MC.GetMaximum(),h_CL.GetMaximum())
    #h_MC.SetMaximum(hmax*100)
    
    h_MC.SetLineColor(1)
    h_CL.SetLineColor(632)
    
    h_MC.SetMarkerColor(1)
    h_CL.SetMarkerColor(632)
    
    h_MC.SetMarkerStyle(4)
    h_CL.SetMarkerStyle(25)

    x_values = getX(h_MC)
    y_MC = getY(h_MC)
    y_CL = getY(h_CL)

    if N_CL == 0 or N_MC == 0:
        # Default params #
        nom   = 2
        slope = 1
        cog   = np.mean(x_values)

        # Plot only upper part (no ratio) #
        C.Divide(1)
        pad1 = C.cd(1)
        pad1.SetLogy()
        pad1.SetTopMargin(0.1)
        pad1.SetBottomMargin(0.1)
        pad1.SetLeftMargin(0.1)
        pad1.SetRightMargin(0.05)
        pad1.Draw()
        
        h_MC.Draw("e1p")
        h_CL.Draw("e1p same")
        
        leg_up = ROOT.TLegend(0.6,0.6,0.9,0.9)
        leg_up.AddEntry(h_MC,"#splitline{{MC Fakes}}{{Integral = {:0.3f}}}".format(N_MC))
        leg_up.AddEntry(h_CL,"#splitline{{MC Closure}}{{Integral = {:0.3f}}}".format(N_CL))
        leg_up.SetBorderSize(0)
        leg_up.SetFillStyle(0)
        leg_up.Draw()

    else:
        cog = h_MC.GetMean() # Center of gravity
        nom = max(0,min(2,N_MC/N_CL))
        h_MC.Scale(1./N_MC)
        h_CL.Scale(1./N_CL)

        
        ratio = h_MC.Clone("ratio")
        ratio.Divide(h_CL)
        ratio.SetTitle("")
        
        C.Divide(1,2)

        pad1 = C.cd(1)
        pad1.SetLogy()
        pad1.SetTopMargin(0.1)
        pad1.SetBottomMargin(0.05)
        pad1.SetLeftMargin(0.1)
        pad1.SetRightMargin(0.05)
        pad1.Draw()
        
        h_MC.Draw("e1p")
        h_CL.Draw("e1p same")
        
        leg_up = ROOT.TLegend(0.6,0.6,0.9,0.9)
        leg_up.AddEntry(h_MC,"#splitline{{MC Fakes}}{{Integral = {:0.3f}}}".format(N_MC))
        leg_up.AddEntry(h_CL,"#splitline{{MC Closure}}{{Integral = {:0.3f}}}".format(N_CL))
        leg_up.SetBorderSize(0)
        leg_up.SetFillStyle(0)
        leg_up.Draw()
        
        pad2 = C.cd(2)
        pad2.SetTopMargin(0.0)
        pad2.SetBottomMargin(0.10)
        pad2.SetLeftMargin(0.1)
        pad2.SetRightMargin(0.05)
        pad2.Draw()
        
        y_values = getY(ratio)
        edges_cog = getEdges(ratio)-cog
        ratio_cog = ROOT.TH1F("ratio_cog","ratio_cog",edges_cog.shape[0]-1,edges_cog)
        for i in range(1,ratio_cog.GetNbinsX()+1):
            ratio_cog.SetBinContent(i,ratio.GetBinContent(i))
            ratio_cog.SetBinError(i,ratio.GetBinError(i))
        ratio_cog.SetTitle("")
        ratio_cog.GetXaxis().SetTitle("x-#bar{x}")
        ratio_cog.Draw("e1p")
#        np.testing.assert_allclose(x_values - cog,getX(ratio_cog))
#        np.testing.assert_allclose(y_values ,getY(ratio_cog))

        #slope, _ = np.polyfit(x_values - cog, y_values, 1) 
        ratio_cog.Fit("pol1","QS","")
        fit = ratio_cog.GetFunction('pol1')
        slope = fit.GetParameter(1)
        k = fit.GetParameter(0)
        func = lambda x : slope*x+k

        print (f'Category {cat} : cog = {cog:5.3f}, nom = {nom:5.3f}, slope = {slope:5.3f}')

        leg_down = ROOT.TLegend(0.55,0.65,0.9,1.0)
        leg_down.AddEntry(ratio,"\splitline{{#frac{{MC Fake}}{{MC Closure}}}}{{Slope = {:.3f} x + {:.3f}}}".format(slope,k))
        leg_down.SetTextSize(0.04)
        leg_down.SetBorderSize(0)
        leg_down.SetFillStyle(0)
        leg_down.Draw()
        line = ROOT.TLine(min(edges_cog),1.,max(edges_cog),1.)
        line.Draw()
   
    # Print canvas and save data #
    C.Print(pdfName,f"Title:{cat}")
    data[cat] = {'cog':cog,'nom':nom,'slope':slope}
#
    
C.Print(pdfName+"]")

jsonPath = f'CompareFakes/factors_{args.name}.json'
with open(jsonPath,'w') as handle:
    json.dump(data,handle,indent=4)
print (f"Saved {jsonPath}")
