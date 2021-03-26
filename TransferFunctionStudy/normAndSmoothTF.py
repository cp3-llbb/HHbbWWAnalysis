#! /usr/bin/env python
"""
Courtesy of Brieux François :
https://github.com/BrieucF/HHTools/blob/MIS_RUNII/histFactory_hh/normAndSmoothTF.py
"""

import os
import sys
import ROOT 
import argparse
import copy
import numpy as np

from rebin import rebin2D

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = 2000
ROOT.gStyle.SetOptStat("rm")
ROOT.gStyle.SetOptFit()

def parseArguments():
    parser = argparse.ArgumentParser(description='Build transfer functions out of 2D histogram created by plotter.')
    parser.add_argument('-i', '--input', type=str, dest='input', required=True, help='Input file')
    parser.add_argument('-o', '--output', type=str, dest='output', required=True, help='Output file')
    parser.add_argument('-p', '--plot', action='store',dest='plot', required=False, type=str, default=None,help='Plot directory in "pdf/" directory')
    return parser.parse_args()


def normalizeDeltaE(hist):
    xAxis = hist.GetXaxis()
    yAxis = hist.GetYaxis()

    for i in range(1, xAxis.GetNbins()+1):
        binWidths = []
        for j in range(1, yAxis.GetNbins()+1):
            width = yAxis.GetBinWidth(j)
            if width not in binWidths:
                binWidths.append(width)

        binWidths.sort()
        minWidth = binWidths[0]

        for j in range(1, yAxis.GetNbins()+1):
            oldContent = hist.GetBinContent(i, j)
            widthRatio = yAxis.GetBinWidth(j)/minWidth
            hist.SetBinContent(i, j, oldContent/widthRatio)

    hist.Smooth(1)
    
    for i in range(1, xAxis.GetNbins()+1):
        
        integral = hist.Integral(i, i, 1, yAxis.GetNbins())

        for j in range(1, yAxis.GetNbins()+1):
            oldContent = hist.GetBinContent(i, j)

            if integral > 0:
                hist.SetBinContent(i, j, oldContent/integral)
            else:
                hist.SetBinContent(i, j, 0)

def plotSlices(hist2D,title):
    xAxis = hist2D.GetXaxis()
    yAxis = hist2D.GetYaxis()
    list_hist1D = []

    graph_mean = ROOT.TGraph(xAxis.GetNbins())
    graph_std = ROOT.TGraph(xAxis.GetNbins())
    graph_res = ROOT.TGraph(xAxis.GetNbins())
    
    for i in range(1, xAxis.GetNbins()+1):
        # Get plot range #
       # Project on Y axis #
        hist = hist2D.ProjectionY("%d_py"%i,i,i)
        if hist.Integral() <= 0.:
            continue

        # plot parameters #
        low_lim = xAxis.GetBinLowEdge(i)
        up_lim = xAxis.GetBinUpEdge(i)
        center = xAxis.GetBinCenter(i)
        hist.SetTitle(title+" (%0.f < E_{parton} < %0.f GeV)"%(low_lim,up_lim))
        hist.GetXaxis().SetTitle("E_{reco}-E_{gen}")
        hist.GetXaxis().SetTitleOffset(1.2)

        # Record mean and std to graph #
        graph_mean.SetPoint(i-1,center,hist.GetMean())
        graph_std.SetPoint(i-1,center,hist.GetStdDev())
        graph_res.SetPoint(i-1,center,hist.GetStdDev()/center)

        # Add to list #
        list_hist1D.append(hist)

    graph_mean.SetTitle("#DeltaE mean;E_{gen} [GeV]];Mean #DeltaE")
    graph_std.SetTitle("#DeltaE standard deviation;E_{gen} [GeV]];Std dev #DeltaE")
    graph_res.SetTitle("Resolution;E_{gen} [GeV]];#frac{#sigma(E)}{E}")
    graph_mean.GetHistogram().GetXaxis().SetRangeUser(0.,up_lim)
    graph_mean.GetHistogram().SetMinimum(0.)
    graph_std.GetHistogram().GetXaxis().SetRangeUser(0.,up_lim)
    graph_std.GetHistogram().SetMinimum(0.)
    graph_res.GetHistogram().GetXaxis().SetRangeUser(0.,up_lim)
    graph_res.GetHistogram().SetMinimum(0.)

    return list_hist1D, graph_mean, graph_std, graph_res

def findAndCorrectArtifacts(h,threshold=10):
    xAxis = h.GetXaxis()
    yAxis = h.GetYaxis()
    Nx = xAxis.GetNbins()
    Ny = yAxis.GetNbins()
    for x in range(1,Nx+1):
        for y in range(1,Ny+1):
            content = h.GetBinContent(x,y)
            c_xlow  = content
            c_xup   = content
            c_ylow  = content
            c_yup   = content
            if x > 1:
                c_xlow = h.GetBinContent(x-1,y)
            if y > 1: 
                c_ylow = h.GetBinContent(x,y-1)
            if x < Nx:
                c_xup = h.GetBinContent(x+1,y)
            if y < Ny:
                c_yup = h.GetBinContent(x,y+1)
            avg = (c_xlow+c_xup+c_ylow+c_yup)/4
            if avg!=0 and abs((content-avg)/avg) > threshold:
                h.SetBinContent(x,y,avg)
            #    print ()
            #    print (x,y)
            #    print (xAxis.GetBinCenter(x),yAxis.GetBinCenter(y))
            #    print (content,c_xlow,c_xup,c_ylow,c_yup)
            #    print (content,avg)
            #    print (content/avg)
            #    print ((content-avg)/avg)

def normAndSmooth(TFset, inFile, outFile, plot):
    inFile = ROOT.TFile.Open(inFile)
    outFile = ROOT.TFile(outFile, "recreate")
   
    for TF in TFset:
        print ("Looking over %s"%(TF['base']))
        inputHists = []
        keys = [key.GetName() for key in inFile.GetListOfKeys()]
        for hist in TF["histNames"]:
            if hist not in keys:
                raise RuntimeError("Histogram {} not available in {}".format(hist,inFile))
            inputHists.append( inFile.Get(hist) )
        
        DeltaEvsE = inputHists[0].Clone( TF["base"] )
        for hist in inputHists[1:]:
            DeltaEvsE.Add(hist)


        DeltaEvsE.SetTitle("Hist 2D from bamboo")
        DeltaEvsE.GetXaxis().SetTitle("E_{gen}")
        DeltaEvsE.GetYaxis().SetTitle("E_{rec} - E_{gen}")

        #findAndCorrectArtifacts(DeltaEvsE,10)

        DeltaEvsE_rebin = DeltaEvsE.Clone("rebin")
        if 'rebinX' in TF.keys() and (isinstance(TF['rebinX'],list) or isinstance(TF['rebinX'],np.ndarray)) and 'rebinY' in TF.keys() and (isinstance(TF['rebinY'],list) or isinstance(TF['rebinY'],np.ndarray)):
            DeltaEvsE_rebin = rebin2D(DeltaEvsE_rebin,TF['rebinX'],TF['rebinY'],TF["base"]+'_rebin')
        if 'rebinX' in TF.keys() and isinstance(TF['rebinX'],int):
            DeltaEvsE.RebinX(TF['rebinX'])
        if 'rebinY' in TF.keys() and isinstance(TF['rebinX'],int):
            DeltaEvsE.RebinY(TF['rebinY'])
        DeltaEvsE_rebin.SetTitle("Hist 2D from bamboo [rebinned]")

        DeltaEvsE_Norm = DeltaEvsE_rebin.Clone( TF["norm"] )
        DeltaEvsE_Norm.SetTitle("Transfer function")

        normalizeDeltaE(DeltaEvsE_Norm)
        
        if plot is not None:
            path_out = "pdf/{}".format(plot)
            if not os.path.exists(path_out):
                os.makedirs(path_out)
            pdfName = os.path.join(path_out,"{}.pdf".format(TF["norm"]))
            list_hist, g_mean, g_std, g_res = plotSlices(DeltaEvsE_Norm,TF["title"])
            C = ROOT.TCanvas("C","C",800,600)
            C.Print(pdfName+"[")

            # Draw TF #
            C.SetLogz(True)
            DeltaEvsE.SetContour(100)
            DeltaEvsE.Draw("colz")
            C.Print(pdfName)
            C.Clear()
            DeltaEvsE_rebin.SetContour(100)
            DeltaEvsE_rebin.Draw("colz")
            C.Print(pdfName)
            C.Clear()
            DeltaEvsE_Norm.SetContour(100)
            DeltaEvsE_Norm.Draw("colz")
            C.Print(pdfName)
            C.SetLogz(False)

            # Draw graphs #
            C.Clear()
            g_mean.Draw()
            C.Print(pdfName)
            C.Clear()
            g_std.Draw()
            #g_std.Fit("pol1")
            #f = g_std.GetListOfFunctions().FindObject("pol1")
            ##f.FixParameter(0,0)
            #print (TF["title"]+" Resolution = %0.2f"%(f.GetParameter(1)*100))
            C.Print(pdfName)
            C.Clear()
            #C.SetLogy(True)
            g_res.Draw()
            C.Print(pdfName)
            #C.SetLogy(False)


            # Draw profiles #
            C.Clear()
            C.SetCanvasSize(2000,900)
            C.Divide(2,1)
            pad1 = C.cd(1)
            pad2 = C.cd(2)
            pad1.SetLogy(True)
            threshold = 0.001
            for hist in list_hist:
                pad1.Clear()
                pad1.cd()
                hist_log = hist.Clone("log")
                hist_log.Draw("hist")

                pad2.Clear()
                pad2.cd()
                bins_above_threshold = [i for i in range(1,hist.GetNbinsX()+1) if hist.GetBinContent(i)>threshold]
                if len(bins_above_threshold) != 0:
                    origin_bin = hist.FindBin(hist.GetMaximum())
                    max_range = max(origin_bin-min(bins_above_threshold),max(bins_above_threshold)-origin_bin)
                    hist.GetXaxis().SetRange(max(1,origin_bin-max_range),min(hist.GetXaxis().GetNbins(),origin_bin+max_range))

                hist.Draw("hist")
                C.Print(pdfName)
            C.Clear()
            C.Print(pdfName+"]")

        
        DeltaEvsE.Write()
        DeltaEvsE_rebin.Write()
        DeltaEvsE_Norm.Write()
    
    outFile.Close()
    inFile.Close()

if __name__ == "__main__":
    options = parseArguments()
    TFset = [
#                {
#                    "histNames": ['ak4b_TF'],
#                    "base" : "b_TF",
#                    "norm" : "b_TF_norm",
#                    "title": "b transfer function",
#                    "rebinX": np.r_[0:100:5, 100:200:10, 200:300:20, 300:500:50, 500:1001:250, 3000 ],
#                    "rebinY": np.r_[-500:501:5]
#                },
#                {
#                    "histNames": ['ak4b_TF_bregCorr'],
#                    "base" : "b_TF_bregCorr",
#                    "norm" : "b_TF_bregCorr_norm",
#                    "title": "b transfer function (with regression correction)",
#                    "rebinX": np.r_[0:100:5, 100:200:10, 200:300:20, 300:500:50, 500:1001:250, 3000 ],
#                    "rebinY": np.r_[-500:501:5]
#                },
#                {
#                    "histNames": ['ak4c_TF'],
#                    "base" : "c_TF",
#                    "norm" : "c_TF_norm",
#                    "title": "c transfer function",
#                    "rebinX": np.r_[0:100:5, 100:200:10, 200:300:20, 300:500:50, 500:1001:250, 3000 ],
#                    "rebinY": np.r_[-500:501:5]
#                },
#                {
#                    "histNames": ['ak4l_TF'],
#                    "base" : "l_TF",
#                    "norm" : "l_TF_norm",
#                    "title": "lightjet transfer function",
#                    "rebinX": np.r_[0:100:5, 100:200:10, 200:300:20, 300:500:50, 500:1001:250, 3000 ],
#                    "rebinY": np.r_[-500:501:5]
#                },
#                {
#                    "histNames": ['e_TF'],
#                    "base" : "e_TF",
#                    "norm" : "e_TF_norm",
#                    "title": "e^{#pm} transfer function",
#                    "rebinX": np.r_[0:100:5, 100:200:10, 200:300:20, 300:500:100, 500:1001:500, 3000 ],
#                    "rebinY": np.r_[-500:501:1]
#                },
                {
                    "histNames": ['m_TF'],
                    "base" : "mu_TF",
                    "norm" : "mu_TF_norm",
                    "title": "#mu^{#pm} transfer function",
                    "rebinX": np.r_[0:100:5, 100:200:10, 200:300:20, 300:500:100, 500:1001:500, 3000 ],
                    "rebinY": np.r_[-500:-300:50, -300:-200:20, -200:-100:10, -100:-40:2 , -40:40:1 , 40:100:2 , 100:200:10, 200:300:20, 300:500:50],
                },
            ]

    normAndSmooth(TFset, options.input, options.output, options.plot)


