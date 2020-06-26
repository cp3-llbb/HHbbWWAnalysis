#! /usr/bin/env python
"""
Courtesy of Brieux François :
https://github.com/BrieucF/HHTools/blob/MIS_RUNII/histFactory_hh/normAndSmoothTF.py
"""

import ROOT as R
import argparse
import copy

R.gROOT.SetBatch(True)
R.gErrorIgnoreLevel = 2000
R.gStyle.SetOptStat("rm")
R.gStyle.SetOptFit()

def parseArguments():
    parser = argparse.ArgumentParser(description='Build transfer functions out of 2D histogram created by plotter.')
    parser.add_argument('-i', '--input', type=str, dest='input', required=True, help='Input file')
    parser.add_argument('-o', '--output', type=str, dest='output', required=True, help='Output file')
    parser.add_argument('-p', '--plot', action='store_true',dest='plot', required=False, default=False,help='Produce plots')
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
    threshold = 0.001
    list_hist1D = []

    graph_mean = R.TGraph(xAxis.GetNbins())
    graph_std = R.TGraph(xAxis.GetNbins())
    graph_res = R.TGraph(xAxis.GetNbins())
    
    for i in range(1, xAxis.GetNbins()+1):
        # Get plot range #
        bins_above_threshold = []
        for j in range(1, yAxis.GetNbins()+1):
            if hist2D.GetBinContent(i,j) > threshold:
                bins_above_threshold.append(j)

        # Project on Y axis #
        hist = hist2D.ProjectionY("%d_py"%i,i,i)
        if hist.GetEntries()==0:
            continue

        # plot parameters #
        origin_bin = hist.FindBin(hist.GetMaximum())
        max_range = max(origin_bin-min(bins_above_threshold),max(bins_above_threshold)-origin_bin)
        hist.GetXaxis().SetRange(origin_bin-max_range,origin_bin+max_range)
        low_lim = xAxis.GetBinLowEdge(i)
        up_lim = xAxis.GetBinUpEdge(i)
        hist.SetTitle(title+" (%0.f < E_{parton} < %0.f GeV)"%(low_lim,up_lim))
        hist.GetXaxis().SetTitle("E_{reco}-E_{gen}")
        hist.GetXaxis().SetTitleOffset(1.2)

        # Record mean and std to graph #
        graph_mean.SetPoint(i-1,(up_lim-low_lim/2),hist.GetMean())
        graph_std.SetPoint(i-1,(up_lim-low_lim/2),hist.GetStdDev())
        graph_res.SetPoint(i-1,(up_lim-low_lim/2),hist.GetStdDev()/(up_lim-low_lim/2))

        # Add to list #
        list_hist1D.append(hist)

    graph_mean.SetTitle("#DeltaE mean;E_{gen} [GeV]];Mean #DeltaE")
    graph_std.SetTitle("#DeltaE standard deviation;E_{gen} [GeV]];Std dev #DeltaE")
    graph_res.SetTitle("Resolution;E_{gen} [GeV]];#frac{#sigma(E)}{E}")
    graph_mean.GetHistogram().GetXaxis().SetRangeUser(0.,(up_lim-low_lim/2))
    graph_mean.GetHistogram().GetYaxis().SetRangeUser(0.,hist.GetMean())
    graph_std.GetHistogram().GetXaxis().SetRangeUser(0.,(up_lim-low_lim/2))
    graph_std.GetHistogram().GetYaxis().SetRangeUser(0.,hist.GetStdDev())
    graph_res.GetHistogram().GetXaxis().SetRangeUser(0.,(up_lim-low_lim/2))

    return list_hist1D, graph_mean, graph_std, graph_res

def normAndSmooth(TFset, inFile, outFile, plot):
    inFile = R.TFile.Open(inFile)
    outFile = R.TFile(outFile, "recreate")
   
    for TF in TFset:
        print ("Looking over %s"%(TF['base']))
        inputHists = []
        for hist in TF["histNames"]:
            inputHists.append( inFile.Get(hist) )
        
        DeltaEvsE = inputHists[0].Clone( TF["base"] )
        for hist in inputHists[1:]:
            DeltaEvsE.Add(hist)

        DeltaEvsE.SetTitle("Hist 2D from bamboo")
        DeltaEvsE.GetXaxis().SetTitle("E_{gen}")
        DeltaEvsE.GetYaxis().SetTitle("E_{rec} - E_{gen}")
        
        DeltaEvsE_Norm = DeltaEvsE.Clone( TF["norm"] )
        DeltaEvsE_Norm.SetTitle("Transfer function")
        normalizeDeltaE(DeltaEvsE_Norm)

        if plot:
            list_hist, g_mean, g_std, g_res = plotSlices(DeltaEvsE_Norm,TF["title"])
            C = R.TCanvas("C","C",800,600)
            C.Print(TF["norm"]+".pdf[")

            # Draw TF #
            DeltaEvsE.Draw("colz")
            C.Print(TF["norm"]+".pdf")
            C.Clear()
            DeltaEvsE_Norm.Draw("colz")
            C.Print(TF["norm"]+".pdf")
            DeltaEvsE_Norm.GetYaxis().SetRangeUser(-30,30)
            C.Clear()
            DeltaEvsE_Norm.Draw("colz")
            C.Print(TF["norm"]+".pdf")

            # Draw graphs #
            C.Clear()
            g_mean.Draw()
            C.Print(TF["norm"]+".pdf")
            C.Clear()
            g_std.Draw()
            #g_std.Fit("pol1")
            #f = g_std.GetListOfFunctions().FindObject("pol1")
            ##f.FixParameter(0,0)
            #print (TF["title"]+" Resolution = %0.2f"%(f.GetParameter(1)*100))
            C.Print(TF["norm"]+".pdf")
            C.Clear()
            g_res.Draw()
            C.Print(TF["norm"]+".pdf")


            # Draw profiles #
            for hist in list_hist:
                C.Clear()
                hist.Draw("hist")
                C.Print(TF["norm"]+".pdf")
            C.Print(TF["norm"]+".pdf]")

        
        DeltaEvsE.Write()
        DeltaEvsE_Norm.Write()
    
    outFile.Close()
    inFile.Close()

if __name__ == "__main__":
    options = parseArguments()
    TFset = [
                {
                    "histNames": ['b_TF'],
                    "base" : "b_TF",
                    "norm" : "b_TF_norm",
                    "title": "b transfer function"
                },
                {
                    "histNames": ['e_TF'],
                    "base" : "e_TF",
                    "norm" : "e_TF_norm",
                    "title": "e^{#pm} transfer function"
                },
                {
                    "histNames": ['mu_TF'],
                    "base" : "mu_TF",
                    "norm" : "mu_TF_norm",
                    "title": "#mu^{#pm} transfer function"
                },
            ]

    normAndSmooth(TFset, options.input, options.output, options.plot)
#    TFset = [
#        {
#            "histNames":
#                [
#                    "tf_bjet1_matchedToAfterFSR_allEta_hh_llmetjj_HWWleptons_nobtag_csv",
#                    "tf_bjet2_matchedToAfterFSR_allEta_hh_llmetjj_HWWleptons_nobtag_csv",
#                ],
#            "base": "ERecMinEGenVSEGen_bjet_matchedToAfterFSR_allEta_hh_llmetjj_HWWleptons_nobtag_csv",
#            "norm": "ERecMinEGenVSEGen_bjet_matchedToAfterFSR_allEta_Norm_hh_llmetjj_HWWleptons_nobtag_csv"
#        },
#        {
#            "histNames":
#                [
#                    "tf_bjet1_matchedToAfterFSR_barrel_hh_llmetjj_HWWleptons_nobtag_csv",
#                    "tf_bjet2_matchedToAfterFSR_barrel_hh_llmetjj_HWWleptons_nobtag_csv",
#                ],
#            "base": "ERecMinEGenVSEGen_bjet_matchedToAfterFSR_barrel_hh_llmetjj_HWWleptons_nobtag_csv",
#            "norm": "ERecMinEGenVSEGen_bjet_matchedToAfterFSR_barrel_Norm_hh_llmetjj_HWWleptons_nobtag_csv"
#        },
#        {
#            "histNames":
#                [
#                    "tf_bjet1_matchedToAfterFSR_endcap_hh_llmetjj_HWWleptons_nobtag_csv",
#                    "tf_bjet2_matchedToAfterFSR_endcap_hh_llmetjj_HWWleptons_nobtag_csv",
#                ],
#            "base": "ERecMinEGenVSEGen_bjet_matchedToAfterFSR_endcap_hh_llmetjj_HWWleptons_nobtag_csv",
#            "norm": "ERecMinEGenVSEGen_bjet_matchedToAfterFSR_endcap_Norm_hh_llmetjj_HWWleptons_nobtag_csv"
#        },
#        {
#            "histNames":
#                [
#                    "tf_bjet1_matchedToAfterFSR_allEta_hh_llmetjj_HWWleptons_btagM_csv",
#                    "tf_bjet2_matchedToAfterFSR_allEta_hh_llmetjj_HWWleptons_btagM_csv",
#                ],
#            "base": "ERecMinEGenVSEGen_bjet_matchedToAfterFSR_allEta_hh_llmetjj_HWWleptons_btagM_csv",
#            "norm": "ERecMinEGenVSEGen_bjet_matchedToAfterFSR_allEta_Norm_hh_llmetjj_HWWleptons_btagM_csv"
#        },
#        {
#            "histNames":
#                [
#                    "tf_bjet1_matchedToAfterFSR_barrel_hh_llmetjj_HWWleptons_btagM_csv",
#                    "tf_bjet2_matchedToAfterFSR_barrel_hh_llmetjj_HWWleptons_btagM_csv",
#                ],
#            "base": "ERecMinEGenVSEGen_bjet_matchedToAfterFSR_barrel_hh_llmetjj_HWWleptons_btagM_csv",
#            "norm": "ERecMinEGenVSEGen_bjet_matchedToAfterFSR_barrel_Norm_hh_llmetjj_HWWleptons_btagM_csv"
#        },
#        {
#            "histNames":
#                [
#                    "tf_bjet1_matchedToAfterFSR_endcap_hh_llmetjj_HWWleptons_btagM_csv",
#                    "tf_bjet2_matchedToAfterFSR_endcap_hh_llmetjj_HWWleptons_btagM_csv",
#                ],
#            "base": "ERecMinEGenVSEGen_bjet_matchedToAfterFSR_endcap_hh_llmetjj_HWWleptons_btagM_csv",
#            "norm": "ERecMinEGenVSEGen_bjet_matchedToAfterFSR_endcap_Norm_hh_llmetjj_HWWleptons_btagM_csv"
#        },
#        #{
#        #    "histNames":
#        #        [
#        #            "tf_bjet1_matchedToBeforeFSR_allEta",
#        #            "tf_bjet2_matchedToBeforeFSR_allEta",
#        #        ],
#        #    "base": "ERecMinEGenVSEGen_bjet_matchedToBeforeFSR_allEta",
#        #    "norm": "ERecMinEGenVSEGen_bjet_matchedToBeforeFSR_allEta_Norm"
#        #},
#        #{
#        #    "histNames":
#        #        [
#        #            "tf_bjet1_matchedToBeforeFSR_barrel",
#        #            "tf_bjet2_matchedToBeforeFSR_barrel",
#        #        ],
#        #    "base": "ERecMinEGenVSEGen_bjet_matchedToBeforeFSR_barrel",
#        #    "norm": "ERecMinEGenVSEGen_bjet_matchedToBeforeFSR_barrel_Norm"
#        #},
#        #{
#        #    "histNames":
#        #        [
#        #            "tf_bjet1_matchedToBeforeFSR_endcap",
#        #            "tf_bjet2_matchedToBeforeFSR_endcap",
#        #        ],
#        #    "base": "ERecMinEGenVSEGen_bjet_matchedToBeforeFSR_endcap",
#        #    "norm": "ERecMinEGenVSEGen_bjet_matchedToBeforeFSR_endcap_Norm"
#        #},
#        {
#            "histNames":
#                [
#                    "tf_mu1_matchedToAfterFSR_allEta_hh_llmetjj_HWWleptons_nobtag_csv",
#                    "tf_mu2_matchedToAfterFSR_allEta_hh_llmetjj_HWWleptons_nobtag_csv",
#                ],
#            "base": "ERecMinEGenVSEGen_mu_matchedToAfterFSR_allEta_hh_llmetjj_HWWleptons_nobtag_csv",
#            "norm": "ERecMinEGenVSEGen_mu_matchedToAfterFSR_allEta_Norm_hh_llmetjj_HWWleptons_nobtag_csv"
#        },
#        {
#            "histNames":
#                [
#                    "tf_mu1_matchedToAfterFSR_barrel_hh_llmetjj_HWWleptons_nobtag_csv",
#                    "tf_mu2_matchedToAfterFSR_barrel_hh_llmetjj_HWWleptons_nobtag_csv",
#                ],
#            "base": "ERecMinEGenVSEGen_mu_matchedToAfterFSR_barrel_hh_llmetjj_HWWleptons_nobtag_csv",
#            "norm": "ERecMinEGenVSEGen_mu_matchedToAfterFSR_barrel_Norm_hh_llmetjj_HWWleptons_nobtag_csv"
#        },
#        {
#            "histNames":
#                [
#                    "tf_mu1_matchedToAfterFSR_endcap_hh_llmetjj_HWWleptons_nobtag_csv",
#                    "tf_mu2_matchedToAfterFSR_endcap_hh_llmetjj_HWWleptons_nobtag_csv",
#                ],
#            "base": "ERecMinEGenVSEGen_mu_matchedToAfterFSR_endcap_hh_llmetjj_HWWleptons_nobtag_csv",
#            "norm": "ERecMinEGenVSEGen_mu_matchedToAfterFSR_endcap_Norm_hh_llmetjj_HWWleptons_nobtag_csv"
#        },
#        {
#            "histNames":
#                [
#                    "tf_mu1_matchedToAfterFSR_allEta_hh_llmetjj_HWWleptons_btagM_csv",
#                    "tf_mu2_matchedToAfterFSR_allEta_hh_llmetjj_HWWleptons_btagM_csv",
#                ],
#            "base": "ERecMinEGenVSEGen_mu_matchedToAfterFSR_allEta_hh_llmetjj_HWWleptons_btagM_csv",
#            "norm": "ERecMinEGenVSEGen_mu_matchedToAfterFSR_allEta_Norm_hh_llmetjj_HWWleptons_btagM_csv"
#        },
#        {
#            "histNames":
#                [
#                    "tf_mu1_matchedToAfterFSR_barrel_hh_llmetjj_HWWleptons_btagM_csv",
#                    "tf_mu2_matchedToAfterFSR_barrel_hh_llmetjj_HWWleptons_btagM_csv",
#                ],
#            "base": "ERecMinEGenVSEGen_mu_matchedToAfterFSR_barrel_hh_llmetjj_HWWleptons_btagM_csv",
#            "norm": "ERecMinEGenVSEGen_mu_matchedToAfterFSR_barrel_Norm_hh_llmetjj_HWWleptons_btagM_csv"
#        },
#        {
#            "histNames":
#                [
#                    "tf_mu1_matchedToAfterFSR_endcap_hh_llmetjj_HWWleptons_btagM_csv",
#                    "tf_mu2_matchedToAfterFSR_endcap_hh_llmetjj_HWWleptons_btagM_csv",
#                ],
#            "base": "ERecMinEGenVSEGen_mu_matchedToAfterFSR_endcap_hh_llmetjj_HWWleptons_btagM_csv",
#            "norm": "ERecMinEGenVSEGen_mu_matchedToAfterFSR_endcap_Norm_hh_llmetjj_HWWleptons_btagM_csv"
#        },
#        #{
#        #    "histNames":
#        #        [
#        #            "tf_mu1_matchedToBeforeFSR_allEta",
#        #            "tf_mu2_matchedToBeforeFSR_allEta",
#        #        ],
#        #    "base": "ERecMinEGenVSEGen_mu_matchedToBeforeFSR_allEta",
#        #    "norm": "ERecMinEGenVSEGen_mu_matchedToBeforeFSR_allEta_Norm"
#        #},
#        #{
#        #    "histNames":
#        #        [
#        #            "tf_mu1_matchedToBeforeFSR_barrel",
#        #            "tf_mu2_matchedToBeforeFSR_barrel",
#        #        ],
#        #    "base": "ERecMinEGenVSEGen_mu_matchedToBeforeFSR_barrel",
#        #    "norm": "ERecMinEGenVSEGen_mu_matchedToBeforeFSR_barrel_Norm"
#        #},
#        #{
#        #    "histNames":
#        #        [
#        #            "tf_mu1_matchedToBeforeFSR_endcap",
#        #            "tf_mu2_matchedToBeforeFSR_endcap",
#        #        ],
#        #    "base": "ERecMinEGenVSEGen_mu_matchedToBeforeFSR_endcap",
#        #    "norm": "ERecMinEGenVSEGen_mu_matchedToBeforeFSR_endcap_Norm"
#        #},
#        {
#            "histNames":
#                [
#                    "tf_el1_matchedToAfterFSR_allEta_hh_llmetjj_HWWleptons_nobtag_csv",
#                    "tf_el2_matchedToAfterFSR_allEta_hh_llmetjj_HWWleptons_nobtag_csv",
#                ],
#            "base": "ERecMinEGenVSEGen_el_matchedToAfterFSR_allEta_hh_llmetjj_HWWleptons_nobtag_csv",
#            "norm": "ERecMinEGenVSEGen_el_matchedToAfterFSR_allEta_Norm_hh_llmetjj_HWWleptons_nobtag_csv"
#        },
#        {
#            "histNames":
#                [
#                    "tf_el1_matchedToAfterFSR_barrel_hh_llmetjj_HWWleptons_nobtag_csv",
#                    "tf_el2_matchedToAfterFSR_barrel_hh_llmetjj_HWWleptons_nobtag_csv",
#                ],
#            "base": "ERecMinEGenVSEGen_el_matchedToAfterFSR_barrel_hh_llmetjj_HWWleptons_nobtag_csv",
#            "norm": "ERecMinEGenVSEGen_el_matchedToAfterFSR_barrel_Norm_hh_llmetjj_HWWleptons_nobtag_csv"
#        },
#        {
#            "histNames":
#                [
#                    "tf_el1_matchedToAfterFSR_endcap_hh_llmetjj_HWWleptons_nobtag_csv",
#                    "tf_el2_matchedToAfterFSR_endcap_hh_llmetjj_HWWleptons_nobtag_csv",
#                ],
#            "base": "ERecMinEGenVSEGen_el_matchedToAfterFSR_endcap_hh_llmetjj_HWWleptons_nobtag_csv",
#            "norm": "ERecMinEGenVSEGen_el_matchedToAfterFSR_endcap_Norm_hh_llmetjj_HWWleptons_nobtag_csv"
#        },
#        {
#            "histNames":
#                [
#                    "tf_el1_matchedToAfterFSR_allEta_hh_llmetjj_HWWleptons_btagM_csv",
#                    "tf_el2_matchedToAfterFSR_allEta_hh_llmetjj_HWWleptons_btagM_csv",
#                ],
#            "base": "ERecMinEGenVSEGen_el_matchedToAfterFSR_allEta_hh_llmetjj_HWWleptons_btagM_csv",
#            "norm": "ERecMinEGenVSEGen_el_matchedToAfterFSR_allEta_Norm_hh_llmetjj_HWWleptons_btagM_csv"
#        },
#        {
#            "histNames":
#                [
#                    "tf_el1_matchedToAfterFSR_barrel_hh_llmetjj_HWWleptons_btagM_csv",
#                    "tf_el2_matchedToAfterFSR_barrel_hh_llmetjj_HWWleptons_btagM_csv",
#                ],
#            "base": "ERecMinEGenVSEGen_el_matchedToAfterFSR_barrel_hh_llmetjj_HWWleptons_btagM_csv",
#            "norm": "ERecMinEGenVSEGen_el_matchedToAfterFSR_barrel_Norm_hh_llmetjj_HWWleptons_btagM_csv"
#        },
#        {
#            "histNames":
#                [
#                    "tf_el1_matchedToAfterFSR_endcap_hh_llmetjj_HWWleptons_btagM_csv",
#                    "tf_el2_matchedToAfterFSR_endcap_hh_llmetjj_HWWleptons_btagM_csv",
#                ],
#            "base": "ERecMinEGenVSEGen_el_matchedToAfterFSR_endcap_hh_llmetjj_HWWleptons_btagM_csv",
#            "norm": "ERecMinEGenVSEGen_el_matchedToAfterFSR_endcap_Norm_hh_llmetjj_HWWleptons_btagM_csv"
#        },
#        #{
#        #    "histNames":
#        #        [
#        #            "tf_el1_matchedToBeforeFSR_allEta",
#        #            "tf_el2_matchedToBeforeFSR_allEta",
#        #        ],
#        #    "base": "ERecMinEGenVSEGen_el_matchedToBeforeFSR_allEta",
#        #    "norm": "ERecMinEGenVSEGen_el_matchedToBeforeFSR_allEta_Norm"
#        #},
#        #{
#        #    "histNames":
#        #        [
#        #            "tf_el1_matchedToBeforeFSR_barrel",
#        #            "tf_el2_matchedToBeforeFSR_barrel",
#        #        ],
#        #    "base": "ERecMinEGenVSEGen_el_matchedToBeforeFSR_barrel",
#        #    "norm": "ERecMinEGenVSEGen_el_matchedToBeforeFSR_barrel_Norm"
#        #},
#        #{
#        #    "histNames":
#        #        [
#        #            "tf_el1_matchedToBeforeFSR_endcap",
#        #            "tf_el2_matchedToBeforeFSR_endcap",
#        #        ],
#        #    "base": "ERecMinEGenVSEGen_el_matchedToBeforeFSR_endcap",
#        #    "norm": "ERecMinEGenVSEGen_el_matchedToBeforeFSR_endcap_Norm"
#        #},
#    ]



