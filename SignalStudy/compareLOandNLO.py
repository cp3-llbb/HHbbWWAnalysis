import os
import sys
from copy import deepcopy
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

#histName = "mHH"
histName = "cosThetaStar"

era = '2016'

C = ROOT.TCanvas()
#C.SetLogy()
#C.SetLogx()
C.Print(f"Comparison_{histName}_{era}.pdf[")

for channel in ['GluGluToHHTo2B2VTo2L2Nu','GluGluToHHTo2B2WToLNu2J','GluGluToHHTo2B2Tau']:
    for node in ['cHHH0','cHHH1','cHHH2p45','cHHH5']:
        path_NLO = f"/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SignalLOHistogramsCorr_{era}/results/{channel}_node_{node}.root".replace('2B2WToLNu2J','2B2VLNu2J')

        paths_LO = {
        'Benchmark SM' : f'/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SignalLOHistogramsCorr_{era}/results/{channel}_node_SM.root',
        'Benchmark 1'  : f'/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SignalLOHistogramsCorr_{era}/results/{channel}_node_10.root',
        'Benchmark 2'  : f'/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SignalLOHistogramsCorr_{era}/results/{channel}_node_11.root',
        'Benchmark 3'  : f'/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SignalLOHistogramsCorr_{era}/results/{channel}_node_12.root',
        'Benchmark 4'  : f'/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SignalLOHistogramsCorr_{era}/results/{channel}_node_1.root',
        'Benchmark 5'  : f'/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SignalLOHistogramsCorr_{era}/results/{channel}_node_2.root',
        'Benchmark 6'  : f'/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SignalLOHistogramsCorr_{era}/results/{channel}_node_4.root',
        'Benchmark 7'  : f'/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SignalLOHistogramsCorr_{era}/results/{channel}_node_5.root',
        'Benchmark 8'  : f'/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SignalLOHistogramsCorr_{era}/results/{channel}_node_6.root',
        'Benchmark 9'  : f'/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SignalLOHistogramsCorr_{era}/results/{channel}_node_7.root',
        'Benchmark 10' : f'/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SignalLOHistogramsCorr_{era}/results/{channel}_node_8.root',
        'Benchmark 11' : f'/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SignalLOHistogramsCorr_{era}/results/{channel}_node_9.root',
        'Benchmark 12' : f'/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SignalLOHistogramsCorr_{era}/results/{channel}_node_SM.root',
        }


        paths_reweightedLO = {
        'Benchmark SM' : f'/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SignalLOHistogramsCorr_{era}/results/{channel}_node_SMBenchmark{node}.root',
        'Benchmark 1'  : f'/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SignalLOHistogramsCorr_{era}/results/{channel}_node_10Benchmark{node}.root',
        'Benchmark 2'  : f'/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SignalLOHistogramsCorr_{era}/results/{channel}_node_11Benchmark{node}.root',
        'Benchmark 3'  : f'/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SignalLOHistogramsCorr_{era}/results/{channel}_node_12Benchmark{node}.root',
        'Benchmark 4'  : f'/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SignalLOHistogramsCorr_{era}/results/{channel}_node_1Benchmark{node}.root',
        'Benchmark 5'  : f'/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SignalLOHistogramsCorr_{era}/results/{channel}_node_2Benchmark{node}.root',
        'Benchmark 6'  : f'/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SignalLOHistogramsCorr_{era}/results/{channel}_node_4Benchmark{node}.root',
        'Benchmark 7'  : f'/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SignalLOHistogramsCorr_{era}/results/{channel}_node_5Benchmark{node}.root',
        'Benchmark 8'  : f'/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SignalLOHistogramsCorr_{era}/results/{channel}_node_6Benchmark{node}.root',
        'Benchmark 9'  : f'/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SignalLOHistogramsCorr_{era}/results/{channel}_node_7Benchmark{node}.root',
        'Benchmark 10' : f'/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SignalLOHistogramsCorr_{era}/results/{channel}_node_8Benchmark{node}.root',
        'Benchmark 11' : f'/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SignalLOHistogramsCorr_{era}/results/{channel}_node_9Benchmark{node}.root',
        'Benchmark 12' : f'/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SignalLOHistogramsCorr_{era}/results/{channel}_node_SMBenchmark{node}.root',
        }

        h_LO = {}
        for name,path in paths_LO.items():
            F_LO = ROOT.TFile(path)
            h_LO[name] = deepcopy(F_LO.Get(histName))
            F_LO.Close()

        F_NLO = ROOT.TFile(path_NLO)
        h_NLO = deepcopy(F_NLO.Get(histName))
        F_NLO.Close()

        h_r = {}
        for name, path in paths_reweightedLO.items():
            F = ROOT.TFile(path)
            h_r[name] = deepcopy(F.Get(histName))
            F.Close()


        min_colors = 51
        max_colors = 100  


# All LO together #
        C.Clear()

        legend = ROOT.TLegend(0.5,0.3,0.88,0.88)
        legend.SetFillStyle(4000)
        legend.SetBorderSize(0)

        h_NLO.Scale(1/h_NLO.Integral())

        h_NLO.SetLineWidth(2)
        h_NLO.SetLineStyle(1)
        h_NLO.SetLineColor(1)
        h_NLO.SetTitle(f"{channel} {era}")
        h_NLO.GetXaxis().SetTitle("m_{HH} [GeV]")
        h_NLO.GetYaxis().SetTitle("#sum w (normalized to unity)")
        h_NLO.SetMaximum(h_NLO.GetMaximum()*1.1)
        h_NLO.SetMinimum(0.)

        h_NLO.Draw("H")
        legend.AddEntry(h_NLO,f"NLO sample {node}")

        for i,(name, h) in enumerate(h_LO.items()):
            h.Scale(1/h.Integral())
            col = round(min_colors+i*(max_colors-min_colors)/len(h_LO))
            h.SetLineStyle(1)
            h.SetLineColor(col)
            legend.AddEntry(h,name + " LO")
            h.Draw("H same")

        legend.Draw()
        C.Print(f"Comparison_{histName}_{era}.pdf")


# NLO + reweighted LO #
        C.Clear()

        legend = ROOT.TLegend(0.5,0.3,0.88,0.88)
        legend.SetFillStyle(4000)
        legend.SetBorderSize(0)

        h_NLO.Draw("H")
        legend.AddEntry(h_NLO,f"NLO sample ({node})")

        for i,(name, h) in enumerate(h_r.items()):
            h.Scale(1/h.Integral())
            col = round(min_colors+i*(max_colors-min_colors)/len(h_r))
            h.SetLineStyle(3)
            h.SetLineColor(col)
            legend.AddEntry(h,name + f" LO#rightarrow NLO {node}")
            h.Draw("H same")
        legend.Draw()
        C.Print(f"Comparison_{histName}_{era}.pdf")


C.Print(f"Comparison_{histName}_{era}.pdf]")




