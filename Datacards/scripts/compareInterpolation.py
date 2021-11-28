import os
import sys
import copy
import ROOT
import argparse
from interpolation import Interpolation

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

p1 = 300
p2 = 400
p3 = 350
branch = 'LowMass'
#p1 = 600
#p2 = 700
#p3 = 650
#branch = 'HighMass'
dim = '2D'
path1 = f'/nfs/scratch/fynu/fbury/Datacards/bbww_dl/TestInterpolation/Exact/datacard_fit_Resonant_{branch}_Graviton_{dim}_noSyst_M_{p1}_testMorphing_2016/'
path2 = f'/nfs/scratch/fynu/fbury/Datacards/bbww_dl/TestInterpolation/Exact/datacard_fit_Resonant_{branch}_Graviton_{dim}_noSyst_M_{p2}_testMorphing_2016/'
path3 = f'/nfs/scratch/fynu/fbury/Datacards/bbww_dl/TestInterpolation/Exact/datacard_fit_Resonant_{branch}_Graviton_{dim}_noSyst_M_{p3}_testMorphing_2016/'
path4 = f'/nfs/scratch/fynu/fbury/Datacards/bbww_dl/TestInterpolation/Interpolation/datacard_fit_Resonant_{branch}_Graviton_{dim}_noSyst_M_{p3}_testMorphing_2016/'

files = {
    'HH resolved 1b': (f'HH_{p1}_resolved1b_GGF_2016.root',f'HH_{p2}_resolved1b_GGF_2016.root',f'HH_{p3}_resolved1b_GGF_2016.root',f'HH_{p3}_resolved1b_GGF_2016.root'),
    'HH resolved 2b': (f'HH_{p1}_resolved2b_GGF_2016.root',f'HH_{p2}_resolved2b_GGF_2016.root',f'HH_{p3}_resolved2b_GGF_2016.root',f'HH_{p3}_resolved2b_GGF_2016.root'),
    'HH boosted': (f'HH_{p1}_boosted_GGF_2016.root',f'HH_{p2}_boosted_GGF_2016.root',f'HH_{p3}_boosted_GGF_2016.root',f'HH_{p3}_boosted_GGF_2016.root')
}

histograms = [
    (f'ggHH_M_{p1}_hbbhwwdl',f'ggHH_M_{p1}_hbbhwwdl'),
    (f'ggHH_M_{p2}_hbbhwwdl',f'ggHH_M_{p2}_hbbhwwdl'),
    (f'ggHH_M_{p3}_hbbhwwdl',f'ggHH_M_{p3}_hbbhwwdl'),
    (f'ggHH_M_{p3}_hbbhwwdl',f'ggHH_M_{p3}_hbbhwwdl'),
]
legends = [
    f'Graviton M = {p1} GeV (exact)',
    f'Graviton M = {p2} GeV (exact)',
    f'Graviton M = {p3} GeV (exact)',
    f'Graviton M = {p3} GeV (interpolated)',
]

pathPdf = f'CompareInterpolation/compareInterpolation_{dim}_{p1}_{p2}_{p3}.pdf'
if dim == '1D':
    C = ROOT.TCanvas('C','C',800,600)
elif dim == '2D':
    C = ROOT.TCanvas('C','C',3000,600)
else:
    raise ValueError

C.Print(pathPdf+'[')

h1_tot = None
h2_tot = None
h3_tot = None

def getHistograms(rootFile,histNames):
    h = None
    F = ROOT.TFile(rootFile)
    for histName in histNames:
        if h is None:
            h = copy.deepcopy(F.Get(histName))
        else:
            h.Add(F.Get(histName))
    F.Close()
    return h

for title,(f1,f2,f3,f4) in files.items():
    h1 = getHistograms(os.path.join(path1,f1),histograms[0])
    h2 = getHistograms(os.path.join(path2,f2),histograms[1])
    h3 = getHistograms(os.path.join(path3,f3),histograms[2])
    h4 = getHistograms(os.path.join(path4,f4),histograms[3])

    if dim == '1D':
        C.Clear()
        C.SetCanvasSize(800,600)
        h1.SetTitle(title)
        h1.SetMaximum(max([h.GetMaximum() for h in [h1,h2,h3,h4]])*1.1)
        h1.SetLineColor(601)
        h2.SetLineColor(634)
        h3.SetLineColor(418)
        h4.SetLineColor(432)
        h4.SetLineStyle(2)
        h1.Draw('hist')
        h2.Draw('hist same')
        h3.Draw('hist same')
        h4.Draw('hist same')
        legend= ROOT.TLegend(0.25,0.60,0.65,0.85)
        for h,leg in zip([h1,h2,h3,h4],legends):
            legend.AddEntry(h,f"#splitline{{{leg}}}{{[yield = {h.Integral():6.2f}]}}")
        legend.Draw()
    if dim == '2D':
        # COLZ plots #
        C.SetCanvasSize(3000,600)
        C.Clear()
        C.Divide(4)
        C.cd(1)
        h1.SetTitle(f'{title} : {legends[0]}')
        h1.Draw("colz")
        C.cd(2)
        h2.SetTitle(f'{title} : {legends[1]}')
        h2.Draw("colz")
        C.cd(3)
        h3.SetTitle(f'{title} : {legends[2]}')
        h3.Draw("colz")
        C.cd(4)
        h4.SetTitle(f'{title} : {legends[3]}')
        h4.Draw("colz")

    C.Print(pathPdf)

    if dim == '2D':
        # projection plots #
        C.Clear()
        C.SetCanvasSize(1200,600)
        C.Divide(2)
        pad1 = C.cd(1)
        h1_x = h1.ProjectionX('h1_x')
        h2_x = h2.ProjectionX('h2_x')
        h3_x = h3.ProjectionX('h3_x')
        h4_x = h4.ProjectionX('h4_x')
        h1_x.SetTitle(f"Projection X")
        h1_x.SetMaximum(max([h.GetMaximum() for h in [h1_x,h2_x,h3_x,h4_x]])*1.1)
        h1_x.SetLineColor(601)
        h2_x.SetLineColor(634)
        h3_x.SetLineColor(418)
        h4_x.SetLineColor(432)
        h4_x.SetLineStyle(2)
        h1_x.Draw('hist')
        h2_x.Draw('hist same')
        h3_x.Draw('hist same')
        h4_x.Draw('hist same')
        legend_x = ROOT.TLegend(0.25,0.60,0.70,0.85)
        for h,leg in zip([h1_x,h2_x,h3_x,h4_x],legends):
            legend_x.AddEntry(h,f"#splitline{{{leg}}}{{[yield = {h.Integral():6.2f}]}}")
        legend_x.Draw()
        pad2 = C.cd(2)
        h1_y = h1.ProjectionY('h1_y')
        h2_y = h2.ProjectionY('h2_y')
        h3_y = h3.ProjectionY('h3_y')
        h4_y = h4.ProjectionY('h4_y')
        h1_y.SetTitle("Projection Y")
        h1_y.SetMaximum(max([h.GetMaximum() for h in [h1_y,h2_y,h3_y,h4_y]])*1.1)
        h1_y.SetLineColor(601)
        h2_y.SetLineColor(634)
        h3_y.SetLineColor(418)
        h4_y.SetLineColor(432)
        h4_y.SetLineStyle(2)
        h1_y.Draw('hist')
        h2_y.Draw('hist same')
        h3_y.Draw('hist same')
        h4_y.Draw('hist same')
        legend_y = ROOT.TLegend(0.40,0.60,0.85,0.85)
        for h,leg in zip([h1_y,h2_y,h3_y,h4_y],legends):
            legend_y.AddEntry(h,f"#splitline{{{leg}}}{{[yield = {h.Integral():6.2f}]}}")
        legend_y.Draw()

        C.Print(pathPdf)

    # Add total histogram #
    if h1_tot is None:
        h1_tot = copy.deepcopy(h1)
        h2_tot = copy.deepcopy(h2)
        h3_tot = copy.deepcopy(h3)
    else:
        h1_tot.Add(h1)
        h2_tot.Add(h2)
        h3_tot.Add(h3)

def returnContours(h):
    # https://root.cern/doc/master/classTHistPainter.html#HP16a
    h.Draw("CONT LIST")
    ROOT.gPad.Update()
    contours = ROOT.gROOT.GetListOfSpecials().FindObject("contours")
    graphs = []
    for i in range(contours.GetSize()):
        g = contours.At(i).First()
        if g:
            g.SetFillColor(0)
            graphs.append(g)
    return copy.deepcopy(graphs)

interpolate = Interpolation(p1,p2,p3)
h4_tot = interpolate(h1_tot,h2_tot,'test')
h1s_tot = interpolate.h1s
h2s_tot = interpolate.h2s

from IPython import embed
#embed()


for h1,h2,l1,l2,c1,c2 in zip([h1_tot,h1s_tot,h3_tot],\
                       [h2_tot,h2s_tot,h4_tot],\
                       [f"M_{{HH}} = {p1} GeV",f"M_{{HH}} = {p1} GeV (shifted at {p3} GeV)",f"M_{{HH}} = {p3} GeV (exact)"],\
                       [f"M_{{HH}} = {p2} GeV",f"M_{{HH}} = {p2} GeV (shifted at {p3} GeV)",f"M_{{HH}} = {p3} GeV (interpolated)"],\
                       [ROOT.kBlue+1,ROOT.kBlue+1,ROOT.kGreen+1],\
                       [ROOT.kRed+1,ROOT.kRed+1,ROOT.kMagenta+1]):
    h1.RebinX(20)
    h1.RebinY(20)
    h2.RebinX(20)
    h2.RebinY(20)
    h1.SetContour(5)
    h2.SetContour(5)

    g1s = returnContours(h1)
    g2s = returnContours(h2)

    C.Clear()
    C.SetCanvasSize(800,600)
    C.SetLeftMargin(0.10)
    C.SetRightMargin(0.05)
    C.SetBottomMargin(0.10)
    C.SetTopMargin(0.05)


    h_dummy = ROOT.TH2D('h_dummy','h_dummy',100,0.,1.,100,0.,600.)
    h_dummy.SetTitle(";DNN score;HME")
    h_dummy.Draw()
    leg = ROOT.TLegend(0.15,0.75,0.8,.93)
    leg.AddEntry(g1s[0],l1)
    leg.AddEntry(g2s[0],l2)
    for g1 in g1s:
        g1.SetLineColor(c1)
        g1.SetLineWidth(2)
        g1.Draw("same")
    for g2 in g2s:
        g2.SetLineColor(c2)
        g2.SetLineWidth(2)
        g2.Draw("same")
    leg.Draw()
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
#leg.SetTitleSize(0.02)
    C.Print(pathPdf)


C.Print(pathPdf+']')

