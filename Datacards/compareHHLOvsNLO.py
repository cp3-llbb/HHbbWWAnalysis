import os
import sys
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

classic_bbww = {
    '2016' : ROOT.TFile("/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/Datacards/datacards/datacard_fit_classic_SM_2016/plotit/root/ggHH_BM_SM_hbbhwwdl.root"),
    '2017' : ROOT.TFile("/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/Datacards/datacards/datacard_fit_classic_SM_2017/plotit/root/ggHH_BM_SM_hbbhwwdl.root"),
    '2018' : ROOT.TFile("/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/Datacards/datacards/datacard_fit_classic_SM_2018/plotit/root/ggHH_BM_SM_hbbhwwdl.root"),
}

classic_bbtt = {
    '2016' : ROOT.TFile("/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/Datacards/datacards/datacard_fit_classic_SM_2016/plotit/root/ggHH_BM_SM_hbbhtt.root"),
    '2017' : ROOT.TFile("/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/Datacards/datacards/datacard_fit_classic_SM_2017/plotit/root/ggHH_BM_SM_hbbhtt.root"),
    '2018' : ROOT.TFile("/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/Datacards/datacards/datacard_fit_classic_SM_2018/plotit/root/ggHH_BM_SM_hbbhtt.root"),
}

LO_bbww = {
    '2016' : ROOT.TFile("/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/Datacards/datacards/datacard_fit_LO_SM_2016/plotit/root/ggHH_BM_SM_hbbhwwdl.root"),
    '2017' : ROOT.TFile("/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/Datacards/datacards/datacard_fit_LO_SM_2017/plotit/root/ggHH_BM_SM_hbbhwwdl.root"),
    '2018' : ROOT.TFile("/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/Datacards/datacards/datacard_fit_LO_SM_2018/plotit/root/ggHH_BM_SM_hbbhwwdl.root"),
}
LO_bbtt = {
    '2016' : ROOT.TFile("/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/Datacards/datacards/datacard_fit_LO_SM_2016/plotit/root/ggHH_BM_SM_hbbhtt.root"),
    '2017' : ROOT.TFile("/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/Datacards/datacards/datacard_fit_LO_SM_2017/plotit/root/ggHH_BM_SM_hbbhtt.root"),
    '2018' : ROOT.TFile("/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/Datacards/datacards/datacard_fit_LO_SM_2018/plotit/root/ggHH_BM_SM_hbbhtt.root"),
}

NLO_bbww = {
    '2016' : ROOT.TFile("/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/Datacards/datacards/datacard_fit_NLO_SM_2016/plotit/root/ggHH_BM_SM_hbbhwwdl.root"),
    '2017' : ROOT.TFile("/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/Datacards/datacards/datacard_fit_NLO_SM_2017/plotit/root/ggHH_BM_SM_hbbhwwdl.root"),
    '2018' : ROOT.TFile("/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/Datacards/datacards/datacard_fit_NLO_SM_2018/plotit/root/ggHH_BM_SM_hbbhwwdl.root"),
}
NLO_bbtt = {
    '2016' : ROOT.TFile("/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/Datacards/datacards/datacard_fit_NLO_SM_2016/plotit/root/ggHH_BM_SM_hbbhtt.root"),
    '2017' : ROOT.TFile("/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/Datacards/datacards/datacard_fit_NLO_SM_2017/plotit/root/ggHH_BM_SM_hbbhtt.root"),
    '2018' : ROOT.TFile("/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/Datacards/datacards/datacard_fit_NLO_SM_2018/plotit/root/ggHH_BM_SM_hbbhtt.root"),
}

plots = [
 'HH_resolved1b_GGF',
 'HH_resolved1b_VBF',
 'HH_resolved1b_H',
 'HH_resolved2b_GGF',
 'HH_resolved2b_VBF',
 'HH_resolved2b_H',
 'HH_boosted_GGF',
 'HH_boosted_VBF',
 'HH_boosted_H',
 'HH_resolved_DY_VVV',
 'HH_boosted_DY_VVV',
 'HH_inclusive_other',
]

yields = {
    '2016' : [0.,0.,0.,0.,0.,0.],
    '2017' : [0.,0.,0.,0.,0.,0.],
    '2018' : [0.,0.,0.,0.,0.,0.],
}

for era in ['2016','2017','2018']:
    C = ROOT.TCanvas("C","C",800,600)
    pdfName = f"SM_LOvsNLO_{era}.pdf"
    C.Print(pdfName+'[')
    for plot in plots:
        leg = ROOT.TLegend(0.6,0.6,0.85,0.85)
        h_LO_bbww = LO_bbww[era].Get(plot)
        h_LO_bbtt = LO_bbtt[era].Get(plot)
        h_NLO_bbww = NLO_bbww[era].Get(plot)
        h_NLO_bbtt = NLO_bbtt[era].Get(plot)
        h_classic_bbww = classic_bbww[era].Get(plot)
        h_classic_bbtt = classic_bbtt[era].Get(plot)

        h_LO_bbww.Rebin(8)
        h_LO_bbtt.Rebin(8)
        h_NLO_bbww.Rebin(8)
        h_NLO_bbtt.Rebin(8)
        h_classic_bbww.Rebin(8)
        h_classic_bbtt.Rebin(8)

        yields[era][0] += h_LO_bbww.Integral()
        yields[era][1] += h_LO_bbtt.Integral()
        yields[era][2] += h_NLO_bbww.Integral()
        yields[era][3] += h_NLO_bbtt.Integral()
        yields[era][4] += h_classic_bbww.Integral()
        yields[era][5] += h_classic_bbtt.Integral()

        leg.AddEntry(h_LO_bbww,"all LO -> SM bbww")
        leg.AddEntry(h_LO_bbtt,"all LO -> SM bbtt")
        leg.AddEntry(h_NLO_bbww,"all NLO -> SM bbww")
        leg.AddEntry(h_NLO_bbtt,"all NLO -> SM bbtt")
        leg.AddEntry(h_classic_bbww,"NLO SM bbww")
        leg.AddEntry(h_classic_bbtt,"NLO SM bbtt")

        h_LO_bbww.SetLineWidth(2)
        h_LO_bbtt.SetLineWidth(2)
        h_NLO_bbww.SetLineWidth(2)
        h_NLO_bbtt.SetLineWidth(2)
        h_classic_bbww.SetLineWidth(2)
        h_classic_bbtt.SetLineWidth(2)

        h_LO_bbww.SetLineColor(ROOT.kRed+2)
        h_LO_bbtt.SetLineColor(ROOT.kGreen+2)
        h_NLO_bbww.SetLineColor(ROOT.kRed+2)
        h_NLO_bbtt.SetLineColor(ROOT.kGreen+2)
        h_classic_bbww.SetLineColor(ROOT.kRed+2)
        h_classic_bbtt.SetLineColor(ROOT.kGreen+2)

        h_LO_bbww.SetLineStyle(1)
        h_LO_bbtt.SetLineStyle(1)
        h_NLO_bbww.SetLineStyle(2)
        h_NLO_bbtt.SetLineStyle(2)
        h_classic_bbww.SetLineStyle(3)
        h_classic_bbtt.SetLineStyle(3)

        h_LO_bbww.SetMaximum(max([h.GetMaximum() for h in [h_LO_bbww,h_LO_bbtt,h_NLO_bbww,h_NLO_bbtt,h_classic_bbww,h_classic_bbtt]])*1.1)
        h_LO_bbww.SetTitle(plot)
        h_LO_bbww.SetMinimum(0.)
        h_LO_bbww.Draw("hist")
        h_LO_bbtt.Draw("hist same")
        h_NLO_bbww.Draw("hist same")
        h_NLO_bbtt.Draw("hist same")
        h_classic_bbww.Draw("hist same")
        h_classic_bbtt.Draw("hist same")
        leg.Draw()
        
        C.Print(pdfName)
        
    print (f'Era {era}') 
    for y,n in zip(yields[era],['all LO -> SM bbww','all LO -> SM bbtt','all NLO -> SM bbww','all NLO -> SM bbtt','NLO SM bbww','NLO SM bbtt']):
        print (f'{n:20s} -> {y:3.2f}')
     
    C.Print(pdfName+']')

