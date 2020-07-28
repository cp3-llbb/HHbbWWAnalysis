import os
import sys
import ROOT

ROOT.gROOT.SetBatch(True)

#filename = os.path.join(os.path.abspath(os.path.dirname(__file__)),'datacards',sys.argv[1])
filename = sys.argv[1]
f = ROOT.TFile(filename)
d = f.Get("HH_2l_0tau/data_obs")
h = f.Get("HH_2l_0tau/TTZ")
h.Add(f.Get("HH_2l_0tau/VH"))
h.Add(f.Get("HH_2l_0tau/TT"))
h.Add(f.Get("HH_2l_0tau/TTW"))
h.Add(f.Get("HH_2l_0tau/WW"))
h.Add(f.Get("HH_2l_0tau/TTH"))
h.Add(f.Get("HH_2l_0tau/ZZ"))
h.Add(f.Get("HH_2l_0tau/DY"))
h.Add(f.Get("HH_2l_0tau/Other"))
h.Add(f.Get("HH_2l_0tau/WZ"))

C = ROOT.TCanvas()
d.SetMaximum(max(d.GetMaximum(),h.GetMaximum())*1.1)
d.SetLineColor(ROOT.kRed)
h.SetLineColor(ROOT.kBlue)
d.SetLineWidth(2)
h.SetLineWidth(2)
d.Draw("H")
h.Draw("H same")
C.Print(filename.replace(".root",".pdf"))

