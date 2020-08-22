import os
import sys
import ROOT

ROOT.gROOT.SetBatch(True)

#filename = os.path.join(os.path.abspath(os.path.dirname(__file__)),'datacards',sys.argv[1])
filename = sys.argv[1]
f = ROOT.TFile(filename)
#d = f.Get("HH_2l_0tau/data_obs")
#h = f.Get("HH_2l_0tau/TTZ")
#h.Add(f.Get("HH_2l_0tau/VH"))
#h.Add(f.Get("HH_2l_0tau/TT"))
#h.Add(f.Get("HH_2l_0tau/TTW"))
#h.Add(f.Get("HH_2l_0tau/WW"))
#h.Add(f.Get("HH_2l_0tau/TTH"))
#h.Add(f.Get("HH_2l_0tau/ZZ"))
#h.Add(f.Get("HH_2l_0tau/DY"))
#h.Add(f.Get("HH_2l_0tau/Other"))
#h.Add(f.Get("HH_2l_0tau/WZ"))

d = f.Get("data_obs")
h = f.Get("DY")
h.Add(f.Get("TT"))
h.Add(f.Get("WJets"))
h.Add(f.Get("VH"))
h.Add(f.Get("TH"))
h.Add(f.Get("TTH"))
h.Add(f.Get("ZZ"))
h.Add(f.Get("WZ"))
h.Add(f.Get("WW"))
h.Add(f.Get("WZ"))
h.Add(f.Get("TTZ"))
h.Add(f.Get("TTW"))
h.Add(f.Get("TTWW"))
h.Add(f.Get("Other"))
h.Add(f.Get("data_fake"))

C = ROOT.TCanvas()
d.SetMaximum(max(d.GetMaximum(),h.GetMaximum())*1.1)
d.SetLineColor(ROOT.kRed)
h.SetLineColor(ROOT.kBlue)
d.SetLineWidth(2)
h.SetLineWidth(2)
d.Draw("H")
h.Draw("H same")
C.Print(filename.replace(".root",".pdf"))

