import ROOT
import numpy as np

f= ROOT.TFile("weight_MuMu_data_2016.root")
h = f.Get("Weight")

xbins = []
max_rebin = 5
current_rebin = 0
cont = []
err = []
for i in range(1,h.GetNbinsX()+1):
    cont.append(h.GetBinContent(i))
    err.append(h.GetBinError(i))
    mean_cont = sum(cont)/len(cont)
    mean_err = sum(err)/len(err)
    print ('i',i,'len',len(cont))
    print (err,cont)
    if mean_cont == 0.:
        xbins.append(h.GetBinLowEdge(i))
        cont = []
        err = []
        continue
    print ('\t',mean_err,mean_cont,mean_err/mean_cont)
    if mean_err/mean_cont <=0.1:
        xbins.append(h.GetBinLowEdge(i))
        cont = []
        err = []
        continue
    else:
        current_rebin += 1
    print (current_rebin,max_rebin)
    if current_rebin == max_rebin:
        xbins.append(h.GetBinLowEdge(i))
        current_rebin = 0
        cont = []
        err = []
        
xbins.append(h.GetBinLowEdge(h.GetNbinsX()+1))

h_rebin = h.Rebin(len(xbins)-1,'rebin',np.array(xbins))
h_rebin.SetLineColor(634)
h_rebin.SetLineWidth(1)

base_width = h.GetBinLowEdge(i+1)-h.GetBinLowEdge(i)
for i in range(1,h_rebin.GetNbinsX()+1):
    width = h_rebin.GetBinLowEdge(i+1)-h_rebin.GetBinLowEdge(i)
    h_rebin.SetBinContent(i,h_rebin.GetBinContent(i)*base_width/width)
    h_rebin.SetBinError(i,h_rebin.GetBinError(i)*base_width/width)

print ('test')
for i in range(1,h_rebin.GetNbinsX()+1):
    if h_rebin.GetBinContent(i)>0:
        print (i,h_rebin.GetBinContent(i),h_rebin.GetBinError(i)/h_rebin.GetBinContent(i),h_rebin.GetBinLowEdge(i),h_rebin.GetBinLowEdge(i+1))

C = ROOT.TCanvas('C','C',1600,1200)
h.Draw()
h_rebin.Draw("same")

C.Print("test_rebin.pdf")
#input("wait")



