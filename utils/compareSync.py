import os
import sys
from copy import deepcopy,copy
from collections import OrderedDict 

from ROOT import TFile, TH1F, TCanvas, TLegend, gROOT, gStyle
    
gROOT.SetBatch(True)
gStyle.SetOptStat(0)

def compareSync(dict_files):
    dict_trees = {}
    list_branches = []
    rootfiles = []

    # Get all trees in dict #
    for k,f in dict_files.items():
        rootfile = TFile(f)
        tree = rootfile.Get("syncTree")
        dict_trees[k] = tree
        list_branches.append([key.GetName() for key in tree.GetListOfBranches()])
        rootfiles.append(rootfile)

    # Find common branches #
    branches = list(set(list_branches[0]).intersection(*list_branches[1:]))
    branches.sort()

    # Loop to get histograms in dict # 
    hist_dict = OrderedDict()

    print ("[INFO] Obtaining the histograms")
    for branch in branches:
        print ("... Branch : %s"%branch,end=' ')
        branch_dict = {}
        for group,tree in dict_trees.items():
            print ("*",end='')
            if group == "UCLA group":
                tree.Draw(branch+">>h"+group+branch+"(100)",branch+"!=-10000","norm")
            else:
                tree.Draw(branch+">>h"+group+branch+"(100)",branch+"!=-9999","norm")
            hist = gROOT.FindObject("h"+group+branch)
            branch_dict[group] = deepcopy(hist)
            del hist
        print()
        hist_dict[branch] = branch_dict


    print ("[INFO] Plotting the histograms")
    # Produce the PDF #
    C = TCanvas("C","C",800,600)
    outName = "compareSync.pdf"
    C.Print(outName+"[")
    # Loop over branches #
    for branch,gHist in hist_dict.items():
        print ("... Branch : %s"%branch,end=' ')
        C.Clear()
        # Find max value for all hist at given branch #
        hmax = 0
        for hist in gHist.values():
            if hist.GetMaximum()>hmax:
                hmax = hist.GetMaximum()
        # Loop over groups and hists #
        legend = TLegend(0.60,0.80,0.99,0.99)
        legend.SetHeader("Legend","C");
        for i,(group, hist) in enumerate(gHist.items()):
            print ("*",end='')
            hist.SetTitle(branch)
            hist.SetLineWidth(2)
            hist.SetLineColor(i+1)
            hist.SetLineStyle(i+1)
            hist.SetMaximum(hmax*1.1)
            legend.AddEntry(hist,group+" : %d entries"%hist.GetEntries(),"l")
            if i==0:
                hist.Draw("hist")
            else:
                hist.Draw("hist same")
        legend.Draw()
        C.Print(outName,'Title:'+branch)
        print()
    
    C.Print(outName+"]")

    for rootfile in rootfiles:
        rootfile.Close()


compareSync({'Louvain group':'/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/Synchronization/results/GluGluToRadionToHHTo2B2VTo2L2Nu_M-750.root',
             'Talinn group':'/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/sync_bbww_Tallinn_2016_v5.root',
             #'UCLA group':'/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/myNanoProdMc2016_NANO_Friend_20191003_v4.root',
            })




