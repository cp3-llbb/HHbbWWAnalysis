import os
import sys
from copy import deepcopy,copy
from collections import OrderedDict 

from ROOT import TFile, TH1F, TCanvas, TLegend, gROOT, gStyle
    
gROOT.SetBatch(True)
gStyle.SetOptStat(0)

def compareSyncHist(dict_files):
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
        aMax = None
        aMin = None
        nBins = None
        for group,tree in dict_trees.items():
            print ("*",end='')
            # Draw arguments #
            if aMax is None and aMin is None and nBins is None:
                branch_str = branch+">>h"+group+branch+"(100)"
            else:
                branch_str = branch+">>h"+group+branch+"("+str(nBins)+","+str(aMin)+","+str(aMax)+")"
            # Draw #
            if group == "UCLA group":
                tree.Draw(branch_str,branch+"!=-10000")
            else:
                tree.Draw(branch_str,branch+"!=-9999")

            # Recover histo #
            hist = gROOT.FindObject("h"+group+branch)
            if aMax is None and aMin is None:
                aMax = hist.GetXaxis().GetBinUpEdge(hist.GetNbinsX())
                aMin = hist.GetXaxis().GetBinLowEdge(1)
                nBins = hist.GetNbinsX()
            branch_dict[group] = deepcopy(hist)
            del hist
        print()
        hist_dict[branch] = branch_dict

    print ("[INFO] Plotting the histograms")
    # Produce the PDF #
    C = TCanvas("C","C",800,600)
    outName = "compareSyncHist.pdf"
    C.Print(outName+"[")
    # Line attributes #
    colors = [1,634,419,601]
    style  = [1,2,3,4]
    width  = [1,2,2,2]
    # Loop over branches #
    for branch,gHist in hist_dict.items():
        print ("... Branch : %s"%branch,end=' ')
        C.Clear()
        # Find max value for all hist at given branch #
        hmax = 0
        hmin = 1e20
        for hist in gHist.values():
            if hist.GetMaximum()>hmax:
                hmax = hist.GetMaximum()
            if hist.GetMinimum()<hmin:
                hmin = hist.GetMinimum()
        # Loop over groups and hists #
        legend = TLegend(0.60,0.80,0.99,0.99)
        legend.SetHeader("Legend","C");
        for i,(group, hist) in enumerate(gHist.items()):
            print ("*",end='')
            hist.SetTitle(branch)
            hist.SetLineWidth(width[i])
            hist.SetLineColor(colors[i])
            hist.SetLineStyle(style[i])
            hist.SetMaximum(hmax*1.1)
            hist.SetMinimum(hmin*0.9)
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


compareSyncHist({
            #'Louvain group (HH sample)':'/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/Synchronization/results/GluGluToRadionToHHTo2B2VTo2L2Nu_M-750.root',
            #'Tallinn group v5 (HH sample)':'/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/utils/sync_bbww_Tallinn_2016_v5.root',
            #'Tallinn group v7 (HH sample)':'/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/utils/sync_bbww_Tallinn_2016_v7.root',
            #'Tallinn group v8 (HH sample)':'/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/utils/sync_bbww_Tallinn_2016_v8.root',
            #'Tallinn group v9 (HH sample)':'/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/utils/sync_bbww_Tallinn_2016_v9.root',
            #'UCLA group (HH sample)':'/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/myNanoProdMc2016_NANO_Friend_20191003_v4.root',
            'Louvain group (TT sample)':'/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/Synchronization/results/TTTo2L2Nu.root',
            'Tallinn group (TT sample)':'/home/users/f/b/fbury/bamboodev/HHbbWWAnalysis/utils/sync_bbww_ttbar_Tallinn_2016_v1.root',
            })




