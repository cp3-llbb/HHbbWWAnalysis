import os
import sys
from copy import deepcopy,copy
from collections import OrderedDict 

from ROOT import TFile, TH1F, TCanvas, TLegend, gROOT, gStyle
    
gROOT.SetBatch(True)
gStyle.SetOptStat(0)

def compareSyncHist(dict_files,name):
    dict_trees = {}
    list_branches = {}
    rootfiles = []

    # Get all trees in dict #
    for k,(f,t) in dict_files.items():
        rootfile = TFile(f)
        tree = rootfile.Get(t)
        dict_trees[k] = tree
        list_branches[k] = [key.GetName() for key in tree.GetListOfBranches()]
        rootfiles.append(rootfile)

    # Find common branches #
    branches = list()
    for ki in list_branches:
        for kj in list_branches:
            if ki==kj:
                continue
            branches += [br for br in list_branches[ki] if br in list_branches[kj]]

    branches = sorted(list(set(branches)))

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
            if not branch in list_branches[group]:
                continue
            # Draw arguments #
            histname = ("h"+group+branch).replace('(','').replace(')','')
            if aMax is None and aMin is None and nBins is None:
                branch_str = branch+">>"+histname+"(100)"
            else:
                branch_str = branch+">>"+histname+"("+str(nBins)+","+str(aMin)+","+str(aMax)+")"
            # Draw #
            if group == "UCLA group":
                tree.Draw(branch_str,branch+"!=-10000")
            else:
                tree.Draw(branch_str,branch+"!=-9999")

            # Recover histo #
            hist = gROOT.FindObject(histname)
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
    outName = "compareSyncHist_%s.pdf"%name
    C.Print(outName+"[")
    # Line attributes #
    colors = [1,634,419,601]
    style  = [1,2,3,4]
    width  = [1,2,2,2]
    group_att = {gr:{'c':colors[i],'s':style[i],'w':width[i]} for i,gr in enumerate(dict_files.keys())}
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
            hist.SetLineWidth(group_att[group]['w'])
            hist.SetLineColor(group_att[group]['c'])
            hist.SetLineStyle(group_att[group]['s'])
            hist.SetMaximum(hmax*1.1)
            #hist.SetMinimum(hmin*0.9)
            hist.SetMinimum(0)
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
            'Louvain group' : ('/home/users/f/b/fbury/bamboodev/HHbbWWAnalysis/Sync/sync_bbww_Louvain_2016_v19.root','syncTree_hhbb2l_SR'),
            'Tallinn group' : ('/home/users/f/b/fbury/bamboodev/HHbbWWAnalysis/Sync/sync_bbww_Tallinn_2016_v18.root','syncTree_hhbb2l_SR'),
            'Aachen group' : ('/home/users/f/b/fbury/bamboodev/HHbbWWAnalysis/Sync/sync_bbww_aachen_2016.root','syncTree_hhbb2l_SR')},
             'SR')

compareSyncHist({
            'Louvain group' : ('/home/users/f/b/fbury/bamboodev/HHbbWWAnalysis/Sync/sync_bbww_Louvain_2016_v19.root','syncTree_hhbb2l_Fake'),
            'Tallinn group' : ('/home/users/f/b/fbury/bamboodev/HHbbWWAnalysis/Sync/sync_bbww_Tallinn_2016_v18.root','syncTree_hhbb2l_Fake'),
            'Aachen group' : ('/home/users/f/b/fbury/bamboodev/HHbbWWAnalysis/Sync/sync_bbww_aachen_2016.root','syncTree_hhbb2l_Fake')},
            'Fake')




