import os
import sys
import copy

from ROOT import TFile

def saveToFile(f,l):
    with open(f, 'w') as handle:
        handle.write("event,Louvain,Tallinn\n")
        for item in l:
            handle.write("%s,%s,%s\n"%(item[0],item[1],item[2]))

def compareSyncTree(file1,file2):
    # Get list of trees #
    f1 = TFile(file1)
    t1 = f1.Get("syncTree")
    f2 = TFile(file2)
    t2 = f2.Get("syncTree")

    # Clean files #
    try:
        os.remove("ele1_isfakeablesel_diff.txt")
    except:
        pass
    try:
        os.remove("mu1_isfakeable_diff.txt")
    except:
        pass
    try:
        os.remove("ele1_ismvasel_diff.txt")
    except:
        pass
    try:
        os.remove("mu1_ismvasel_diff.txt")
    except:
        pass

    ele1_ismvasel_diff = []
    mu1_ismvasel_diff = []
    ele1_isfakeablesel_diff = []
    mu1_isfakeablesel_diff = []

    for e1 in t1:
        print ("Event ",e1.event)
        for e2 in t2:
            if e2.event > e1.event+1000:
                break
            if e1.event == e2.event:
                # Fakeable #
                if e1.ele1_isfakeablesel != -9999 and e2.ele1_isfakeablesel != -9999 and e1.ele1_isfakeablesel != e2.ele1_isfakeablesel:
                    ele1_isfakeablesel_diff.append([e1.event,e1.ele1_isfakeablesel,e2.ele1_isfakeablesel])
                    saveToFile("ele1_isfakeablesel_diff.txt",ele1_isfakeablesel_diff)
                    print ("\tele1_isfakeablesel_diff")
                if e1.mu1_isfakeablesel != -9999 and e2.mu1_isfakeablesel != -9999 and e1.mu1_isfakeablesel != e2.mu1_isfakeablesel:
                    mu1_isfakeablesel_diff.append([e1.event,e1.mu1_isfakeablesel,e2.mu1_isfakeablesel])
                    saveToFile("mu1_isfakeable_diff.txt",mu1_isfakeablesel_diff)
                    print ("\tmu1_isfakeablesel_diff")

                # mvasel #
                if e1.ele1_ismvasel != -9999 and e2.ele1_ismvasel != -9999 and e1.ele1_ismvasel != e2.ele1_ismvasel:
                    ele1_ismvasel_diff.append([e1.event,e1.ele1_ismvasel,e2.ele1_ismvasel])
                    saveToFile("ele1_ismvasel_diff.txt",ele1_ismvasel_diff)
                    print ("\tele1_ismvasel_diff")
                if e1.mu1_ismvasel != -9999 and e2.mu1_ismvasel != -9999 and e1.mu1_ismvasel != e2.mu1_ismvasel:
                    mu1_ismvasel_diff.append([e1.event,e1.mu1_ismvasel,e2.mu1_ismvasel])
                    saveToFile("mu1_ismvasel_diff.txt",mu1_ismvasel_diff)
                    print ("\tmu1_ismvasel_diff")

# Execution #
compareSyncTree(file1 = '/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/Synchronization/results/GluGluToRadionToHHTo2B2VTo2L2Nu_M-750.root',
                file2 = '/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/sync_bbww_Tallinn_2016_v7.root')




