import ROOT
import sys

HLTList =  [   
            "HLT_IsoMu22",
            "HLT_IsoTkMu22",
            "HLT_IsoMu22_eta2p1",
            "HLT_IsoTkMu22_eta2p1",
            "HLT_IsoMu24",
            "HLT_IsoTkMu24",
            "HLT_Ele27_WPTight_Gsf",
            "HLT_Ele25_eta2p1_WPTight_Gsf",
            "HLT_Ele27_eta2p1_WPLoose_Gsf",
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
            "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL",
            "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",
            "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
            "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
            "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL",
            "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ",
            ]

def checkHLTpath(filepath,substr=None):
    f = ROOT.TFile(filepath,"READ")
    t = f.Get("Events")

    check = {HLT:False for HLT in HLTList}

    for HLT in HLTList:
                print()
        for branch in t.GetListOfBranches():
            if branch.GetName() == HLT:
                print ("Found HLT path : %s"%branch.GetName())
                check[HLT] = True
    print ()
    if substr is not None:
        for branch in t.GetListOfBranches():
            name = branch.GetName()
            if name.find(substr) != -1:
                print ("\t%s"%name)

    if not all(check.values()):
        print ("Missing HLT paths in file")
        for HLT,found in check.items():
            if not found:
                print('...',HLT)
    else:
        print ("All HLT paths have been found")

if __name__ == '__main__':
    if len(sys.argv)>2:
        checkHLTpath(sys.argv[1],sys.argv[2])
    else:
        checkHLTpath(sys.argv[1])
