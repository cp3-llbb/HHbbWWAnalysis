import ROOT
import sys

HLTList =  [   
              "HLT_IsoMu24",
              "HLT_Ele27_WPTight_Gsf",
              "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",
              "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
              "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
              "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
              "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ"
              "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
              "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL"
            ]

def checkHLTpath(filepath,substr=None):
    f = ROOT.TFile(filepath,"READ")
    t = f.Get("Events")

    for HLT in HLTList:
        for branch in t.GetListOfBranches():
            if branch.GetName().find(HLT) != -1:
                print ("Found HLT path : %s"%branch.GetName())
    print ()
    if substr is not None:
        for branch in t.GetListOfBranches():
            name = branch.GetName()
            if name.find(substr) != -1:
                print ("\t%s"%name)
        


if __name__ == '__main__':
    if len(sys.argv)>2:
        checkHLTpath(sys.argv[1],sys.argv[2])
    else:
        checkHLTpath(sys.argv[1])
