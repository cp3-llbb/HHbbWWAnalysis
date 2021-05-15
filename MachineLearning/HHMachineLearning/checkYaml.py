import os
import sys
import yaml
import ROOT

with open(sys.argv[-1],"r") as handle:
    d = yaml.load(handle)

sampleDir = d['sampleDir']

for era in d['sampleDict'].keys():
    print ("Checking era",era)
    for cat,sample_list in d['sampleDict'][era].items():
        print ("Checking cat",cat)
        for sample in sample_list:
            path = os.path.join(sampleDir,sample)
            if not os.path.exists(path):
                print("Could not find :",path)
            F = ROOT.TFile(path)
            if not F:
                print ("Could not open root file :",path) 
            F.Close()
                
