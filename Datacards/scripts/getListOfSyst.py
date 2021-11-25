import os
import glob
import sys
import ROOT
import enlighten


systematics = []
list_dir = sys.argv[1:]
for directory in list_dir:
    print ("Looking at {}".format(directory))
    files = glob.glob(os.path.join(directory,"*.root"))
    pbar = enlighten.Counter(total=len(files), desc='Progress', unit='files')
    for f in files:
        pbar.update()
        F = ROOT.TFile(f)
        for key in F.GetListOfKeys():
            if not 'TH1' in key.GetClassName():
                continue
            if not '__' in key.GetName():
                continue
            syst = key.GetName().split('__')[1].replace('up','').replace('down','')
            if syst not in systematics:
                systematics.append(syst)
        F.Close()

print ("List of systematics is below :")
for syst in sorted(systematics):
    print ("  - {}".format(syst))
