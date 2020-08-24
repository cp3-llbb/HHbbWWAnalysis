import os
import sys
import glob
import json

from ROOT import TCanvas, gROOT, TGraphAsymmErrors
gROOT.SetBatch(True)

# sys.argv[1] can either be a json of a directory containg jsons

path = os.path.abspath(sys.argv[1])

# Get dict of path #
dict_files = {}
if os.path.isdir(path):
    for f in glob.glob(os.path.join(path,"*.json")):
         dict_files[os.path.basename(f).replace(".json","")] = f
elif os.path.isfile(path) and path.endswith(".json"):
    dict_files[os.path.basename(path).replace(".json","")] = path

C = TCanvas()
C.Print("PUWeights.pdf[")
# Loop over json content #
for name,path in dict_files.items():
    print ("Looking at %s"%name)
    # Load json content #
    with open(path,"r") as handle:
        data = json.load(handle)
    # Get binning content #
    bins = data['data'] # List of dict (1 dict per bin)
    g = TGraphAsymmErrors(len(bins))
    # Loop over bins #
    for i,b in enumerate(bins):
        val = b['value']
        up = b['error_high']
        down = b['error_low']
        bincenter = (b['bin'][0]+b['bin'][1])/2
        g.SetPoint(i,bincenter,val)
        g.SetPointError(i,0,0,up,down)

    g.GetHistogram().SetMinimum(0.)
    g.GetHistogram().SetMaximum(2.)
    g.SetTitle("Sample : %s"%name)
    g.GetXaxis().SetTitle("nTrueInt")
    g.GetYaxis().SetTitle("PU weight")
    C.Clear()
    g.Draw()
    C.Print("PUWeights.pdf","Title:%s"%name)

C.Print("PUWeights.pdf]")




