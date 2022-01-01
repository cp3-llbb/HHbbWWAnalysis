import os
import sys
import glob
import copy
import json 
import math
import argparse
from array import array
import numpy as np
import ROOT
from IPython import embed

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)

parser = argparse.ArgumentParser(description='Comparing DY')
parser.add_argument('--dir', action='store', required=True, type=str,
                    help='Path to shapes')
parser.add_argument('--suffix', action='store', required=True, type=str,
                    help='Suffix to file')
args = parser.parse_args()


def findEra(cat):
    for era in ['2016','2017','_2018']:
        if era in cat:
            return era
    raise ValueError(f'Category {cat}')

plotDict = {}

for f in glob.glob(os.path.join(args.dir,'*root')):
    category = os.path.basename(f).replace('.root','')
    era = findEra(category)
    if category not in plotDict.keys():
        plotDict[category] = {}

    F = ROOT.TFile(f)

    for key in F.GetListOfKeys():
        key = key.GetName()
        sample,shape = key.split('_UE_')
        if sample not in plotDict[category].keys():
            plotDict[category][sample] = {}
        plotDict[category][sample][shape] = copy.deepcopy(F.Get(key))

    F.Close()

results = {}
for cat in plotDict.keys():
    results[cat] = {}
    for sample in plotDict[cat].keys():
        results[cat][sample] = []
        h_nom = plotDict[cat][sample]['nom']
        h_down = plotDict[cat][sample]['down']
        h_up = plotDict[cat][sample]['up']

        assert h_nom.GetNbinsX() == h_up.GetNbinsX()
        assert h_nom.GetNbinsX() == h_down.GetNbinsX()

        for i in range(1,h_nom.GetNbinsX()+1):
            edges = [round(h_nom.GetXaxis().GetBinLowEdge(i),6),
                     round(h_nom.GetXaxis().GetBinUpEdge(i),6)]
            y_nom = h_nom.GetBinContent(i)
            y_down = h_down.GetBinContent(i)
            y_up = h_up.GetBinContent(i)
            results[cat][sample].append({'bin': edges,
                                         'up': y_up / y_nom if y_nom > 1e-5 else 1.,
                                         'down': y_down / y_nom if y_nom > 1e-5 else 1.})
            # Datacard script put min at 1e-5

with open(f'CompareUE/UE_{args.suffix}.json','w') as handle:
    json.dump(results,handle,indent=4)
print (f'Saved numerical factors to CompareUE/UE_{args.suffix}.json')
