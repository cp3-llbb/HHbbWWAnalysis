import os
import yaml
from copy import deepcopy
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

from IPython import embed

baseDir = '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/'
era = '2016'
paths_LO = [
    f'full{era}NanoV7_TTHIDLoose_Fake_DNN11_Boosted1B_GGFLO_BenchmarkcHHH0',
    f'full{era}NanoV7_TTHIDLoose_Fake_DNN11_Boosted1B_GGFLO_BenchmarkcHHH1',
    f'full{era}NanoV7_TTHIDLoose_Fake_DNN11_Boosted1B_GGFLO_BenchmarkcHHH2p45',
    f'full{era}NanoV7_TTHIDLoose_Fake_DNN11_Boosted1B_GGFLO_BenchmarkcHHH5',
    f'full{era}NanoV7_TTHIDLoose_Fake_DNN11_Resolved1B_GGFLO_BenchmarkcHHH0',
    f'full{era}NanoV7_TTHIDLoose_Fake_DNN11_Resolved1B_GGFLO_BenchmarkcHHH1',
    f'full{era}NanoV7_TTHIDLoose_Fake_DNN11_Resolved1B_GGFLO_BenchmarkcHHH2p45',
    f'full{era}NanoV7_TTHIDLoose_Fake_DNN11_Resolved1B_GGFLO_BenchmarkcHHH5',
    f'full{era}NanoV7_TTHIDLoose_Fake_DNN11_Resolved2B_GGFLO_BenchmarkcHHH0',
    f'full{era}NanoV7_TTHIDLoose_Fake_DNN11_Resolved2B_GGFLO_BenchmarkcHHH1',
    f'full{era}NanoV7_TTHIDLoose_Fake_DNN11_Resolved2B_GGFLO_BenchmarkcHHH2p45',
    f'full{era}NanoV7_TTHIDLoose_Fake_DNN11_Resolved2B_GGFLO_BenchmarkcHHH5',
]

paths_NLO = [
    f'full{era}NanoV7_TTHIDLoose_Fake_DNN11_Resolved1B_GGFNLO',
    f'full{era}NanoV7_TTHIDLoose_Fake_DNN11_Resolved2B_GGFNLO',
    f'full{era}NanoV7_TTHIDLoose_Fake_DNN11_Boosted1B_GGFNLO',
]

def RunDirectory(path):
    # Load yaml #
    with open(os.path.join(path,'plots.yml')) as handle:
        d = yaml.load(handle)
    
    plotsPerSamples = {}
    lumi = d['configuration']['luminosity'][era]
    for sampleName,sampleCfg in d['files'].items():
        F = ROOT.TFile(os.path.join(path,'results',sampleName))
        plots = {}
        factor = lumi
        #if 'cross-section' in sampleCfg.keys(): # Xsec = 1
        #    factor *= sampleCfg['cross-section']
        if 'branching-ratio' in sampleCfg.keys():
            factor *= sampleCfg['branching-ratio']
        if 'generated-events' in sampleCfg.keys():
            factor /= sampleCfg['generated-events']

        listKeys = [key.GetName() for key in F.GetListOfKeys()]
        for plotName,plotCfg in d['plots'].items():
            if plotName not in listKeys:
                continue
            h = F.Get(plotName) 
            h.Scale(factor)
            plots[plotName] = deepcopy(h)
            plots[plotName].GetXaxis().SetTitle(plotCfg['x-axis'])
            channel = plotCfg['labels'][0]['text']
            if "ResolvedOneBtag" in plotName:
                cat = "Resolved 1b"
            elif "ResolvedTwoBtags" in plotName:
                cat = "Resolved 2b"
            elif "Boosted" in plotName:
                cat = "Boosted 1b"
            else:
                cat = ""
            plots[plotName].GetYaxis().SetTitle('yield')
            plots[plotName].SetTitle(f"{cat} ({channel} channel)")

        F.Close()

        plotsPerSamples[sampleName.replace('.root','')] = plots
    return plotsPerSamples

plots_LO = {}
for path_LO in paths_LO:
    print (f'Looking at {path_LO}')
    plotsPerSamples = RunDirectory(os.path.join(baseDir,path_LO))
    for sample,plots in plotsPerSamples.items():
        if sample in plots_LO.keys():
            plots_LO[sample].update(plots)
        else:
            plots_LO[sample] = plots

plots_NLO = {}
for path_NLO in paths_NLO:
    print (f'Looking at {path_NLO}')
    plotsPerSamples = RunDirectory(os.path.join(baseDir,path_NLO))
    for sample,plots in plotsPerSamples.items():
        if sample in plots_NLO.keys():
            plots_NLO[sample].update(plots)
        else:
            plots_NLO[sample] = plots

def findKey(aDict,parts):
    for k in aDict.keys():
        if all([part in k for part in parts]):
            return k
    return None
          
LONames = list(plots_LO.keys())
NLONames = list(plots_NLO.keys())
path_plots = "ValidationPlots/"
if not os.path.exists(path_plots):
    os.makedirs(path_plots)

for channel in ['GluGluToHHTo2B2Tau','GluGluToHHTo2B2VTo2L2Nu']:
    for benchmark in ["cHHH0","cHHH1","cHHH2p45","cHHH5"]:
        text = ROOT.TPaveText(0.6,0.5,0.88,0.6,"NDC")
        text.SetBorderSize(0)
        text.SetFillStyle(4000)
        text.AddText(f'Channel : {channel}')
        text.AddText(f'Node    : {benchmark}')

        C = ROOT.TCanvas()
        path_pdf = os.path.join(path_plots,f'{channel}_{benchmark}_{era}.pdf')
        C.Print(path_pdf+'[')
        print (benchmark,channel)
        LOKey = findKey(plots_LO,[channel,benchmark])
        NLOKey = findKey(plots_NLO,[channel,benchmark])
        plotNames_LO = set(plots_LO[LOKey].keys())
        plotNames_NLO = set(plots_NLO[NLOKey].keys())

        for plotName in sorted(plotNames_LO.intersection(plotNames_NLO)):
            pLO = plots_LO[LOKey][plotName]
            pNLO = plots_NLO[NLOKey][plotName]
            
            pLO.SetLineColor(ROOT.kGreen+2)
            pLO.SetMarkerColor(ROOT.kGreen+2)
            pLO.SetMarkerStyle(8)
            pLO.SetMarkerSize(0.5)
            pNLO.SetLineColor(ROOT.kRed+1)
            pNLO.SetMarkerColor(ROOT.kRed+1)
            pNLO.SetMarkerStyle(4)
            pNLO.SetMarkerSize(0.5)

            if pLO.GetNbinsX() > 800:
                pLO.Rebin(10)
                pNLO.Rebin(10)
            elif pLO.GetNbinsX() > 300:
                pLO.Rebin(5)
                pNLO.Rebin(5)


            pLO.SetMaximum(max([pLO.GetMaximum(),pNLO.GetMaximum()])*1.2)

            pLO.Draw("e1p")
            pNLO.Draw("e1p same")

            legend = ROOT.TLegend(0.6,0.65,0.88,0.88)
            legend.AddEntry(pLO,'#splitline{Reweighted LO}{Integral = %0.5f}'%pLO.Integral())
            legend.AddEntry(pNLO,'#splitline{NLO}{Integral = %0.5f}'%pNLO.Integral())
            legend.Draw()
            text.Draw()

            C.Print(path_pdf)
        C.Print(path_pdf+']')



