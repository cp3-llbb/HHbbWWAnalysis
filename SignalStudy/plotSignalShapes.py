import os
import sys
import ROOT
import copy
import glob

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

path = '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SignalLOHistogramsCorr_2018/results/'
era = '2018'

benchmarks = [
    'BenchmarkSM',
    'Benchmark1',
    'Benchmark2',
    'Benchmark3',
    'Benchmark4',
    'Benchmark5',
    'Benchmark6',
    'Benchmark7',
    'Benchmark8',
    'Benchmark9',
    'Benchmark10',
    'Benchmark11',
    'Benchmark12',
]


samples = ['GluGluToHHTo2B2WToLNu2J','GluGluToHHTo2B2VTo2L2Nu','GluGluToHHTo2B2Tau']
variables = ['mHH','cosThetaStar']

min_colors = 51
max_colors = 100   

for sample in samples:
    print (sample)
    Ndict = {var:{} for var in variables}
    for benchmark in benchmarks:
        print ('\t',benchmark)
        hists = {var:{} for var in variables}
        files = glob.glob(os.path.join(path,'*{}*{}.root'.format(sample,benchmark)))
        if len(files) == 0:
            continue
        for f in files:
            print (f)
            FNLO = ROOT.TFile(f,"READ")
            FLO = ROOT.TFile(f.replace(benchmark,''))
            for variable in variables:
                key = os.path.basename(f).replace(benchmark,'').replace('.root','')
                hNLO = copy.deepcopy(FNLO.Get(variable))
                hLO = copy.deepcopy(FLO.Get(variable))
                hists[variable][key] = hNLO
                Ndict[variable][key] = (hNLO.Integral(),hLO.Integral())
            FLO.Close()
            FNLO.Close()
        C = ROOT.TCanvas()
        pdfName = f'pdf/{sample}_{benchmark}_{era}.pdf'
        C.Print(pdfName+'[')
        C.Clear()
        C.SetLogy()
        for variable in variables:
            C.Clear()
            legend = ROOT.TLegend(0.5,0.6,0.98,0.98)
            option = 'hist'
            hmax = max([hist.GetMaximum() for hist in hists[variable].values()])
            for i,name in enumerate(sorted(hists[variable].keys())):
                hist = hists[variable][name]
                color = round(min_colors+i*(max_colors-min_colors)/len(hists[variable]))
                legName = "{} reweighted to NLO (LO->NLO yield change = {:3.5f}%)".format(name,abs(Ndict[variable][name][0]-Ndict[variable][name][1])/Ndict[variable][name][1]*100)
                legend.AddEntry(hist,legName)
                hist.SetLineColor(color)
                hist.SetTitle(benchmark)
                hist.SetMaximum(hmax*1.1)
                hist.Draw(option)
                if 'same' not in option:
                    option += ' same'
            legend.Draw()
            C.Print(pdfName)
        C.Print(pdfName+']')


