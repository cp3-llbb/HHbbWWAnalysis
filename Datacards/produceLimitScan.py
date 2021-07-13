import os
import sys
import json
import shutil
import argparse
import subprocess
import ROOT
import numpy as np
import multiprocess as mp
from IPython import embed

from yamlLoader import YMLIncludeLoader

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(55)


def runSingleLimit(config):
    limit_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              'datacards',
                              os.path.basename(config).replace('.yml',''),
                              'limits',
                              'limits.json')
    if not os.path.exists(limit_path):
        print (f'Limits for config {config} not found -> will compute them')
        #process = subprocess.Popen(['python','-u','produceDataCards.py','--worker','--combine','limits','--yaml',config],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
        process = subprocess.Popen(['python','-u','produceDataCards.py','--combine','--yaml',config],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
        while True:
            nextline = process.stdout.readline()
            if nextline == '' and process.poll() is not None:
                break
            sys.stdout.write(nextline)
            sys.stdout.flush()
        process.communicate()
    if not os.path.exists(limit_path):
        #raise RuntimeError(f'Cannot find {limit_path}')
        limits = {}
    else:
        print (f'Limits for config {config} found')
        with open(limit_path,'r') as handle:
            limits = json.load(handle)
        limits = {float(key):float(val) for key,val in limits.items()}
    return (config,limits)
            

class LimitScan:
    def __init__(self,suffix,configs,labels,jobs,unblind,path_out):
        self.configs = configs

        if all([isinstance(config,list) for config in self.configs]):
            self.multiscan = True
        elif all([isinstance(config,str) for config in self.configs]):
            self.multiscan = False
        else:
            raise RuntimeError("Not correct configs")

        self.labels = labels
        if all([isinstance(label,list) for label in self.labels]):
            self.multilabel = True
        elif all([isinstance(label,str) for label in self.labels]):
            self.multilabel = False
        else:
            raise RuntimeError("Not correct labels")

        if len(labels) == 0:
            labels = [str(i) for i in range(1,len(self.configs)+1)]

        if self.multiscan:
            N = 0
            for config in self.configs:
                if self.multilabel:
                    for labels,configs in zip(self.labels,self.configs):
                        if len(labels) != len(config):
                            raise RuntimeError(f"You provided {len(config)} configs but {len(labels)} labels")
                else:
                    if len(labels) != len(config):
                        raise RuntimeError(f"You provided {len(config)} configs but {len(labels)} labels")
                N += len(config)

        else:
            if len(labels) != len(self.configs):
                raise RuntimeError(f"You provided {len(self.configs)} configs but {len(labels)} labels")
            N = len(self.configs)

        self.path_out = path_out

        self.jobs = jobs
        if self.jobs == -1:
            self.jobs = N
        else:
            self.jobs = min(N,self.jobs)

        self.unblind = unblind
            
        if self.multiscan:
            results = self.getLimits([config for configs in self.configs for config in configs])
            if self.multilabel:
                limitsPerPoint = [{labels[i]:results[config] for i,config in enumerate(configs)} for labels,configs in zip(self.labels,self.configs)] 
            else:
                limitsPerPoint = [{self.labels[i]:results[config] for i,config in enumerate(configs)} for configs in self.configs] 
            graphsList = []
            for i in range(len(limitsPerPoint)):
                graphs,labelConv = self.produceLimitGraphs(limitsPerPoint[i])
                self.plotLimits(graphs,labelConv,suffix[i])
                graphsList.append(graphs[1])
            self.plotMultipleLimits(graphsList,labelConv,suffix)
        else:
            results = self.getLimits(self.configs)
            limitsPerPoint = {self.labels[i]:results[config] for i,config in enumerate(self.configs)}
            graphs,labelConv = self.produceLimitGraphs(limitsPerPoint)
            self.plotLimits(graphs,labelConv)


    def getLimits(self,configs):
        pool = mp.Pool(min(mp.cpu_count(),self.jobs))
        results = {result[0]:result[1] for result in pool.map(runSingleLimit, configs)}
        return results

    def produceLimitGraphs(self,limitsPerPoint):
        # Init #
        onesigma_up = []
        twosigma_up = []
        onesigma_low = []
        twosigma_low = []
        central = []
        data = []
        xpoints = []

        # Check if labels are all float or string labels #
        labels = list(limitsPerPoint.keys())
        floatLabels = True
        for label in labels:
            try:
                float(label)
            except:
                floatLabels = False

        if floatLabels:
            labelConv = {label:float(label) for label in labels}
        else:
            labelConv = {label:float(i) for i,label in enumerate(labels)}
            
        # Loop over limits #    
        for label,limits in limitsPerPoint.items():
            if len(limits.keys()) < 6:
                raise RuntimeError(f'Not all limits found for label {label}')
            central.append(limits[50.0])
            twosigma_low.append(limits[2.5])
            onesigma_low.append(limits[16.0])
            onesigma_up.append(limits[84.0])
            twosigma_up.append(limits[97.5])
            data.append(limits[-1.])
            xpoints.append(labelConv[label])

        # Bands #
        onesigma_low.reverse()
        twosigma_low.reverse()
        onesigma_all = onesigma_up + onesigma_low
        twosigma_all = twosigma_up + twosigma_low

        xpoints_all = xpoints + list(reversed(xpoints))
        xpoints_f =  np.array(xpoints) 

        # Graphs #
        g_data = ROOT.TGraph(len(xpoints), np.array(xpoints),  np.array(data))
        g_central = ROOT.TGraph(len(xpoints),  np.array(xpoints),  np.array(central))
        g_onesigma = ROOT.TGraph(len(xpoints)*2,  np.array(xpoints_all),  np.array(onesigma_all))
        g_twosigma = ROOT.TGraph(len(xpoints)*2,  np.array(xpoints_all),  np.array(twosigma_all))

        g_twosigma.SetFillColorAlpha(ROOT.kOrange, 0.5)
        g_twosigma.SetLineStyle(2)
        g_onesigma.SetFillColorAlpha(ROOT.kGreen+2, 0.5)
        g_onesigma.SetLineStyle(2)
        g_central.SetLineWidth(2)
        g_central.SetLineStyle(7)
        g_central.SetLineColor(9)
        g_central.SetMarkerStyle(20)
        g_central.SetMarkerSize(1)
        g_central.SetMarkerColor(9)
        g_data.SetLineWidth(2)
        g_data.SetLineStyle(1)
        g_data.SetLineColor(ROOT.kBlack)
        g_data.SetMarkerStyle(20)
        g_data.SetMarkerSize(1)
        g_data.SetMarkerColor(ROOT.kBlack)

        if floatLabels:
            labelConv = None

        return (g_data,g_central,g_onesigma,g_twosigma),labelConv

    def plotLimits(self,graphs,labelConv,suffix=""):
        suffix = suffix.replace(' ','_')
        g_data,g_central,g_onesigma,g_twosigma = graphs

        rootFile = ROOT.TFile(os.path.join(self.path_out,f'limits{suffix}.root'),"UPDATE")
        
        # Canvas #
        c1 = ROOT.TCanvas("c1","c1",600, 600)
        c1.SetLeftMargin(0.15)
        c1.SetRightMargin(0.1)
        c1.SetTopMargin(0.05)
        c1.SetBottomMargin(0.10)
        if labelConv is not None:
            c1.SetBottomMargin(0.20)
        c1.SetLogy()
        c1.SetGridx()  
        c1.SetGridy()  
        c1.SetTickx()  
        c1.SetTicky()  

        # Plot #
        xpoints = list(g_central.GetX())
        minx = min(xpoints)*0.9
        maxx = max(xpoints)*1.1
        b1 = ROOT.TH1F("b2","b2", len(xpoints)*2, minx, maxx)
        b1.SetTitle(" ")
        b1.GetYaxis().SetTitle("95% C.L. limits on production rate (pb)")
        ylow = g_onesigma.GetHistogram().GetMinimum()
        yhigh = g_twosigma.GetHistogram().GetMaximum()
        b1.GetYaxis().SetRangeUser(ylow*0.9, yhigh*1.2)
        b1.SetStats(0)
        #b1.GetXaxis().SetNDivisions(len(xpoints))
        if labelConv is not None:
            for label,x in labelConv.items():
                b1.GetXaxis().SetBinLabel(b1.GetXaxis().FindBin(x),label)
        b1.LabelsOption("v")
        b1.LabelsDeflate("X")
        b1.LabelsDeflate("Y")

        b1.Draw()
        g_twosigma.Draw("fe3same")
        g_onesigma.Draw("fe3same")
        g_central.Draw("lpsame")
        if self.unblind:
            g_data.Draw("lpsame")

        g_data.SetName("data")
        g_central.SetName("central")
        g_onesigma.SetName("onesigma")
        g_twosigma.SetName("twosigma")
        g_central.Write()
        g_data.Write()
        g_onesigma.Write()
        g_twosigma.Write()

        
        # Legend #
        leg = ROOT.TLegend(0.60,0.75,0.85,0.92)
        leg.SetFillColor(ROOT.kWhite)
        leg.SetTextFont(42)
        if self.unblind:
            leg.AddEntry(g_data,"Observed","pl")
        leg.AddEntry(g_central,"Expected 95% upper limit","l")
        leg.AddEntry(g_onesigma,"1 std. deviation, Expected","f")
        leg.AddEntry(g_twosigma,"2 std. deviation, Expected","f")
#        tex0 = ROOT.TLatex(0.08,0.91, "#scale[1.4]{#font[61]{CMS}} Internal"+" "*5+"35.87 fb^{-1}(13 TeV),2016")
#        tex0.SetNDC()
#        tex0.SetTextSize(.04)
#        tex0.SetTextFont(42)
#        tex0.Draw("same")
#        tex1 = ROOT.TLatex(0.2,0.21, text)
#        tex1.SetNDC(); tex1.SetTextSize(.04); tex1.SetTextFont(42)
#        tex1.Draw("same")
        
        leg.Draw("same")
        c1.Print(os.path.join(self.path_out,f'limits{suffix}.pdf'))

    def plotMultipleLimits(self,graphs,labelConv,suffixes):
        # Canvas #
        c1 = ROOT.TCanvas("c1","c1",600, 600)
        c1.SetLeftMargin(0.15)
        c1.SetRightMargin(0.1)
        c1.SetTopMargin(0.05)
        c1.SetBottomMargin(0.10)
        c1.SetLogy()
        c1.SetGridx()  
        c1.SetGridy()  
        c1.SetTickx()  
        c1.SetTicky()  

        # Plot #
        xpoints = sorted(list(set([val for graph in graphs for val in list(graph.GetX())])))
        minx = min(xpoints)*0.9
        maxx = max(xpoints)*1.1
        b1 = ROOT.TH1F("b2","b2", len(xpoints)*2, minx, maxx)
        b1.SetTitle(" ")
        b1.GetYaxis().SetTitle("95% C.L. limits on production rate (pb)")
        ylow  = min([g.GetHistogram().GetMinimum() for g in graphs])
        yhigh = max([g.GetHistogram().GetMaximum() for g in graphs])
        b1.GetYaxis().SetRangeUser(ylow, yhigh)
        b1.SetStats(0)
        #b1.GetXaxis().SetNDivisions(len(xpoints))
        if labelConv is not None:
            for label,x in labelConv.items():
                b1.GetXaxis().SetBinLabel(b1.GetXaxis().FindBin(x),label)
        b1.LabelsDeflate("X")
        b1.LabelsDeflate("Y")

        b1.Draw()
        for g in graphs:
            g.Draw("lpsame PLC")
        
        # Legend #
        leg = ROOT.TLegend(0.50,0.55,0.85,0.90)
        leg.SetFillColor(ROOT.kWhite)
        leg.SetTextFont(42)
        for g,suffix in zip(graphs,suffixes):
            leg.AddEntry(g,f"#splitline{{Expected 95% upper limit}}{{{suffix}}}","l")
#        tex0 = ROOT.TLatex(0.08,0.91, "#scale[1.4]{#font[61]{CMS}} Internal"+" "*5+"35.87 fb^{-1}(13 TeV),2016")
#        tex0.SetNDC()
#        tex0.SetTextSize(.04)
#        tex0.SetTextFont(42)
#        tex0.Draw("same")
#        tex1 = ROOT.TLatex(0.2,0.21, text)
#        tex1.SetNDC(); tex1.SetTextSize(.04); tex1.SetTextFont(42)
#        tex1.Draw("same")
        
        leg.Draw("same")
        c1.Print(os.path.join(self.path_out,f'limits.pdf'))




if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Limit scans')
    parser.add_argument('--configs', required=False, type=str, default=None,
                        help='Yaml config files for single scan with uncertainty bands, separated by ":"')
    parser.add_argument('--multiconfigs', nargs='+', required=False, default=None,
                        help='Yaml config files for multiple scans without uncertainty bands')
    parser.add_argument('--labels', required=False, type=str, default=None,
                        help='Labels for each config')
    parser.add_argument('--multilabels', nargs='+', required=False, default=None,
                        help='Labels for multiconfigs (can either be one string if labels are the same, or one string per configs)')
    parser.add_argument('--suffix', required=True, type=str,
                        help='Suffix for the subdirectory (one per sets for multiconfigs, separated by ":"')
    parser.add_argument('-j','--jobs', action='store', required=False, type=int, default=1,
                        help='Number of commands to run in parallel (default = 1), using -1 will spawn all the commands')
    parser.add_argument('--unblind', action='store_true', required=False, default=False,
                        help='Unblind data')
    args = parser.parse_args()
    configs = None
    if args.configs is not None:
        if ':' not in args.configs:
            raise RuntimeError("Not correct format for --configs")
        configs = args.configs.split(':')
        subdirName = args.suffix
        suffixes = args.suffix

        if args.labels is None:
            raise RuntimeError('--configs requires --labels')
        if ':' not in args.labels:  
            raise RuntimeError("Not correct format for --labels")
        labels = args.labels.split(':')


    if args.multiconfigs is not None:
        if not all([':' in config for config in args.multiconfigs]):
            raise RuntimeError("Not correct format for --multiconfigs")
        configs = [config.split(':') for config in args.multiconfigs]
        suffixes = args.suffix.split(':')
        if len(suffixes) != len(args.multiconfigs) +1:
            raise RuntimeError("For multiconfigs the format for thr suffix should be : subdir:label1:label2:...")
        subdirName = suffixes[0]
        suffixes = suffixes[1:]

        if args.multilabels is None:
            raise RuntimeError('--multiconfigs requires --multilabels')
        if not all([':' in label for label in args.multilabels]):
            raise RuntimeError("Not correct format for --multilabels")
        labels = [label.split(':') for label in args.multilabels]
            
    if configs is not None:
        # Create output directory #
        path_out = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                'Limits',
                                subdirName)
        if os.path.exists(path_out):
            shutil.rmtree(path_out)
        os.makedirs(path_out)

        # Save command #
        cmd = "python produceLimitScan.py "
        for key,val in args.__dict__.items():
            if isinstance(val,list):
                val = ' '.join([str(v) for v in val])
                cmd += f"--{key} {val} "
            elif isinstance(val,bool):
                if not val:
                    continue
                else:
                    val = ""
                cmd += f"--{key} {val} "
            elif isinstance(val,str):
                cmd += f'--{key} "{val}" '
            elif val is None:
                continue
        with open(os.path.join(path_out,'command.txt'),'w') as handle:
            handle.write(cmd)

        # Run #
        instance = LimitScan(suffixes,configs,labels,args.jobs,args.unblind,path_out)
