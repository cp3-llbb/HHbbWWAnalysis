import os
import sys
import json
import yaml
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


def runCommand(command,subdir,combine_args=None,debug=False):
    # Modify command #
    command = command.split(' ')
    if '--combine' in command:
        raise RuntimeError('Please do not use --combine in your commands')
    if '-u' not in command:
        command.insert(1,'-u')
    if combine_args is not None and isinstance(combine_args,list):
        command += ['--combine'] + combine_args
    if debug:
        command.append('--debug')

    #command += ['--plotIt']

    #print (' '.join(command))

    # Run command #
    process = subprocess.Popen(command,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT,
                               universal_newlines=True)

    # Browse output #
    output_path = None
    while True:
        nextline = process.stdout.readline().strip()
        if "Output path" in nextline:
            output_path = nextline.split(' ')[-1]
        if nextline == '' and process.poll() is not None:
            break
        sys.stdout.write(nextline+'\n')
        sys.stdout.flush()
    process.communicate()

    # Finalize #
    if output_path is None:
        print ('Could not find output path from command')
    elif not os.path.isdir(output_path):
        print (output_path)
        print (f'Output path {output_path} not valid')
    else:
        limit_json = os.path.join(output_path,subdir,'limits.json')
        if not os.path.exists(limit_json):
            print (f'Limit file {limit_json} does not exist')
        else:
            with open(limit_json,'r') as handle:
                limits = json.load(handle)
            limits = {float(key):float(val) for key,val in limits.items()}
            return limits

    return None
            

            


def runSingleLimit(config,combine_args=None,custom_args=None,additional_args=None,debug=False):
    content = parseYaml(config,custom_args)

    path_single = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               content['outputDir'],
                               'limits',
                               'limits.json')

    path_combine = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                content['outputDir'],
                                'limits',
                                'combination',
                                'limits.json')

    print (f"Running {content['outputDir']}")
    cmd = ['python','-u','produceDataCards.py','--yaml',config]
    if combine_args is not None and isinstance(combine_args,list):
        cmd += ['--combine'] + combine_args
    if custom_args is not None and isinstance(custom_args,list):
        cmd += ['--custom'] + custom_args
    if additional_args is not None and isinstance(additional_args,str):
        cmd += [arg for arg in additional_args.split(' ') if len(arg)>0]
    process = subprocess.Popen(cmd,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT,
                               universal_newlines=True)
    while True:
        nextline = process.stdout.readline()
        if nextline == '' and process.poll() is not None:
            break
        sys.stdout.write(nextline)
        sys.stdout.flush()
    process.communicate()

    path_found = None
    for path in [path_single,path_combine]:
        if os.path.exists(path):
            path_found = path
    if path_found is None:
        limits = {}
    else:
        print (f'Limits for config {config} found in {path_found}')
        with open(path_found,'r') as handle:
            limits = json.load(handle)
        limits = {float(key):float(val) for key,val in limits.items()}
    return limits
            

class LimitScan:
    def __init__(self,curves,legend,combine_args,path_out,jobs,force,unblind,debug):
        self.curves         = curves
        self.legend         = legend
        self.combine_args   = combine_args
        self.path_out       = path_out
        self.jobs           = jobs
        self.force          = force
        self.unblind        = unblind
        self.debug          = debug
        #self.configs = {name:config['configs'] for name,config in configs.items()}
        #self.attributes = {name:{key:val for key,val in config.items() if key != 'configs'} for name,config in configs.items()}
        #self.additional_args = additional_args


        # Make directory #
        if not os.path.isdir(self.path_out):
            os.makedirs(self.path_out)

        path_save = os.path.join(self.path_out,'graphs.root')
        is_save = os.path.exists(path_save) and not self.force

#        # Multi scan or not #
#        if isinstance(self.configs,dict) and \
#           all([isinstance(config,list) for config in self.configs.values()]) and \
#           all([isinstance(innerCfg,dict)  for outerCfg in self.configs.values() for innerCfg in outerCfg]):
#            self.multiscan = True
#        elif isinstance(self.configs,list) and \
#             all([isinstance(config,dict) for config in self.configs]):
#            self.multiscan = False
#        else:
#            raise RuntimeError("Not correct configs")

        if not is_save:
            if self.force:
                print ('Root save found, but will produce it anyway')
            else:
                print ('Root save not found, will produce it')
#            # Run the processes #
#            if self.multiscan:
#                configs_to_run = [{'config': innerCfg['config'],
#                                   'custom': innerCfg['custom'] if 'custom' in innerCfg.keys() else None} 
#                                        for outerCfg in self.configs.values() for innerCfg in outerCfg if 'config' in innerCfg.keys()]
#            else:
#                configs_to_run = [{'config': config['config'],
#                                   'custom': config['custom'] if 'custom' in config.keys() else None}
#                                        for config in self.configs if 'config' in config.keys()]
#            N = len(configs_to_run) 
#            if jobs == -1:
#                jobs = N
#            self.jobs = min(N,jobs,mp.cpu_count())
#            results = self.getLimits(configs_to_run)
            self.getLimits()

#            # Get harcoded values #
#            if self.multiscan:
#                configs_values = {suffix:[innerCfg for innerCfg in outerCfg if 'values' in innerCfg.keys()] for suffix,outerCfg in self.configs.items()}
#            else:
#                configs_values = [config for config in self.configs if 'values' in config.keys()]

        else:
            print ('Root save found, will load it')
            graphs_save = {}
            F = ROOT.TFile(path_save,'READ')
            for key in F.GetListOfKeys():
                if key.GetClassName() == "TDirectoryFile":
                    curveName = key.GetName()
                    self.curves[curveName]['graphs'] = [F.Get(f'{curveName}/{name}') for name in ['data','central','onesigma','twosigma']]
            

        # Computing graphs #
        if not is_save:
            for curveName in self.curves.keys():
                limitsPerPoint = {point['label']:point['limit'] for point in self.curves[curveName]['points'] if point['limit'] is not None}
                graphs,labelConv = self.produceLimitGraphs(limitsPerPoint)
                self.curves[curveName]['graphs'] = graphs 
                self.curves[curveName]['labelConv'] = labelConv
                

        # Per curve limit plot #
        for curveName in self.curves.keys():
            if 'labelConv' in self.curves[curveName].keys():
                labelConv = self.curves[curveName]['labelConv']
            else:
                labelConv = None
            self.plotLimits(graphs=self.curves[curveName]['graphs'],suffix=curveName,labelConv=labelConv)

        # Multi limit plot #
        multiGraphs = []
        curveNames = []
        attributes = []
        for curveName in self.curves.keys():
            if self.curves[curveName]['graphs'] is not None:
                curveNames.append(curveName)
                multiGraphs.append(self.curves[curveName]['graphs'][1])
                attributes.append(self.curves[curveName]['attributes'])
        self.plotMultipleLimits(graphs=multiGraphs,suffixes=curveNames,attributes=attributes,labelConv=labelConv)

#        if self.multiscan:
#            if not is_save:
#                idx = 0
#                multiLimitsPerPoint = {}
#                for suffix,outerCfg in self.configs.items():
#                    multiLimitsPerPoint[suffix] = {}
#                    for innerCfg in outerCfg:
#                        if 'values' in innerCfg:
#                            continue
#                        multiLimitsPerPoint[suffix][innerCfg['label']] = results[idx]
#                        idx += 1
#                for suffix, values in configs_values.items():
#                    if len(values) == 0:
#                        continue
#                    for value in values:
#                        multiLimitsPerPoint[suffix][value['label']] = value['values']
#                # Single scan plots #
#                multiGraphs = []
#                suffixes = []
#                attributes = []
#                for suffix,limitsPerPoint in multiLimitsPerPoint.items():
#                    print (f'Running suffix {suffix}')
#                    suffixes.append(suffix)
#                    attributes.append(self.attributes[suffix])
#                    graphs,labelConv = self.produceLimitGraphs(limitsPerPoint)
#                    graphs_save[suffix] = graphs
#                    self.plotLimits(graphs,suffix,labelConv)
#                    multiGraphs.append(graphs[1])
#                self.plotMultipleLimits(graphs=multiGraphs,suffixes=suffixes,attributes=attributes,labelConv=labelConv)
#            else:
#                multiGraphs = []
#                suffixes = []
#                attributes = []
#                for suffix, graphs in graphs_save.items():
#                    suffixes.append(suffix)
#                    attributes.append(self.attributes[suffix])
#                    self.plotLimits(graphs=graphs,suffix=suffix)
#                    multiGraphs.append(graphs[1])
#                self.plotMultipleLimits(graphs=multiGraphs,suffixes=suffixes,attributes=attributes)
#        else:
#            if not is_save:
#                limitsPerPoint = {config['label']:results[i] 
#                                        for i,config in enumerate(self.configs) if 'values' not in config.keys()} 
#                for values in configs_values:
#                    for value in values:
#                        limitsPerPoint[value['label']] = value['values']
#                graphs,labelConv = self.produceLimitGraphs(limitsPerPoint)
#                graphs_save['limits'] = graphs
#                self.plotLimits(graphs,labelConv)
#            else:
#                graphs = graphs_save['limits']
#                self.plotLimits(graphs)

        # Save root files #
        if not is_save:
            print ('Root save created')
            F = ROOT.TFile(path_save,'RECREATE')
            for curveName in self.curves.keys():
                if self.curves[curveName]['graphs'] is not None:
                    d = F.mkdir(curveName,curveName)
                    d.cd()
                    for graph in self.curves[curveName]['graphs']:
                        graph.Write(graph.GetName())
            F.Close()

 
    def getLimits(self):
        # Run the pool #
        pool_cmds = [(point['command'],point['subdir'],self.combine_args,self.debug) 
                            for values in self.curves.values() for point in values['points'] 
                            if 'command' in point.keys() ]
        with mp.Pool(processes=min(mp.cpu_count(),self.jobs)) as pool:
            results = pool.starmap(runCommand, pool_cmds)
        # Add values to curves dict #
        idx = 0
        for icurve,curveName in enumerate(self.curves.keys()):
            for ipoint in range(len(self.curves[curveName]['points'])):
                if 'command' in self.curves[curveName]['points'][ipoint].keys():
                    self.curves[curveName]['points'][ipoint]['limit'] = \
                            results[idx]
                            #results[icurve*len(self.curves[curveName]['points'])+ipoint]
                    idx += 1

    def produceLimitGraphs(self,limitsPerPoint):
        if len(limitsPerPoint) == 0:
            return None,None

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
            if limits is None:
                continue
            if len(limits.keys()) < 6:
                print(f'[WARNING] Not all limits found for label {label}')
                continue
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
        if len(xpoints) == 0:
            g_data = ROOT.TGraph()
            g_central = ROOT.TGraph()
            g_onesigma = ROOT.TGraph()
            g_twosigma = ROOT.TGraph()
        else:
            g_data = ROOT.TGraph(len(xpoints), np.array(xpoints),  np.array(data))
            g_central = ROOT.TGraph(len(xpoints),  np.array(xpoints),  np.array(central))
            g_onesigma = ROOT.TGraph(len(xpoints)*2,  np.array(xpoints_all),  np.array(onesigma_all))
            g_twosigma = ROOT.TGraph(len(xpoints)*2,  np.array(xpoints_all),  np.array(twosigma_all))

        if floatLabels:
            labelConv = None

        return (g_data,g_central,g_onesigma,g_twosigma),labelConv

    def plotLimits(self,graphs,suffix="",labelConv=None):
        if graphs is None:
            print (f'No graph found for label {suffix}')
            return 

        suffix = str(suffix)
        suffix = suffix.replace(' ','_')
        g_data,g_central,g_onesigma,g_twosigma = graphs

        g_twosigma.SetFillColorAlpha(ROOT.kOrange, 0.5)
        g_twosigma.SetLineStyle(2)
        g_onesigma.SetFillColorAlpha(ROOT.kGreen+2, 0.5)
        g_onesigma.SetLineStyle(2)
        g_central.SetLineWidth(2)
        g_central.SetLineStyle(7)
        g_central.SetLineColor(9)
        g_central.SetMarkerStyle(20)
        g_central.SetMarkerSize(0.6)
        g_central.SetMarkerColor(9)
        g_data.SetLineWidth(2)
        g_data.SetLineStyle(1)
        g_data.SetLineColor(ROOT.kBlack)
        g_data.SetMarkerStyle(20)
        g_data.SetMarkerSize(0.6)
        g_data.SetMarkerColor(ROOT.kBlack)


        if g_central.GetN() == 0:
            return

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

        pdfPath = os.path.join(self.path_out,f'limits{suffix}.pdf')
        c1.Print(pdfPath+'[')

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
        leg = ROOT.TLegend(*self.legend['position'])
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
        c1.Print(pdfPath)

        c1.SetLogy(False)
        c1.Print(pdfPath)

        c1.Print(pdfPath+']')

    def plotMultipleLimits(self,graphs,suffixes,attributes,labelConv=None):
        if len(graphs) == 0:
            print (f'No multi graph found for label {suffix}')
            return 

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

        pdfPath = os.path.join(self.path_out,f'limits.pdf')
        c1.Print(pdfPath+'[')

        # Plot #
        xpoints = sorted(list(set([val for graph in graphs for val in list(graph.GetX()) if graph is not None])))
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

        # Check if color set anywhere, otherwise use Palette #
        custom_color = False
        for attribute in attributes:
            for key in attribute.keys():
                if 'Color' in key:
                    custom_color = True

        # Legend #
        N = len(graphs)
        leg = ROOT.TLegend(*self.legend['position'])
        leg.SetFillColor(ROOT.kWhite)
        leg.SetTextFont(42)
        leg.SetTextSize(self.legend['textsize'])
        b1.Draw()
        for g in graphs:
            g.SetLineWidth(2)
            g.SetLineStyle(7)
            g.SetLineColor(9)
            g.SetMarkerStyle(20)
            g.SetMarkerSize(0.6)
            g.SetMarkerColor(9)
            if custom_color:
                g.Draw("lpsame")
            else:
                g.Draw("lpsame PLC")

        # Attributes #
        for g,suffix,attribute in zip(graphs,suffixes,attributes):
            for funcName,val in attribute.items():
                if not hasattr(g,funcName):
                    raise RuntimeError(f'Graph {g.GetName()} of type {g.__class__.__name__} has not method called {funcName}')
                if isinstance(val,str) and val.startswith('#'):
                    val = ROOT.TColor.GetColor(val)
                getattr(g,funcName)(val)
            leg.AddEntry(g,f"#splitline{{Expected 95% upper limit}}{{{str(suffix)}}}","l")

       
#        tex0 = ROOT.TLatex(0.08,0.91, "#scale[1.4]{#font[61]{CMS}} Internal"+" "*5+"35.87 fb^{-1}(13 TeV),2016")
#        tex0.SetNDC()
#        tex0.SetTextSize(.04)
#        tex0.SetTextFont(42)
#        tex0.Draw("same")
#        tex1 = ROOT.TLatex(0.2,0.21, text)
#        tex1.SetNDC(); tex1.SetTextSize(.04); tex1.SetTextFont(42)
#        tex1.Draw("same")
        
        leg.Draw("same")
        c1.Print(pdfPath)

        c1.SetLogy(False)
        c1.Print(pdfPath)

        c1.Print(pdfPath+']')


if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Limit scans')
    parser.add_argument('--yaml', action='store', required=True, type=str,
                        help='Yaml containing parameters')
    parser.add_argument('-j','--jobs', action='store', required=False, type=int, default=1,
                        help='Number of commands to run in parallel (default = 1), using -1 will spawn all the commands')
    parser.add_argument('--unblind', action='store_true', required=False, default=False,
                        help='Unblind data')
    parser.add_argument('--force', action='store_true', required=False, default=False,
                        help='Force recreation of save')
    parser.add_argument('--debug', action='store_true', required=False, default=False,
                        help='Avoids sending jobs')
    parser.add_argument('-a','--additional', action='store', required=False, default=None, type=str,
                        help='Additional arguments to pass to the commands [as a string, beware to include a space right before]')
    args = parser.parse_args()

    with open(args.yaml,'r') as handle:
        content = yaml.load(handle,Loader=YMLIncludeLoader)

   
    instance = LimitScan(curves          = content['curves'],
                         legend          = content['legend'],
                         combine_args    = content['combine'],
                         path_out        = content['outputDir'],
                         jobs            = args.jobs,
                         force           = args.force,
                         debug           = args.debug,
                         unblind         = args.unblind)


