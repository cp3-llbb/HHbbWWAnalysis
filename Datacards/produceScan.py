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


def runCommand(command,combine_args=None,subdirs={},debug=False):
    # Modify command #
    command = command.split(' ')
    if '--combine' in command:
        raise RuntimeError('Please do not use --combine in your commands')
    if '-u' not in command:
        command.insert(1,'-u')
    if combine_args is not None and isinstance(combine_args,list) and len(combine_args) > 0:
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
        print (f'Output path {output_path} not valid')
    else:
        # Get results based on subdirs #
        results = {}
        for combineType,subdirNames in subdirs.items(): 
            results[combineType] = {}
            for subdirName in subdirNames:
                json_file = os.path.join(output_path,combineType,subdirName,f'{combineType}.json')
                subdirName = os.path.basename(subdirName)
                if os.path.exists(json_file):
                    with open(json_file,'r') as handle:
                        result = json.load(handle)
                    results[combineType][subdirName] = result
                else:
                    #print (f'Combine type {combineType} file {json_file} does not exist')
                    results[combineType][subdirName] = None

    return results
            
            

class Scan:
    def __init__(self,curves,legend,combine_args,path_out,jobs,force,unblind,debug):
        self.curves         = curves
        self.legend         = legend
        self.combine_args   = combine_args
        self.path_out       = path_out
        self.jobs           = jobs
        self.force          = force
        self.unblind        = unblind
        self.debug          = debug

        # Make directory #
        if not os.path.isdir(self.path_out):
            os.makedirs(self.path_out)

        path_save = os.path.join(self.path_out,'graphs.root')
        is_save = os.path.exists(path_save) and not self.force

        if self.combine_args is None:
            self.combine_args = []

        if not is_save:
            if self.force:
                print ('Root save found, but will produce it anyway')
            else:
                print ('Root save not found, will produce it')
            self.runTasks()

        else:
            print ('Root save found, will load it')
            F = ROOT.TFile(path_save,'READ')
            for curveName in self.curves.keys():
                # Get curve directory #
                if not hasattr(F,curveName):
                    raise RuntimeError(f'Curve name `{curveName}` not found in {path_save}')
                curveDir = getattr(F,curveName)
                self.curves[curveName]['graphs'] = {}
                for combineType in self.combine_args:
                    # Get combine type directory #
                    if not hasattr(curveDir,combineType):
                        raise RuntimeError(f'Curve name `{curveName}` combine type {combineType} not found in {path_save}')
                    combineTypeDir = getattr(curveDir,combineType)
                    self.curves[curveName]['graphs'][combineType] = {}
                    for key in combineTypeDir.GetListOfKeys():
                        # Get subname directory #
                        if not hasattr(combineTypeDir,key.GetName()):
                            raise RuntimeError(f'Curve name `{curveName}` combine type {combineType} entry {key.GetName()} not found in {path_save}')
                        subnameDir = getattr(combineTypeDir,key.GetName())
                        self.curves[curveName]['graphs'][combineType][key.GetName()] = []
                        for g in subnameDir.GetListOfKeys():
                            self.curves[curveName]['graphs'][combineType][key.GetName()].append(subnameDir.Get(g.GetName()))

        # Producing graphs #
        if not is_save:
            for curveName in self.curves.keys():
                self.curves[curveName]['graphs'] = {}
                self.curves[curveName]['labelConv'] = {}
                resultsPerPoint = {}
                for point in self.curves[curveName]['points']:
                    for combineType, results in point['results'].items():
                        if combineType not in resultsPerPoint.keys():
                            self.curves[curveName]['graphs'][combineType] = {}
                            self.curves[curveName]['labelConv'][combineType] = {}
                            resultsPerPoint[combineType] = {}
                        for resultName,resultVals in results.items():
                            if resultName not in resultsPerPoint[combineType].keys():
                                resultsPerPoint[combineType][resultName] = {}
                            resultsPerPoint[combineType][resultName][point['label']] = resultVals
                                
                for combineType in resultsPerPoint.keys():
                    if combineType == 'limits':
                        for resultName,resultVals in resultsPerPoint[combineType].items():
                            graphs,labelConv = self.produceLimitGraphs(resultsPerPoint[combineType][resultName]) 
                            self.curves[curveName]['graphs'][combineType][resultName] = graphs
                            self.curves[curveName]['labelConv'][combineType][resultName] = labelConv
                    elif combineType == 'gof':
                        for resultName,resultVals in resultsPerPoint[combineType].items():
                            self.curves[curveName]['graphs'][combineType][resultName] = self.produceGofGraph(resultsPerPoint[combineType][resultName])
                    else:
                        raise RuntimeError(f'Combine type {combineType} not understood')
                    
        embed()
                
        ###   LIMIT   ###
        if 'limits' in combine_args:
            # Per curve limit plot #
            subnames = []
            for curveName in self.curves.keys():
                if 'labelConv' in self.curves[curveName].keys():
                    labelConv = self.curves[curveName]['labelConv']
                else:
                    labelConv = None
                for subname in self.curves[curveName]['graphs']['limits'].keys():
                    self.plotLimits(graphs=self.curves[curveName]['graphs']['limits'][subname],suffix=f'{curveName}_{subname}',labelConv=labelConv)
                    if subname not in subnames:
                        subnames.append(subname)

            # Multi limit plot #
            for subname in subnames:
                multiGraphs = []
                curveNames = []
                attributes = []
                for curveName in self.curves.keys():
                    if self.curves[curveName]['graphs']['limits'][subname] is not None:
                        curveNames.append(curveName)
                        multiGraphs.append(self.curves[curveName]['graphs']['limits'][subname][1])
                        attributes.append(self.curves[curveName]['attributes'])
                self.plotMultipleGraphs(graphs=multiGraphs,suffixes=curveNames,attributes=attributes,title=f'limits_{subname}',labelConv=labelConv)


        ###   GOF   ### 
        if 'gof' in combine_args:
            # get subnames #
            subnames = []
            for curveName in self.curves.keys():
                subnames += [key for key in self.curves[curveName]['graphs']['gof'].keys() if key not in subnames]
            # Produce GOF plots #
            for subname in subnames:
                multiGraphs = []
                curveNames = []
                for curveName in self.curves.keys():
                    graph = self.curves[curveName]['graphs']['gof'][subname]
                    if graph is not None:
                        curveNames.append(curveName)
                        multiGraphs.append(graph)
                self.plotMultipleGraphs(graphs=multiGraphs,suffixes=curveNames,attributes=None,title=f'gof_{subname}',labelConv=None)


        ###   SAVE   ###
        if not is_save:
            print (f'Creating save file {path_save}')
            F = ROOT.TFile(path_save,'RECREATE')
            for curveName in self.curves.keys():
                curveDir = F.mkdir(curveName,curveName)
                curveDir.cd()
                for combineType in combine_args:
                    combineTypeDir = curveDir.mkdir(combineType,combineType)
                    combineTypeDir.cd()
                    for subname in self.curves[curveName]['graphs'][combineType]:
                        subnameDir = combineTypeDir.mkdir(subname,subname)
                        subnameDir.cd()
                        if self.curves[curveName]['graphs'][combineType][subname] is not None:
                            graphs = self.curves[curveName]['graphs'][combineType][subname]
                            if not isinstance(graphs,tuple) and not isinstance(graphs,list):
                                graphs = [graphs]
                            for graph in graphs:
                                graph.Write(graph.GetName(),ROOT.TObject.kOverwrite)
                    combineTypeDir.cd()
                F.Close()
                print ('... done')

 
    def runTasks(self):
        # Create list of args #
        pool_cmds = []
        for curveName,values in self.curves.items():
            for ipoint,point in enumerate(values['points']):
                if not 'command' in point.keys():
                    continue
                subdirs = {}
                for combine_arg in self.combine_args:
                    if not combine_arg in point.keys():
                        print(f'Key {combine_arg} missing in entry {ipoint} of curve {curveName} : will run it but no plotting used')
                    else:
                        subdirs[combine_arg] = point[combine_arg]
                pool_cmds.append((point['command'],     # command 
                                  self.combine_args,    # combine args 
                                  subdirs,              # subdir
                                  self.debug))          # debug
        # Run the pool #
        with mp.Pool(processes=min(mp.cpu_count(),self.jobs)) as pool:
            results = pool.starmap(runCommand, pool_cmds)
        # Add values to curves dict #
        idx = 0
        for icurve,curveName in enumerate(self.curves.keys()):
            for ipoint in range(len(self.curves[curveName]['points'])):
                if 'command' in self.curves[curveName]['points'][ipoint].keys():
                    result = results[idx]
                    for combineType in result.keys():
                        for subname in result[combineType].keys():
                            if result[combineType][subname] is None:
                                continue
                            if combineType == 'limits':
                                result[combineType][subname] = {float(k):float(v) for k,v in result[combineType][subname].items()}
                            if combineType == 'gof':
                                result[combineType][subname] = self.produceGofValue(result[combineType][subname])
                    self.curves[curveName]['points'][ipoint]['results'] = result
                    idx += 1

    def produceGofValue(self,gofResult):
        data = gofResult['data']
        toys = np.sort(np.array(gofResult['toys']))
        if toys.shape[0] > 0:
            return round((1-(np.abs(toys - data)).argmin()/toys.shape[0]) * 100,6)
        else:
            return -1.

    def produceGofGraph(self,gof):
        g = ROOT.TGraph(len(gof))
        for i,(x,y) in enumerate(gof.items()):
            if y is not None:
                g.SetPoint(i,x,y)
        return g
            
        

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

        g_data.SetName("data")
        g_central.SetName("central")
        g_onesigma.SetName("onesigma")
        g_twosigma.SetName("twosigma")
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

#        g_central.Write()
#        g_data.Write()
#        g_onesigma.Write()
#        g_twosigma.Write()

        
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

    def plotMultipleGraphs(self,graphs,suffixes,attributes,title,labelConv=None):
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

        pdfPath = os.path.join(self.path_out,f'{title}.pdf')
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

   
    instance = Scan(curves          = content['curves'],
                    legend          = content['legend'],
                    combine_args    = content['combine'],
                    path_out        = content['outputDir'],
                    jobs            = args.jobs,
                    force           = args.force,
                    debug           = args.debug,
                    unblind         = args.unblind)


