import os
import sys
import yaml
import copy
import glob
import math
import json
import ROOT
import numpy as np
from root_numpy import hist2array   
from array import array

from CDFShift import CDFShift

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

class WeightDY:
    def __init__(self, channel, config, title, outputname, mode, era, xaxis, yaxis, rebin_1D=None, rebin_2D=None):
        self.channel    = channel
        self.config     = config
        self.mode       = mode
        self.era        = era
        self.title      = title
        self.xaxis      = xaxis
        self.yaxis      = yaxis
        self.rebin_1D   = rebin_1D
        self.rebin_2D   = rebin_2D
        path_out = os.path.join(os.path.abspath(os.path.dirname(__file__)),'weights')
        if not os.path.exists(path_out):
            os.makedirs(path_out)
        self.outputname = os.path.join(path_out,outputname)
        self.hist_type  = None # TH1 or TH2

        self.histograms = {}
        for cat, inp in self.config.items():
            print ('Recovering %s'%cat)
            dirInfo = self.extractDirInformation(inp['path'],inp['histname'])
            self.histograms[cat] = self.produceDistribution(cat,dirInfo)

        self.produceShape(self.histograms)
        self.produceCDFShift(self.histograms)

        if self.hist_type == 'TH1':
            self.plotWeights1D()
        elif self.hist_type == 'TH2':
            self.plotWeights2D()
        else:
            raise RuntimeError('Hist type not initialized')

        self.saveToJson()
        self.saveToRoot()

    def extractDirInformation(self,directory,histogram):
        """ 
        Returns dict with 
          ->keys : 'luminosity', 'plot_options', 'samples', 'histograms'
              -> 'luminosity' : dict
                  -> key : era
                  -> val : lumi
              -> 'plot_options': dict
                  -> key : option name 
                  -> val : option value
              -> 'samples': dict
                  -> key : name rootfile
                  -> val : dict information (cross section, group, ...)
              -> 'histograms : dict
                  -> key : name rootfile 
                  -> val : histogram 
        """

        # Get YAML info #
        yaml_dict = self.loadYaml(os.path.join(directory,'plots_FakeExtrapolation.yml'),histogram)

        # Loop over file to recover histograms #
        hist_dict = {}
        path_file = os.path.abspath(os.path.join(directory, 'results','*root'))
        for f in glob.glob(path_file):
            name = os.path.basename(f)
            #print ("Looking at %s"%name)
            # Check if in YAML #
            if name not in yaml_dict['samples'].keys():
                print ("[WARNING] File %s will be ignored"%(name))
                continue
            # Save hist in dict #
            h = self.getHistogram(f,histogram)
            hist_dict[name] =h 

        yaml_dict.update({'histograms':hist_dict})
        return yaml_dict

    def loadYaml(self,yaml_path,histogram):
        # Parse YAML #
        with open(yaml_path,"r") as handle:
            full_dict = yaml.load(handle,Loader=yaml.FullLoader)
        # Get Lumi per era # 
        lumi_dict = full_dict["configuration"]["luminosity"]

        # Get plot options #
        opt_to_keep = ['x-axis','y-axis']
        try:
            options_dict = {k:full_dict["plots"][histogram][k] for k in full_dict["plots"][histogram].keys() & opt_to_keep}
        except KeyError:
            print ("Could not find hist %s in YAML %s, will proceed without"%(histogram,yaml_path))
            options_dict = {k:'' for k in opt_to_keep}

        # Get data per sample #
        sample_dict = {}
        info_to_keep = ['cross-section','generated-events','group','type','era']
        for sample,data in full_dict['files'].items():
            sample_dict[sample] = {k:data[k] for k in data.keys() & info_to_keep}

        return {'luminosity':lumi_dict,'plot_options':options_dict,'samples':sample_dict}

    def getHistogram(self,rootfile,histname):
        f = ROOT.TFile(rootfile)
        if f.GetListOfKeys().Contains(histname):
            h = copy.deepcopy(f.Get(histname))
            if self.hist_type is None:
                if 'TH1' in h.ClassName():
                    self.hist_type = 'TH1'
                elif 'TH2' in h.ClassName():
                    self.hist_type = 'TH2'
        else:
            print ("Could not find hist %s in %s"%(histname,rootfile))
            h = None
        f.Close()
        return h

    def produceDistribution(self,category,info_dict):
        hist_dict = info_dict['histograms']
        lumi_dict = info_dict['luminosity']
        samp_dict = info_dict['samples']

        if category in ['ZVeto_1b','ZVeto_2b']: # Cannot compare with Data-MC(except DY) here
            mode = "mc" 
        else:
            mode = self.mode
        #mode = self.mode

        dist = None
        # Loop over hist per sample #
        for sample, data_dict in samp_dict.items():
            # Get hist #
            h = hist_dict[sample]
            if h is None:
                continue
            if dist is None:
                if self.hist_type == 'TH1':
                    dist = ROOT.TH1D(category,category,h.GetNbinsX(),h.GetBinLowEdge(1),h.GetBinLowEdge(h.GetNbinsX())+h.GetBinWidth(h.GetNbinsX()))
                elif self.hist_type == 'TH2':
                    dist = ROOT.TH2D(category,category,
                                     h.GetNbinsX(),h.GetXaxis().GetBinLowEdge(1),h.GetXaxis().GetBinUpEdge(h.GetNbinsX()),
                                     h.GetNbinsY(),h.GetYaxis().GetBinLowEdge(1),h.GetYaxis().GetBinUpEdge(h.GetNbinsY()))
                    dist.Sumw2()
                else:
                    raise RuntimeError('Hist type not initialized')
            # Add histograms depending on type #
            if 'cross-section' in data_dict.keys():
                factor = data_dict["cross-section"]*lumi_dict[self.era]/data_dict["generated-events"]
            else:
                factor = 1.
            if 'branching-ratio' in data_dict.keys():
                factor *= data_dict["branching-ratio"]

            if mode == 'data':
                if data_dict['type'] == 'data':
                    dist.Add(h)
                if data_dict['type'] == 'mc':
                    if data_dict['group'] != 'DY':
                        dist.Add(h,-factor)     
            elif mode == 'mc':
                if data_dict['type'] == 'mc' and data_dict['group'] == 'DY':
                    dist.Add(h,factor)  
            else:
                raise RuntimeError("Unknown mode")
        dist.GetXaxis().SetTitle(self.xaxis)
        if self.hist_type == 'TH2':
            dist.GetYaxis().SetTitle(self.yaxis)
        return dist

    def produceShape(self,hist_dict):
        hist_dict = copy.deepcopy(hist_dict)

        self.N_ZPeak0b = hist_dict['ZPeak_0b'].Integral()
        self.N_ZPeak1b = hist_dict['ZPeak_1b'].Integral()
        self.N_ZPeak2b = hist_dict['ZPeak_2b'].Integral()
        self.N_ZVeto0b = hist_dict['ZVeto_0b'].Integral()
        self.N_ZVeto1b = hist_dict['ZVeto_1b'].Integral()
        self.N_ZVeto2b = hist_dict['ZVeto_2b'].Integral()

        # Rebin if requested #
        if self.hist_type == 'TH1' and self.rebin_1D is not None:
            hist_dict['ZPeak_0b'] = self.rebinWeights1D(hist_dict['ZPeak_0b'])
            hist_dict['ZPeak_1b'] = self.rebinWeights1D(hist_dict['ZPeak_1b'])
            hist_dict['ZPeak_2b'] = self.rebinWeights1D(hist_dict['ZPeak_2b'])
            hist_dict['ZVeto_0b'] = self.rebinWeights1D(hist_dict['ZVeto_0b'])
            hist_dict['ZVeto_1b'] = self.rebinWeights1D(hist_dict['ZVeto_1b'])
            hist_dict['ZVeto_2b'] = self.rebinWeights1D(hist_dict['ZVeto_2b'])
        if self.hist_type == 'TH2' and self.rebin_2D is not None:
            hist_dict['ZPeak_0b'] = self.rebinWeights2D(hist_dict['ZPeak_0b'])
            hist_dict['ZPeak_1b'] = self.rebinWeights2D(hist_dict['ZPeak_1b'])
            hist_dict['ZPeak_2b'] = self.rebinWeights2D(hist_dict['ZPeak_2b'])
            hist_dict['ZVeto_0b'] = self.rebinWeights2D(hist_dict['ZVeto_0b'])
            hist_dict['ZVeto_1b'] = self.rebinWeights2D(hist_dict['ZVeto_1b'])
            hist_dict['ZVeto_2b'] = self.rebinWeights2D(hist_dict['ZVeto_2b'])

        if self.mode == "data":
            # mm0b 1257451.0613
            # mm1b 161402.5533
            # mm2b 27150.4709
            # ee0b 274633.8001
            # ee1b 35004.4788
            # ee2b 6098.4302
            if self.channel == 'MuMu':
                self.factor_1b = 161402.5533 / 1257451.0613
                self.factor_2b = 27150.4709  / 1257451.0613
                self.factor_ZVeto = 180623.891 / 258598.634
            if self.channel == 'ElEl':
                self.factor_1b = 35004.4788 / 274633.8001
                self.factor_2b = 6098.4302  / 274633.8001
                self.factor_ZVeto = 32679.407 / 51517.366 
        elif self.mode == "mc":
            self.factor_1b = self.N_ZPeak1b/self.N_ZPeak0b
            self.factor_2b = self.N_ZPeak2b/self.N_ZPeak0b
            self.factor_ZVeto = 1.

        print ("Sum weight ZVeto_0b : ",hist_dict['ZVeto_0b'].Integral())
        print ("Sum weight ZVeto_1b : ",hist_dict['ZVeto_1b'].Integral())
        print ("Sum weight ZVeto_2b : ",hist_dict['ZVeto_2b'].Integral())
        print ("Sum weight ZPeak_0b : ",hist_dict['ZPeak_0b'].Integral())
        print ("Sum weight ZPeak_1b : ",hist_dict['ZPeak_1b'].Integral())
        print ("Sum weight ZPeak_2b : ",hist_dict['ZPeak_2b'].Integral())

        hist_dict['ZPeak_0b'].Scale(1./hist_dict['ZPeak_0b'].Integral())
        hist_dict['ZPeak_1b'].Scale(1./hist_dict['ZPeak_1b'].Integral())
        hist_dict['ZPeak_2b'].Scale(1./hist_dict['ZPeak_2b'].Integral())

        self.weight_1b = hist_dict['ZPeak_1b'].Clone("Weight1B")
        self.weight_1b.Divide(hist_dict['ZPeak_0b'])
        self.weight_1b.Scale(self.factor_1b*self.factor_ZVeto)
        self.weight_2b = hist_dict['ZPeak_2b'].Clone("Weight2B")
        self.weight_2b.Divide(hist_dict['ZPeak_0b'])
        self.weight_2b.Scale(self.factor_2b*self.factor_ZVeto)
            # self.weight = ZPeak_2b / ZPeak_0b (normalized histogram ratio) * N_2b/N_0b (stats factor)
        if self.hist_type == 'TH2':
            hist_dict['ZPeak_0b'].Scale(self.N_ZPeak0b)
            hist_dict['ZPeak_1b'].Scale(self.N_ZPeak1b)
            hist_dict['ZPeak_2b'].Scale(self.N_ZPeak2b)

        self.shape_1b = hist_dict['ZVeto_0b'].Clone("Shape1B")
        self.shape_1b.Multiply(self.weight_1b)
        self.shape_2b = hist_dict['ZVeto_0b'].Clone("Shape2B")
        self.shape_2b.Multiply(self.weight_2b)
            # self.shape = (ZPeak_2b / ZPeak_0b) * ZVeto_0b (unnormalized histogram) = shape of the data DY in SR

        if self.hist_type == 'TH1':
            hist_dict['ZVeto_0b'].Scale(1./hist_dict['ZVeto_0b'].Integral())

        # Rebin weights # (after computing shape because N bins will be different 
        #self.weight_1b = self.rebinWeights1D(self.weight_1b,0.10,5)
        #self.weight_2b = self.rebinWeights1D(self.weight_2b,0.15,5)

        print ("Factor (1b)         : ",self.factor_1b)
        print ("Shape (1b)          : ",self.shape_1b.Integral())
        print ("Factor (2b)         : ",self.factor_2b)
        print ("Shape (2b)          : ",self.shape_2b.Integral())

    def produceCDFShift(self,hist_dict):
        hist_dict = copy.deepcopy(hist_dict)

        inst1b = CDFShift(hist_dict['ZPeak_0b'],
                          hist_dict['ZPeak_1b'],
                          hist_dict['ZVeto_0b'],
                          norm = hist_dict['ZVeto_0b'].Integral()*self.factor_1b*self.factor_ZVeto,
                          #norm = hist_dict['ZVeto_1b'].Integral(),
                          additional_hist={'ZVeto_1b':hist_dict['ZVeto_1b']})
        inst2b = CDFShift(hist_dict['ZPeak_0b'],
                          hist_dict['ZPeak_2b'],
                          hist_dict['ZVeto_0b'],
                          norm = hist_dict['ZVeto_0b'].Integral()*self.factor_2b*self.factor_ZVeto,
                          #norm = hist_dict['ZVeto_2b'].Integral(),
                          additional_hist={'ZVeto_2b':hist_dict['ZVeto_2b']})
        name1b = self.outputname+'_CDFShift1b.root'
        name2b = self.outputname+'_CDFShift2b.root'
        inst1b.saveToRoot(name1b)
        inst2b.saveToRoot(name2b)

        hist_dict1b = hist_dict
        hist_dict1b.update({'ZPeak_0b_cdf'              : inst1b.cdf1,
                            'ZPeak_1b_cdf'              : inst1b.cdf2,
                            'ZVeto_0b_cdf'              : inst1b.cdf3,
                            'ZVeto_1b_cdf'              : inst1b.additional_cdf['cdf_ZVeto_1b'],
                            'ZVeto_1b_extrapolated'     : inst1b.returnHist(),
                            'ZVeto_1b_extrapolated_cdf' : inst1b.ncdf3,
                            'cdf_graph'                 : inst1b.gncdf3,
                            'shift'                     : inst1b.gshift,
                            'shift_inv'                 : inst1b.gshift_inv})
        hist_dict2b = hist_dict
        hist_dict2b.update({'ZPeak_0b_cdf'              : inst2b.cdf1,
                            'ZPeak_2b_cdf'              : inst2b.cdf2,
                            'ZVeto_0b_cdf'              : inst2b.cdf3,
                            'ZVeto_2b_cdf'              : inst2b.additional_cdf['cdf_ZVeto_2b'],
                            'ZVeto_2b_extrapolated'     : inst2b.returnHist(),
                            'ZVeto_2b_extrapolated_cdf' : inst2b.ncdf3,
                            'cdf_graph'                 : inst1b.gncdf3,
                            'shift'                     : inst1b.gshift,
                            'shift_inv'                 : inst1b.gshift_inv})

        self.plotCDFShift(hist_dict1b,name1b.replace('.root','.pdf'),cat='1b',mode=self.mode) 
        self.plotCDFShift(hist_dict2b,name2b.replace('.root','.pdf'),cat='2b',mode=self.mode)

    @staticmethod
    def plotInCanvas(hist_dict,legend_pos,pdfname,logy=False,opt='',norm=False,text=None,textpos=None,maxy=None):
#        for hist in hist_dict.values():
#            for i in range(1,hist.GetNbinsX()+1):
#                if hist.GetBinContent(i)>0:
#                    hist.SetBinError(i,math.sqrt(hist.GetBinContent(i)))
#                else:
#                    hist.SetBinError(i,math.sqrt(-hist.GetBinContent(i)))

        hist_dict = copy.deepcopy(hist_dict)
        if norm:
            for hist in hist_dict.values():
                hist.Scale(1./hist.Integral())
        c = ROOT.TCanvas('c1','c1',800,600)
        if logy:
            c.SetLogy()
        leg = ROOT.TLegend(*legend_pos)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        maxh = max([h.GetMaximum() for h in hist_dict.values()])*1.1
        if maxy is not None:
            maxh = min(maxh,maxy)
        for name,hist in hist_dict.items():
            leg.AddEntry(hist,name)
        for hist in hist_dict.values():
            hist.GetYaxis().SetRangeUser(0,maxh)
            hist.SetTitle('')
            hist.Draw(opt)
            if 'same' not in opt: opt += " same"
        leg.Draw()
        if text is not None and textpos is not None:
            pt = ROOT.TPaveText(*textpos,"NDC")
            pt.SetBorderSize(0)
            pt.SetFillStyle(0)
            for t in text:
                pt.AddText(t)
            pt.Draw()
        c.Print(pdfname)
            

    def plotCDFShift(self,hist_dict,pdfname,cat,mode):
        if mode == 'mc': mode = 'DY MC'
        if mode == 'data': mode = 'DY data'

        h1 = hist_dict['ZPeak_0b']
        h2 = hist_dict['ZPeak_%s'%(cat)]
        h3 = hist_dict['ZVeto_0b']
        h4 = hist_dict['ZVeto_%s'%(cat)]
        nh4 = hist_dict['ZVeto_%s_extrapolated'%(cat)]

        cdf1 = hist_dict['ZPeak_0b_cdf']
        cdf2 = hist_dict['ZPeak_%s_cdf'%(cat)]
        cdf3 = hist_dict['ZVeto_0b_cdf']
        cdf4 = hist_dict['ZVeto_%s_cdf'%(cat)]
        ncdf4 = hist_dict['ZVeto_%s_extrapolated_cdf'%(cat)]

        h1.SetLineColor(ROOT.kRed+1)
        h2.SetLineColor(ROOT.kBlue+1)
        h3.SetLineColor(ROOT.kGreen+1)
        h4.SetLineColor(ROOT.kOrange+2)
        nh4.SetLineColor(ROOT.kMagenta+1)
               
        cdf1.SetLineColor(ROOT.kRed+1)
        cdf2.SetLineColor(ROOT.kBlue+1)
        cdf3.SetLineColor(ROOT.kGreen+1)
        cdf4.SetLineColor(ROOT.kOrange+2)
        ncdf4.SetLineColor(ROOT.kMagenta+1)

        c = ROOT.TCanvas('c','c',800,600)
        c.Print(pdfname+'[')

        d = {'#splitline{Z peak 0b [%s]}{Integral = %0.3f}'%(mode,h1.Integral())                    : h1,
             '#splitline{Z peak %s [%s]}{Integral = %0.3f}'%(cat,mode,h2.Integral())                : h2,
             '#splitline{Z veto 0b [%s]}{Integral = %0.3f}'%(mode,h3.Integral())                    : h3,
             '#splitline{Z veto %s [DY MC]}{Integral = %0.3f}'%(cat,h4.Integral())                  : h4,
             '#splitline{Z veto %s [%s extrapolated]}{Integral = %0.3f}'%(cat,mode,nh4.Integral())  : nh4}
        self.plotInCanvas(d,[0.4,0.4,0.85,0.85],pdfname,norm=True,opt="hist")
        d = {'Z peak 0b [%s]'%(mode)                    : cdf1,
             'Z peak %s [%s]'%(cat,mode)                : cdf2,
             'Z veto 0b [%s]'%(mode)                    : cdf3,
             'Z veto %s [DY MC]'%(cat)                  : cdf4,
             'Z veto %s [%s extrapolated]'%(cat,mode)   : ncdf4}
        self.plotInCanvas(d,[0.4,0.15,0.85,0.5],pdfname,opt="hist",maxy=1.)

        d = {'#splitline{Z peak 0b [%s]}{Integral = %0.3f}'%(mode,h1.Integral())      : h1,
             '#splitline{Z peak %s [%s]}{Integral = %0.3f}'%(cat,mode,h2.Integral())  : h2,
             '#splitline{Z veto 0b [%s]}{Integral = %0.3f}'%(mode,h3.Integral())      : h3}
        self.plotInCanvas(d,[0.5,0.5,0.85,0.85],pdfname,opt="hist")

        d = {'#splitline{Z veto %s [DY MC]}{Integral = %0.3f}'%(cat,h4.Integral()):                 h4,
             '#splitline{Z veto %s [%s extrapolated]}{Integral = %0.3f}'%(cat,mode,nh4.Integral()): nh4}  
        text = ["KS test (no option) : %0.3e"%(nh4.KolmogorovTest(h4)),
                "KS test (option 'M') : %0.3e"%(nh4.KolmogorovTest(h4,"M")),
                "KS test (option 'X') : %0.3e"%(nh4.KolmogorovTest(h4,"X")),
                "Z Veto 0b correction factor : %0.5f"%(self.factor_ZVeto),
                "Factor Z Peak %s / Z Peak 0b : %0.5f"%(cat,self.factor_1b if cat == '1b' else self.factor_2b)]
        self.plotInCanvas(d,[0.5,0.6,0.85,0.85],pdfname,text=text,textpos=[0.5,0.35,0.85,0.55],opt='hist')

        d = {'Z peak 0b [%s]'%(mode):                     cdf1,
             'Z peak %s [%s]'%(cat,mode):                 cdf2,
             'Z veto 0b [%s]'%(mode):                     cdf3}
        self.plotInCanvas(d,[0.4,0.15,0.85,0.5],pdfname,opt="hist",maxy=1.)
        
        c.Clear()
        hist_dict['cdf_graph'].SetTitle('Non corrected cdf')
        hist_dict['cdf_graph'].Draw()

        c.Print(pdfname)
        c.Clear()

        hist_dict['shift'].SetTitle('Shift;Shift;CDF y')
        hist_dict['shift'].Draw()

        c.Print(pdfname)
        c.Clear()
        
        hist_dict['shift_inv'].SetTitle('Inverted Shift;CDF y;Shift')
        hist_dict['shift_inv'].Draw()

        c.Print(pdfname)
        c.Clear()

        c.Print(pdfname+']')

    def rebinWeights1D(self,h):
        xn = self.rebin_1D
        hn = h.Rebin(len(xn)-1,'rebin',array('d',xn))
        return hn

    def rebinWeights2D(self,h):
        x = [h.GetXaxis().GetBinLowEdge(i) for i in range(1,h.GetNbinsX()+2)]
        y = [h.GetYaxis().GetBinLowEdge(i) for i in range(1,h.GetNbinsY()+2)]
        xn = self.rebin_2D[0]
        yn = self.rebin_2D[1]

        assert all([ix in x for ix in xn])
        assert all([iy in y for iy in yn])
        
        hn = ROOT.TH2F("new","new",len(xn)-1,array('f',xn),len(yn)-1,array('f',yn))
        content, edges = hist2array(h,return_edges=True)
        xedges = edges[0]
        yedges = edges[1]
        
        binCont = {(i+1,j+1):content[i][j] for i in range(len(xedges)-1) for j in range(len(yedges)-1)}
        binErr  = {(i+1,j+1):0. for i in range(len(xedges)-1) for j in range(len(yedges)-1)}
        binNewCont = {}
        binNewErr  = {}
        
        convDict = {}
        for i in range(1,h.GetNbinsX()+1):
            for j in range(1,h.GetNbinsY()+1):
                xLow = h.GetXaxis().GetBinLowEdge(i)
                xUp  = h.GetXaxis().GetBinUpEdge(i)
                yLow = h.GetYaxis().GetBinLowEdge(j)
                yUp  = h.GetYaxis().GetBinUpEdge(j)
                assert h.GetBinContent(i,j) == binCont[(i,j)]
                binErr[(i,j)] = h.GetBinError(i,j)
                ii = hn.GetXaxis().FindBin(xLow)
                jj = hn.GetYaxis().FindBin(yLow)
                convDict[(i,j)] = (ii,jj)
                binNewCont[(ii,jj)] = 0.
                binNewErr[(ii,jj)] = []
        for oldidx, newidx in convDict.items():
            binNewCont[newidx] += binCont[oldidx]
            binNewErr[newidx].append(binCont[oldidx])
        for key,val in binNewErr.items():
            binNewErr[key] = math.sqrt(sum([v**2 for v in val]))
        
        for ix,iy in binNewCont.keys():
            hn.SetBinContent(ix,iy,binNewCont[(ix,iy)])
            hn.SetBinError(ix,iy,binNewErr[(ix,iy)])

        assert abs(h.Integral()-hn.Integral())<0.1

        return hn

    def plotWeights1D(self):
        histograms = copy.deepcopy(self.histograms)

        for name in ['ZVeto_0b','ZPeak_0b','ZPeak_1b','ZPeak_2b']:
            histograms[name].Scale(1./ histograms[name].Integral()) 

        amax = max([histograms[k].GetMaximum() for k in ['ZVeto_0b','ZPeak_0b','ZPeak_1b','ZPeak_2b']])

        C = ROOT.TCanvas("c1", "c1", 600, 900)

        pad1 = ROOT.TPad("pad1", "pad1", 0., 0.0, 1., 1.0)
        pad1.SetTopMargin(0.08)
        pad1.SetBottomMargin(0.61)
        pad1.SetLeftMargin(0.12)
        pad1.SetRightMargin(0.1)
        pad1.SetGridx()
        pad1.SetGridy()
        
        pad2 = ROOT.TPad("pad2", "pad2", 0, 0.0, 1, 0.6)
        pad2.SetTopMargin(0.)
        pad2.SetBottomMargin(0.56)
        pad2.SetLeftMargin(0.12)
        pad2.SetRightMargin(0.1)
        pad2.SetGridx()
        pad2.SetGridy()

        pad3 = ROOT.TPad("pad3", "pad3", 0, 0.0, 1, 0.33)
        pad3.SetTopMargin(0.02)
        pad3.SetBottomMargin(0.18)
        pad3.SetLeftMargin(0.12)
        pad3.SetRightMargin(0.1)
        pad3.SetGridx()
        pad3.SetGridy()

        pad1.Draw()
        pad2.Draw()
        pad3.Draw()

        ##########  PAD 1 ##########
        amax = max([histograms[k].GetMaximum() for k in ['ZVeto_0b','ZPeak_0b','ZPeak_2b','ZPeak_1b','ZPeak_1b']])

        pad1.cd()

        histograms['ZVeto_0b'].SetMaximum(amax*1.2)
        histograms['ZVeto_0b'].SetMinimum(0)
        histograms['ZVeto_0b'].GetYaxis().SetTitle("Normalized events")
        histograms['ZVeto_0b'].GetYaxis().SetTitleSize(0.035)
        histograms['ZVeto_0b'].GetYaxis().SetTitleOffset(1.30)
        histograms['ZVeto_0b'].GetYaxis().SetLabelSize(0.02)
        histograms['ZVeto_0b'].GetXaxis().SetTitleSize(0.)
        histograms['ZVeto_0b'].GetXaxis().SetLabelSize(0.)
        histograms['ZVeto_0b'].SetTitle(self.title)

        histograms['ZVeto_0b'].SetLineWidth(2) 
        histograms['ZVeto_1b'].SetLineWidth(2) 
        histograms['ZVeto_2b'].SetLineWidth(2) 
        histograms['ZPeak_0b'].SetLineWidth(2) 
        histograms['ZPeak_1b'].SetLineWidth(2) 
        histograms['ZPeak_2b'].SetLineWidth(2) 
        histograms['ZVeto_0b'].SetLineColor(418) 
        histograms['ZPeak_0b'].SetLineColor(602)
        histograms['ZPeak_1b'].SetLineColor(619)
        histograms['ZPeak_2b'].SetLineColor(634)


        histograms['ZVeto_0b'].Draw("hist")
        histograms['ZPeak_0b'].Draw("same hist")
        histograms['ZPeak_1b'].Draw("same hist")
        histograms['ZPeak_2b'].Draw("same hist")

        leg1 = ROOT.TLegend(0.50,0.73,0.89,0.91)
        leg1.SetHeader("Legend","C")
        leg1.SetTextSize(0.02)
        mode = "DY data estimation" if self.mode == "data" else "DY MC"
        leg1.AddEntry(histograms['ZVeto_0b'],"#splitline{Z Veto 0 btag (%s)}{Integral = %0.2f}"%(mode,self.N_ZVeto0b))
        leg1.AddEntry(histograms['ZPeak_0b'],"#splitline{Z Peak 0 btag (%s)}{Integral = %0.2f}"%(mode,self.N_ZPeak0b))
        leg1.AddEntry(histograms['ZPeak_1b'],"#splitline{Z Peak 1 btag (%s)}{Integral = %0.2f}"%(mode,self.N_ZPeak1b))
        leg1.AddEntry(histograms['ZPeak_2b'],"#splitline{Z Peak 2 btag (%s)}{Integral = %0.2f}"%(mode,self.N_ZPeak2b))
        leg1.Draw()

        text = ROOT.TPaveText(0.6,0.65,0.85,0.72,"NB NDC")
        text.SetFillStyle(1001)
        text.AddText("N_{1b}/N_{0b} = %0.5f"%self.factor_1b)
        text.AddText("N_{2b}/N_{0b} = %0.5f"%self.factor_2b)
        text.AddText("N_{0b}^{corr} = %0.5f"%self.factor_ZVeto)
        text.Draw()


        ##########  PAD 2 ##########
        pad2.cd()

        #self.weight_1b.SetMaximum(max(self.weight_1b.GetMaximum(),self.weight_2b.GetMaximum())*1.1)
        #self.weight_1b.SetMaximum(self.weight_1b.GetMaximum()*1.8)

        content_w1b = np.ravel(hist2array(self.weight_1b))
        self.weight_1b.SetMaximum(np.percentile(content_w1b,95))
        self.weight_1b.SetMinimum(0.)
        self.weight_1b.SetLineWidth(2)
        self.weight_2b.SetLineWidth(2)
        self.weight_1b.SetLineColor(602)
        self.weight_2b.SetLineColor(600)
        #self.weight.SetLineWidth(2)
        #self.weight.SetLineColor(634)
        self.weight_1b.SetTitle("")

        self.weight_1b.GetYaxis().SetTitle("Weight")
        self.weight_1b.GetYaxis().SetTitleSize(0.035)
        self.weight_1b.GetYaxis().SetTitleOffset(1.30)
        self.weight_1b.GetYaxis().SetLabelSize(0.02)

        self.weight_1b.GetXaxis().SetTitleSize(0.)
        self.weight_1b.GetXaxis().SetLabelSize(0.)

        self.weight_1b.Draw("ep")
        self.weight_2b.Draw("ep same")

        leg2 = ROOT.TLegend(0.65,0.85,0.89,0.98)
        leg2.SetHeader("Legend","C")
        leg2.SetTextSize(0.03)
        leg2.AddEntry(self.weight_1b,"Weight (1b)")
        leg2.AddEntry(self.weight_2b,"Weight (2b)")
        leg2.Draw()

        ##########  PAD 3 ##########
        pad3.cd()

        max_shape = max([h.GetMaximum() for h in [histograms['ZVeto_1b'],histograms['ZVeto_2b'],self.shape_1b,self.shape_2b]])

        self.shape_1b.SetLineWidth(2) 
        self.shape_2b.SetLineWidth(2) 
        self.shape_1b.SetLineColor(602) 
        self.shape_2b.SetLineColor(600) 

        histograms['ZVeto_1b'].SetLineWidth(1)
        histograms['ZVeto_2b'].SetLineWidth(1)
        histograms['ZVeto_1b'].SetLineColor(602)
        histograms['ZVeto_2b'].SetLineColor(600)

        self.shape_1b.SetMinimum(0)
        self.shape_1b.SetMaximum(max_shape*1.1)
        self.shape_1b.SetTitle("")
        self.shape_1b.GetYaxis().SetTitle("DY shape")
        self.shape_1b.GetYaxis().SetTitleSize(0.08)
        self.shape_1b.GetYaxis().SetTitleOffset(0.7)
        self.shape_1b.GetYaxis().SetLabelSize(0.05)

        self.shape_1b.GetXaxis().SetTitleSize(0.08)
        self.shape_1b.GetXaxis().SetTitleOffset(1.15)
        self.shape_1b.GetXaxis().SetLabelSize(0.05)

        self.shape_1b.Draw("H")
        self.shape_2b.Draw("H same")
        histograms['ZVeto_1b'].Draw("H same")
        histograms['ZVeto_2b'].Draw("H same")

        leg3 = ROOT.TLegend(0.5,0.40,0.89,0.93)
        leg3.SetHeader("Legend","C")
        leg3.SetTextSize(0.04)
        leg3.AddEntry(self.shape_1b,"#splitline{DY shape from data (1b)}{Integral = %0.2f}"%self.shape_1b.Integral() if self.mode == 'data' else "#splitline{DY shape from MC (1b)}{Integral = %0.2f}"%self.shape_1b.Integral())
        leg3.AddEntry(histograms['ZVeto_1b'],"#splitline{Z Veto (1b) (DY MC)}{Integral = %0.2f}"%self.N_ZVeto1b)
        leg3.AddEntry(self.shape_2b,"#splitline{DY shape from data (2b)}{Integral = %0.2f}"%self.shape_2b.Integral() if self.mode == 'data' else "#splitline{DY shape from MC (2b)}{Integral = %0.2f}"%self.shape_2b.Integral())
        leg3.AddEntry(histograms['ZVeto_2b'],"#splitline{Z Veto (2b) (DY MC)}{Integral = %0.2f}"%self.N_ZVeto2b)
        leg3.Draw()



        C.Print(self.outputname+'_%s.pdf'%self.era)


    def plotWeights2D(self):
        histograms = copy.deepcopy(self.histograms)

        ROOT.gStyle.SetTitleW(1.)
        C1 = ROOT.TCanvas("c1", "c1", 4000, 4000)
        C1.Divide(3,3)

        mode = "DY data estimation" if self.mode == "data" else "DY MC"

#        print ('-'*80)
#        print ("Examining ZVeto_0b")
#        self.findNegativeBinsInTh2(histograms['ZVeto_0b'])
#        print ('-'*80)
#        print ("Examining ZPeak_0b")
#        self.findNegativeBinsInTh2(histograms['ZPeak_0b'])
#        print ('-'*80)
#        print ("Examining ZPeak_2b")
#        self.findNegativeBinsInTh2(histograms['ZPeak_2b'])
#        print ('-'*80)
#        print ("Examining ZVeto_0b")
#        self.findNegativeBinsInTh2(histograms['ZPeak_1b'])
#        print ('-'*80)



        amax = max([histograms[k].GetMaximum() for k in ['ZVeto_0b','ZPeak_0b','ZPeak_2b','ZPeak_1b']])

        #----- Zpeak and Zveto 0b -----#
        C1.cd(1)
        ROOT.gPad.SetRightMargin(0.15)
        histograms['ZPeak_1b'].SetTitle('Z Peak 1 btag (%s) : Integral = %0.2f'%(mode,self.N_ZPeak1b))
        histograms['ZPeak_1b'].SetTitleSize(0.9,"t")
        #histograms['ZPeak_1b'].GetZaxis().SetRangeUser(0,amax)
        histograms['ZPeak_1b'].GetXaxis().SetRangeUser(0,histograms['ZPeak_1b'].GetXaxis().GetBinUpEdge(histograms['ZPeak_1b'].GetNbinsX()))
        histograms['ZPeak_1b'].GetYaxis().SetRangeUser(0,histograms['ZPeak_1b'].GetYaxis().GetBinUpEdge(histograms['ZPeak_1b'].GetNbinsY()))
        histograms['ZPeak_1b'].Draw("colz")
        C1.cd(7)
        ROOT.gPad.SetRightMargin(0.15)
        histograms['ZPeak_2b'].SetTitle('Z Peak 2 btag (%s) : Integral = %0.2f'%(mode,self.N_ZPeak2b))
        #histograms['ZPeak_2b'].GetZaxis().SetRangeUser(0,amax)
        histograms['ZPeak_2b'].GetXaxis().SetRangeUser(0,histograms['ZPeak_2b'].GetXaxis().GetBinUpEdge(histograms['ZPeak_2b'].GetNbinsX()))
        histograms['ZPeak_2b'].GetYaxis().SetRangeUser(0,histograms['ZPeak_2b'].GetYaxis().GetBinUpEdge(histograms['ZPeak_2b'].GetNbinsY()))
        histograms['ZPeak_2b'].Draw("colz")
        C1.cd(5)
        ROOT.gPad.SetRightMargin(0.15)
        histograms['ZVeto_0b'].Draw("colz")
        #histograms['ZVeto_0b'].GetZaxis().SetRangeUser(0,amax)
        histograms['ZVeto_0b'].GetXaxis().SetRangeUser(0,histograms['ZVeto_0b'].GetXaxis().GetBinUpEdge(histograms['ZVeto_0b'].GetNbinsX()))
        histograms['ZVeto_0b'].GetYaxis().SetRangeUser(0,histograms['ZVeto_0b'].GetYaxis().GetBinUpEdge(histograms['ZVeto_0b'].GetNbinsY()))
        histograms['ZVeto_0b'].SetTitle('Z Veto 0 btag (%s) : Integral = %0.2f'%(mode,self.N_ZVeto0b))
        C1.cd(4)
        ROOT.gPad.SetRightMargin(0.15)
        histograms['ZPeak_0b'].SetTitle('Z Peak 0 btag (%s) : Integral = %0.2f'%(mode,self.N_ZPeak0b))
        histograms['ZPeak_0b'].GetXaxis().SetRangeUser(0,histograms['ZPeak_0b'].GetXaxis().GetBinUpEdge(histograms['ZPeak_0b'].GetNbinsX()))
        histograms['ZPeak_0b'].GetYaxis().SetRangeUser(0,histograms['ZPeak_0b'].GetYaxis().GetBinUpEdge(histograms['ZPeak_0b'].GetNbinsY()))
        #histograms['ZPeak_0b'].GetZaxis().SetRangeUser(0,amax)
        histograms['ZPeak_0b'].Draw("colz")

        #----- Weights -----#
        C1.cd(2)
        ROOT.gPad.SetRightMargin(0.15)
        #ROOT.gPad.SetLogz()
        self.weight_1b.SetTitle('Weight (1b);Lead lepton #eta;Lead lepton P_{T}')
        #content_w1b = np.ravel(hist2array(self.weight_1b))
        #self.weight_1b.GetZaxis().SetRangeUser(0,np.percentile(content_w1b,95))
        self.weight_1b.GetXaxis().SetRangeUser(25,self.weight_1b.GetXaxis().GetBinUpEdge(self.weight_1b.GetNbinsX()))
        self.weight_1b.GetYaxis().SetRangeUser(25,self.weight_1b.GetYaxis().GetBinUpEdge(self.weight_1b.GetNbinsY()))
        self.weight_1b.Draw("colz")
        C1.cd(8)
        ROOT.gPad.SetRightMargin(0.15)
        #ROOT.gPad.SetLogz()
        self.weight_2b.SetTitle('Weight (2b);Lead lepton #eta;Lead lepton P_{T};')
        #content_w2b = np.ravel(hist2array(self.weight_2b))
        #self.weight_2b.GetZaxis().SetRangeUser(0,np.percentile(content_w2b,95))
        self.weight_2b.GetXaxis().SetRangeUser(25,self.weight_2b.GetXaxis().GetBinUpEdge(self.weight_2b.GetNbinsX()))
        self.weight_2b.GetYaxis().SetRangeUser(25,self.weight_2b.GetYaxis().GetBinUpEdge(self.weight_2b.GetNbinsY()))
        self.weight_2b.Draw("colz")

        #----- DY shape and MC DY in Zveto -----#
        #max_shape = max([h.GetMaximum() for h in [histograms['ZVeto_1b'],histograms['ZVeto_2b'],self.shape_1b,self.shape_2b]])
        mode = "DY shape from data" if self.mode == "data" else "DY shape from MC"

        C1.cd(3)
        ROOT.gPad.SetRightMargin(0.15)
        self.shape_1b.SetTitle(mode+' (Z Veto 1b) : Integral = %0.2f'%self.shape_1b.Integral())
        self.shape_1b.GetXaxis().SetRangeUser(0,self.shape_1b.GetXaxis().GetBinUpEdge(self.shape_1b.GetNbinsX()))
        self.shape_1b.GetYaxis().SetRangeUser(0,self.shape_1b.GetYaxis().GetBinUpEdge(self.shape_1b.GetNbinsY()))
        self.shape_1b.Draw("colz")
        C1.cd(9)
        ROOT.gPad.SetRightMargin(0.15)
        self.shape_2b.SetTitle(mode+' (Z Veto 2b) : Integral = %0.2f'%self.shape_2b.Integral())
        self.shape_2b.GetXaxis().SetRangeUser(0,self.shape_2b.GetXaxis().GetBinUpEdge(self.shape_2b.GetNbinsX()))
        self.shape_2b.GetYaxis().SetRangeUser(0,self.shape_2b.GetYaxis().GetBinUpEdge(self.shape_2b.GetNbinsY()))
        self.shape_2b.Draw("colz")


        #----- Summary -----#
        C1.cd(6)
        pt = ROOT.TPaveText(.1,.1,.9,.9)
        pt.AddText("Weight_{1b/2b} = #frac{N_{1b/2b}^{Z Peak}}{N_{0b}^{Z Peak}} x #frac{h(Z Peak)_{1b/2b}^{norm}}{h(Z Peak)_{0b}^{norm}}")
        pt.AddText("DY shape_{1b/2b} = Weight_{1b/2b} x h(Z Veto)_{0b}")
        pt.AddText("#frac{N_{1b}^{Z Peak}}{N_{0b}^{Z Peak}} = %0.5f"%self.factor_1b)
        pt.AddText("#frac{N_{2b}^{Z Peak}}{N_{0b}^{Z Peak}} = %0.5f"%self.factor_2b)
        pt.Draw()


        #----- Arrows -----#
        C1.cd(0)
        div1b = ROOT.TText(0.32,0.655,"Divide")
        div1b.SetTextSize(0.02)
        div1b.Draw()

        arrdiv1b_1 = ROOT.TArrow(0.29,0.62,0.32,0.65,0.01)
        arrdiv1b_1.Draw()
        arrdiv1b_2 = ROOT.TArrow(0.29,0.71,0.32,0.67,0.01)
        arrdiv1b_2.Draw()
        arrdiv1b_3 = ROOT.TArrow(0.37,0.67,0.40,0.70,0.01)
        arrdiv1b_3.Draw()

        div2b = ROOT.TText(0.32,0.325,"Divide")
        div2b.SetTextSize(0.02)
        div2b.Draw()

        arrdiv2b_1 = ROOT.TArrow(0.29,0.29,0.32,0.32,0.01)
        arrdiv2b_1.Draw()
        arrdiv2b_2 = ROOT.TArrow(0.29,0.38,0.32,0.34,0.01)
        arrdiv2b_2.Draw()
        arrdiv2b_3 = ROOT.TArrow(0.37,0.32,0.40,0.30,0.01)
        arrdiv2b_3.Draw()

        mul1b = ROOT.TText(0.65,0.655,"Multiply")
        mul1b.SetTextSize(0.02)
        mul1b.Draw()

        arrmul1b_1 = ROOT.TArrow(0.62,0.62,0.65,0.65,0.01)
        arrmul1b_1.Draw()
        arrmul1b_2 = ROOT.TArrow(0.62,0.71,0.65,0.67,0.01)
        arrmul1b_2.Draw()
        arrmul1b_3 = ROOT.TArrow(0.72,0.67,0.74,0.70,0.01)
        arrmul1b_3.Draw()

        mul2b = ROOT.TText(0.65,0.325,"Multiply")
        mul2b.SetTextSize(0.02)
        mul2b.Draw()

        arrmul2b_1 = ROOT.TArrow(0.62,0.29,0.65,0.32,0.01)
        arrmul2b_1.Draw()
        arrmul2b_2 = ROOT.TArrow(0.62,0.38,0.65,0.34,0.01)
        arrmul2b_2.Draw()
        arrmul2b_3 = ROOT.TArrow(0.71,0.32,0.74,0.30,0.01)
        arrmul2b_3.Draw()


        C1.Print(self.outputname+'1_%s.pdf'%self.era)

        #----- Shape and DY MC comparison -----#
        max_1b = max([histograms['ZVeto_1b'].GetMaximum(),self.shape_1b.GetMaximum()])
        max_2b = max([histograms['ZVeto_2b'].GetMaximum(),self.shape_2b.GetMaximum()])

        C2 = ROOT.TCanvas("c2", "c2", 4000, 4000)
        C2.Divide(2,2)

        C2.cd(1)
        ROOT.gPad.SetRightMargin(0.15)
        self.shape_1b.GetZaxis().SetRangeUser(0,max_1b)
        self.shape_1b.Draw("colz")
        C2.cd(2)
        ROOT.gPad.SetRightMargin(0.15)
        histograms['ZVeto_1b'].SetTitle('DY MC (Z Veto 1b) : Integral = %0.2f'%self.N_ZVeto1b)
        histograms['ZVeto_1b'].GetZaxis().SetRangeUser(0,max_1b)
        histograms['ZVeto_1b'].GetXaxis().SetRangeUser(0,histograms['ZVeto_1b'].GetXaxis().GetBinUpEdge(histograms['ZVeto_1b'].GetNbinsX()))
        histograms['ZVeto_1b'].GetYaxis().SetRangeUser(0,histograms['ZVeto_1b'].GetYaxis().GetBinUpEdge(histograms['ZVeto_1b'].GetNbinsY()))
        histograms['ZVeto_1b'].Draw("colz")
        C2.cd(3)
        ROOT.gPad.SetRightMargin(0.15)
        self.shape_2b.GetZaxis().SetRangeUser(0,max_2b)
        self.shape_2b.Draw("colz")
        C2.cd(4)
        ROOT.gPad.SetRightMargin(0.15)
        histograms['ZVeto_2b'].SetTitle('DY MC (Z Veto 2b) : Integral = %0.2f'%self.N_ZVeto2b)
        histograms['ZVeto_2b'].GetZaxis().SetRangeUser(0,max_2b)
        histograms['ZVeto_2b'].GetXaxis().SetRangeUser(0,histograms['ZVeto_2b'].GetXaxis().GetBinUpEdge(histograms['ZVeto_2b'].GetNbinsX()))
        histograms['ZVeto_2b'].GetYaxis().SetRangeUser(0,histograms['ZVeto_2b'].GetYaxis().GetBinUpEdge(histograms['ZVeto_2b'].GetNbinsY()))
        histograms['ZVeto_2b'].Draw("colz")
        C2.Print(self.outputname+'2_%s.pdf'%self.era)
    


    def saveToJson(self):
        json_dict = {'dimension': 2, 'variables': ['Eta','Pt'], 'error_type': 'absolute'}
        if self.hist_type == 'TH1':
            for name,weight in zip(['weight_1b','weight_2b'],[self.weight_1b, self.weight_2b]):
                # Fill data #
                binning = []
                data = []
                xAxis = weight.GetXaxis()
                for i in range(1,weight.GetNbinsX()+1):
                    # Get content #
                    cont = weight.GetBinContent(i)
                    if cont <= 0:
                        continue
                    # Get error #
                    err = weight.GetBinError(i)

                    # Get binning #
                    xLow  = xAxis.GetBinLowEdge(i)
                    xHigh = xAxis.GetBinUpEdge(i)
                    
                    if len(binning) == 0:
                        binning.append(xLow)
                    binning.append(xHigh)

                    # Save data #
                    data.append({'bin':[xLow,xHigh],'value':cont,'error_low':err,'error_high':err}) 
               
                json_dict['binning'] = {'x':[-2.5,2.5],'y':binning}
                json_dict['data'] = [{'bin':[-2.5,2.5],'values':data}]

                with open(self.outputname+'_%s_%s.json'%(name,self.era),'w') as f:
                    json.dump(json_dict,f,indent=4)
                print ("Saved json to",self.outputname+'_%s_%s.json'%(name,self.era))
        if self.hist_type == 'TH2':
            for name,weight in zip(['weight_1b','weight_2b'],[self.weight_1b, self.weight_2b]):
                # Fill data #
                eta_binning = []
                pt_binning = []
                data = []
                xAxis = weight.GetXaxis()
                yAxis = weight.GetYaxis()
                # Loop over eta bins in Y #
                for j in range(1,weight.GetNbinsX()+1):
                    bin_dict = {'values':[]}
                    etaLow = xAxis.GetBinLowEdge(j)
                    etaHigh = xAxis.GetBinUpEdge(j)
                    bin_dict['bin'] = [etaLow,etaHigh]
                    if len(eta_binning) == 0:
                        eta_binning.append(etaLow)
                    eta_binning.append(etaHigh)
                    # Loop over PT bins in Y #
                    for i in range(1,weight.GetNbinsY()+1):
                        #if weight.GetBinContent(i,j) == 0:
                        #    continue
                        PtLow  = yAxis.GetBinLowEdge(i)
                        PtHigh = yAxis.GetBinUpEdge(i)
                        bin_dict['values'].append({'bin' : [PtLow,PtHigh],
                                                   'value' : max(0.,weight.GetBinContent(i,j)),
                                                   'error_low' : weight.GetBinError(i,j),
                                                   'error_high' : weight.GetBinError(i,j)})
                        if len(pt_binning) == 0:
                            if PtLow not in pt_binning:
                                pt_binning.append(PtLow)
                        if PtHigh not in pt_binning:
                            pt_binning.append(PtHigh)
                    data.append(bin_dict)

                json_dict['binning'] = {'x':eta_binning,'y':pt_binning}
                json_dict['data'] = data

                # Save to json #
                with open(self.outputname+'_%s_%s.json'%(name,self.era),'w') as f:
                    json.dump(json_dict,f,indent=2)
                print ("Saved json to",self.outputname+'_%s_%s.json'%(name,self.era))
                    

            
    def saveToRoot(self):
        root_file = ROOT.TFile(self.outputname+'_%s.root'%self.era,"recreate")
        self.weight_1b.Write("Weight1B")
        self.weight_2b.Write("Weight2B")
        root_file.Write()
        root_file.Close()

        
        
if __name__ == "__main__":
    with open (sys.argv[1],'r') as handle:
        d = yaml.load(handle,Loader=yaml.SafeLoader)

    instance = WeightDY(channel     = 'ElEl',
                        config      = d['ElElChannel'],
                        title       = d['title']+' (e^{+}e^{-} channel)',
                        outputname  = d['filename'].format(**{'channel':'ElEl','type':'data'}),
                        mode        = 'data',
                        era         = d['era'],
                        xaxis       = d['xaxis'] if 'xaxis' in d.keys() else None,
                        yaxis       = d['yaxis'] if 'yaxis' in d.keys() else None,
                        rebin_1D    = d['rebin_1D'] if 'rebin_1D' in d.keys() else None,
                        rebin_2D    = d['rebin_2D'] if 'rebin_2D' in d.keys() else None)
    instance = WeightDY(channel     = 'ElEl',
                        config      = d['ElElChannel'],
                        title       = d['title']+' (e^{+}e^{-} channel)',
                        outputname  = d['filename'].format(**{'channel':'ElEl','type':'mc'}),
                        mode        = 'mc',
                        era         = d['era'],
                        xaxis       = d['xaxis'] if 'xaxis' in d.keys() else None,
                        yaxis       = d['yaxis'] if 'yaxis' in d.keys() else None,
                        rebin_1D    = d['rebin_1D'] if 'rebin_1D' in d.keys() else None,
                        rebin_2D    = d['rebin_2D'] if 'rebin_2D' in d.keys() else None)
    instance = WeightDY(channel     = 'MuMu',
                        config      = d['MuMuChannel'],
                        title       = d['title']+' (#mu^{+}#mu^{-} channel)',
                        outputname  = d['filename'].format(**{'channel':'MuMu','type':'data'}),
                        mode        = 'data',
                        era         = d['era'],
                        xaxis       = d['xaxis'] if 'xaxis' in d.keys() else None,
                        yaxis       = d['yaxis'] if 'yaxis' in d.keys() else None,
                        rebin_1D    = d['rebin_1D'] if 'rebin_1D' in d.keys() else None,
                        rebin_2D    = d['rebin_2D'] if 'rebin_2D' in d.keys() else None)
    instance = WeightDY(channel     = 'MuMu',
                        config      = d['MuMuChannel'],
                        title       = d['title']+' (#mu^{+}#mu^{-} channel)',
                        outputname  = d['filename'].format(**{'channel':'MuMu','type':'mc'}),
                        mode        = 'mc',
                        era         = d['era'],
                        xaxis       = d['xaxis'] if 'xaxis' in d.keys() else None,
                        yaxis       = d['yaxis'] if 'yaxis' in d.keys() else None,
                        rebin_1D    = d['rebin_1D'] if 'rebin_1D' in d.keys() else None,
                        rebin_2D    = d['rebin_2D'] if 'rebin_2D' in d.keys() else None)
