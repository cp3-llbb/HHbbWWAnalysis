import os
import sys
import yaml
import copy
import glob
import json
import ROOT
import numpy as np

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

class WeightDY:
    def __init__(self, config, title, outputname, mode, era):
        self.config     = config
        self.mode       = mode
        self.era        = era
        self.title      = title
        self.outputname = outputname
        self.hist_type  = None # TH1 or TH2

        self.histograms = {}
        for cat, inp in self.config.items():
            dirInfo = self.extractDirInformation(inp['path'],inp['histname'])
            self.histograms[cat] = self.produceDistribution(cat,dirInfo)

        self.produceShape(self.histograms)

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
        yaml_dict = self.loadYaml(os.path.join(directory,'plots.yml'),histogram)

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
        plot_dict = info_dict['plot_options']
        samp_dict = info_dict['samples']

        if category in ['ZVeto_1b','ZVeto_2b']: # Cannot compare with Data-MC(except DY) here
            mode = "mc" 
        else:
            mode = self.mode

        dist = None
        # Loop over hist per sample #
        for sample, data_dict in samp_dict.items():
            # Get hist #
            h = hist_dict[sample]
            if self.hist_type == 'TH2':
                h.Rebin2D(2,3)
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
            if mode == 'data':
                if data_dict['type'] == 'data':
                    dist.Add(h)
                if data_dict['type'] == 'mc':
                    if data_dict['group'] != 'DY':
                        factor = data_dict["cross-section"]*lumi_dict[self.era]/data_dict["generated-events"]    
                        #if data_dict['group'] == 'ttbar':
                        #    factor *= 0.7
                        dist.Add(h,-factor)     
            elif mode == 'mc':
                if data_dict['type'] == 'mc' and data_dict['group'] == 'DY':
                    factor = data_dict["cross-section"]*lumi_dict[self.era]/data_dict["generated-events"]    
                    if data_dict['group'] == 'ttbar':
                        factor *= 0.7
                    dist.Add(h,factor)  
            else:
                raise RuntimeError("Unknown mode")
        dist.GetXaxis().SetTitle(plot_dict['x-axis']) 
        dist.GetYaxis().SetTitle(plot_dict['y-axis']) 
        #dist.Sumw2()
        return dist

    def produceShape(self,hist_dict):

        hist_dict['ZPeak_0b'] = self.padNegativeBinsWithZeros(hist_dict['ZPeak_0b'])
        hist_dict['ZPeak_1b'] = self.padNegativeBinsWithZeros(hist_dict['ZPeak_1b'])
        hist_dict['ZPeak_2b'] = self.padNegativeBinsWithZeros(hist_dict['ZPeak_2b'])
        hist_dict['ZVeto_0b'] = self.padNegativeBinsWithZeros(hist_dict['ZVeto_0b'])
        hist_dict['ZVeto_1b'] = self.padNegativeBinsWithZeros(hist_dict['ZVeto_1b'])
        hist_dict['ZVeto_2b'] = self.padNegativeBinsWithZeros(hist_dict['ZVeto_2b'])


        self.N_ZPeak0b = hist_dict['ZPeak_0b'].Integral()
        self.N_ZPeak1b = hist_dict['ZPeak_1b'].Integral()
        self.N_ZPeak2b = hist_dict['ZPeak_2b'].Integral()
        self.N_ZVeto0b= hist_dict['ZVeto_0b'].Integral()
        self.N_ZVeto1b= hist_dict['ZVeto_1b'].Integral()
        self.N_ZVeto2b= hist_dict['ZVeto_2b'].Integral()

        self.factor_1b = self.N_ZPeak1b/self.N_ZPeak0b
        self.factor_2b = self.N_ZPeak2b/self.N_ZPeak0b

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
        self.weight_1b.Scale(self.factor_1b)
        self.weight_2b = hist_dict['ZPeak_2b'].Clone("Weight2B")
        self.weight_2b.Divide(hist_dict['ZPeak_0b'])
        self.weight_2b.Scale(self.factor_2b)
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
        #self.weight_1b = self.rebinWeights(self.weight_1b,0.10,5)
        #self.weight_2b = self.rebinWeights(self.weight_2b,0.15,5)

        print ("Factor (1b)         : ",self.factor_1b)
        print ("Shape (1b)          : ",self.shape_1b.Integral())
        print ("Factor (2b)         : ",self.factor_2b)
        print ("Shape (2b)          : ",self.shape_2b.Integral())

    @staticmethod
    def findNegativeBinsInTh2(h):
        from ctypes import c_int
        x = c_int()
        y = c_int()
        z = c_int()
        for i in range(h.GetNbinsX()*h.GetNbinsY()):
            cont = h.GetBinContent(i)
            if cont<0:
                h.GetBinXYZ(i,x,y,z)
                print ("\tGlobal bin %d : X bin %d, Y bin %d, Z bin %d"%(i,x.value,y.value,z.value))
                print ("\tBin content :",cont)
                print ("\tBin low edge : X %0.2f, Y %0.2f"%(h.GetXaxis().GetBinLowEdge(x.value),h.GetYaxis().GetBinLowEdge(y.value)))

    @staticmethod
    def padNegativeBinsWithZeros(h):
        if 'TH1' in h.ClassName():
            for i in range(h.GetNbinsX()):
                if h.GetBinContent(i)<0:
                    h.SetBinContent(i,0)
        elif 'TH2' in h.ClassName():
            for i in range(h.GetNbinsX()*h.GetNbinsY()):
                if h.GetBinContent(i)<0:
                    h.SetBinContent(i,0)
        return h 


    @staticmethod
    def rebinWeights(h,threshold,max_rebin):
        # Produce bins for rebin operation #
        xbins = []
        current_rebin = 0
        cont = []
        err = []
        for i in range(1,h.GetNbinsX()+1):
            cont.append(h.GetBinContent(i))
            err.append(h.GetBinError(i))
            mean_cont = sum(cont)/len(cont)
            mean_err = sum(err)/len(err)
            if mean_cont == 0.:
                xbins.append(h.GetBinLowEdge(i))
                cont = []
                err = []
                continue
            if mean_err/mean_cont <= threshold:
                xbins.append(h.GetBinLowEdge(i))
                cont = []
                err = []
                continue
            else:
                current_rebin += 1
            if current_rebin == max_rebin:
                xbins.append(h.GetBinLowEdge(i))
                current_rebin = 0
                cont = []
                err = []
        if xbins[-1] != h.GetBinLowEdge(h.GetNbinsX()+1):
            xbins[-1] = h.GetBinLowEdge(h.GetNbinsX()+1)

        # Rebin #
        h_rebin = h.Rebin(len(xbins)-1,'rebin',np.array(xbins))

        # Rescale bins (rebin adds the bin while we need the average) #
        base_width = h.GetBinLowEdge(i+1)-h.GetBinLowEdge(i)
        for i in range(1,h_rebin.GetNbinsX()+1):
            width = h_rebin.GetBinLowEdge(i+1)-h_rebin.GetBinLowEdge(i)
            h_rebin.SetBinContent(i,h_rebin.GetBinContent(i)*base_width/width)
            h_rebin.SetBinError(i,h_rebin.GetBinError(i)*base_width/width)

        return h_rebin


    def plotWeights1D(self):
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
        amax = max([self.histograms[k].GetMaximum() for k in ['ZVeto_0b','ZPeak_0b','ZPeak_2b','ZPeak_1b','ZPeak_1b']])

        pad1.cd()

        self.histograms['ZVeto_0b'].SetMaximum(amax*1.1)
        self.histograms['ZVeto_0b'].SetMinimum(0)
        self.histograms['ZVeto_0b'].GetYaxis().SetTitle("Normalized events")
        self.histograms['ZVeto_0b'].GetYaxis().SetTitleSize(0.035)
        self.histograms['ZVeto_0b'].GetYaxis().SetTitleOffset(1.30)
        self.histograms['ZVeto_0b'].GetYaxis().SetLabelSize(0.02)
        self.histograms['ZVeto_0b'].GetXaxis().SetTitleSize(0.)
        self.histograms['ZVeto_0b'].GetXaxis().SetLabelSize(0.)
        self.histograms['ZVeto_0b'].SetTitle(self.title)

        self.histograms['ZVeto_0b'].SetLineWidth(2) 
        self.histograms['ZVeto_1b'].SetLineWidth(2) 
        self.histograms['ZVeto_2b'].SetLineWidth(2) 
        self.histograms['ZPeak_0b'].SetLineWidth(2) 
        self.histograms['ZPeak_1b'].SetLineWidth(2) 
        self.histograms['ZPeak_2b'].SetLineWidth(2) 
        self.histograms['ZVeto_0b'].SetLineColor(418) 
        self.histograms['ZPeak_0b'].SetLineColor(602)
        self.histograms['ZPeak_1b'].SetLineColor(619)
        self.histograms['ZPeak_2b'].SetLineColor(634)


        self.histograms['ZVeto_0b'].Draw("hist")
        self.histograms['ZPeak_0b'].Draw("same hist")
        self.histograms['ZPeak_1b'].Draw("same hist")
        self.histograms['ZPeak_2b'].Draw("same hist")

        leg1 = ROOT.TLegend(0.50,0.73,0.89,0.91)
        leg1.SetHeader("Legend","C")
        leg1.SetTextSize(0.02)
        mode = "DY data estimation" if self.mode == "data" else "DY MC"
        leg1.AddEntry(self.histograms['ZVeto_0b'],"#splitline{Z Veto 0 btag (%s)}{Integral = %0.2f}"%(mode,self.N_ZVeto0b))
        leg1.AddEntry(self.histograms['ZPeak_0b'],"#splitline{Z Peak 0 btag (%s)}{Integral = %0.2f}"%(mode,self.N_ZPeak0b))
        leg1.AddEntry(self.histograms['ZPeak_1b'],"#splitline{Z Peak 1 btag (%s)}{Integral = %0.2f}"%(mode,self.N_ZPeak1b))
        leg1.AddEntry(self.histograms['ZPeak_2b'],"#splitline{Z Peak 2 btag (%s)}{Integral = %0.2f}"%(mode,self.N_ZPeak2b))
        leg1.Draw()

        text = ROOT.TPaveText(0.6,0.65,0.85,0.72,"NB NDC")
        text.SetFillStyle(1001)
        text.AddText("N_{1b}/N_{0b} = %0.5f"%self.factor_1b)
        text.AddText("N_{2b}/N_{0b} = %0.5f"%self.factor_2b)
        text.Draw()


        ##########  PAD 2 ##########
        pad2.cd()

        #self.weight_1b.SetMaximum(max(self.weight_1b.GetMaximum(),self.weight_2b.GetMaximum())*1.1)
        self.weight_1b.SetMaximum(0.2)
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

        max_shape = max([h.GetMaximum() for h in [self.histograms['ZVeto_1b'],self.histograms['ZVeto_2b'],self.shape_1b,self.shape_2b]])

        self.shape_1b.SetLineWidth(2) 
        self.shape_2b.SetLineWidth(2) 
        self.shape_1b.SetLineColor(602) 
        self.shape_2b.SetLineColor(600) 

        self.histograms['ZVeto_1b'].SetLineWidth(1)
        self.histograms['ZVeto_2b'].SetLineWidth(1)
        self.histograms['ZVeto_1b'].SetLineColor(602)
        self.histograms['ZVeto_2b'].SetLineColor(600)

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
        self.histograms['ZVeto_1b'].Draw("H same")
        self.histograms['ZVeto_2b'].Draw("H same")

        leg3 = ROOT.TLegend(0.5,0.40,0.89,0.93)
        leg3.SetHeader("Legend","C")
        leg3.SetTextSize(0.04)
        leg3.AddEntry(self.shape_1b,"#splitline{DY shape from data (1b)}{Integral = %0.2f}"%self.shape_1b.Integral() if self.mode == 'data' else "#splitline{DY shape from MC (1b)}{Integral = %0.2f}"%self.shape_1b.Integral())
        leg3.AddEntry(self.histograms['ZVeto_1b'],"#splitline{Z Veto (1b) (DY MC)}{Integral = %0.2f}"%self.N_ZVeto1b)
        leg3.AddEntry(self.shape_2b,"#splitline{DY shape from data (2b)}{Integral = %0.2f}"%self.shape_2b.Integral() if self.mode == 'data' else "#splitline{DY shape from MC (2b)}{Integral = %0.2f}"%self.shape_2b.Integral())
        leg3.AddEntry(self.histograms['ZVeto_2b'],"#splitline{Z Veto (2b) (DY MC)}{Integral = %0.2f}"%self.N_ZVeto2b)
        leg3.Draw()



        C.Print(self.outputname+'_%s.pdf'%self.era)


    def plotWeights2D(self):
        ROOT.gStyle.SetTitleW(1.)
        C1 = ROOT.TCanvas("c1", "c1", 4000, 4000)
        C1.Divide(3,3)

        mode = "DY data estimation" if self.mode == "data" else "DY MC"

#        print ('-'*80)
#        print ("Examining ZVeto_0b")
#        self.findNegativeBinsInTh2(self.histograms['ZVeto_0b'])
#        print ('-'*80)
#        print ("Examining ZPeak_0b")
#        self.findNegativeBinsInTh2(self.histograms['ZPeak_0b'])
#        print ('-'*80)
#        print ("Examining ZPeak_2b")
#        self.findNegativeBinsInTh2(self.histograms['ZPeak_2b'])
#        print ('-'*80)
#        print ("Examining ZVeto_0b")
#        self.findNegativeBinsInTh2(self.histograms['ZPeak_1b'])
#        print ('-'*80)



        amax = max([self.histograms[k].GetMaximum() for k in ['ZVeto_0b','ZPeak_0b','ZPeak_2b','ZPeak_1b']])

        #----- Zpeak and Zveto 0b -----#
        C1.cd(1)
        ROOT.gPad.SetRightMargin(0.15)
        self.histograms['ZPeak_1b'].SetTitle('Z Peak 1 btag (%s) : Integral = %0.2f;Lead lepton #eta;Lead lepton P_{T}'%(mode,self.N_ZPeak1b))
        self.histograms['ZPeak_1b'].SetTitleSize(0.9,"t")
        #self.histograms['ZPeak_1b'].GetZaxis().SetRangeUser(0,amax)
        self.histograms['ZPeak_1b'].Draw("colz")
        C1.cd(7)
        ROOT.gPad.SetRightMargin(0.15)
        self.histograms['ZPeak_2b'].SetTitle('Z Peak 2 btag (%s) : Integral = %0.2f;Lead lepton #eta;Lead lepton P_{T}'%(mode,self.N_ZPeak2b))
        #self.histograms['ZPeak_2b'].GetZaxis().SetRangeUser(0,amax)
        self.histograms['ZPeak_2b'].Draw("colz")
        C1.cd(5)
        ROOT.gPad.SetRightMargin(0.15)
        self.histograms['ZVeto_0b'].Draw("colz")
        #self.histograms['ZVeto_0b'].GetZaxis().SetRangeUser(0,amax)
        self.histograms['ZVeto_0b'].SetTitle('Z Veto 0 btag (%s) : Integral = %0.2f;Lead lepton #eta;Lead lepton P_{T}'%(mode,self.N_ZVeto0b))
        C1.cd(4)
        ROOT.gPad.SetRightMargin(0.15)
        self.histograms['ZPeak_0b'].SetTitle('Z Peak 0 btag (%s) : Integral = %0.2fLead lepton #eta;Lead lepton P_{T};'%(mode,self.N_ZPeak0b))
        #self.histograms['ZPeak_0b'].GetZaxis().SetRangeUser(0,amax)
        self.histograms['ZPeak_0b'].Draw("colz")

        #----- Weights -----#
        C1.cd(2)
        ROOT.gPad.SetRightMargin(0.15)
        #ROOT.gPad.SetLogz()
        self.weight_1b.SetTitle('Weight (1b);Lead lepton #eta;Lead lepton P_{T}')
        self.weight_1b.GetZaxis().SetRangeUser(0,0.3)
        self.weight_1b.Draw("colz")
        C1.cd(8)
        ROOT.gPad.SetRightMargin(0.15)
        #ROOT.gPad.SetLogz()
        self.weight_2b.SetTitle('Weight (2b)Lead lepton #eta;Lead lepton P_{T};')
        self.weight_2b.GetZaxis().SetRangeUser(0,0.1)
        self.weight_2b.Draw("colz")

        #----- DY shape and MC DY in Zveto -----#
        #max_shape = max([h.GetMaximum() for h in [self.histograms['ZVeto_1b'],self.histograms['ZVeto_2b'],self.shape_1b,self.shape_2b]])
        mode = "DY shape from data" if self.mode == "data" else "DY shape from MC"

        C1.cd(3)
        ROOT.gPad.SetRightMargin(0.15)
        self.shape_1b.SetTitle(mode+' (Z Veto 1b) : Integral = %0.2f'%self.shape_1b.Integral()+';Lead lepton P_{T};Lead lepton #eta')
        self.shape_1b.Draw("colz")
        C1.cd(9)
        ROOT.gPad.SetRightMargin(0.15)
        self.shape_2b.SetTitle(mode+' (Z Veto 2b) : Integral = %0.2f'%self.shape_2b.Integral()+';Lead lepton P_{T};Lead lepton #eta')
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
        max_1b = max([self.histograms['ZVeto_1b'].GetMaximum(),self.shape_1b.GetMaximum()])
        max_2b = max([self.histograms['ZVeto_2b'].GetMaximum(),self.shape_2b.GetMaximum()])

        C2 = ROOT.TCanvas("c2", "c2", 4000, 4000)
        C2.Divide(2,2)

        C2.cd(1)
        ROOT.gPad.SetRightMargin(0.15)
        self.shape_1b.GetZaxis().SetRangeUser(0,max_1b)
        self.shape_1b.Draw("colz")
        C2.cd(2)
        ROOT.gPad.SetRightMargin(0.15)
        self.histograms['ZVeto_1b'].SetTitle('DY MC (Z Veto 1b) : Integral = %0.2f'%self.N_ZVeto1b)
        self.histograms['ZVeto_1b'].GetZaxis().SetRangeUser(0,max_1b)
        self.histograms['ZVeto_1b'].Draw("colz")
        C2.cd(3)
        ROOT.gPad.SetRightMargin(0.15)
        self.shape_2b.GetZaxis().SetRangeUser(0,max_2b)
        self.shape_2b.Draw("colz")
        C2.cd(4)
        ROOT.gPad.SetRightMargin(0.15)
        self.histograms['ZVeto_2b'].SetTitle('DY MC (Z Veto 2b) : Integral = %0.2f'%self.N_ZVeto2b)
        self.histograms['ZVeto_2b'].GetZaxis().SetRangeUser(0,max_2b)
        self.histograms['ZVeto_2b'].Draw("colz")
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
                for j in range(1,weight.GetNbinsY()+1):
                    bin_dict = {'values':[]}
                    etaLow = yAxis.GetBinLowEdge(j)
                    etaHigh = yAxis.GetBinUpEdge(j)
                    bin_dict['bin'] = [etaLow,etaHigh]
                    if len(eta_binning) == 0:
                        eta_binning.append(etaLow)
                    eta_binning.append(etaHigh)
                    # Loop over PT bins in X #
                    for i in range(1,weight.GetNbinsX()+1):
                        if weight.GetBinContent(i,j) == 0:
                            continue
                        PtLow  = xAxis.GetBinLowEdge(i)
                        PtHigh = xAxis.GetBinUpEdge(i)
                        bin_dict['values'].append({'bin' : [PtLow,PtHigh],
                                                   'value' : weight.GetBinContent(i,j),
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
    #----- Path -----#
    path_ZVeto = '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/full2016_AutoPULeptonTriggerJetMETBtagSF_tight_withRatio_noSyst_v2'
    path_ZPeak = '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/full2016_AutoPULeptonTriggerJetMETBtagSF_tight_withRatio_ZPeak_noSyst_v2'

    #----- 1D weight -----#
    ElElChannel = {'ZVeto_0b' : {'path': path_ZVeto,
                                 'histname': 'ElEl_HasElElTightTwoAk4JetsExclusiveResolvedNoBtag_firstlepton_pt'},
                   'ZVeto_1b' : {'path': path_ZVeto,
                                 'histname': 'ElEl_HasElElTightTwoAk4JetsExclusiveResolvedOneBtag_firstlepton_pt'},
                   'ZVeto_2b' : {'path': path_ZVeto,
                                 'histname': 'ElEl_HasElElTightTwoAk4JetsExclusiveResolvedTwoBtags_firstlepton_pt'},
                   'ZPeak_0b' : {'path': path_ZPeak,
                                 'histname': 'ElEl_HasElElTightZPeakTwoAk4JetsExclusiveResolvedNoBtag_firstlepton_pt'},
                   'ZPeak_1b' : {'path': path_ZPeak,
                                 'histname': 'ElEl_HasElElTightZPeakTwoAk4JetsExclusiveResolvedOneBtag_firstlepton_pt'},
                   'ZPeak_2b' : {'path': path_ZPeak,
                                 'histname': 'ElEl_HasElElTightZPeakTwoAk4JetsExclusiveResolvedTwoBtags_firstlepton_pt'}}

    MuMuChannel = {'ZVeto_0b' : {'path': path_ZVeto,
                                 'histname': 'MuMu_HasMuMuTightTwoAk4JetsExclusiveResolvedNoBtag_firstlepton_pt'},
                   'ZVeto_1b' : {'path': path_ZVeto,
                                 'histname': 'MuMu_HasMuMuTightTwoAk4JetsExclusiveResolvedOneBtag_firstlepton_pt'},
                   'ZVeto_2b' : {'path': path_ZVeto,
                                 'histname': 'MuMu_HasMuMuTightTwoAk4JetsExclusiveResolvedTwoBtags_firstlepton_pt'},
                   'ZPeak_0b' : {'path': path_ZPeak,
                                 'histname': 'MuMu_HasMuMuTightZPeakTwoAk4JetsExclusiveResolvedNoBtag_firstlepton_pt'},
                   'ZPeak_1b' : {'path': path_ZPeak,
                                 'histname': 'MuMu_HasMuMuTightZPeakTwoAk4JetsExclusiveResolvedOneBtag_firstlepton_pt'},
                   'ZPeak_2b' : {'path': path_ZPeak,
                                 'histname': 'MuMu_HasMuMuTightZPeakTwoAk4JetsExclusiveResolvedTwoBtags_firstlepton_pt'}}


#    instance = WeightDY(ElElChannel,'First lepton P_{T} (e^{+}e^{-} channel)','weight_ElEl_data_1D','data','2016')
#    instance = WeightDY(MuMuChannel,'First lepton P_{T} (#mu^{+}#mu^{-} channel)','weight_MuMu_data_1D','data','2016')
#    instance = WeightDY(ElElChannel,'First lepton P_{T} (e^{+}e^{-} channel)','weight_ElEl_mc_1D','mc','2016')
#    instance = WeightDY(MuMuChannel,'First lepton P_{T} (#mu^{+}#mu^{-} channel)','weight_MuMu_mc_1D','mc','2016')

    #----- 2D weight -----#
    ElElChannel = {'ZVeto_0b' : {'path': path_ZVeto,
                                 'histname': 'ElEl_HasElElTightTwoAk4JetsExclusiveResolvedNoBtag_firstlepton_ptVSeta'},
                   'ZVeto_1b' : {'path': path_ZVeto,
                                 'histname': 'ElEl_HasElElTightTwoAk4JetsExclusiveResolvedOneBtag_firstlepton_ptVSeta'},
                   'ZVeto_2b' : {'path': path_ZVeto,
                                 'histname': 'ElEl_HasElElTightTwoAk4JetsExclusiveResolvedTwoBtags_firstlepton_ptVSeta'},
                   'ZPeak_0b' : {'path': path_ZPeak,
                                 'histname': 'ElEl_HasElElTightZPeakTwoAk4JetsExclusiveResolvedNoBtag_firstlepton_ptVSeta'},
                   'ZPeak_1b' : {'path': path_ZPeak,
                                 'histname': 'ElEl_HasElElTightZPeakTwoAk4JetsExclusiveResolvedOneBtag_firstlepton_ptVSeta'},
                   'ZPeak_2b' : {'path': path_ZPeak,
                                 'histname': 'ElEl_HasElElTightZPeakTwoAk4JetsExclusiveResolvedTwoBtags_firstlepton_ptVSeta'}}
    MuMuChannel = {'ZVeto_0b' : {'path': path_ZVeto,
                                 'histname': 'MuMu_HasMuMuTightTwoAk4JetsExclusiveResolvedNoBtag_firstlepton_ptVSeta'},
                   'ZVeto_1b' : {'path': path_ZVeto,
                                 'histname': 'MuMu_HasMuMuTightTwoAk4JetsExclusiveResolvedOneBtag_firstlepton_ptVSeta'},
                   'ZVeto_2b' : {'path': path_ZVeto,
                                 'histname': 'MuMu_HasMuMuTightTwoAk4JetsExclusiveResolvedTwoBtags_firstlepton_ptVSeta'},
                   'ZPeak_0b' : {'path': path_ZPeak,
                                 'histname': 'MuMu_HasMuMuTightZPeakTwoAk4JetsExclusiveResolvedNoBtag_firstlepton_ptVSeta'},
                   'ZPeak_1b' : {'path': path_ZPeak,
                                 'histname': 'MuMu_HasMuMuTightZPeakTwoAk4JetsExclusiveResolvedOneBtag_firstlepton_ptVSeta'},
                   'ZPeak_2b' : {'path': path_ZPeak,
                                 'histname': 'MuMu_HasMuMuTightZPeakTwoAk4JetsExclusiveResolvedTwoBtags_firstlepton_ptVSeta'}}

    instance = WeightDY(ElElChannel,'First lepton P_{T} (e^{+}e^{-} channel)','weight_ElEl_data_2D','data','2016')
    instance = WeightDY(MuMuChannel,'First lepton P_{T} (#mu^{+}#mu^{-} channel)','weight_MuMu_data_2D','data','2016')
    instance = WeightDY(ElElChannel,'First lepton P_{T} (e^{+}e^{-} channel)','weight_ElEl_mc_2D','mc','2016')
    instance = WeightDY(MuMuChannel,'First lepton P_{T} (#mu^{+}#mu^{-} channel)','weight_MuMu_mc_2D','mc','2016')
