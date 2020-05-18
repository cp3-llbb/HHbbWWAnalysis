import os
import sys
import yaml
import copy
import glob
import json
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

class WeightDY:
    def __init__(self, config, title, outputname, mode, era):
        self.config     = config
        self.mode       = mode
        self.era        = era
        self.title      = title
        self.outputname = outputname

        self.histograms = {}
        for cat, inp in self.config.items():
            dirInfo = self.extractDirInformation(inp['path'],inp['histname'])
            self.histograms[cat] = self.produceDistribution(cat,dirInfo)

        self.produceShape(self.histograms)
        self.plotWeights()
        self.saveToJson()

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

        mode = "mc" if category == "ZVeto_2b" else self.mode

        dist = None
        # Loop over hist per sample #
        for sample, data_dict in samp_dict.items():
            # Get hist #
            h = hist_dict[sample]
            if h is None:
                continue
            if dist is None:
                dist = ROOT.TH1F(category,category,h.GetNbinsX(),h.GetBinLowEdge(1),h.GetBinLowEdge(h.GetNbinsX())+h.GetBinWidth(h.GetNbinsX()))
            # Add histograms depending on type #
            if mode == 'data':
                if data_dict['type'] == 'data':
                    dist.Add(h)
                if data_dict['type'] == 'mc':
                    if data_dict['group'] != 'DY':
                        factor = data_dict["cross-section"]*lumi_dict[self.era]/data_dict["generated-events"]    
                        if data_dict['group'] == 'ttbar':
                            factor *= 0.7
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
        dist.Sumw2()
        return dist

    def produceShape(self,hist_dict):
        N_0b = hist_dict['ZPeak_0b'].Integral()
        N_2b = hist_dict['ZPeak_2b'].Integral()
        self.factor = N_2b/N_0b

        print ("Sum weight ZVeto_0b : ",hist_dict['ZVeto_0b'].Integral())
        print ("Sum weight ZVeto_2b : ",hist_dict['ZVeto_2b'].Integral())
        print ("Sum weight ZPeak_0b : ",hist_dict['ZPeak_0b'].Integral())
        print ("Sum weight ZPeak_2b : ",hist_dict['ZPeak_2b'].Integral())

        hist_dict['ZPeak_2b'].Scale(1./hist_dict['ZPeak_2b'].Integral())
        hist_dict['ZPeak_0b'].Scale(1./hist_dict['ZPeak_0b'].Integral())

        self.ratio = hist_dict['ZPeak_2b'].Clone("Ratio")
        self.ratio.Divide(hist_dict['ZPeak_0b'])
        self.ratio.Scale(self.factor)
            # self.ratio = ZPeak_2b / ZPeak_0b (normalized histogram ratio) * N_2b/N_0b (stats factor)

        
        self.shape = hist_dict['ZVeto_0b'].Clone("Weights")
        self.shape.Multiply(self.ratio)
            # self.shape = (ZPeak_2b / ZPeak_0b) * ZVeto_0b (unnormalized histogram) = shape of the data DY in SR

        hist_dict['ZVeto_0b'].Scale(1./hist_dict['ZVeto_0b'].Integral())

        print ("Factor              : ",self.factor)
        print ("Weights             : ",self.shape.Integral())

    def plotWeights(self):
        C = ROOT.TCanvas("c1", "c1", 600, 900)

        pad1 = ROOT.TPad("pad1", "pad1", 0., 0.0, 1., 1.0)
        pad1.SetTopMargin(0.08)
        pad1.SetBottomMargin(0.61)
        pad1.SetLeftMargin(0.12)
        pad1.SetRightMargin(0.1)
        pad1.SetGridx()
        pad1.SetGridy()
        
        pad2 = ROOT.TPad("pad2", "pad2", 0, 0.00, 1, 0.6)
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

        amax = max([self.histograms[k].GetMaximum() for k in ['ZVeto_0b','ZPeak_0b','ZPeak_2b']])

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
        self.histograms['ZVeto_2b'].SetLineWidth(2) 
        self.histograms['ZPeak_0b'].SetLineWidth(2) 
        self.histograms['ZPeak_2b'].SetLineWidth(2) 
        self.histograms['ZVeto_0b'].SetLineColor(418) 
        self.histograms['ZPeak_0b'].SetLineColor(602)
        self.histograms['ZPeak_2b'].SetLineColor(634)


        self.histograms['ZVeto_0b'].Draw("hist")
        self.histograms['ZPeak_0b'].Draw("same hist")
        self.histograms['ZPeak_2b'].Draw("same hist")

        leg1 = ROOT.TLegend(0.50,0.80,0.89,0.91)
        leg1.SetHeader("Legend","C")
        leg1.SetTextSize(0.02)
        mode = "DY data estimation" if self.mode == "data" else "DY MC"
        leg1.AddEntry(self.histograms['ZVeto_0b'],"Z Veto 0 btag (%s)"%mode)
        leg1.AddEntry(self.histograms['ZPeak_0b'],"Z Peak 0 btag (%s)"%mode)
        leg1.AddEntry(self.histograms['ZPeak_2b'],"Z Peak 2 btag (%s)"%mode)
        leg1.Draw()

        pad2.cd()

        #self.ratio.SetMaximum(3.)
        #self.ratio.SetMinimum(0.)
        self.ratio.SetLineWidth(2)
        self.ratio.SetTitle("")

        self.ratio.GetYaxis().SetTitle("Weight")
        self.ratio.GetYaxis().SetTitleSize(0.035)
        self.ratio.GetYaxis().SetTitleOffset(1.30)
        self.ratio.GetYaxis().SetLabelSize(0.02)

        self.ratio.GetXaxis().SetTitleSize(0.)
        self.ratio.GetXaxis().SetLabelSize(0.)

        self.ratio.Draw("ep")


        pad3.cd()

        self.shape.SetLineWidth(2) 
        self.shape.SetLineColor(1) 
        self.shape.SetMinimum(0)
        self.shape.SetMaximum(max(self.histograms['ZVeto_2b'].GetMaximum(),self.shape.GetMaximum())*1.1)
        self.shape.SetTitle("")
        self.shape.GetYaxis().SetTitle("DY shape")
        self.shape.GetYaxis().SetTitleSize(0.08)
        self.shape.GetYaxis().SetTitleOffset(0.7)
        self.shape.GetYaxis().SetLabelSize(0.05)

        self.shape.GetXaxis().SetTitleSize(0.08)
        self.shape.GetXaxis().SetTitleOffset(1.15)
        self.shape.GetXaxis().SetLabelSize(0.05)

        self.shape.Draw("ep")
        self.histograms['ZVeto_2b'].Draw("hist same")

        leg3 = ROOT.TLegend(0.5,0.75,0.89,0.89)
        leg3.SetHeader("Legend","C")
        leg3.AddEntry(self.shape,"DY shape from data" if self.mode == 'data' else "DY shape from MC")
        leg3.AddEntry(self.histograms['ZVeto_2b'],"Z Veto 2 btag (DY MC)")
        leg3.Draw()

        text = ROOT.TPaveText(0.6,0.65,0.85,0.70,"NB NDC")
        text.SetFillStyle(1001)
        text.AddText("N_{2b}/N_{0b} = %0.5f"%self.factor)
        text.Draw()



        C.Print(self.outputname+'_%s.pdf'%self.era)

    def saveToJson(self):
        json_dict = {'dimension': 2, 'variables': ['AbsEta','Pt'], 'error_type': 'absolute'}
        # Fill data #
        binning = []
        data = []
        xAxis = self.ratio.GetXaxis()
        for i in range(1,self.ratio.GetNbinsX()+1):
            # Get content #
            cont = self.ratio.GetBinContent(i)
            if cont <= 0:
                continue
            # Get error #
            err = self.ratio.GetBinError(i)

            # Get binning #
            xLow  = xAxis.GetBinLowEdge(i)
            xHigh = xAxis.GetBinUpEdge(i)
            
            if len(binning) == 0:
                binning.append(xLow)
            binning.append(xHigh)

            # Save data #
            data.append({'bin':[xLow,xHigh],'value':cont,'error_low':err,'error_high':err}) 
       
        json_dict['binning'] = {'x':[0,2.5],'y':binning}
        json_dict['data'] = [{'bin':[0,2.5],'values':data}]

        with open(self.outputname+'_%s.json'%self.era,'w') as f:
            json.dump(json_dict,f,indent=4)
    
        
        
        

ElElChannel = {'ZVeto_0b' : {'path': '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/full2016_AutoPULeptonTriggerJetMETBtagSF_DYStudy0And2Btag_ZVeto/',
                             'histname': 'ElEl_HasElElTightTwoAk4JetsExclusiveResolvedNoBtag_firstlepton_pt'},
               'ZVeto_2b' : {'path': '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/full2016_AutoPULeptonTriggerJetMETBtagSF_DYStudy0And2Btag_ZVeto/',
                             'histname': 'ElEl_HasElElTightTwoAk4JetsExclusiveResolvedTwoBtags_firstlepton_pt'},
               'ZPeak_0b' : {'path': '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/full2016_AutoPULeptonTriggerJetMETBtagSF_DYStudy0And2Btag_inZpeak/',
                             'histname': 'ElEl_HasElElTightZPeakTwoAk4JetsExclusiveResolvedNoBtag_firstlepton_pt'},
               'ZPeak_2b' : {'path': '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/full2016_AutoPULeptonTriggerJetMETBtagSF_DYStudy0And2Btag_inZpeak/',
                             'histname': 'ElEl_HasElElTightZPeakTwoAk4JetsExclusiveResolvedTwoBtags_firstlepton_pt'}}

MuMuChannel = {'ZVeto_0b' : {'path': '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/full2016_AutoPULeptonTriggerJetMETBtagSF_DYStudy0And2Btag_ZVeto/',
                             'histname': 'MuMu_HasMuMuTightTwoAk4JetsExclusiveResolvedNoBtag_firstlepton_pt'},
               'ZVeto_2b' : {'path': '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/full2016_AutoPULeptonTriggerJetMETBtagSF_DYStudy0And2Btag_ZVeto/',
                             'histname': 'MuMu_HasMuMuTightTwoAk4JetsExclusiveResolvedTwoBtags_firstlepton_pt'},
               'ZPeak_0b' : {'path': '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/full2016_AutoPULeptonTriggerJetMETBtagSF_DYStudy0And2Btag_inZpeak/',
                             'histname': 'MuMu_HasMuMuTightZPeakTwoAk4JetsExclusiveResolvedNoBtag_firstlepton_pt'},
               'ZPeak_2b' : {'path': '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/full2016_AutoPULeptonTriggerJetMETBtagSF_DYStudy0And2Btag_inZpeak/',
                             'histname': 'MuMu_HasMuMuTightZPeakTwoAk4JetsExclusiveResolvedTwoBtags_firstlepton_pt'}}


instance = WeightDY(ElElChannel,'First lepton P_{T} (e^{+}e^{-} channel)','weight_ElEl_data','data','2016')
instance = WeightDY(MuMuChannel,'First lepton P_{T} (#mu^{+}#mu^{-} channel)','weight_MuMu_data','data','2016')
instance = WeightDY(ElElChannel,'First lepton P_{T} (e^{+}e^{-} channel)','weight_ElEl_mc','mc','2016')
instance = WeightDY(MuMuChannel,'First lepton P_{T} (#mu^{+}#mu^{-} channel)','weight_MuMu_mc','mc','2016')
