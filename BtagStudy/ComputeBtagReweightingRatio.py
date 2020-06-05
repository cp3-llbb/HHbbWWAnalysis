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

class BtagReweightingRatio:
    def __init__(self, config, title, outputname, era, select_group=None):
        self.config     = config
        self.era        = era
        self.title      = title
        self.outputname = outputname
        self.select_group = select_group
        if self.select_group:
            self.outputname += '_'+self.select_group

        self.histograms = {}
        for cat, inp in self.config.items():
            dirInfo = self.extractDirInformation(inp['path'],inp['histname'])
            self.histograms[cat] = self.produceDistribution(cat,dirInfo)

        self.plot()
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

        dist = None
        # Loop over hist per sample #
        for sample, data_dict in samp_dict.items():
            # Get hist #
            h = hist_dict[sample]
            if h is None:
                continue
            if dist is None:
                dist = ROOT.TH1D(category,category,h.GetNbinsX(),h.GetBinLowEdge(1),h.GetBinLowEdge(h.GetNbinsX())+h.GetBinWidth(h.GetNbinsX()))
                dist.Sumw2()
            # Add histograms (only mc)# 
            if data_dict['type'] == 'mc':
                if self.select_group:
                    if data_dict['group'] != self.select_group:
                        continue
                #factor = data_dict["cross-section"]*lumi_dict[self.era]/data_dict["generated-events"]     
                factor = 1 # We want event weight sum
                dist.Add(h,factor)  
        dist.GetXaxis().SetTitle(plot_dict['x-axis']) 
        dist.GetYaxis().SetTitle(plot_dict['y-axis']) 
        return dist

    def plot(self):
        C = ROOT.TCanvas("c1", "c1", 600, 600)
        C.Divide(1,2)
        colors = [418,602]
        leg = ROOT.TLegend(0.65,0.70,0.95,0.90)
        leg.SetBorderSize(0)
        leg.SetTextSize(0.05)

        C.cd(1)
        ROOT.gPad.SetRightMargin(0.05)
        ROOT.gPad.SetLeftMargin(0.15)
        ROOT.gPad.SetTopMargin(0.05)
        ROOT.gPad.SetBottomMargin(0.01)
        ROOT.gPad.SetLogy()

        ROOT.gStyle.SetTitleFontSize(0.1)

        BtagAfter = self.histograms['After btag weight ElEl'].Clone("BtagAfter")
        BtagAfter.Add(self.histograms['After btag weight MuMu'])
        BtagAfter.Add(self.histograms['After btag weight ElMu'])

        BtagBefore = self.histograms['Before btag weight ElEl'].Clone("BtagBefore")
        BtagBefore.Add(self.histograms['Before btag weight MuMu'])
        BtagBefore.Add(self.histograms['Before btag weight ElMu'])

        BtagAfter.SetTitle(self.title)
        BtagAfter.GetYaxis().SetTitle("Event weight sum")
        BtagAfter.GetYaxis().SetTitleOffset(1.)
        BtagAfter.GetYaxis().SetLabelSize(0.05)
        BtagAfter.GetYaxis().SetTitleSize(0.06)
        BtagAfter.GetXaxis().SetTitle("")
        BtagAfter.GetXaxis().SetLabelSize(0.)

        BtagAfter.SetLineWidth(2)
        BtagAfter.SetLineColor(418)
        BtagBefore.SetLineWidth(2)
        BtagBefore.SetLineColor(602)

        hmax = max([BtagAfter.GetMaximum(),BtagBefore.GetMaximum()])
        BtagAfter.SetMaximum(hmax*1.5)

        leg.AddEntry(BtagAfter,'After btag weight')
        leg.AddEntry(BtagBefore,'Before btag weight')

        BtagAfter.Draw("H")
        BtagBefore.Draw("H same")
    
        leg.Draw()
    
        C.cd(2)
        ROOT.gPad.SetRightMargin(0.05)
        ROOT.gPad.SetLeftMargin(0.15)
        ROOT.gPad.SetTopMargin(0.01)
        ROOT.gPad.SetBottomMargin(0.15)
        ROOT.gPad.SetGridy()

        self.ratio = BtagBefore.Clone("Ratio")
        self.ratio.Divide(BtagAfter)

        self.ratio.SetTitle("")
        self.ratio.SetLineWidth(2)
        self.ratio.SetLineColor(1)

        self.ratio.GetYaxis().SetTitle("#frac{#sum w(before)}{#sum w(after)}")
        self.ratio.GetYaxis().SetLabelSize(0.05)
        self.ratio.GetYaxis().SetTitleSize(0.06)
        self.ratio.GetYaxis().SetRangeUser(0.5,1.5)
        self.ratio.GetYaxis().SetTitleOffset(1.)

        #self.ratio.GetXaxis().SetTitle("N selected Ak4 jets")
        self.ratio.GetXaxis().SetLabelSize(0.05)
        self.ratio.GetXaxis().SetTitleSize(0.06)

        self.ratio.Draw("")

        C.Print(self.outputname+'_%s.pdf'%(self.era))


    def saveToJson(self):
        json_dict = {'dimension': 1, 'variables': ['NumJets'], 'error_type': 'absolute'}
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
       
        json_dict['binning'] = {'x':binning}
        json_dict['data'] = data

        with open(self.outputname+'_%s.json'%(self.era),'w') as f:
            json.dump(json_dict,f,indent=4)
        print ("Saved json to",self.outputname+'_%s.json'%(self.era))

            
    def saveToRoot(self):
        root_file = ROOT.TFile(self.outputname+'_%s.root'%self.era,"recreate")
        self.ratio.Write("ratio")
        root_file.Write()
        root_file.Close()
        
        
        
if __name__ == "__main__":
    #----- Path -----#
    path_reweighting_on = '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/full2016_BtagReweighting_On_v2/'
    path_reweighting_off = '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/full2016_BtagReweighting_Off_v2/'

    #----- Instantiantion -----#
    #ElElChannel = {'Before btag weight' : {'path': path_reweighting_off,
    #                                   'histname': 'ElEl_HasElElTightTwoAk4Jets_PreSelAk4Jets_N'},
    #               'After btag weight' : {'path': path_reweighting_on,
    #                                   'histname': 'ElEl_HasElElTightTwoAk4Jets_PreSelAk4Jets_N'}}
    #MuMuChannel = {'Before btag weight' : {'path': path_reweighting_off,
    #                                   'histname': 'MuMu_HasMuMuTightTwoAk4Jets_PreSelAk4Jets_N'},
    #               'After btag weight' : {'path': path_reweighting_on,
    #                                   'histname': 'MuMu_HasMuMuTightTwoAk4Jets_PreSelAk4Jets_N'}}
    #ElMuChannel = {'Before btag weight' : {'path': path_reweighting_off,
    #                                   'histname': 'ElMu_HasElMuTightTwoAk4Jets_PreSelAk4Jets_N'},
    #               'After btag weight' : {'path': path_reweighting_on,
    #                                   'histname': 'ElMu_HasElMuTightTwoAk4Jets_PreSelAk4Jets_N'}}
    AllChannels = {'Before btag weight ElEl' :  {'path': path_reweighting_off,
                                                 'histname': 'ElEl_HasElElTightTwoAk4Jets_PreSelAk4Jets_N'},
                   'After btag weight ElEl' : {'path': path_reweighting_on,
                                               'histname': 'ElEl_HasElElTightTwoAk4Jets_PreSelAk4Jets_N'},
                   'Before btag weight MuMu' :  {'path': path_reweighting_off,
                                                 'histname': 'MuMu_HasMuMuTightTwoAk4Jets_PreSelAk4Jets_N'},
                   'After btag weight MuMu' : {'path': path_reweighting_on,
                                               'histname': 'MuMu_HasMuMuTightTwoAk4Jets_PreSelAk4Jets_N'},
                   'Before btag weight ElMu' :  {'path': path_reweighting_off,
                                                 'histname': 'ElMu_HasElMuTightTwoAk4Jets_PreSelAk4Jets_N'},
                   'After btag weight ElMu' : {'path': path_reweighting_on,
                                               'histname': 'ElMu_HasElMuTightTwoAk4Jets_PreSelAk4Jets_N'}}


    #instance = BtagReweightingRatio(ElElChannel,'e^{+}e^{-} channel','BtagReweightingRatio_ElEl','2016')
    #instance = BtagReweightingRatio(MuMuChannel,'#mu^{+}#mu^{-} channel','BtagReweightingRatio_MuMu','2016')
    #instance = BtagReweightingRatio(ElMuChannel,'e^{#pm}#mu^{#mp} channel','BtagReweightingRatio_ElMu','2016')
    instance = BtagReweightingRatio(AllChannels,'','BtagReweightingRatio','2016')
    instance = BtagReweightingRatio(AllChannels,'','BtagReweightingRatio','2016',select_group='ttbar')
    instance = BtagReweightingRatio(AllChannels,'','BtagReweightingRatio','2016',select_group='DY')

