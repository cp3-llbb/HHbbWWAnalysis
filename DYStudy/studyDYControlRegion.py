import os
import sys
import glob
import argparse
import copy
import yaml
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

class DYControlRegion():
    def __init__(self,variable,variable_name,output,plot_data=False,plot_MC=False,exclude_DY=False):
        self.outname = "DYPlots/"+output
        self.era     = '2016'
        self.variable = variable
        self.variable_name = variable_name
        self.data_hist_dict = {}
        self.MC_hist_dict = {}
        self.plot_data = plot_data
        self.plot_MC = plot_MC
        self.exclude_DY = exclude_DY

        self.extractAllInformation()
        self.processHistograms()
        self.saveHistograms()

    def extractAllInformation(self):
        self.info_dict = {}
        for (directory,legend),histogram in self.variable.items():
            directory = os.path.abspath(directory)
            self.info_dict[legend] = self.extractDirInformation(directory,histogram)

            # self.info_dict : dict
            #   -> key : directory basename
            #   -> val : dict:
            #       ->keys : 'luminosity', 'plot_options', 'samples', 'histograms'
            #           -> 'luminosity' : dict
            #               -> key : era
            #               -> val : lumi
            #           -> 'plot_options': dict
            #               -> key : option name 
            #               -> val : option value
            #           -> 'samples': dict
            #               -> key : name rootfile
            #               -> val : dict information (cross section, group, ...)
            #           -> 'histograms : dict
            #               -> key : name rootfile
            #               -> val : histogram 


    def extractDirInformation(self,directory,histogram):
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
        else:
            print ("Could not find hist %s in %s"%(histname,rootfile))
            h = None
        # Overbins #
        h.SetBinContent(1,h.GetBinContent(0)+h.GetBinContent(1))
        h.SetBinContent(h.GetNbinsX(),h.GetBinContent(h.GetNbinsX())+h.GetBinContent(h.GetNbinsX()+1))
        f.Close()
        return h

    def processHistograms(self):
        colors = [601,634,418,808,402,618,874,835] # Must be extended if more required
        for i,(key, val) in enumerate(self.info_dict.items()):
            # Parse #
            data_hist = None
            MC_hist = None
            hist_dict = val['histograms'] 
            lumi_dict = val['luminosity'] 
            plot_dict = val['plot_options'] 
            samp_dict = val['samples']
            # Loop over hist per sample #
            for sample, data_dict in samp_dict.items():
                # Get hist #
                h = hist_dict[sample]
                if h is None:
                    continue
                if data_hist is None:
                    data_hist = ROOT.TH1F(self.variable_name+"DataMinusMCDY",self.variable_name,h.GetNbinsX(),h.GetBinLowEdge(1),h.GetBinLowEdge(h.GetNbinsX())+h.GetBinWidth(h.GetNbinsX()))
                if MC_hist is None:
                    MC_hist = ROOT.TH1F(self.variable_name+"MCDY",self.variable_name,h.GetNbinsX(),h.GetBinLowEdge(1),h.GetBinLowEdge(h.GetNbinsX())+h.GetBinWidth(h.GetNbinsX()))

                # Add histograms depending on type #
                if 'cross-section' in data_dict.keys():
                    factor = data_dict["cross-section"]*lumi_dict[self.era]/data_dict["generated-events"]
                else:
                    factor = 1.
                if 'branching-ratio' in data_dict.keys():
                    factor *= data_dict["branching-ratio"]
                    
                if self.exclude_DY:
                    if data_dict['type'] == 'data':
                        data_hist.Add(h)
                    elif data_dict['type'] == 'mc':
                        if data_dict['group'] != 'DY':
                            data_hist.Add(h,-factor)
                        else:  
                            MC_hist.Add(h,factor)
                else:
                    if data_dict['type'] == 'data':
                        data_hist.Add(h)
                    elif data_dict['type'] == 'mc':
                        MC_hist.Add(h,factor)
                        

                    
            # Esthetics #
            if self.plot_MC and self.plot_data:
                i *= 2
            data_hist.SetLineWidth(2)
            data_hist.SetLineColor(colors[i])
            data_hist.GetXaxis().SetTitle(plot_dict['x-axis'])
            data_hist.GetYaxis().SetTitle(plot_dict['y-axis'])
            data_hist.SetTitle(self.variable_name)
            MC_hist.SetLineWidth(2)
            if self.plot_MC and self.plot_data:
                MC_hist.SetLineColor(colors[i+1])
            else:
                MC_hist.SetLineColor(colors[i])
            MC_hist.GetXaxis().SetTitle(plot_dict['x-axis'])
            MC_hist.GetYaxis().SetTitle(plot_dict['y-axis'])
            MC_hist.SetTitle(self.variable_name)
            # Normalize to unity #
            if data_hist.Integral() != 0:
                data_hist.Scale(1/data_hist.Integral())
            if MC_hist.Integral() != 0:
                MC_hist.Scale(1/MC_hist.Integral())
            # Save #
            self.data_hist_dict[key] = data_hist
            self.MC_hist_dict[key] = MC_hist

    def saveHistograms(self):
        num_plots = len(self.data_hist_dict.keys())*self.plot_data + len(self.MC_hist_dict.keys())*self.plot_MC
        plot_ratio = num_plots==2
        # Legend #
        if num_plots>=4:
            legend = ROOT.TLegend(0.2,0.8,0.89,0.89)
            legend.SetNColumns(2)
            legend.SetTextSize(0.012)
        else:
            legend = ROOT.TLegend(0.5,0.8,0.89,0.89)
            legend.SetTextSize(0.015)
        legend.SetHeader("Legend","C")

        # Canvas and pad #
        C = ROOT.TCanvas("c1", "c1", 600, 600)
        pad1 = ROOT.TPad("pad1", "pad1", 0, 0.0, 1, 1.0)
        pad1.SetBottomMargin(0.15)
        pad1.SetLeftMargin(0.15)
        pad1.SetRightMargin(0.1)
        if plot_ratio:
            pad1.SetBottomMargin(0.32)
        pad1.SetGridx()
        pad1.SetGridy()
        pad1.Draw()
        pad1.cd()

        # Get Max values #
        max_data = max([h.GetMaximum() for h in self.data_hist_dict.values()])
        max_MC = max([h.GetMaximum() for h in self.MC_hist_dict.values()])
        min_data = min([0]+[h.GetMinimum() for h in self.data_hist_dict.values()])
        min_MC = min([0]+[h.GetMinimum() for h in self.MC_hist_dict.values()])
        if self.plot_data and self.plot_MC:
            amax = max(max_data,max_MC)
            amin = min(min_data,min_MC)
        elif self.plot_data:
            amax = max_data
            amin = min_data
        elif self.plot_MC:
            amax = max_MC
            amin = min_MC

        # Plot and save #
        opt = "hist"
        for key in self.data_hist_dict.keys():
            if self.plot_data:
                hist_data = self.data_hist_dict[key]
                hist_data.SetMaximum(amax*1.2)
                hist_data.SetMinimum(amin*1.1)
                hist_data.Draw(opt)
                if opt.find("same") == -1:
                    opt += " same"
            if self.plot_MC:
                hist_MC = self.MC_hist_dict[key]
                hist_MC.SetMaximum(amax*1.2)
                hist_MC.SetMinimum(amin*1.1)
                hist_MC.Draw(opt)
                if opt.find("same") == -1:
                    opt += " same"
            if self.exclude_DY:
                if self.plot_data:
                    legend.AddEntry(hist_data,"%s : %s"%(key,'Data - MC except DY'),'l')
                if self.plot_MC:
                    legend.AddEntry(hist_MC,"%s : %s"%(key,'MC DY'),'l')
            else:
                if self.plot_data:
                    legend.AddEntry(hist_data,"%s : %s"%(key,'Data'),'l')
                if self.plot_MC:
                    legend.AddEntry(hist_MC,"%s : %s"%(key,'MC'),'l')

        legend.Draw()
        # Ratio #
        if plot_ratio:
            keys = list(self.data_hist_dict.keys())
            if self.plot_data and not self.plot_MC:
                hist1 = self.data_hist_dict[keys[0]]
                hist2 = self.data_hist_dict[keys[1]]
            elif self.plot_MC and not self.plot_data:
                hist1 = self.MC_hist_dict[keys[0]]
                hist2 = self.MC_hist_dict[keys[1]]
            elif self.plot_data and self.plot_MC:
                hist1 = self.data_hist_dict[keys[0]]
                hist2 = self.MC_hist_dict[keys[0]]
            ratio = hist1.Clone()
            ratio.Sumw2()
            ratio.Divide(hist2)

            # Redraw axis to avoid clipping 0
            hist1.GetXaxis().SetLabelSize(0.)
            hist1.GetXaxis().SetTitle('')
            
            pad2 = ROOT.TPad("pad2", "pad2", 0, 0.0, 1, 0.3)
            pad2.SetTopMargin(0)
            pad2.SetBottomMargin(0.4)
            pad2.SetLeftMargin(0.15)
            pad1.SetRightMargin(0.1)
            pad2.SetGridx()
            pad2.SetGridy()
            pad2.Draw()
            pad2.cd()

            ratio.SetLineColor(ROOT.kBlack)
            ratio.SetMinimum(0.0)
            ratio.SetMaximum(2.0)
            ratio.SetStats(0)
            ratio.SetMarkerStyle(21)
            ratio.Draw("ep")

            ratio.SetTitle("")
            ratio.GetYaxis().SetTitle("Ratio")
            ratio.GetYaxis().SetNdivisions(505)
            ratio.GetYaxis().SetTitleSize(20)
            ratio.GetYaxis().SetTitleFont(43)
            ratio.GetYaxis().SetTitleOffset(1.8)
            ratio.GetYaxis().SetLabelFont(43)
            ratio.GetYaxis().SetLabelSize(15)
            
            ratio.GetXaxis().SetNdivisions(510)
            ratio.GetXaxis().SetTitleSize(20)
            ratio.GetXaxis().SetTitleFont(43)
            ratio.GetXaxis().SetTitleOffset(4.)
            ratio.GetXaxis().SetLabelFont(43)
            ratio.GetXaxis().SetLabelSize(15)

        if self.plot_data:
            self.outname += "_data"
        if self.plot_MC:
            self.outname += "_MC"
        if self.exclude_DY:
            self.outname += "_DY"
        else:
            self.outname += "_All"
        C.Print(self.outname+".pdf")
            


if __name__ == "__main__":
    path_Zveto = '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/full2016NanoV6_Tight_Fake/'
    path_Zpeak = '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/full2016NanoV6_Tight_ZPeak_Fake/'

    # ElEl First lepton pt #
    var_0btagZveto_dict = {(path_Zveto,'Z Veto 0 btag'): 'ElEl_HasElElTightTwoAk4JetsExclusiveResolvedNoBtag_firstlepton_pt'}
    var_1btagZveto_dict = {(path_Zveto,'Z Veto 1 btag'): 'ElEl_HasElElTightTwoAk4JetsExclusiveResolvedOneBtag_firstlepton_pt'}
    var_2btagZveto_dict = {(path_Zveto,'Z Veto 2 btag'): 'ElEl_HasElElTightTwoAk4JetsExclusiveResolvedTwoBtags_firstlepton_pt'}

    var_0btagZpeak_dict = {(path_Zpeak,'Z Peak 0 btag'): 'ElEl_HasElElTightZPeakTwoAk4JetsExclusiveResolvedNoBtag_firstlepton_pt'}
    var_1btagZpeak_dict = {(path_Zpeak,'Z Peak 1 btag'): 'ElEl_HasElElTightZPeakTwoAk4JetsExclusiveResolvedOneBtag_firstlepton_pt'}
    var_2btagZpeak_dict = {(path_Zpeak,'Z Peak 2 btag'): 'ElEl_HasElElTightZPeakTwoAk4JetsExclusiveResolvedTwoBtags_firstlepton_pt'}

    variable_name = 'e^{+}e^{-} First lepton P_{T}'
    instance = DYControlRegion(var_0btagZveto_dict,variable_name,"ElEl_firstlepton_PT_0btagZVeto",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_0btagZveto_dict,variable_name,"ElEl_firstlepton_PT_0btagZVeto",plot_data=True,plot_MC=True,exclude_DY=False)
    instance = DYControlRegion(var_1btagZveto_dict,variable_name,"ElEl_firstlepton_PT_1btagZVeto",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_1btagZveto_dict,variable_name,"ElEl_firstlepton_PT_1btagZVeto",plot_data=True,plot_MC=True,exclude_DY=False)
    instance = DYControlRegion(var_2btagZveto_dict,variable_name,"ElEl_firstlepton_PT_2btagZVeto",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_2btagZveto_dict,variable_name,"ElEl_firstlepton_PT_2btagZVeto",plot_data=True,plot_MC=True,exclude_DY=False)

    instance = DYControlRegion(var_0btagZpeak_dict,variable_name,"ElEl_firstlepton_PT_0btagZpeak",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_0btagZpeak_dict,variable_name,"ElEl_firstlepton_PT_0btagZpeak",plot_data=True,plot_MC=True,exclude_DY=False)
    instance = DYControlRegion(var_1btagZpeak_dict,variable_name,"ElEl_firstlepton_PT_1btagZpeak",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_1btagZpeak_dict,variable_name,"ElEl_firstlepton_PT_1btagZpeak",plot_data=True,plot_MC=True,exclude_DY=False)
    instance = DYControlRegion(var_2btagZpeak_dict,variable_name,"ElEl_firstlepton_PT_2btagZpeak",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_2btagZpeak_dict,variable_name,"ElEl_firstlepton_PT_2btagZpeak",plot_data=True,plot_MC=True,exclude_DY=False)

    instance = DYControlRegion({**var_0btagZveto_dict,**var_0btagZpeak_dict},variable_name,"ElEl_firstlepton_PT_0btag",plot_data=True,plot_MC=False,exclude_DY=True)
    instance = DYControlRegion({**var_0btagZveto_dict,**var_0btagZpeak_dict},variable_name,"ElEl_firstlepton_PT_0btag",plot_data=False,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion({**var_1btagZveto_dict,**var_1btagZpeak_dict},variable_name,"ElEl_firstlepton_PT_1btag",plot_data=True,plot_MC=False,exclude_DY=True)
    instance = DYControlRegion({**var_1btagZveto_dict,**var_1btagZpeak_dict},variable_name,"ElEl_firstlepton_PT_1btag",plot_data=False,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion({**var_2btagZveto_dict,**var_2btagZpeak_dict},variable_name,"ElEl_firstlepton_PT_2btag",plot_data=True,plot_MC=False,exclude_DY=True)
    instance = DYControlRegion({**var_2btagZveto_dict,**var_2btagZpeak_dict},variable_name,"ElEl_firstlepton_PT_2btag",plot_data=False,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion({**var_0btagZveto_dict,**var_1btagZveto_dict,**var_2btagZveto_dict},variable_name,"ElEl_firstlepton_PT_ZVeto",plot_data=False,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion({**var_0btagZveto_dict,**var_1btagZveto_dict,**var_2btagZveto_dict},variable_name,"ElEl_firstlepton_PT_ZVeto",plot_data=True,plot_MC=False,exclude_DY=True)
    instance = DYControlRegion({**var_0btagZpeak_dict,**var_1btagZpeak_dict,**var_2btagZpeak_dict},variable_name,"ElEl_firstlepton_PT_ZPeak",plot_data=False,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion({**var_0btagZpeak_dict,**var_1btagZpeak_dict,**var_2btagZpeak_dict},variable_name,"ElEl_firstlepton_PT_ZPeak",plot_data=True,plot_MC=False,exclude_DY=True)

    # MuMu First lepton pt #
    var_0btagZveto_dict = {(path_Zveto,'Z Veto 0 btag'): 'MuMu_HasMuMuTightTwoAk4JetsExclusiveResolvedNoBtag_firstlepton_pt'}
    var_1btagZveto_dict = {(path_Zveto,'Z Veto 1 btag'): 'MuMu_HasMuMuTightTwoAk4JetsExclusiveResolvedOneBtag_firstlepton_pt'}
    var_2btagZveto_dict = {(path_Zveto,'Z Veto 2 btag'): 'MuMu_HasMuMuTightTwoAk4JetsExclusiveResolvedTwoBtags_firstlepton_pt'}

    var_0btagZpeak_dict = {(path_Zpeak,'Z Peak 0 btag'): 'MuMu_HasMuMuTightZPeakTwoAk4JetsExclusiveResolvedNoBtag_firstlepton_pt'}
    var_1btagZpeak_dict = {(path_Zpeak,'Z Peak 1 btag'): 'MuMu_HasMuMuTightZPeakTwoAk4JetsExclusiveResolvedOneBtag_firstlepton_pt'}
    var_2btagZpeak_dict = {(path_Zpeak,'Z Peak 2 btag'): 'MuMu_HasMuMuTightZPeakTwoAk4JetsExclusiveResolvedTwoBtags_firstlepton_pt'}

    variable_name = '#mu^{+}#mu^{-} First lepton P_{T}'
    instance = DYControlRegion(var_0btagZveto_dict,variable_name,"MuMu_firstlepton_PT_0btagZVeto",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_0btagZveto_dict,variable_name,"MuMu_firstlepton_PT_0btagZVeto",plot_data=True,plot_MC=True,exclude_DY=False)
    instance = DYControlRegion(var_1btagZveto_dict,variable_name,"MuMu_firstlepton_PT_1btagZVeto",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_1btagZveto_dict,variable_name,"MuMu_firstlepton_PT_1btagZVeto",plot_data=True,plot_MC=True,exclude_DY=False)
    instance = DYControlRegion(var_2btagZveto_dict,variable_name,"MuMu_firstlepton_PT_2btagZVeto",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_2btagZveto_dict,variable_name,"MuMu_firstlepton_PT_2btagZVeto",plot_data=True,plot_MC=True,exclude_DY=False)

    instance = DYControlRegion(var_0btagZpeak_dict,variable_name,"MuMu_firstlepton_PT_0btagZpeak",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_0btagZpeak_dict,variable_name,"MuMu_firstlepton_PT_0btagZpeak",plot_data=True,plot_MC=True,exclude_DY=False)
    instance = DYControlRegion(var_1btagZpeak_dict,variable_name,"MuMu_firstlepton_PT_1btagZpeak",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_1btagZpeak_dict,variable_name,"MuMu_firstlepton_PT_1btagZpeak",plot_data=True,plot_MC=True,exclude_DY=False)
    instance = DYControlRegion(var_2btagZpeak_dict,variable_name,"MuMu_firstlepton_PT_2btagZpeak",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_2btagZpeak_dict,variable_name,"MuMu_firstlepton_PT_2btagZpeak",plot_data=True,plot_MC=True,exclude_DY=False)

    instance = DYControlRegion({**var_0btagZveto_dict,**var_0btagZpeak_dict},variable_name,"MuMu_firstlepton_PT_0btag",plot_data=True,plot_MC=False,exclude_DY=True)
    instance = DYControlRegion({**var_0btagZveto_dict,**var_0btagZpeak_dict},variable_name,"MuMu_firstlepton_PT_0btag",plot_data=False,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion({**var_1btagZveto_dict,**var_1btagZpeak_dict},variable_name,"MuMu_firstlepton_PT_1btag",plot_data=True,plot_MC=False,exclude_DY=True)
    instance = DYControlRegion({**var_1btagZveto_dict,**var_1btagZpeak_dict},variable_name,"MuMu_firstlepton_PT_1btag",plot_data=False,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion({**var_2btagZveto_dict,**var_2btagZpeak_dict},variable_name,"MuMu_firstlepton_PT_2btag",plot_data=True,plot_MC=False,exclude_DY=True)
    instance = DYControlRegion({**var_2btagZveto_dict,**var_2btagZpeak_dict},variable_name,"MuMu_firstlepton_PT_2btag",plot_data=False,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion({**var_0btagZveto_dict,**var_1btagZveto_dict,**var_2btagZveto_dict},variable_name,"MuMu_firstlepton_PT_ZVeto",plot_data=False,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion({**var_0btagZveto_dict,**var_1btagZveto_dict,**var_2btagZveto_dict},variable_name,"MuMu_firstlepton_PT_ZVeto",plot_data=True,plot_MC=False,exclude_DY=True)
    instance = DYControlRegion({**var_0btagZpeak_dict,**var_1btagZpeak_dict,**var_2btagZpeak_dict},variable_name,"MuMu_firstlepton_PT_ZPeak",plot_data=False,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion({**var_0btagZpeak_dict,**var_1btagZpeak_dict,**var_2btagZpeak_dict},variable_name,"MuMu_firstlepton_PT_ZPeak",plot_data=True,plot_MC=False,exclude_DY=True)

    # ElEl leadjet Pt #
    var_0btagZveto_dict = {(path_Zveto,'Z Veto 0 btag'): 'ElEl_HasElElTightTwoAk4JetsExclusiveResolvedNoBtag_leadjet_pt'}
    var_1btagZveto_dict = {(path_Zveto,'Z Veto 1 btag'): 'ElEl_HasElElTightTwoAk4JetsExclusiveResolvedOneBtag_leadbjet_pt'}
    var_2btagZveto_dict = {(path_Zveto,'Z Veto 2 btag'): 'ElEl_HasElElTightTwoAk4JetsExclusiveResolvedTwoBtags_leadbjet_pt'}

    var_0btagZpeak_dict = {(path_Zpeak,'Z Peak 0 btag'): 'ElEl_HasElElTightZPeakTwoAk4JetsExclusiveResolvedNoBtag_leadjet_pt'}
    var_1btagZpeak_dict = {(path_Zpeak,'Z Peak 1 btag'): 'ElEl_HasElElTightZPeakTwoAk4JetsExclusiveResolvedOneBtag_leadbjet_pt'}
    var_2btagZpeak_dict = {(path_Zpeak,'Z Peak 2 btag'): 'ElEl_HasElElTightZPeakTwoAk4JetsExclusiveResolvedTwoBtags_leadbjet_pt'}

    variable_name = 'e^{+}e^{-} Lead jet P_{T}'
    instance = DYControlRegion(var_0btagZveto_dict,variable_name,"ElEl_leadjet_PT_0btagZVeto",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_0btagZveto_dict,variable_name,"ElEl_leadjet_PT_0btagZVeto",plot_data=True,plot_MC=True,exclude_DY=False)
    instance = DYControlRegion(var_1btagZveto_dict,variable_name,"ElEl_leadjet_PT_1btagZVeto",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_1btagZveto_dict,variable_name,"ElEl_leadjet_PT_1btagZVeto",plot_data=True,plot_MC=True,exclude_DY=False)
    instance = DYControlRegion(var_2btagZveto_dict,variable_name,"ElEl_leadjet_PT_2btagZVeto",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_2btagZveto_dict,variable_name,"ElEl_leadjet_PT_2btagZVeto",plot_data=True,plot_MC=True,exclude_DY=False)

    instance = DYControlRegion(var_0btagZpeak_dict,variable_name,"ElEl_leadjet_PT_0btagZpeak",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_0btagZpeak_dict,variable_name,"ElEl_leadjet_PT_0btagZpeak",plot_data=True,plot_MC=True,exclude_DY=False)
    instance = DYControlRegion(var_1btagZpeak_dict,variable_name,"ElEl_leadjet_PT_1btagZpeak",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_1btagZpeak_dict,variable_name,"ElEl_leadjet_PT_1btagZpeak",plot_data=True,plot_MC=True,exclude_DY=False)
    instance = DYControlRegion(var_2btagZpeak_dict,variable_name,"ElEl_leadjet_PT_2btagZpeak",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_2btagZpeak_dict,variable_name,"ElEl_leadjet_PT_2btagZpeak",plot_data=True,plot_MC=True,exclude_DY=False)

    instance = DYControlRegion({**var_0btagZveto_dict,**var_0btagZpeak_dict},variable_name,"ElEl_leadjet_PT_0btag",plot_data=True,plot_MC=False,exclude_DY=True)
    instance = DYControlRegion({**var_0btagZveto_dict,**var_0btagZpeak_dict},variable_name,"ElEl_leadjet_PT_0btag",plot_data=False,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion({**var_1btagZveto_dict,**var_1btagZpeak_dict},variable_name,"ElEl_leadjet_PT_1btag",plot_data=True,plot_MC=False,exclude_DY=True)
    instance = DYControlRegion({**var_1btagZveto_dict,**var_1btagZpeak_dict},variable_name,"ElEl_leadjet_PT_1btag",plot_data=False,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion({**var_2btagZveto_dict,**var_2btagZpeak_dict},variable_name,"ElEl_leadjet_PT_2btag",plot_data=True,plot_MC=False,exclude_DY=True)
    instance = DYControlRegion({**var_2btagZveto_dict,**var_2btagZpeak_dict},variable_name,"ElEl_leadjet_PT_2btag",plot_data=False,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion({**var_0btagZveto_dict,**var_1btagZveto_dict,**var_2btagZveto_dict},variable_name,"ElEl_leadjet_PT_ZVeto",plot_data=False,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion({**var_0btagZveto_dict,**var_1btagZveto_dict,**var_2btagZveto_dict},variable_name,"ElEl_leadjet_PT_ZVeto",plot_data=True,plot_MC=False,exclude_DY=True)
    instance = DYControlRegion({**var_0btagZpeak_dict,**var_1btagZpeak_dict,**var_2btagZpeak_dict},variable_name,"ElEl_leadjet_PT_ZPeak",plot_data=False,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion({**var_0btagZpeak_dict,**var_1btagZpeak_dict,**var_2btagZpeak_dict},variable_name,"ElEl_leadjet_PT_ZPeak",plot_data=True,plot_MC=False,exclude_DY=True)

    # MuMu leadjet Pt #
    var_0btagZveto_dict = {(path_Zveto,'Z Veto 0 btag'): 'MuMu_HasMuMuTightTwoAk4JetsExclusiveResolvedNoBtag_leadjet_pt'}
    var_1btagZveto_dict = {(path_Zveto,'Z Veto 1 btag'): 'MuMu_HasMuMuTightTwoAk4JetsExclusiveResolvedOneBtag_leadbjet_pt'}
    var_2btagZveto_dict = {(path_Zveto,'Z Veto 2 btag'): 'MuMu_HasMuMuTightTwoAk4JetsExclusiveResolvedTwoBtags_leadbjet_pt'}

    var_0btagZpeak_dict = {(path_Zpeak,'Z Peak 0 btag'): 'MuMu_HasMuMuTightZPeakTwoAk4JetsExclusiveResolvedNoBtag_leadjet_pt'}
    var_1btagZpeak_dict = {(path_Zpeak,'Z Peak 1 btag'): 'MuMu_HasMuMuTightZPeakTwoAk4JetsExclusiveResolvedOneBtag_leadbjet_pt'}
    var_2btagZpeak_dict = {(path_Zpeak,'Z Peak 2 btag'): 'MuMu_HasMuMuTightZPeakTwoAk4JetsExclusiveResolvedTwoBtags_leadbjet_pt'}

    variable_name = '#mu^{+}#mu^{-} Lead jet P_{T}'
    instance = DYControlRegion(var_0btagZveto_dict,variable_name,"MuMu_leadjet_PT_0btagZVeto",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_0btagZveto_dict,variable_name,"MuMu_leadjet_PT_0btagZVeto",plot_data=True,plot_MC=True,exclude_DY=False)
    instance = DYControlRegion(var_1btagZveto_dict,variable_name,"MuMu_leadjet_PT_1btagZVeto",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_1btagZveto_dict,variable_name,"MuMu_leadjet_PT_1btagZVeto",plot_data=True,plot_MC=True,exclude_DY=False)
    instance = DYControlRegion(var_2btagZveto_dict,variable_name,"MuMu_leadjet_PT_2btagZVeto",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_2btagZveto_dict,variable_name,"MuMu_leadjet_PT_2btagZVeto",plot_data=True,plot_MC=True,exclude_DY=False)

    instance = DYControlRegion(var_0btagZpeak_dict,variable_name,"MuMu_leadjet_PT_0btagZpeak",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_0btagZpeak_dict,variable_name,"MuMu_leadjet_PT_0btagZpeak",plot_data=True,plot_MC=True,exclude_DY=False)
    instance = DYControlRegion(var_1btagZpeak_dict,variable_name,"MuMu_leadjet_PT_1btagZpeak",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_1btagZpeak_dict,variable_name,"MuMu_leadjet_PT_1btagZpeak",plot_data=True,plot_MC=True,exclude_DY=False)
    instance = DYControlRegion(var_2btagZpeak_dict,variable_name,"MuMu_leadjet_PT_2btagZpeak",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_2btagZpeak_dict,variable_name,"MuMu_leadjet_PT_2btagZpeak",plot_data=True,plot_MC=True,exclude_DY=False)

    instance = DYControlRegion({**var_0btagZveto_dict,**var_0btagZpeak_dict},variable_name,"MuMu_leadjet_PT_0btag",plot_data=True,plot_MC=False,exclude_DY=True)
    instance = DYControlRegion({**var_0btagZveto_dict,**var_0btagZpeak_dict},variable_name,"MuMu_leadjet_PT_0btag",plot_data=False,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion({**var_1btagZveto_dict,**var_1btagZpeak_dict},variable_name,"MuMu_leadjet_PT_1btag",plot_data=True,plot_MC=False,exclude_DY=True)
    instance = DYControlRegion({**var_1btagZveto_dict,**var_1btagZpeak_dict},variable_name,"MuMu_leadjet_PT_1btag",plot_data=False,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion({**var_2btagZveto_dict,**var_2btagZpeak_dict},variable_name,"MuMu_leadjet_PT_2btag",plot_data=True,plot_MC=False,exclude_DY=True)
    instance = DYControlRegion({**var_2btagZveto_dict,**var_2btagZpeak_dict},variable_name,"MuMu_leadjet_PT_2btag",plot_data=False,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion({**var_0btagZveto_dict,**var_1btagZveto_dict,**var_2btagZveto_dict},variable_name,"MuMu_leadjet_PT_ZVeto",plot_data=False,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion({**var_0btagZveto_dict,**var_1btagZveto_dict,**var_2btagZveto_dict},variable_name,"MuMu_leadjet_PT_ZVeto",plot_data=True,plot_MC=False,exclude_DY=True)
    instance = DYControlRegion({**var_0btagZpeak_dict,**var_1btagZpeak_dict,**var_2btagZpeak_dict},variable_name,"MuMu_leadjet_PT_ZPeak",plot_data=False,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion({**var_0btagZpeak_dict,**var_1btagZpeak_dict,**var_2btagZpeak_dict},variable_name,"MuMu_leadjet_PT_ZPeak",plot_data=True,plot_MC=False,exclude_DY=True)

    # ElEl N jets #
    var_0btagZveto_dict = {(path_Zveto,'Z Veto 0 btag'): 'ElEl_HasElElTightTwoAk4JetsExclusiveResolvedNoBtag_Ak4Jets_N'}
    var_1btagZveto_dict = {(path_Zveto,'Z Veto 1 btag'): 'ElEl_HasElElTightTwoAk4JetsExclusiveResolvedOneBtag_Ak4Jets_N'}
    var_2btagZveto_dict = {(path_Zveto,'Z Veto 2 btag'): 'ElEl_HasElElTightTwoAk4JetsExclusiveResolvedTwoBtags_Ak4Jets_N'}

    var_0btagZpeak_dict = {(path_Zpeak,'Z Peak 0 btag'): 'ElEl_HasElElTightZPeakTwoAk4JetsExclusiveResolvedNoBtag_Ak4Jets_N'}
    var_1btagZpeak_dict = {(path_Zpeak,'Z Peak 1 btag'): 'ElEl_HasElElTightZPeakTwoAk4JetsExclusiveResolvedOneBtag_Ak4Jets_N'}
    var_2btagZpeak_dict = {(path_Zpeak,'Z Peak 2 btag'): 'ElEl_HasElElTightZPeakTwoAk4JetsExclusiveResolvedTwoBtags_Ak4Jets_N'}

    variable_name = 'e^{+}e^{-} Number Ak4 jets'
    instance = DYControlRegion(var_0btagZveto_dict,variable_name,"ElEl_NAk4Jets_0btagZVeto",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_0btagZveto_dict,variable_name,"ElEl_NAk4Jets_0btagZVeto",plot_data=True,plot_MC=True,exclude_DY=False)
    instance = DYControlRegion(var_1btagZveto_dict,variable_name,"ElEl_NAk4Jets_1btagZVeto",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_1btagZveto_dict,variable_name,"ElEl_NAk4Jets_1btagZVeto",plot_data=True,plot_MC=True,exclude_DY=False)
    instance = DYControlRegion(var_2btagZveto_dict,variable_name,"ElEl_NAk4Jets_2btagZVeto",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_2btagZveto_dict,variable_name,"ElEl_NAk4Jets_2btagZVeto",plot_data=True,plot_MC=True,exclude_DY=False)

    instance = DYControlRegion(var_0btagZpeak_dict,variable_name,"ElEl_NAk4Jets_0btagZpeak",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_0btagZpeak_dict,variable_name,"ElEl_NAk4Jets_0btagZpeak",plot_data=True,plot_MC=True,exclude_DY=False)
    instance = DYControlRegion(var_1btagZpeak_dict,variable_name,"ElEl_NAk4Jets_1btagZpeak",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_1btagZpeak_dict,variable_name,"ElEl_NAk4Jets_1btagZpeak",plot_data=True,plot_MC=True,exclude_DY=False)
    instance = DYControlRegion(var_2btagZpeak_dict,variable_name,"ElEl_NAk4Jets_2btagZpeak",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_2btagZpeak_dict,variable_name,"ElEl_NAk4Jets_2btagZpeak",plot_data=True,plot_MC=True,exclude_DY=False)

    instance = DYControlRegion({**var_0btagZveto_dict,**var_0btagZpeak_dict},variable_name,"ElEl_leadjet_PT_0btag",plot_data=True,plot_MC=False,exclude_DY=True)
    instance = DYControlRegion({**var_0btagZveto_dict,**var_0btagZpeak_dict},variable_name,"ElEl_leadjet_PT_0btag",plot_data=False,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion({**var_1btagZveto_dict,**var_1btagZpeak_dict},variable_name,"ElEl_leadjet_PT_1btag",plot_data=True,plot_MC=False,exclude_DY=True)
    instance = DYControlRegion({**var_1btagZveto_dict,**var_1btagZpeak_dict},variable_name,"ElEl_leadjet_PT_1btag",plot_data=False,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion({**var_2btagZveto_dict,**var_2btagZpeak_dict},variable_name,"ElEl_leadjet_PT_2btag",plot_data=True,plot_MC=False,exclude_DY=True)
    instance = DYControlRegion({**var_2btagZveto_dict,**var_2btagZpeak_dict},variable_name,"ElEl_leadjet_PT_2btag",plot_data=False,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion({**var_0btagZveto_dict,**var_1btagZveto_dict,**var_2btagZveto_dict},variable_name,"ElEl_leadjet_PT_ZVeto",plot_data=False,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion({**var_0btagZveto_dict,**var_1btagZveto_dict,**var_2btagZveto_dict},variable_name,"ElEl_leadjet_PT_ZVeto",plot_data=True,plot_MC=False,exclude_DY=True)
    instance = DYControlRegion({**var_0btagZpeak_dict,**var_1btagZpeak_dict,**var_2btagZpeak_dict},variable_name,"ElEl_leadjet_PT_ZPeak",plot_data=False,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion({**var_0btagZpeak_dict,**var_1btagZpeak_dict,**var_2btagZpeak_dict},variable_name,"ElEl_leadjet_PT_ZPeak",plot_data=True,plot_MC=False,exclude_DY=True)

    # MuMu N jets #
    var_0btagZveto_dict = {(path_Zveto,'Z Veto 0 btag'): 'MuMu_HasMuMuTightTwoAk4JetsExclusiveResolvedNoBtag_Ak4Jets_N'}
    var_1btagZveto_dict = {(path_Zveto,'Z Veto 1 btag'): 'MuMu_HasMuMuTightTwoAk4JetsExclusiveResolvedOneBtag_Ak4Jets_N'}
    var_2btagZveto_dict = {(path_Zveto,'Z Veto 2 btag'): 'MuMu_HasMuMuTightTwoAk4JetsExclusiveResolvedTwoBtags_Ak4Jets_N'}

    var_0btagZpeak_dict = {(path_Zpeak,'Z Peak 0 btag'): 'MuMu_HasMuMuTightZPeakTwoAk4JetsExclusiveResolvedNoBtag_Ak4Jets_N'}
    var_1btagZpeak_dict = {(path_Zpeak,'Z Peak 1 btag'): 'MuMu_HasMuMuTightZPeakTwoAk4JetsExclusiveResolvedOneBtag_Ak4Jets_N'}
    var_2btagZpeak_dict = {(path_Zpeak,'Z Peak 2 btag'): 'MuMu_HasMuMuTightZPeakTwoAk4JetsExclusiveResolvedTwoBtags_Ak4Jets_N'}

    variable_name = '#mu^{+}#mu^{-} Number Ak4 jets'
    instance = DYControlRegion(var_0btagZveto_dict,variable_name,"MuMu_NAk4Jets_0btagZVeto",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_0btagZveto_dict,variable_name,"MuMu_NAk4Jets_0btagZVeto",plot_data=True,plot_MC=True,exclude_DY=False)
    instance = DYControlRegion(var_1btagZveto_dict,variable_name,"MuMu_NAk4Jets_1btagZVeto",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_1btagZveto_dict,variable_name,"MuMu_NAk4Jets_1btagZVeto",plot_data=True,plot_MC=True,exclude_DY=False)
    instance = DYControlRegion(var_2btagZveto_dict,variable_name,"MuMu_NAk4Jets_2btagZVeto",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_2btagZveto_dict,variable_name,"MuMu_NAk4Jets_2btagZVeto",plot_data=True,plot_MC=True,exclude_DY=False)

    instance = DYControlRegion(var_0btagZpeak_dict,variable_name,"MuMu_NAk4Jets_0btagZpeak",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_0btagZpeak_dict,variable_name,"MuMu_NAk4Jets_0btagZpeak",plot_data=True,plot_MC=True,exclude_DY=False)
    instance = DYControlRegion(var_1btagZpeak_dict,variable_name,"MuMu_NAk4Jets_1btagZpeak",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_1btagZpeak_dict,variable_name,"MuMu_NAk4Jets_1btagZpeak",plot_data=True,plot_MC=True,exclude_DY=False)
    instance = DYControlRegion(var_2btagZpeak_dict,variable_name,"MuMu_NAk4Jets_2btagZpeak",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_2btagZpeak_dict,variable_name,"MuMu_NAk4Jets_2btagZpeak",plot_data=True,plot_MC=True,exclude_DY=False)

    instance = DYControlRegion({**var_0btagZveto_dict,**var_0btagZpeak_dict},variable_name,"MuMu_leadjet_PT_0btag",plot_data=True,plot_MC=False,exclude_DY=True)
    instance = DYControlRegion({**var_0btagZveto_dict,**var_0btagZpeak_dict},variable_name,"MuMu_leadjet_PT_0btag",plot_data=False,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion({**var_1btagZveto_dict,**var_1btagZpeak_dict},variable_name,"MuMu_leadjet_PT_1btag",plot_data=True,plot_MC=False,exclude_DY=True)
    instance = DYControlRegion({**var_1btagZveto_dict,**var_1btagZpeak_dict},variable_name,"MuMu_leadjet_PT_1btag",plot_data=False,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion({**var_2btagZveto_dict,**var_2btagZpeak_dict},variable_name,"MuMu_leadjet_PT_2btag",plot_data=True,plot_MC=False,exclude_DY=True)
    instance = DYControlRegion({**var_2btagZveto_dict,**var_2btagZpeak_dict},variable_name,"MuMu_leadjet_PT_2btag",plot_data=False,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion({**var_0btagZveto_dict,**var_1btagZveto_dict,**var_2btagZveto_dict},variable_name,"MuMu_leadjet_PT_ZVeto",plot_data=False,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion({**var_0btagZveto_dict,**var_1btagZveto_dict,**var_2btagZveto_dict},variable_name,"MuMu_leadjet_PT_ZVeto",plot_data=True,plot_MC=False,exclude_DY=True)
    instance = DYControlRegion({**var_0btagZpeak_dict,**var_1btagZpeak_dict,**var_2btagZpeak_dict},variable_name,"MuMu_leadjet_PT_ZPeak",plot_data=False,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion({**var_0btagZpeak_dict,**var_1btagZpeak_dict,**var_2btagZpeak_dict},variable_name,"MuMu_leadjet_PT_ZPeak",plot_data=True,plot_MC=False,exclude_DY=True)


    # ElEl N bjets #
    var_0btagZveto_dict = {(path_Zveto,'Z Veto 0 btag'): 'ElEl_HasElElTightTwoAk4JetsExclusiveResolvedNoBtag_Ak4BJets_N'}
    var_1btagZveto_dict = {(path_Zveto,'Z Veto 1 btag'): 'ElEl_HasElElTightTwoAk4JetsExclusiveResolvedOneBtag_Ak4BJets_N'}
    var_2btagZveto_dict = {(path_Zveto,'Z Veto 2 btag'): 'ElEl_HasElElTightTwoAk4JetsExclusiveResolvedTwoBtags_Ak4BJets_N'}

    var_0btagZpeak_dict = {(path_Zpeak,'Z Peak 0 btag'): 'ElEl_HasElElTightZPeakTwoAk4JetsExclusiveResolvedNoBtag_Ak4BJets_N'}
    var_1btagZpeak_dict = {(path_Zpeak,'Z Peak 1 btag'): 'ElEl_HasElElTightZPeakTwoAk4JetsExclusiveResolvedOneBtag_Ak4BJets_N'}
    var_2btagZpeak_dict = {(path_Zpeak,'Z Peak 2 btag'): 'ElEl_HasElElTightZPeakTwoAk4JetsExclusiveResolvedTwoBtags_Ak4BJets_N'}

    variable_name = 'e^{+}e^{-} Number Ak4 bjets'
    instance = DYControlRegion(var_0btagZveto_dict,variable_name,"ElEl_NAk4BJets_0btagZVeto",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_0btagZveto_dict,variable_name,"ElEl_NAk4BJets_0btagZVeto",plot_data=True,plot_MC=True,exclude_DY=False)
    instance = DYControlRegion(var_1btagZveto_dict,variable_name,"ElEl_NAk4BJets_1btagZVeto",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_1btagZveto_dict,variable_name,"ElEl_NAk4BJets_1btagZVeto",plot_data=True,plot_MC=True,exclude_DY=False)
    instance = DYControlRegion(var_2btagZveto_dict,variable_name,"ElEl_NAk4BJets_2btagZVeto",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_2btagZveto_dict,variable_name,"ElEl_NAk4BJets_2btagZVeto",plot_data=True,plot_MC=True,exclude_DY=False)

    instance = DYControlRegion(var_0btagZpeak_dict,variable_name,"ElEl_NAk4BJets_0btagZpeak",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_0btagZpeak_dict,variable_name,"ElEl_NAk4BJets_0btagZpeak",plot_data=True,plot_MC=True,exclude_DY=False)
    instance = DYControlRegion(var_1btagZpeak_dict,variable_name,"ElEl_NAk4BJets_1btagZpeak",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_1btagZpeak_dict,variable_name,"ElEl_NAk4BJets_1btagZpeak",plot_data=True,plot_MC=True,exclude_DY=False)
    instance = DYControlRegion(var_2btagZpeak_dict,variable_name,"ElEl_NAk4BJets_2btagZpeak",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_2btagZpeak_dict,variable_name,"ElEl_NAk4BJets_2btagZpeak",plot_data=True,plot_MC=True,exclude_DY=False)

    instance = DYControlRegion({**var_0btagZveto_dict,**var_0btagZpeak_dict},variable_name,"ElEl_NAk4BJets_0btag",plot_data=True,plot_MC=False,exclude_DY=True)
    instance = DYControlRegion({**var_0btagZveto_dict,**var_0btagZpeak_dict},variable_name,"ElEl_NAk4BJets_0btag",plot_data=False,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion({**var_1btagZveto_dict,**var_1btagZpeak_dict},variable_name,"ElEl_NAk4BJets_1btag",plot_data=True,plot_MC=False,exclude_DY=True)
    instance = DYControlRegion({**var_1btagZveto_dict,**var_1btagZpeak_dict},variable_name,"ElEl_NAk4BJets_1btag",plot_data=False,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion({**var_2btagZveto_dict,**var_2btagZpeak_dict},variable_name,"ElEl_NAk4BJets_2btag",plot_data=True,plot_MC=False,exclude_DY=True)
    instance = DYControlRegion({**var_2btagZveto_dict,**var_2btagZpeak_dict},variable_name,"ElEl_NAk4BJets_2btag",plot_data=False,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion({**var_0btagZveto_dict,**var_1btagZveto_dict,**var_2btagZveto_dict},variable_name,"ElEl_NAk4BJets_ZVeto",plot_data=False,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion({**var_0btagZveto_dict,**var_1btagZveto_dict,**var_2btagZveto_dict},variable_name,"ElEl_NAk4BJets_ZVeto",plot_data=True,plot_MC=False,exclude_DY=True)
    instance = DYControlRegion({**var_0btagZpeak_dict,**var_1btagZpeak_dict,**var_2btagZpeak_dict},variable_name,"ElEl_NAk4BJets_ZPeak",plot_data=False,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion({**var_0btagZpeak_dict,**var_1btagZpeak_dict,**var_2btagZpeak_dict},variable_name,"ElEl_NAk4BJets_ZPeak",plot_data=True,plot_MC=False,exclude_DY=True)

    # MuMu N bjets #
    var_0btagZveto_dict = {(path_Zveto,'Z Veto 0 btag'): 'MuMu_HasMuMuTightTwoAk4JetsExclusiveResolvedNoBtag_Ak4BJets_N'}
    var_1btagZveto_dict = {(path_Zveto,'Z Veto 1 btag'): 'MuMu_HasMuMuTightTwoAk4JetsExclusiveResolvedOneBtag_Ak4BJets_N'}
    var_2btagZveto_dict = {(path_Zveto,'Z Veto 2 btag'): 'MuMu_HasMuMuTightTwoAk4JetsExclusiveResolvedTwoBtags_Ak4BJets_N'}

    var_0btagZpeak_dict = {(path_Zpeak,'Z Peak 0 btag'): 'MuMu_HasMuMuTightZPeakTwoAk4JetsExclusiveResolvedNoBtag_Ak4BJets_N'}
    var_1btagZpeak_dict = {(path_Zpeak,'Z Peak 1 btag'): 'MuMu_HasMuMuTightZPeakTwoAk4JetsExclusiveResolvedOneBtag_Ak4BJets_N'}
    var_2btagZpeak_dict = {(path_Zpeak,'Z Peak 2 btag'): 'MuMu_HasMuMuTightZPeakTwoAk4JetsExclusiveResolvedTwoBtags_Ak4BJets_N'}

    variable_name = '#mu^{+}#mu^{-} Number Ak4 bjets'
    instance = DYControlRegion(var_0btagZveto_dict,variable_name,"MuMu_NAk4BJets_0btagZVeto",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_0btagZveto_dict,variable_name,"MuMu_NAk4BJets_0btagZVeto",plot_data=True,plot_MC=True,exclude_DY=False)
    instance = DYControlRegion(var_1btagZveto_dict,variable_name,"MuMu_NAk4BJets_1btagZVeto",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_1btagZveto_dict,variable_name,"MuMu_NAk4BJets_1btagZVeto",plot_data=True,plot_MC=True,exclude_DY=False)
    instance = DYControlRegion(var_2btagZveto_dict,variable_name,"MuMu_NAk4BJets_2btagZVeto",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_2btagZveto_dict,variable_name,"MuMu_NAk4BJets_2btagZVeto",plot_data=True,plot_MC=True,exclude_DY=False)

    instance = DYControlRegion(var_0btagZpeak_dict,variable_name,"MuMu_NAk4BJets_0btagZpeak",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_0btagZpeak_dict,variable_name,"MuMu_NAk4BJets_0btagZpeak",plot_data=True,plot_MC=True,exclude_DY=False)
    instance = DYControlRegion(var_1btagZpeak_dict,variable_name,"MuMu_NAk4BJets_1btagZpeak",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_1btagZpeak_dict,variable_name,"MuMu_NAk4BJets_1btagZpeak",plot_data=True,plot_MC=True,exclude_DY=False)
    instance = DYControlRegion(var_2btagZpeak_dict,variable_name,"MuMu_NAk4BJets_2btagZpeak",plot_data=True,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion(var_2btagZpeak_dict,variable_name,"MuMu_NAk4BJets_2btagZpeak",plot_data=True,plot_MC=True,exclude_DY=False)

    instance = DYControlRegion({**var_0btagZveto_dict,**var_0btagZpeak_dict},variable_name,"MuMu_NAk4BJets_0btag",plot_data=True,plot_MC=False,exclude_DY=True)
    instance = DYControlRegion({**var_0btagZveto_dict,**var_0btagZpeak_dict},variable_name,"MuMu_NAk4BJets_0btag",plot_data=False,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion({**var_1btagZveto_dict,**var_1btagZpeak_dict},variable_name,"MuMu_NAk4BJets_1btag",plot_data=True,plot_MC=False,exclude_DY=True)
    instance = DYControlRegion({**var_1btagZveto_dict,**var_1btagZpeak_dict},variable_name,"MuMu_NAk4BJets_1btag",plot_data=False,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion({**var_2btagZveto_dict,**var_2btagZpeak_dict},variable_name,"MuMu_NAk4BJets_2btag",plot_data=True,plot_MC=False,exclude_DY=True)
    instance = DYControlRegion({**var_2btagZveto_dict,**var_2btagZpeak_dict},variable_name,"MuMu_NAk4BJets_2btag",plot_data=False,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion({**var_0btagZveto_dict,**var_1btagZveto_dict,**var_2btagZveto_dict},variable_name,"MuMu_NAk4BJets_ZVeto",plot_data=False,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion({**var_0btagZveto_dict,**var_1btagZveto_dict,**var_2btagZveto_dict},variable_name,"MuMu_NAk4BJets_ZVeto",plot_data=True,plot_MC=False,exclude_DY=True)
    instance = DYControlRegion({**var_0btagZpeak_dict,**var_1btagZpeak_dict,**var_2btagZpeak_dict},variable_name,"MuMu_NAk4BJets_ZPeak",plot_data=False,plot_MC=True,exclude_DY=True)
    instance = DYControlRegion({**var_0btagZpeak_dict,**var_1btagZpeak_dict,**var_2btagZpeak_dict},variable_name,"MuMu_NAk4BJets_ZPeak",plot_data=True,plot_MC=False,exclude_DY=True)

import subprocess
import glob

list_pdf = sorted(glob.glob("DYPlots/*firstlepton*pdf"))
subprocess.Popen(['pdfunite']+list_pdf+["DYPlots_firstlepton.pdf"])

list_pdf = sorted(glob.glob("DYPlots/*leadjet*pdf"))
subprocess.Popen(['pdfunite']+list_pdf+["DYPlots_leadjet.pdf"])

list_pdf = sorted(glob.glob("DYPlots/*NAk4Jets*pdf"))
subprocess.Popen(['pdfunite']+list_pdf+["DYPlots_NAk4Jets.pdf"])

list_pdf = sorted(glob.glob("DYPlots/*NAk4BJets*pdf"))
subprocess.Popen(['pdfunite']+list_pdf+["DYPlots_NAk4BJets.pdf"])
