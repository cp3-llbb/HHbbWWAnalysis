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
    def __init__(self,variable,variable_name,output,plot_MC=False):
        self.outname = output
        self.era     = '2016'
        self.variable = variable
        self.variable_name = variable_name
        self.data_hist_dict = {}
        self.MC_hist_dict = {}
        self.plot_MC = plot_MC
#        self.directories = directories
#        self.variables   = variables 
#        self.era         = era
#        self.outname     = outname

#        print ("Directories to look at :")
#        for d in self.directories:
#            print ("... %s"%d)
#        print ("variables to look at :")
#        for v in self.variables:
#            print ("... %s"%v)
        # inputs : dict
        #   -> key : dir path
        #   -> val : dict
        #       -> key : variable name
        #       -> val : histogram 

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
        options_dict = {k:full_dict["plots"][histogram][k] for k in full_dict["plots"][histogram].keys() & opt_to_keep}

        # Get data per sample #
        sample_dict = {}
        info_to_keep = ['cross-section','generated-events','group','type','era']
        for sample,data in full_dict['files'].items():
            sample_dict[sample] = {k:data[k] for k in data.keys() & info_to_keep}

        return {'luminosity':lumi_dict,'plot_options':options_dict,'samples':sample_dict}
            
    def getHistogram(self,rootfile,histname):
        f = ROOT.TFile(rootfile)
        h = copy.deepcopy(f.Get(histname))
        f.Close()
        return h

    def processHistograms(self):
        colors = [601,634,418,808,402,618] # Must be extended if more required
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
                if data_hist is None:
                    data_hist = ROOT.TH1F(self.variable_name+"DataMinusMCDY",self.variable_name,h.GetNbinsX(),h.GetBinLowEdge(1),h.GetBinLowEdge(h.GetNbinsX())+h.GetBinWidth(h.GetNbinsX()))
                if MC_hist is None:
                    MC_hist = ROOT.TH1F(self.variable_name+"MCDY",self.variable_name,h.GetNbinsX(),h.GetBinLowEdge(1),h.GetBinLowEdge(h.GetNbinsX())+h.GetBinWidth(h.GetNbinsX()))

                # Add histograms depending on type #
                if data_dict['type'] == 'data':
                    data_hist.Add(h)
                elif data_dict['type'] == 'mc':
                    factor = data_dict["cross-section"]*lumi_dict[self.era]/data_dict["generated-events"]
                    if data_dict['group'] != 'DY':
                        data_hist.Add(h,-factor)
                    else:  
                        MC_hist.Add(h,factor)

                    
            # Esthetics #
            data_hist.SetLineWidth(2)
            data_hist.SetLineColor(colors[2*i])
            #data_hist.SetFillColor(colors[2*i])
            data_hist.GetXaxis().SetTitle(plot_dict['x-axis'])
            data_hist.GetYaxis().SetTitle(plot_dict['y-axis'])
            data_hist.SetTitle(self.variable_name)
            MC_hist.SetLineWidth(2)
            MC_hist.SetLineColor(colors[2*i+1])
            #MC_hist.SetFillColor(colors[2*i+1])
            MC_hist.GetXaxis().SetTitle(plot_dict['x-axis'])
            MC_hist.GetYaxis().SetTitle(plot_dict['y-axis'])
            MC_hist.SetTitle(self.variable_name)
            # Normalize to unity #
            data_hist.Scale(1/data_hist.Integral())
            MC_hist.Scale(1/MC_hist.Integral())
            # Save #
            self.data_hist_dict[key] = data_hist
            self.MC_hist_dict[key] = MC_hist

    def saveHistograms(self):
        # Legend #
        legend = ROOT.TLegend(0.5,0.6,0.87,0.87)
        legend.SetHeader("Legend","C")

        plot_ratio = (len(self.data_hist_dict.keys()) == 2 and not self.plot_MC ) or (self.plot_MC and len(self.data_hist_dict.keys()) ==1 )

        # Canvas and pad #
        C = ROOT.TCanvas("c1", "c1", 600, 600)
        pad1 = ROOT.TPad("pad1", "pad1", 0, 0.0, 1, 1.0)
        pad1.SetBottomMargin(0.15)
        pad1.SetLeftMargin(0.15)
        if plot_ratio:
            pad1.SetBottomMargin(0.32)
        pad1.SetGridx()
        pad1.SetGridy()
        pad1.Draw()
        pad1.cd()

        # Get Max values #
        amax = max([h.GetMaximum() for h in self.data_hist_dict.values()]+[h.GetMaximum() for h in self.MC_hist_dict.values()])

        # Plot and save #
        opt = "hist"
        for i,(key,hist) in enumerate(self.data_hist_dict.items()):
            hist.SetMaximum(amax*1.1)
            hist.SetMinimum(0)
            hist.Draw(opt)
            legend.AddEntry(hist,"#splitline{%s}{%s}"%(key,'Data - MC except DY'),'l')
            if i == 0: opt += " same"
            if self.plot_MC:
                self.MC_hist_dict[key].Draw(opt)
                legend.AddEntry(self.MC_hist_dict[key],"#splitline{%s}{%s}"%(key,'MC DY'),'l')

        legend.Draw()
        # Ratio #
        if plot_ratio:
            keys = list(self.data_hist_dict.keys())
            if not self.plot_MC:
                hist1 = self.data_hist_dict[keys[0]]
                hist2 = self.data_hist_dict[keys[1]]
            else: 
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
            pad2.SetGridx()
            pad2.SetGridy()
            pad2.Draw()
            pad2.cd()

            ratio.SetLineColor(ROOT.kBlack)
            ratio.SetMinimum(0.5)
            ratio.SetMaximum(1.5)
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

        C.Print(self.outname+".pdf")
            


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Integration with Cuba')
    #parser.add_argument('--directories', nargs='+',action='store', required=True, type=str,
    #                    help='Directory from bamboo')
    #parser.add_argument('--variables', nargs='+', action='store', required=True, type=str,
    #                    help='Variables (histograms) names')
    #parser.add_argument('--era', action='store', required=True, type=str,
    #                    help='Era : 2016, 2017, 2018')
    #parser.add_argument('--output', action='store', required=True, type=str, 
    #                    help='Output name for pdf (without extension)')
    args = parser.parse_args()


    # ElEl Dilepton pt #
    var_0btagZpeak_dict = {('../full2016_AutoPULeptonTriggerJetMETBtagSF_DYStudy0btag_inZpeak/','Resolved 0 btag (Z peak)') : 'ElEl_HasElElTightZPeakTwoAk4JetsExclusiveResolvedNoBtag_dilepton_pt'}
    var_0btagnoZveto_dict = {('../full2016_AutoPULeptonTriggerJetMETBtagSF_DYStudy0btag_noZveto/','Resolved 0 btag (no Z veto)') : 'ElEl_HasElElTightTwoAk4JetsExclusiveResolvedNoBtag_dilepton_pt'}
    var_2btagZpeak_dict = {('../full2016_AutoPULeptonTriggerJetMETBtagSF_DYStudy2btag_inZpeak/','Resolved 2 btag (Zpeak)') : 'ElEl_HasElElTightZPeakTwoAk4JetsExclusiveResolvedTwoBtags_dilepton_pt'}
    variable_name = 'e^{+}e^{-} Dilepton P_{T}'
    instance = DYControlRegion(var_0btagZpeak_dict,variable_name,"ElEl_dilepton_PT_0btagInZpeak",True)
    instance = DYControlRegion(var_0btagnoZveto_dict,variable_name,"ElEl_dilepton_PT_0btagNoZVeto",True)
    instance = DYControlRegion(var_2btagZpeak_dict,variable_name,"ElEl_dilepton_PT_2btagInZpeak",True)
    instance = DYControlRegion({**var_0btagZpeak_dict,**var_0btagnoZveto_dict,**var_2btagZpeak_dict},variable_name,"ElEl_dilepton_PT",False)

    # ElEl Dilepton pt #
    var_0btagZpeak_dict = {('../full2016_AutoPULeptonTriggerJetMETBtagSF_DYStudy0btag_inZpeak/','Resolved 0 btag (Z peak)') : 'ElEl_HasElElTightZPeakTwoAk4JetsExclusiveResolvedNoBtag_firstlepton_pt'}
    var_0btagnoZveto_dict = {('../full2016_AutoPULeptonTriggerJetMETBtagSF_DYStudy0btag_noZveto/','Resolved 0 btag (no Z veto)') : 'ElEl_HasElElTightTwoAk4JetsExclusiveResolvedNoBtag_firstlepton_pt'}
    var_2btagZpeak_dict = {('../full2016_AutoPULeptonTriggerJetMETBtagSF_DYStudy2btag_inZpeak/','Resolved 2 btag (Zpeak)') : 'ElEl_HasElElTightZPeakTwoAk4JetsExclusiveResolvedTwoBtags_firstlepton_pt'}
    variable_name = 'e^{+}e^{-} first lepton P_{T}'
    instance = DYControlRegion(var_0btagZpeak_dict,variable_name,"ElEl_firstlepton_PT_0btagInZpeak",True)
    instance = DYControlRegion(var_0btagnoZveto_dict,variable_name,"ElEl_firstlepton_PT_0btagNoZVeto",True)
    instance = DYControlRegion(var_2btagZpeak_dict,variable_name,"ElEl_firstlepton_PT_2btagInZpeak",True)
    instance = DYControlRegion({**var_0btagZpeak_dict,**var_0btagnoZveto_dict,**var_2btagZpeak_dict},variable_name,"ElEl_firstlepton_PT",False)


    # MuMu Dilepton pt #
    var_0btagZpeak_dict = {('../full2016_AutoPULeptonTriggerJetMETBtagSF_DYStudy0btag_inZpeak/','Resolved 0 btag (Z peak)') : 'MuMu_HasMuMuTightZPeakTwoAk4JetsExclusiveResolvedNoBtag_dilepton_pt'}
    var_0btagnoZveto_dict = {('../full2016_AutoPULeptonTriggerJetMETBtagSF_DYStudy0btag_noZveto/','Resolved 0 btag (no Z veto)') : 'MuMu_HasMuMuTightTwoAk4JetsExclusiveResolvedNoBtag_dilepton_pt'}
    var_2btagZpeak_dict = {('../full2016_AutoPULeptonTriggerJetMETBtagSF_DYStudy2btag_inZpeak/','Resolved 2 btag (Zpeak)') : 'MuMu_HasMuMuTightZPeakTwoAk4JetsExclusiveResolvedTwoBtags_dilepton_pt'}
    variable_name = '#mu^{+}#mu^{-} Dilepton P_{T}'
    instance = DYControlRegion(var_0btagZpeak_dict,variable_name,"MuMu_dilepton_PT_0btagInZpeak",True)
    instance = DYControlRegion(var_0btagnoZveto_dict,variable_name,"MuMu_dilepton_PT_0btagNoZVeto",True)
    instance = DYControlRegion(var_2btagZpeak_dict,variable_name,"MuMu_dilepton_PT_2btagInZpeak",True)
    instance = DYControlRegion({**var_0btagZpeak_dict,**var_0btagnoZveto_dict,**var_2btagZpeak_dict},variable_name,"MuMu_dilepton_PT",False)

    # MuMu Dilepton pt #
    var_0btagZpeak_dict = {('../full2016_AutoPULeptonTriggerJetMETBtagSF_DYStudy0btag_inZpeak/','Resolved 0 btag (Z peak)') : 'MuMu_HasMuMuTightZPeakTwoAk4JetsExclusiveResolvedNoBtag_firstlepton_pt'}
    var_0btagnoZveto_dict = {('../full2016_AutoPULeptonTriggerJetMETBtagSF_DYStudy0btag_noZveto/','Resolved 0 btag (no Z veto)') : 'MuMu_HasMuMuTightTwoAk4JetsExclusiveResolvedNoBtag_firstlepton_pt'}
    var_2btagZpeak_dict = {('../full2016_AutoPULeptonTriggerJetMETBtagSF_DYStudy2btag_inZpeak/','Resolved 2 btag (Zpeak)') : 'MuMu_HasMuMuTightZPeakTwoAk4JetsExclusiveResolvedTwoBtags_firstlepton_pt'}
    variable_name = '#mu^{+}#mu^{-} first lepton P_{T}'
    instance = DYControlRegion(var_0btagZpeak_dict,variable_name,"MuMu_firstlepton_PT_0btagInZpeak",True)
    instance = DYControlRegion(var_0btagnoZveto_dict,variable_name,"MuMu_firstlepton_PT_0btagNoZVeto",True)
    instance = DYControlRegion(var_2btagZpeak_dict,variable_name,"MuMu_firstlepton_PT_2btagInZpeak",True)
    instance = DYControlRegion({**var_0btagZpeak_dict,**var_0btagnoZveto_dict,**var_2btagZpeak_dict},variable_name,"MuMu_firstlepton_PT",False)










    sys.exit()
    # ElEl Dilepton pt #
    variable = {('../full2016_AutoPULeptonTriggerJetMETBtagSF_DYStudy0btag_inZpeak/','Resolved 0 btag') : 'ElEl_HasElElTightZPeakTwoAk4JetsExclusiveResolvedNoBtag_dilepton_pt',
                ('../full2016_AutoPULeptonTriggerJetMETBtagSF_DYStudy2btag_inZpeak/','Resolved 2 btag') : 'ElEl_HasElElTightZPeakTwoAk4JetsExclusiveResolvedTwoBtags_dilepton_pt'}
    outputname = 'ElEl_dilepton_pt.pdf'
    variable_name = 'e^{+}e^{-} Dilepton P_{T}'
    instance = DYControlRegion(variable,variable_name,outputname)

    # MuMu first lepton pt #
    variable = {('../full2016_AutoPULeptonTriggerJetMETBtagSF_DYStudy0btag_inZpeak/','Resolved 0 btag') : 'MuMu_HasMuMuTightZPeakTwoAk4JetsExclusiveResolvedNoBtag_firstlepton_pt',
                ('../full2016_AutoPULeptonTriggerJetMETBtagSF_DYStudy2btag_inZpeak/','Resolved 2 btag') : 'MuMu_HasMuMuTightZPeakTwoAk4JetsExclusiveResolvedTwoBtags_firstlepton_pt'}
    outputname = 'MuMu_firstlepton_pt.pdf'
    variable_name = '#mu^{+}#mu^{-} First lepton P_{T}'
    instance = DYControlRegion(variable,variable_name,outputname)

    # ElEl first lepton pt #
    variable = {('../full2016_AutoPULeptonTriggerJetMETBtagSF_DYStudy0btag_inZpeak/','Resolved 0 btag') : 'ElEl_HasElElTightZPeakTwoAk4JetsExclusiveResolvedNoBtag_firstlepton_pt',
                ('../full2016_AutoPULeptonTriggerJetMETBtagSF_DYStudy2btag_inZpeak/','Resolved 2 btag') : 'ElEl_HasElElTightZPeakTwoAk4JetsExclusiveResolvedTwoBtags_firstlepton_pt'}
    outputname = 'ElEl_firstlepton_pt.pdf'
    variable_name = 'e^{+}e^{-} First lepton P_{T}'
    instance = DYControlRegion(variable,variable_name,outputname)

