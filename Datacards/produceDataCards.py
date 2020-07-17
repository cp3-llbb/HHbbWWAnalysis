import os
import sys
import glob
import copy
import yaml
import root_numpy
import math
import numpy as np
from pprint import pprint

import ROOT


class DataCard:
    def __init__(self,path,groups,hist_conv,era):
        self.path = path
        self.groups = groups
        self.hist_conv = hist_conv
        self.era = era
        self.content = {k:{g:None for g in self.groups.keys()} for k in self.hist_conv.keys()}

        self.yaml_dict = self.loadYaml(os.path.join(self.path,'plots.yml'))
        self.loopOverFiles()
        pprint (self.content)
        self.saveDatacard()

    def loopOverFiles(self):
        for f in glob.glob(os.path.join(self.path,'results','*.root')):
            sample = os.path.basename(f)
            # Check if in the group list #
            group = self.findGroup(sample)
            if group is None:
                print ("[WARNING] Could not find sample %s in group list"%sample)
                continue

            hist_dict = self.getHistograms(f)
            self.addSampleToGroup(hist_dict,group)

    def addSampleToGroup(self,hist_dict,group):
        for histname,hist in hist_dict.items():
            if self.content[histname][group] is None:
                self.content[histname][group] = hist
            else:
                self.content[histname][group].Add(hist)

    def findGroup(self,sample):
        gr = None
        for key,val in self.groups.items():
            if not isinstance(val,list):
                raise RuntimeError("Unvalid group dict")
            if sample in val:
                gr = key
                break
        return gr
                
    def getHistograms(self,rootfile):
        f = ROOT.TFile(rootfile)
        sample = os.path.basename(rootfile)
        # Get config info #
        lumi = self.yaml_dict["luminosity"][self.era]
        sample_type = self.yaml_dict["samples"][sample]['type']
        if sample_type == "mc" or sample_type == "signal":
            xsec = self.yaml_dict["samples"][sample]['cross-section']
            sumweight = self.yaml_dict["samples"][sample]['generated-events']

        # Get list of hist names #
        list_histnames = self.getHistList(f)

        # Loop through hists #
        hist_dict = {}
        for datacardname, histname in self.hist_conv.items():
            # Check #
            if not histname in list_histnames:
                print ("Could not find hist %s in %s"%(histname,rootfile))
                continue
            listsyst = [hn for hn in list_histnames if histname in hn and '__' in hn]
            h = self.getHistogram(f,histname,listsyst)
            # Normalize hist to data #
            if sample_type == "mc" or sample_type == "signal":
                h.Scale(lumi*xsec/sumweight)
            # Save in dict #
            hist_dict[datacardname] = h
        f.Close()
        return hist_dict

    @staticmethod
    def getHistogram(f,histnom,listsyst):
        # Get hist #
        h = copy.deepcopy(f.Get(histnom))
        if len(listsyst) != 0:
            dict_syst = {}
            nom_cont = root_numpy.hist2array(h,include_overflow=True,copy=True)
            for syst in listsyst:
                s = syst.replace(histnom,'').replace('__','')
                hs = f.Get(syst)
                dict_syst[s] = np.abs(root_numpy.hist2array(hs,include_overflow=True,copy=True)-nom_cont)
            # Get maximum var for each sys and add quadraticaly #
            syst_names = [n[:-2] for n in dict_syst.keys() if n.endswith('up')]
            syst_tot_squared = np.zeros(nom_cont.shape[0])
            for syst_name in syst_names:
                maxvar = np.maximum(dict_syst[syst_name+'down'],dict_syst[syst_name+'up'])
                syst_tot_squared += np.power(maxvar,2)
            for i in range(h.GetNbinsX()):
                err = h.GetBinError(i)
                new_err = math.sqrt(err**2+syst_tot_squared[i])
                h.SetBinError(i,new_err)
        return h
             
    
    @staticmethod
    def getHistList(f):
        l = []
        for key in f.GetListOfKeys():
            obj = key.ReadObj()
            if "TH1" in obj.ClassName():
                l.append(obj.GetName())
        return l
        

    def loadYaml(self,yaml_path):
        # Parse YAML #
        with open(yaml_path,"r") as handle:
            full_dict = yaml.load(handle,Loader=yaml.FullLoader)
        # Get Lumi per era #  
        lumi_dict = full_dict["configuration"]["luminosity"]

        # Get data per sample #
        sample_dict = {}
        info_to_keep = ['cross-section','generated-events','group','type','era']
        for sample,data in full_dict['files'].items():
            sample_dict[sample] = {k:data[k] for k in data.keys() & info_to_keep}

        return {'luminosity':lumi_dict,'samples':sample_dict}


    def saveDatacard(self):
        path_datacard = os.path.join(os.path.abspath(os.path.dirname(__file__)),'datacards')
        if not os.path.exists(path_datacard):
            os.makedirs(path_datacard)
        for key, gdict in self.content.items():
            # Create file #
            filename = os.path.join(path_datacard,key+'_'+self.era+'.root')
            f = ROOT.TFile(filename,'recreate')
            d = f.mkdir("HH_2l_0tau","HH_2l_0tau")
            d.cd()
            for group,hist in gdict.items():
                if hist is None:
                    continue
                hist.SetTitle(group)
                hist.SetName(group)
                hist.Write(group)
            f.Write()
            f.Close()
            print ("Saved file %s"%filename)




if __name__=="__main__":
    path = '/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/BambooOutputHHtobbWW/full2016_tight_noSyst/'
    groups = {'TTZ':   ['TTZToQQ.root','TTZToLLNuNu.root'],
              'TH' :   [],
              'VH' :   ['HZJ_HToWW.root',
                        'ggZH_HToBB_ZToLL.root',
                        'ggZH_HToBB_ZToNuNu.root'],
              'TT' :   ['TTTo2L2Nu.root',
                        'TTToSemiLeptonic.root',
                        'TTToHadronic.root'],
              'TTW':   ['TTWJetsToQQ.root',
                        'TTWJetsToLNu.root'],
              'WW' :   ['WWToLNuQQ.root',
                        'WWTo2L2Nu.root',
                        'WJetsToLNu.root'],
              'TTH':   ['ttHTobb.root',
                        'ttHToNonbb.root'],
              'ZZ' :   ['ZZTo2L2Nu.root',
                        'ZZTo2L2Q.root',
                        'ZZTo4L.root'],
              'DY' :   ['DYJetsToLL_M-10to50.root',
                        'DYToLL_0J.root',
                        'DYToLL_1J.root',
                        'DYToLL_2J.root'],
              'Other': ['ST_tW_top_5f.root',
                        'ST_tW_antitop_5f.root',
                        'ST_tchannel_antitop_4f.root',
                        'ST_tchannel_top_4f.root',
                        'WWW.root',
                        'WWZ.root',
                        'WZZ.root',
                        'ZZZ.root'],
              'WZ' :   ['WZTo2L2Q.root',
                        'WZTo1L3Nu.root',
                        'WZ1L1Nu2Q.root',
                        'WZ1L1Nu2Q.root',
                        'WZTo3LNu.root'],
              'TTWW' : [],
              'data_fake': [],
              'data_obs' : ['DoubleEGamma_2016B.root',
                            'DoubleEGamma_2016C.root',
                            'DoubleEGamma_2016D.root',
                            'DoubleEGamma_2016E.root',
                            'DoubleEGamma_2016F.root',
                            'DoubleEGamma_2016G.root',
                            'DoubleEGamma_2016H.root',
                            'DoubleMuon_2016B.root',
                            'DoubleMuon_2016C.root',
                            'DoubleMuon_2016D.root',
                            'DoubleMuon_2016E.root',
                            'DoubleMuon_2016F.root',
                            'DoubleMuon_2016G.root',
                            'DoubleMuon_2016H.root',
                            'MuonEG_2016B.root',
                            'MuonEG_2016C.root',
                            'MuonEG_2016D.root',
                            'MuonEG_2016E.root',
                            'MuonEG_2016F.root',
                            'MuonEG_2016G.root',
                            'MuonEG_2016H.root',
                            'SingleElectron_2016B.root',
                            'SingleElectron_2016C.root',
                            'SingleElectron_2016D.root',
                            'SingleElectron_2016E.root',
                            'SingleElectron_2016F.root',
                            'SingleElectron_2016G.root',
                            'SingleElectron_2016H.root',
                            'SingleMuon_2016B.root',
                            'SingleMuon_2016C.root',
                            'SingleMuon_2016D.root',
                            'SingleMuon_2016E.root',
                            'SingleMuon_2016F.root',
                            'SingleMuon_2016G.root',
                            'SingleMuon_2016H.root']}    

    hists = {'HH_cat_ee_1b_jet1_pt' : 'ElEl_HasElElTightTwoAk4JetsExclusiveResolvedOneBtag_leadbjet_pt',
             'HH_cat_mm_1b_jet1_pt' : 'MuMu_HasMuMuTightTwoAk4JetsExclusiveResolvedOneBtag_leadbjet_pt',
             'HH_cat_em_1b_jet1_pt' : 'ElMu_HasElMuTightTwoAk4JetsExclusiveResolvedOneBtag_leadbjet_pt',
             'HH_cat_ee_2b_jet1_pt' : 'ElEl_HasElElTightTwoAk4JetsExclusiveResolvedTwoBtags_leadbjet_pt',
             'HH_cat_mm_2b_jet1_pt' : 'MuMu_HasMuMuTightTwoAk4JetsExclusiveResolvedTwoBtags_leadbjet_pt',
             'HH_cat_em_2b_jet1_pt' : 'ElMu_HasElMuTightTwoAk4JetsExclusiveResolvedTwoBtags_leadbjet_pt',
             'HH_cat_ee_1b_lep1_pt' : 'ElEl_HasElElTightTwoAk4JetsExclusiveResolvedOneBtag_firstlepton_pt',
             'HH_cat_mm_1b_lep1_pt' : 'MuMu_HasMuMuTightTwoAk4JetsExclusiveResolvedOneBtag_firstlepton_pt',
             'HH_cat_em_1b_lep1_pt' : 'ElMu_HasElMuTightTwoAk4JetsExclusiveResolvedOneBtag_firstlepton_pt',
             'HH_cat_ee_2b_lep1_pt' : 'ElEl_HasElElTightTwoAk4JetsExclusiveResolvedTwoBtags_firstlepton_pt',
             'HH_cat_mm_2b_lep1_pt' : 'MuMu_HasMuMuTightTwoAk4JetsExclusiveResolvedTwoBtags_firstlepton_pt',
             'HH_cat_em_2b_lep1_pt' : 'ElMu_HasElMuTightTwoAk4JetsExclusiveResolvedTwoBtags_firstlepton_pt',
             'HH_cat_ee_1b_jet1_eta' : 'ElEl_HasElElTightTwoAk4JetsExclusiveResolvedOneBtag_leadbjet_eta',
             'HH_cat_mm_1b_jet1_eta' : 'MuMu_HasMuMuTightTwoAk4JetsExclusiveResolvedOneBtag_leadbjet_eta',
             'HH_cat_em_1b_jet1_eta' : 'ElMu_HasElMuTightTwoAk4JetsExclusiveResolvedOneBtag_leadbjet_eta',
             'HH_cat_ee_2b_jet1_eta' : 'ElEl_HasElElTightTwoAk4JetsExclusiveResolvedTwoBtags_leadbjet_eta',
             'HH_cat_mm_2b_jet1_eta' : 'MuMu_HasMuMuTightTwoAk4JetsExclusiveResolvedTwoBtags_leadbjet_eta',
             'HH_cat_em_2b_jet1_eta' : 'ElMu_HasElMuTightTwoAk4JetsExclusiveResolvedTwoBtags_leadbjet_eta',
             'HH_cat_ee_1b_lep1_eta' : 'ElEl_HasElElTightTwoAk4JetsExclusiveResolvedOneBtag_firstlepton_eta',
             'HH_cat_mm_1b_lep1_eta' : 'MuMu_HasMuMuTightTwoAk4JetsExclusiveResolvedOneBtag_firstlepton_eta',
             'HH_cat_em_1b_lep1_eta' : 'ElMu_HasElMuTightTwoAk4JetsExclusiveResolvedOneBtag_firstlepton_eta',
             'HH_cat_ee_2b_lep1_eta' : 'ElEl_HasElElTightTwoAk4JetsExclusiveResolvedTwoBtags_firstlepton_eta',
             'HH_cat_mm_2b_lep1_eta' : 'MuMu_HasMuMuTightTwoAk4JetsExclusiveResolvedTwoBtags_firstlepton_eta',
             'HH_cat_em_2b_lep1_eta' : 'ElMu_HasElMuTightTwoAk4JetsExclusiveResolvedTwoBtags_firstlepton_eta',
            }

    instance = DataCard(path,groups,hists,'2016')
