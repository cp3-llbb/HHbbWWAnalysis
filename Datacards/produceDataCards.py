import os
import sys
import glob
import copy
from pprint import pprint

import ROOT


class DataCard:
    def __init__(self,path,groups,hist_conv,era):
        self.path = path
        self.groups = groups
        self.hist_conv = hist_conv
        self.era = era
        self.content = {k:{g:None for g in self.groups.keys()} for k in self.hist_conv.keys()}

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
        hist_dict = {}
        for datacardname, histname in self.hist_conv.items():
            if f.GetListOfKeys().Contains(histname):
                h = copy.deepcopy(f.Get(histname))
            else:
                print ("Could not find hist %s in %s"%(histname,rootfile))
                h = None
            hist_dict[datacardname] = h
        f.Close()
        return hist_dict

    def saveDatacard(self):
        path_datacard = os.path.join(os.path.abspath(os.path.dirname(__file__)),'datacards')
        if not os.path.exists(path_datacard):
            os.makedirs(path_datacard)
        for key, gdict in self.content.items():
            # Create file #
            filename = os.path.join(path_datacard,key+'_'+self.era+'.root')
            f = ROOT.TFile(filename,'recreate')
            d = ROOT.TDirectory("HH_2l_0tau","HH_2l_0tau")
            d.cd()
            for group,hist in gdict.items():
                if hist is None:
                    continue
                hist.SetTitle(group)
                hist.Write(group)
            f.Write()
            f.Close()
            print ("Saved file %s"%filename)




if __name__=="__main__":
    path = '/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/BambooOutputHHtobbWW/full2016_tight_noSyst/'
    groups = {'TTZ':   ['TTZToQQ.root','TTZToLLNuNu.root'],
              'TH' :   [],
              'VH' :   ['HZJ_HToWW.root','ggZH_HToBB_ZToLL.root','ggZH_HToBB_ZToNuNu.root'],
              'TT' :   ['TTTo2L2Nu.root','TTToSemiLeptonic.root','TTToHadronic.root'],
              'TTW':   ['TTWJetsToQQ.root','TTWJetsToLNu.root'],
              'WW' :   ['WWToLNuQQ.root','WWTo2L2Nu.root','WJetsToLNu.root'],
              'TTH':   ['ttHTobb.root','ttHToNonbb.root'],
              'ZZ' :   ['ZZTo2L2Nu.root','ZZTo2L2Q.root','ZZTo4L.root'],
              'DY' :   ['DYJetsToLL_M-10to50.root','DYToLL_0J.root','DYToLL_1J.root','DYToLL_2J.root'],
              'Other': ['ST_tW_top_5f.root','ST_tW_antitop_5f.root','ST_tchannel_antitop_4f.root','ST_tchannel_top_4f.root','WWW.root','WWZ.root','WZZ.root','ZZZ.root'],
              'WZ' :   ['WZTo2L2Q.root','WZTo1L3Nu.root','WZ1L1Nu2Q.root','WZ1L1Nu2Q.root','WZTo3LNu.root'],
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
