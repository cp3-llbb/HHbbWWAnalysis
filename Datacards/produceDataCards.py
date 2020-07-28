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

from datasetConfig import * 


class DataCard:
    def __init__(self,path,yamlName,groups,hist_conv,era,use_syst=False):
        self.path = path
        self.groups = groups
        self.hist_conv = hist_conv
        self.era = era
        self.use_syst = use_syst
        self.content = {k:{g:None for g in self.groups.keys()} for k in self.hist_conv.keys()}

        self.yaml_dict = self.loadYaml(os.path.join(self.path,yamlName))
        self.loopOverFiles()
        #pprint (self.content)
        self.saveDatacard()

    def loopOverFiles(self):
        for f in glob.glob(os.path.join(self.path,'results','*.root')):
            if '__skeleton__' in f:
                continue
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
            br = self.yaml_dict["samples"][sample]["branching-ratio"] if "branching-ratio" in self.yaml_dict["samples"][sample].keys() else 1

        # Get list of hist names #
        list_histnames = self.getHistList(f)

        # Loop through hists #
        hist_dict = {}
        for datacardname, histname in self.hist_conv.items():
            # Check #
            if not histname in list_histnames:
                print ("Could not find hist %s in %s"%(histname,rootfile))
                continue
            listsyst = [hn for hn in list_histnames if histname in hn and '__' in hn] if self.use_syst else []
            h = self.getHistogram(f,histname,listsyst)
            # Normalize hist to data #
            if sample_type == "mc" or sample_type == "signal":
                h.Scale(lumi*xsec*br/sumweight)
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
        path_datacard = os.path.join(os.path.abspath(os.path.dirname(__file__)),'datacards_%s'%self.era)
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
    #instance = DataCard(*returnConfig('2016'),'2016')
    #instance = DataCard(*returnConfig('2017'),'2017')
    #instance = DataCard(*returnConfig('2018'),'2018')
    instance = DataCard(*returnZPeakConfig('2016'),'2016')
