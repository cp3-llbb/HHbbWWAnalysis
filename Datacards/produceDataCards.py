import os
import sys
import glob
import copy
import yaml
import argparse
import numpy as np
import math
from itertools import chain
from pprint import pprint
import ROOT


class DataCard:
    def __init__(self,datacardName,path,yamlName,groups,hist_conv,era,use_syst=False,root_subdir=None,pseudodata=False):
        self.datacardName = datacardName
        self.path = path
        self.groups = groups
        self.hist_conv = hist_conv
        self.era = str(era)
        self.use_syst = use_syst
        self.root_subdir = root_subdir
        self.pseudodata = pseudodata

        if self.pseudodata:
            self.groups = self.generatePseudoData(self.groups)
        self.content = {k:{g:None for g in self.groups.keys()} for k in self.hist_conv.keys()}

        self.yaml_dict = self.loadYaml(os.path.join(self.path,yamlName))
        self.loopOverFiles()
        if self.pseudodata:
            self.roundFakeData()
        self.saveDatacard()

    @staticmethod
    def generatePseudoData(groups):
        newg = {k:v for k,v in groups.items() if k!='data_obs'} 
        newg['data_obs'] = list(chain.from_iterable([v for v in newg.values()]))
        newg['data_real'] = groups['data_obs']
        return newg

    def loopOverFiles(self):
        for f in glob.glob(os.path.join(self.path,'results','*.root')):
            if '__skeleton__' in f:
                continue
            sample = os.path.basename(f)
            # Check if in the group list #
            groups = self.findGroup(sample)
            if self.pseudodata and 'data_real' in groups:
                continue
            if len(groups) == 0:
                raise RuntimeError("[WARNING] Could not find sample %s in group list"%sample)

            hist_dict = self.getHistograms(f)
            for group in groups:
                self.addSampleToGroup(copy.deepcopy(hist_dict),group)
                # deepcopy is needed so that if the histogram is in two groups
                # acting on one version will not change the other

    def addSampleToGroup(self,hist_dict,group):
        for histname,hist in hist_dict.items():
            if self.content[histname][group] is None:
                self.content[histname][group] = hist
            else:
                self.content[histname][group].Add(hist)

    def findGroup(self,sample):
        gr = []
        for key,val in self.groups.items():
            if not isinstance(val,list):
                raise RuntimeError("Group %s does not consist in a list"%key)
            if sample in val:
                gr.append(key)
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
        else:
            xsec = None
            sumweight = None
            br = None

        # Get list of hist names #
        list_histnames = self.getHistList(f)

        # Loop through hists #
        hist_dict = {}
        for datacardname, histnames in self.hist_conv.items():
            if not isinstance(histnames,list):
                histnames = [histnames]
            hist_dict[datacardname] = None
            for histname in histnames:
                # Check #
                if not histname in list_histnames:
                    print ("Could not find hist %s in %s"%(histname,rootfile))
                    continue
                listsyst = [hn for hn in list_histnames if histname in hn and '__' in hn] if self.use_syst else []
                # Nominal histogram #
                h = self.getHistogram(f,histname,lumi,br,xsec,sumweight)
                if hist_dict[datacardname] is None:
                    hist_dict[datacardname] = copy.deepcopy(h)
                else:
                    hist_dict[datacardname].Add(h)
                # Systematic histograms #
                for syst in listsyst:
                    systName = syst.split('__')[-1]
                    systName.replace('up','Up')
                    systName.replace('down','Down')
                    h = self.getHistogram(f,syst,lumi,br,xsec,sumweight)
                    hist_dict[datacardname+'_'+systName] = h
                # Save in dict #
        f.Close()
        return hist_dict

    @staticmethod
    def getHistogram(f,histnom,lumi=None,xsec=None,br=None,sumweight=None):
        # Get hist #
        h = copy.copy(f.Get(histnom))
        # Normalize hist to data #
        if lumi is not None and xsec is not None and br is not None and sumweight is not None:
            h.Scale(lumi*xsec*br/sumweight)
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

    def roundFakeData(self):
        for hist_dict in self.content.values():
            for group, hist in hist_dict.items():
                if group == 'data_obs':
                    for i in range(0,hist.GetNbinsX()+2):
                        if hist.GetBinContent(i) > 0:
                            hist.SetBinContent(i,round(hist.GetBinContent(i)))
                            hist.SetBinError(i,math.sqrt(hist.GetBinContent(i)))
                        else:
                            hist.SetBinContent(i,0)
                            hist.SetBinError(i,0)

    def saveDatacard(self):
        path_datacard = os.path.join(os.path.abspath(os.path.dirname(__file__)),self.datacardName)
        if self.pseudodata:
            path_datacard += "_pseudodata"

        if not os.path.exists(path_datacard):
            os.makedirs(path_datacard)
        for key, gdict in self.content.items():
            # Create file #
            filename = os.path.join(path_datacard,key+'_'+self.era+'.root')
            f = ROOT.TFile(filename,'recreate')
            if self.root_subdir is not None:
                d = f.mkdir(self.root_subdir,self.root_subdir)
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
    parser = argparse.ArgumentParser(description='Produce datacards')
    parser.add_argument('--yaml', action='store', required=True, type=str,
                        help='Yaml containing parameters')
    parser.add_argument('--pseudodata', action='store_true', required=False, default=False,
                        help='Whether to use pseudo data (data = sum of MC)')
    args = parser.parse_args()

    if args.yaml is None:
        raise RuntimeError("Must provide the YAML file")
    with open(args.yaml,'r') as handle:
        f = yaml.load(handle)
    instance = DataCard(datacardName=args.yaml.replace('.yml',''),pseudodata=args.pseudodata,**f)
