import os
import sys
import glob
import copy
import yaml
import argparse
import subprocess
import numpy as np
import math
from itertools import chain
from pprint import pprint
import numpy as np
import enlighten
import ROOT

from IPython import embed
ROOT.gROOT.SetBatch(True)

class DataCard:
    def __init__(self,datacardName=None,path=None,yamlName=None,groups=None,hist_conv=None,era=None,use_syst=False,root_subdir=None,DYnonclosure=None,pseudodata=False,rebin=None,produce_plots=False,**kwargs):
        self.datacardName   = datacardName
        self.path           = path
        self.groups         = groups
        self.hist_conv      = hist_conv
        self.era            = str(era)
        self.use_syst       = use_syst
        self.root_subdir    = root_subdir
        self.pseudodata     = pseudodata
        self.DYnonclosure   = DYnonclosure
        self.rebin          = rebin
        self.produce_plots  = produce_plots
        
        self.yaml_dict = self.loadYaml(self.path,yamlName)

        with open('systematics.yml','r') as handle:
            self.systConvention = yaml.load(handle,Loader=yaml.FullLoader)

        if self.pseudodata:
            print ('Will use pseudodata')
            self.groups = self.generatePseudoData(self.groups)
            if self.datacardName is not None:
                self.datacardName += "_pseudodata"
        self.content = {histName:{g:{} for g in self.groups.keys()} for histName in self.hist_conv.keys()}
        self.systPresent = {histName:{group:[] for group in self.groups.keys()} for histName in self.hist_conv.keys()}

        self.loopOverFiles()
        if self.pseudodata:
            self.roundFakeData()
        if self.DYnonclosure is not None:
            self.applyDYNonClosure()
        if self.rebin is not None:
            self.applyRebinning()
        if self.datacardName is not None:
            self.saveDatacard()
        if self.produce_plots:
            self.preparePlotIt()

    def generatePseudoData(self,groups):
        back_samples = [sample for sample,sampleCfg in self.yaml_dict['samples'].items() if sampleCfg['type']=='mc']
        newg = {k:v for k,v in groups.items() if k!='data_obs'} 
        newg['data_obs'] = {k:v for k,v in groups['data_obs'].items() if k != 'files'}
        newg['data_obs']['files'] = [sample for sample,sampleCfg in self.yaml_dict['samples'].items() if sampleCfg['type']=='mc']
        newg['data_obs']['legend'] = 'pseudo-data'
        newg['data_real'] = copy.deepcopy(groups['data_obs'])
        return newg

    def loopOverFiles(self):
        self.fileHistList = []
        for val in self.hist_conv.values():
            if isinstance(val,list):
                self.fileHistList.extend(val)
            else:
                self.fileHistList.append(val)
        if isinstance(self.path,list):
            files = set()
            for path in self.path:
                files = files.union(set([f for f in glob.glob(os.path.join(path,'results','*.root'))]))
            files = list(files)
        else:
            files = glob.glob(os.path.join(self.path,'results','*.root'))
        pbar = enlighten.Counter(total=len(files), desc='Progress', unit='files')
        for f in sorted(files):
            pbar.update()
            if '__skeleton__' in f:
                continue
            sample = os.path.basename(f)
            # Check if in the group list #
            groups = self.findGroup(sample)
            if self.pseudodata and 'data_real' in groups:
                continue
            if len(groups) == 0:
                print("[WARNING] Could not find sample %s in group list"%sample)
                continue
            hist_dict = self.getHistograms(f)
            for histName in hist_dict.keys():
                if len(hist_dict[histName]) == 0:
                    continue
                for group in groups:
                    self.systPresent[histName][group].append({'systematics': [systname for systname in hist_dict[histName].keys() if systname != 'nominal'],
                                                              'nominal': copy.deepcopy(hist_dict[histName]['nominal'])})

            for group in groups:
                self.addSampleToGroup(copy.deepcopy(hist_dict),group)
                # deepcopy is needed so that if the histogram is in two groups
                # acting on one version will not change the other

        for histName in self.content.keys(): # renaming to avoid overwritting ROOT warnings
            for group in self.content[histName].keys():
                for systName,hist in self.content[histName][group].items():
                    if self.pseudodata and group == 'data_real':
                        continue
                    if hist is None:
                        continue
                    hist.SetName(hist.GetName()+'_'+group+'__'+systName)

        self.correctMissingSystematics()

    def addSampleToGroup(self,hist_dict,group):
        for histname,hists in hist_dict.items():
            if len(hists) == 0:
                continue
            nominal = hists['nominal']
            if not 'nominal' in self.content[histname][group].keys():
                self.content[histname][group]['nominal'] = copy.deepcopy(nominal)
            else:
                self.content[histname][group]['nominal'].Add(nominal)
            if self.use_syst:
                for systName in hists.keys():
                    if systName == 'nominal':
                        continue
                    hist = hists[systName]
                    if systName not in self.content[histname][group].keys():
                        self.content[histname][group][systName] = copy.deepcopy(hist)
                    else:
                        self.content[histname][group][systName].Add(hist)

    def correctMissingSystematics(self):
        for histName in self.systPresent.keys():
            for group in self.systPresent[histName].keys():
                for sampleDict in self.systPresent[histName][group]:
                    for systName in self.content[histName][group].keys():
                        if systName == "nominal":
                            continue
                        if systName not in sampleDict['systematics']:
                            self.content[histName][group][systName].Add(sampleDict['nominal'])

    def findGroup(self,sample):
        group_of_sample = []
        for group in self.groups.keys():
            if 'files' not in self.groups[group]:
                raise RuntimeError("No 'files' item in group {}".format(group))
            files = self.groups[group]['files']
            if not isinstance(files,list):
                raise RuntimeError("Group %s does not consist in a list"%group)
            if sample in files:
                group_of_sample.append(group)
        return group_of_sample
                
    def getHistograms(self,rootfile):
        f = ROOT.TFile(rootfile)
        sample = os.path.basename(rootfile)
        # Get config info #
        lumi = self.yaml_dict["luminosity"][self.era]
        sample_type = self.yaml_dict["samples"][sample]['type']
        if sample_type == "mc" or sample_type == "signal":
            xsec = self.yaml_dict["samples"][sample]['cross-section']
            sumweight = self.yaml_dict["samples"][sample]['generated-events']
            br = self.yaml_dict["samples"][sample]["branching-ratio"] if "branching-ratio" in self.yaml_dict["samples"][sample].keys() else 1.
        else:
            xsec = None
            sumweight = None
            br = None
        # Get list of hist names #
        list_histnames = []
        for key in f.GetListOfKeys():
            keyName = key.GetName()
            if not self.use_syst and '__' in keyName:
                continue
            if any([histname in keyName for histname in self.fileHistList]):
                list_histnames.append(keyName)

        # Loop through hists #
        hist_dict = {}
        for datacardname, histnames in self.hist_conv.items():
            hist_dict[datacardname] = {}
            if not isinstance(histnames,list):
                histnames = [histnames]
            # Get all defined systematics in the file histogram #
            systPresentInFile = []
            if self.use_syst:
                for histname in list_histnames:
                    if '__' not in histname:
                        continue
                    hName, sName = histname.split('__')
                    if hName not in histnames:
                        continue
                    if sName not in systPresentInFile:
                        systPresentInFile.append(sName)
                # this is because if several histograms are to be summed in the conversion
                # and one has a systematic that the other does not have, one needs to add 
                # the nominal to fill up this missing systematic

            for histname in histnames:
                # Check #
                if not histname in list_histnames:
                    #print ("Could not find hist %s in %s"%(histname,rootfile))
                    continue
                # Nominal histogram #
                hnom = self.getHistogram(f,histname,lumi,br,xsec,sumweight)
                if not 'nominal' in hist_dict[datacardname].keys():
                    hist_dict[datacardname]['nominal'] = copy.deepcopy(hnom)
                else:
                    hist_dict[datacardname]['nominal'].Add(hnom)
                # Systematic histograms #
                listsyst = [histname + '__' + systName for systName in systPresentInFile]
                for syst in listsyst:
                    if syst in list_histnames: # systematic is there
                        h = self.getHistogram(f,syst,lumi,br,xsec,sumweight) 
                    else: # systematic for histogram is not present, use the nominal
                        h = copy.deepcopy(hnom) 
                    systName = syst.split('__')[-1]
                    if systName.endswith('up'):
                        systName = systName[:-2]+"Up"
                    elif systName.endswith('down'):
                        systName = systName[:-4]+"Down"
                    else:
                        raise RuntimeError("Could not understand systematics {}".format(systName))
                    if not systName in hist_dict[datacardname].keys():
                        hist_dict[datacardname][systName] = copy.deepcopy(h)
                    else:
                        hist_dict[datacardname][systName].Add(h)
        f.Close()
        return hist_dict

    @staticmethod
    def getHistogram(f,histnom,lumi=None,xsec=None,br=None,sumweight=None):
        # Get hist #
        h = copy.deepcopy(f.Get(histnom))
        # Normalize hist to data #
        if lumi is not None and xsec is not None and br is not None and sumweight is not None:
            h.Scale(lumi*xsec*br/sumweight)
        return h
             
    def loadYaml(self,paths,yamlName):
        if not isinstance(paths,list):
            paths = [paths]
        yamlDict = {}
        for path in paths:
            yamlPath = os.path.join(path,yamlName)
            if not os.path.exists(yamlPath):
                print ("{} -> not found, skipped".format(yamlPath))
                continue
            # Parse YAML #
            with open(yamlPath,"r") as handle:
                full_dict = yaml.load(handle,Loader=yaml.FullLoader)
            # Get Lumi per era #  
            lumi_dict = full_dict["configuration"]["luminosity"]
            if 'luminosity' not in yamlDict.keys():
                yamlDict['luminosity'] = lumi_dict
            else:
                yamlDict['luminosity'].update(lumi_dict)

            # Get data per sample #
            info_to_keep = ['cross-section','generated-events','group','type','era','branching-ratio']
            if 'samples' not in yamlDict.keys():
                sample_dict = {}
                for sample,data in full_dict['files'].items():
                    sample_dict[sample] = {k:data[k] for k in data.keys() & info_to_keep}
                yamlDict['samples'] = sample_dict
            else:
                for sample,data in full_dict['files'].items():
                    if sample not in yamlDict['samples'].keys():
                        yamlDict['samples'].update({sample:data})
                

            # Get plot options #
            if 'plots' not in yamlDict.keys(): 
                yamlDict['plots'] = full_dict['plots']
            else:
                yamlDict['plots'].update(full_dict['plots'])

        return yamlDict

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
        if self.datacardName is None:
            raise RuntimeError("Datacard name is not set")
        path_datacard = os.path.join(os.path.abspath(os.path.dirname(__file__)),self.datacardName)

        if not os.path.exists(path_datacard):
            os.makedirs(path_datacard)
        for histName in self.content.keys():
            filename = os.path.join(path_datacard,histName+'_'+self.era+'.root')
            f = ROOT.TFile(filename,'recreate')
            if self.root_subdir is not None:
                d = f.mkdir(self.root_subdir,self.root_subdir)
                d.cd()
            for group in self.content[histName].keys():
                for systName,hist in self.content[histName][group].items():
                    if self.pseudodata and group == 'data_real':
                        continue
                    if hist is None:
                        continue
                    for i in range(1,hist.GetNbinsX()+1):
                        if hist.GetBinContent(i) < 0.:
                            hist.SetBinContent(i,0.) 
                            hist.SetBinError(i,0.) 
                    
                    if systName ==  "nominal":
                        save_name = group
                    else:
                        if systName.endswith('Up'):
                            baseSyst = systName[:-2]
                        elif systName.endswith('Down'):
                            baseSyst = systName[:-4]
                        else:
                            raise RuntimeError("Problem with syst {} : cannot find if Up or Down".format(systName))
                        if baseSyst not in self.systConvention.keys():
                            raise RuntimeError("Could not find {} in systematic dict".format(baseSyst))
                        CMSName = self.systConvention[baseSyst]
                        systName = systName.replace(baseSyst,CMSName)
                        if '{era}' in systName:
                            systName = systName.format(era=self.era)
                        save_name = '{}_{}'.format(group,systName)
                        
                    hist.SetTitle(save_name)
                    hist.SetName(save_name)
                    hist.Write(save_name)
            f.Write()
            f.Close()
            print ("Saved file %s"%filename)


    def applyDYNonClosure(self):
        print ('Applying DY non closure shift')
        from DYNonClosure import NonClosureSystematic
        dyObj = NonClosureSystematic.load_from_json(self.DYnonclosure['path'])
        for histName in self.content.keys():
            if histName not in self.DYnonclosure['hist']:
                continue
            catName = '{}_{}'.format(histName,self.era)
            dyObj.SetCategory(catName)
            group = self.DYnonclosure['group']

            # Add non closure systematics #
            if self.use_syst:
                h_up, h_down = dyObj.GetNonClosureShape(self.content[histName][group]['nominal'])
                self.content[histName][group]['dy_nonclosureUp'] = h_up
                self.content[histName][group]['dy_nonclosureDown'] = h_down

            # Scale all DY histograms based on non closure #
            for systName in self.content[histName][group].keys():
                if 'dy_nonclosure' in systName:
                    continue # avoid applying twice
                dyObj(self.content[histName][group][systName])
                

    def applyRebinning(self):
        print ('Applying rebinning')
        for histName, histContent in self.content.items():
            if histName in self.rebin.keys():
                method = self.rebin[histName]['method']
                params = self.rebin[histName]['params']
                print ('... {:20s} -> rebinning scheme : {}'.format(histName,self.rebin[histName]['method']))
                if method == 'classic':
                    rebinFunc = self.rebinClassic
                elif method == 'quantile':
                    rebinFunc = self.rebinInQuantile
                elif method == 'threshold':
                    rebinFunc = self.rebinThreshold
                else:
                    raise RuntimeError("Could not understand rebinning method {} for histogram {}".format(method,histName))
                rebinFunc(histName,params)
            else:
                print ('... {:20s} -> rebinning scheme not requested'.format(histName))
        
    def rebinClassic(self,histName,params):
        """
            Using the classic ROOT rebinning
            histName: name of histogram in keys of self.content
            params : (int) number of rebinning groups
        """
        assert isinstance(params,list)
        assert len(params) == 1
        assert isinstance(params[0],int)
        for group,histDict in self.content[histName].items():
            for systName,hist in histDict.items():
                self.content[histName][group][systName].Rebin(params[0])

    def rebinInQuantile(self,histName,params):
        """
            Using the quantile rebinning
            histName: name of histogram in keys of self.content
            params : (list) [quantiles (list), groups (list)]
        """
        from Rebinning import Quantile
        assert isinstance(params,list)
        assert len(params) == 2
        assert isinstance(params[0],list)
        assert isinstance(params[1],list)
        quantiles = params[0]
        groups_for_binning = params[1]
        if not set(groups_for_binning).issubset(set(self.groups.keys())):
            raise RuntimeError('Groups {'+','.join([g for g in groups_for_binning if g not in self.groups.keys()])+'} for quantile are not in group dict')

        hists_for_binning = [histDict['nominal'] for group,histDict in self.content[histName].items() if group in groups_for_binning if histDict['nominal'] is not None]
        qObj = Quantile(hists_for_binning,quantiles)
        for group in self.content[histName].keys():
            for systName,hist in self.content[histName][group].items():
                self.content[histName][group][systName] = qObj(hist) 

    def rebinThreshold(self,histName,params):
        """
            Using the threshold rebinning
            histName: name of histogram in keys of self.content
            params : (list) [thresholds (list), backgrounds (list), extra contributions (list)]
        """
        from Rebinning import Threshold
        assert isinstance(params,list)
        assert len(params) == 3
        assert isinstance(params[0],list)
        assert isinstance(params[1],list)
        assert isinstance(params[2],list)
        thresholds  = params[0]
        backgrounds = params[1]
        extra       = params[2]
        if not set(backgrounds).issubset(set(self.groups.keys())):
            raise RuntimeError('Groups {'+','.join([g for g in backgrounds if g not in self.groups.keys()])+'} for threshold are not in group dict')
        if not set(extra).issubset(set(self.groups.keys())):
            raise RuntimeError('Groups {'+','.join([g for g in extra if g not in self.groups.keys()])+'} for extra contributions are not in group dict')
        hists_for_binning = [histDict['nominal'] for group,histDict in self.content[histName].items() if group in backgrounds]
        hists_extra       = [histDict['nominal'] for group,histDict in self.content[histName].items() if group in extra]
        tObj = Threshold(thresholds,hists_for_binning,True,hists_extra)
        for group in self.content[histName].keys():
            for systName,hist in self.content[histName][group].items():
                self.content[histName][group][systName] = tObj(hist) 
                
    def preparePlotIt(self,suffix=''):
        print ("Preparing plotIt root files")
        # Make new directory with root files for plotIt #
        path_plotIt = os.path.join(os.path.abspath(os.path.dirname(__file__)),self.datacardName,'plotit')
        path_rootfiles = os.path.join(path_plotIt,'root')
        systematics = []
        if not os.path.exists(path_rootfiles):
            os.makedirs(path_rootfiles)
        for group in self.groups.keys():    
            if self.pseudodata and group == 'data_real':
                continue
            rootfile = ROOT.TFile(os.path.join(path_rootfiles,group+'.root'),'recreate')
            for histName in self.content.keys():
                for systName,hist in self.content[histName][group].items():
                    if hist is None:
                        continue
                    if systName == "nominal":
                        hist.Write(histName)
                    else:
                        baseSyst = systName.replace('Up','').replace('Down','')
                        if baseSyst not in systematics:
                            systematics.append(baseSyst)
                        hist.Write(histName+'__'+systName.replace('Up','up').replace('Down','down'))
            rootfile.Close()

        # Create yaml file #
        lumi = self.yaml_dict["luminosity"][self.era]
        config = {}
        config['configuration'] = {'eras'                     : [self.era],
                                   'experiment'               : 'CMS',
                                   'extra-label'              : '%s Datacard'%self.era,
                                   'luminosity-label'         : '%1$.2f fb^{-1} (13 TeV)',
                                   'luminosity'               : {self.era:lumi},
                                   'luminosity-error'         : 0.025,
                                   'blinded-range-fill-style' : 4050,
                                   'blinded-range-fill-color' : "#FDFBFB",
                                   'margin-bottom'            : 0.13,
                                   'margin-left'              : 0.15,
                                   'margin-right'             : 0.03,
                                   'margin-top'               : 0.05,
                                   'height'                   : 599,
                                   'width'                    : 800,
                                   'root'                     : 'root',
                                   'show-overflow'            : 'true'}
        
        config['files'] = {}
        for group,gconfig in self.groups.items():
            if self.pseudodata and group == 'data_real':
                continue
            config['files'][group+'.root'] = {'cross-section'   : 1./lumi,
                                              'era'             : str(self.era),
                                              'generated-events': 1.,
                                              'type'            : gconfig['type']}
            if gconfig['type'] != 'signal':
                config['files'][group+'.root'].update({'group':group})
            else:
                config['files'][group+'.root'].update({k:v for k,v in gconfig.items() if k not in ['files','type']})
        config['groups'] = {}
        for group,gconfig in self.groups.items():
            if self.pseudodata and group == 'data_real':
                continue
            if gconfig['type'] != 'signal':
                config['groups'][group] = {k:v for k,v in gconfig.items() if k not in ['files','type']}


        config['plots'] = {}
        for h1,h2 in self.hist_conv.items():
            if isinstance(h2,list):
                config['plots'][h1] = self.yaml_dict['plots'][h2[0]]
            else:
                config['plots'][h1] = self.yaml_dict['plots'][h2]
            config['plots'][h1]['legend-columns'] = 2
            if 'VBF' in h1 or 'GGF' in h1:
                config['plots'][h1]['blinded-range'] = [0.5,1.]
            else:
                config['plots'][h1].pop('blinded-range',None)

        if len(systematics) > 0:
            config['systematics'] = systematics

        # Write yaml file #
        path_yaml = os.path.join(path_plotIt,'plots.yml')
        with open(path_yaml,'w') as handle:
            yaml.dump(config,handle)
        print ("New yaml file for plotIt : %s"%path_yaml)

        # PlotIt command #
        path_pdf = os.path.join(path_plotIt,'plots_%s'%self.era)
        if not os.path.exists(path_pdf):
            os.makedirs(path_pdf)
        print ("plotIt command :")
        cmd = "plotIt -i {input} -o {output} -e {era} {yaml}".format(**{'input': path_plotIt,
                                                                        'output': path_pdf,   
                                                                        'era': self.era,
                                                                        'yaml': path_yaml})
        print (cmd)
        print ("Calling plotIt")
        rc = self.run_command(cmd.split(" "),print_output=True)
        if rc == 0:
            return ('... success')
        else:
            return ('... failure')

    @staticmethod
    def run_command(command,return_output=False,print_output=False):
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
        outputs = []
        while True:
            out = process.stdout.readline()
            if len(out) == 0 and process.poll() is not None:
                break
            if out and print_output:
                print (out.strip())
            if return_output:
                outputs.append(out)
        rc = process.poll()
        if return_output:
            return rc,outputs
        else:
            return rc
 


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
        f = yaml.load(handle,Loader=yaml.FullLoader)
    instance = DataCard(datacardName=args.yaml.replace('.yml',''),pseudodata=args.pseudodata,**f)
