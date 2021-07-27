import os
import io
import sys
import re
import json
import glob
import copy
import shlex
import yaml
import random
import shutil
import string
import argparse
import logging
import subprocess
import importlib
import numpy as np
import math
import multiprocessing as mp
import numpy as np
from functools import partial
import enlighten
import ROOT

from yamlLoader import YMLIncludeLoader

from IPython import embed

#ROOT.gROOT.SetBatch(True) 

INFERENCE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             'inference')
SETUP_PATH = os.path.join(INFERENCE_PATH,'setup.sh')

COMBINE_ENV = os.path.join(INFERENCE_PATH,
                          'data',
                          'software',
                          'HiggsAnalysis',
                          'CombinedLimit',
                          'env_standalone.sh')

class DataCard:
    def __init__(self,outputDir=None,configPath=None,path=None,yamlName=None,worker=None,groups=None,shapeSyst=None,normSyst=None,hist_conv=None,era=None,use_syst=False,root_subdir=None,histCorrections=None,correctionBeforeAggregation=None,pseudodata=False,rebin=None,regroup=None,histEdit=None,textfiles=None,legend=None,combineConfigs=None,logName=None,**kwargs):
        self.outputDir      = outputDir
        self.configPath     = configPath
        self.path           = path
        self.worker         = worker
        self.era            = str(era)
        self.use_syst       = use_syst
        self.root_subdir    = root_subdir
        self.pseudodata     = pseudodata
        self.hist_conv      = hist_conv
        self.rebin          = rebin
        self.histEdit       = histEdit
        self.textfiles      = textfiles
        self.legend         = legend
        self.combineConfigs = combineConfigs
        self.histCorrections = histCorrections

        # Get entries that can either be dicts or list of dicts (eg with !include in yaml) #
        self.groups = self.includeEntry(groups,'Groups')
        self.normSyst = self.includeEntry(normSyst,'Norm syst')
        self.shapeSyst = self.includeEntry(shapeSyst,'Shape syst')
        self.regroup = self.includeEntry(regroup,'Regrouping')
               
        # Load Yaml configs from bamboo #
        self.yaml_dict = self.loadYaml(self.path,yamlName)

        # Apply pseudo-data #
        if self.pseudodata:
            logging.info('Will use pseudodata')
            self.groups = self.generatePseudoData(self.groups)
            if self.outputDir is not None:
                self.outputDir += "_pseudodata"

        # Create output directory #
        self.datacardPath = os.path.join(os.path.abspath(os.path.dirname(__file__)),self.outputDir)
        if not os.path.exists(self.datacardPath):
            os.makedirs(self.datacardPath)

        # Add logging output to log file #
        handler = logging.FileHandler(os.path.join(self.datacardPath,logName),mode='w')
        logger = logging.getLogger()
        logger.addHandler(handler)

        # Initialise containers #
        self.content = {histName:{g:{} for g in self.groups.keys()} for histName in self.hist_conv.keys()}
        self.systPresent = {histName:{group:[] for group in self.groups.keys()} for histName in self.hist_conv.keys()}

    def run_production(self):
        logging.info('Running production over following categories')
        for cat in self.hist_conv.keys():
            logging.info(f'... {cat}')

        self.loopOverFiles()
        if self.pseudodata:
            self.roundFakeData()
        if self.histCorrections is not None:
            self.applyCorrectionAfterAggregation()
        if self.regroup is not None:
            self.applyRegrouping()
        if self.histEdit is not None:
            self.applyEditing()
        if self.rebin is not None:
            self.applyRebinning()
        if self.outputDir is not None:
            self.saveDatacard()

    @staticmethod
    def includeEntry(entries,name):
        if isinstance(entries,dict):
            return entries
        elif isinstance(entries,list) or isinstance(entries,tuple):
            combined = {}
            for ientries,subentries in enumerate(entries):
                if not isinstance(subentries,dict):
                    raise RuntimeError(f'Subentry index {ientries} from {name} is not a dict')
                combined.update(subentries)
            return combined
        elif entries is None:
            return None
        else: 
            raise RuntimeError(f"{name} format {type(entries)} not understood")

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
                logging.warning("Could not find sample %s in group list"%sample)
                continue
            hist_dict = self.getHistograms(f)
            for histName in hist_dict.keys():
                if len(hist_dict[histName]) == 0:
                    continue
                for group in groups:
                    self.systPresent[histName][group].append({'systematics': [systname for systname in hist_dict[histName].keys() if systname != 'nominal'],
                                                              'nominal': copy.deepcopy(hist_dict[histName]['nominal'])})
            for group in groups:
                # if several groups, need to fetch for each group 
                # -> need to have different memory adress
                self.addSampleToGroup(copy.deepcopy(hist_dict),group)
            del hist_dict

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
             
    @staticmethod
    def loadYaml(paths,yamlNames):
        if not isinstance(paths,list):
            paths = [paths]
        if not isinstance(yamlNames,list):
            yamlNames = [yamlNames]
        yamlDict = {}
        for path in paths:
            for yamlName in yamlNames:
                yamlPath = os.path.join(path,yamlName)
                if not os.path.exists(yamlPath):
                    logging.warning("{} -> not found, skipped".format(yamlPath))
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
        from txtwriter import Writer
        if self.outputDir is None:
            raise RuntimeError("Datacard output path is not set")

        shapes = {histName:os.path.join(self.datacardPath,f"{histName}_{self.era}.root") for histName in self.content.keys()}

        # Save root file #
        for histName in self.content.keys():
            f = ROOT.TFile(shapes[histName],'recreate')
            if self.root_subdir is not None:
                d = f.mkdir(self.root_subdir,self.root_subdir)
                d.cd()
            for group in self.content[histName].keys():
                # Loop over systematics first to fix and save #
                for systName in self.content[histName][group].keys():
                    if self.pseudodata and group == 'data_real': # avoid writing data when using pseudodata
                        continue
                    if systName ==  "nominal": # will do afterwards
                        continue
                    if systName.endswith('Down'): # wait for the Up and do both at the same time
                        continue

                    if systName.endswith('Up'): 
                        # Get correct name #
                        systName = systName[:-2]
                        if systName not in self.shapeSyst.keys():
                            logging.warning("Could not find {} in systematic dict".format(systName))
                            continue
                        CMSName = self.shapeSyst[systName]
                        if '{era}' in CMSName:
                            CMSName = CMSName.format(era=self.era)

                        if CMSName == 'discard':
                            continue
        
                        systNameUp   = systName + "Up"
                        systNameDown = systName + "Down"
                         
                        if systNameUp not in self.content[histName][group].keys():
                            raise RuntimeError(f"Could not find syst named {systNameUp} in group {group} for histogram {histName}")
                        if systNameDown not in self.content[histName][group].keys():
                            raise RuntimeError(f"Could not find syst named {systNameDown} in group {group} for histogram {histName}")

                        # Fix shape in case needed #
                        self.fixHistograms(hnom  = self.content[histName][group]['nominal'],
                                            hup   = self.content[histName][group][systNameUp],
                                            hdown = self.content[histName][group][systNameDown])

                        # Write to file #
                        CMSNameUp   = f"{group}__{CMSName}Up"
                        CMSNameDown = f"{group}__{CMSName}Down"

                        self.content[histName][group][systNameUp].SetTitle(CMSNameUp)
                        self.content[histName][group][systNameUp].SetName(CMSNameUp)
                        self.content[histName][group][systNameDown].SetTitle(CMSNameDown)
                        self.content[histName][group][systNameDown].SetName(CMSNameDown)
        
                        self.content[histName][group][systNameUp].Write(CMSNameUp) 
                        self.content[histName][group][systNameDown].Write(CMSNameDown) 
                    else:
                        raise RuntimeError(f"Something wrong happened with {systName} in group {group} for histogram {histName}")

                # Save nominal (done after because can be correct when fixing systematics) #
                if 'nominal' not in self.content[histName][group].keys():
                    raise RuntimeError(f"Group {group} nominal histogram {histName} was not found")
                self.fixHistograms(hnom=self.content[histName][group]['nominal'])

                self.content[histName][group]['nominal'].SetTitle(group)
                self.content[histName][group]['nominal'].SetName(group)
                self.content[histName][group]['nominal'].Write(group)
                        
            f.Write()
            f.Close()
            logging.info(f"Saved file {shapes[histName]}")

        # Save txt file #
        if self.textfiles is None:
            self.textfiles = '{}.txt'
        if '{}' not in self.textfiles:
            logging.info('Will create a single datacard txt file for all the categories')
            writer = Writer([f'{histName}_{self.era}' for histName in self.content.keys()])
        else:
            logging.info('Will create one datacard txt file per category')
        for histName in self.content.keys():
            binName = f'{histName}_{self.era}'
            if '{}' in self.textfiles:
                writer = Writer(binName)
            # Add processes #
            for group in self.content[histName]:
                writer.addProcess(binName       = binName,
                                  processName   = group,
                                  rate          = self.content[histName][group]['nominal'].Integral(),
                                  processType   = self.groups[group]['type'])
                # Add shape systematics #
                if self.use_syst:
                    for systName in self.content[histName][group]:
                        if systName == 'nominal':
                            continue
                        if systName.endswith('Up'):
                            systName = systName[:-2]
                            # Check integral of systematics
                            normUp = self.content[histName][group][f'{systName}Up'].Integral()
                            normDown = self.content[histName][group][f'{systName}Down'].Integral()
                            if normUp <= 0. or normDown <=0.:
                                continue

                            # Get convention naming #
                            if systName not in self.shapeSyst.keys():
                                raise RuntimeError(f"Could not find {systName} in systematic dict")
                            CMSName = self.shapeSyst[systName]
                            if CMSName == 'discard':
                                logging.debug(f'Discarding systematic {systName}')
                                continue
                            if '{era}' in CMSName:
                                CMSName = CMSName.format(era=self.era)
                            # Record #
                            writer.addShapeSystematic(binName,group,CMSName)
                # Add norm systematics #
                if self.use_syst:
                    for systName,systList in self.normSyst.items():
                        if not isinstance(systList,list):
                            systList = [systList]
                        if systName == 'footer':
                            continue
                        for systContent in systList:
                            if systContent is None:
                                raise RuntimeError(f"Problem with lnL systematic {systName}")
                            processes = list(self.groups.keys())
                            if 'era' in systContent.keys():
                                if str(systContent['era']) != str(self.era):
                                    continue
                            if 'proc' in systContent.keys():
                                if not isinstance(systContent['proc'],list):
                                    selectProcess = [systContent['proc']]
                                else:
                                    selectProcess = systContent['proc']
                                if all([proc.startswith("\\") for proc in selectProcess]):
                                    for proc in selectProcess:
                                        proc = proc.replace('\\','')
                                        if proc in processes:
                                            processes.remove(proc)
                                else:
                                    processes = selectProcess
                            if 'hist' in systContent.keys():
                                if histName != systContent['hist']:
                                    continue
                            for process in processes:
                                writer.addLnNSystematic(binName,process,systName,systContent['val'])
            # Add potential footer #
            if self.use_syst and "footer" in self.normSyst.keys():
                for footer in self.normSyst['footer']:
                    writer.addFooter(footer)
                
            if '{}' in self.textfiles:
                textPath = os.path.join(self.datacardPath,self.textfiles.format(f"{histName}_{self.era}"))
                writer.dump(textPath,os.path.basename(shapes[histName]))
                logging.info(f"Saved file {textPath}")
        if '{}' not in self.textfiles:
            textPath = os.path.join(self.datacardPath,self.textfiles)
            writer.dump(textPath,[os.path.basename(shape) for shape in shapes.values()])
            logging.info(f"Saved file {textPath}")

    @staticmethod
    def fixHistograms(hnom,hdown=None,hup=None):
        val     = 1e-5 * min(1.,hnom.Integral()) 
        valup   = 1e-5 * min(1.,hup.Integral()) if hdown is not None else None
        valdown = 1e-5 * min(1.,hdown.Integral()) if hdown is not None else None
        # Loop over bins #
        for i in range(1,hnom.GetNbinsX()+1):
            # First clip to 0 all negative bin content #
            if hnom.GetBinContent(i) <= 0.:
                hnom.SetBinContent(i,val)
                hnom.SetBinError(i,val) 
            if hup is not None and hup.GetBinContent(i) <= 0.:
                hup.SetBinContent(i,valup)
                hup.SetBinError(i,valup)
            if hdown is not None and hdown.GetBinContent(i) <= 0.:
                hdown.SetBinContent(i,valdown)
                hdown.SetBinError(i,valdown)

            # Second, check the up and down compared to nominal #
            if hnom.GetBinContent(i) > 0 and hdown is not None and hup is not None:
                # Nominal bin not zero #
                # Check if zero bin -> apply nominal**2 / up or down
                if hdown.GetBinContent(i) == 0. and abs(hup.GetBinContent(i)) > 0:
                    hdown.SetBinContent(i, hnom.GetBinContent(i)**2 / hup.GetBinContent(i)) 
                if hup.GetBinContent(i) == 0. and abs(hdown.GetBinContent(i)) > 0:
                    hup.SetBinContent(i, hnom.GetBinContent(i)**2 / hdown.GetBinContent(i))
                # Check if too big, deflate in case #
                if hdown.GetBinContent(i)/hnom.GetBinContent(i) > 100:
                    hdown.SetBinContent(i, 100 * hnom.GetBinContent(i))
                if hup.GetBinContent(i)/hnom.GetBinContent(i) > 100:
                    hup.SetBinContent(i, 100 * hnom.GetBinContent(i))
                # Check if too small, inlate in case #
                if hdown.GetBinContent(i)/hnom.GetBinContent(i) < 1./100:
                    hdown.SetBinContent(i, 1./100 * hnom.GetBinContent(i))
                if hup.GetBinContent(i)/hnom.GetBinContent(i) < 1./100:
                    hup.SetBinContent(i, 1./100 * hnom.GetBinContent(i))
            else:
                # Nominal is == 0 #
                if hup is not None and hdown is not None:
                    if abs(hup.GetBinContent(i)) > 0 or abs(hdown.GetBinContent(i)) > 0:
                        # zero nominal but non zero systematics -> set all at 0.00001 #
                        hnom.SetBinContent(i,0.00001)
                        hup.SetBinContent(i,0.00001)
                        hdown.SetBinContent(i,0.00001)

    def applyCorrectionAfterAggregation(self):
        for corrConfig in self.histCorrections:
            if ":" not in corrConfig['module']:
                raise RuntimeError(f"`:` needs to be in the module arg {corrConfig['module']}")
            lib, clsName = corrConfig['module'].split(':')
            if not os.path.isfile(lib):
                raise RuntimeError(f'File {lib} does not exist')
            mod = importlib.import_module(lib.replace('.py',''))
            cls = getattr(mod, clsName, None)(**corrConfig['init'])
            logging.info(f"Applying {corrConfig['module']}")
            for cat in self.content.keys():
                if cat not in corrConfig['categories'].keys():
                    continue
                catCfg = corrConfig['categories'][cat]
                logging.info(f'... Histogram {cat:20s} in {cls.group} group')
                # Correct nominal and all syst hists #
                if hasattr(cls,'modify'):
                    logging.info(f'\t\t-> applying corrections')
                    for systName,hist in self.content[cat][cls.group].items():
                        cls.modify(hist,**catCfg)
                # Add potential additional shapes #
                if hasattr(cls,'additional') and self.use_syst:
                    add = cls.additional(self.content[cat][cls.group]['nominal'],**catCfg)
                    for key in add.keys():
                        logging.info(f"\t\t-> Adding systematic shape {catCfg['systName']}")
                    self.content[cat][cls.group].update(add)

    def applyRegrouping(self):
        logging.info('Applying Regrouping of categories')
        self.checkForMissingNominals()

        # Initialize and checks #
        regroup_conv = {v:key for key,val in self.regroup.items() for v in val}
        regrouped = {regroup_conv[histName]:{g:{} for g in self.groups.keys()} for histName in self.content.keys()}
        innerCat = [v for val in self.regroup.values() for v in val]
        if len(set(self.content.keys())-set(innerCat)) > 0:
            raise RuntimeError("Regrouping has missing categories : "+','.join([cat for cat in self.content.keys() if cat not in innerCat]))
        if len(set(self.content.keys())-set(innerCat)) < 0:
            raise RuntimeError("Regrouping has excess categories : "+','.join([cat for cat in innerCat if cat not in self.content.keys()]))
        #for outCat,inCats in self.regroup.items():
        for outCat in regrouped.keys():
            inCats = self.regroup[outCat]
            logging.info(f"New category {outCat} will aggregate :")
            for inCat in inCats:
                logging.info(f'... {inCat}')
        # Neeed to compensate missing systematics in categories to aggregate 
        # Eg if someone wants to merge electron and muon channels, there will 
        # be different systematics and we must use the nominal hist when one is missing
        list_syst = {}
        for group in self.groups.keys():
            for outCat in regrouped.keys():
                inCats = self.regroup[outCat]
                list_syst = list(set([systName for cat in inCats for systName in self.content[cat][group].keys()]))
                for cat in inCats:
                    for systName in list_syst:
                        if systName != 'nominal' and systName not in self.content[cat][group].keys():
                            self.content[cat][group][systName] = copy.deepcopy(self.content[cat][group]['nominal'])
                            self.content[cat][group][systName].SetName(self.content[cat][group][systName].GetName()+f'_{systName}')
    
        # Aggregate and save as new content #
        for cat in self.content.keys():
            recat = regroup_conv[cat] 
            for group in self.content[cat].keys():
                for systName,h in self.content[cat][group].items():
                    if systName not in regrouped[recat][group].keys():
                        regrouped[recat][group][systName] = self.content[cat][group][systName]
                    else:
                        regrouped[recat][group][systName].Add(self.content[cat][group][systName])
        self.content = copy.deepcopy(regrouped)
        

    def checkForMissingNominals(self):
        # Check at least one nominal per group #
        missings = []
        for histName in self.content.keys():
            for group in self.content[histName]:
                if 'nominal' not in self.content[histName][group].keys():
                    missings.append([histName,group])
        if len(missings) != 0 :
            logging.error('Following histograms are missing :')
            for histName,group in missings:
                logging.error(f'... {group:20s} -> {histName}')
            raise RuntimeError('Some nominal files are missing')

                
    def applyEditing(self):
        logging.info('Applying editing')
        self.checkForMissingNominals()

        # Apply editing schemes #
        for histName in self.content.keys():
            # Check name and type #
            if histName not in self.histEdit.keys():
                continue
            histType = self.content[histName][list(self.groups.keys())[0]]['nominal'].__class__.__name__ 
            if histType.startswith('TH1'):
                valid_args = ['x','i','val']
            elif histType.startswith('TH2'):
                valid_args = ['x','y','i','j','val']
            else:   
                raise RuntimeError(f'Histogram {histName} type not understood : {histType}')

            logging.info(f'Editing histogram {histName}')
            
            # Get the functions with partial #
            editFuncs = []
            for iedit,editCfg in enumerate(self.histEdit[histName]):
                if 'val' not in editCfg.keys():
                    logging.warning(f'Histogram editing for {histName} is missing the `val` key, will assume 0.')
                    editCfg['val'] = 0.
                if any([key not in valid_args for key in editCfg.keys()]):
                    raise RuntimeError('Keys not understood for config {iedit} in histogram {histName} : '+\
                                    ','.join([key not in valid_args for key in editCfg.keys()]))
                if histType.startswith('TH1'):
                    editFuncs.append(partial(self.editTH1,**editCfg))
                if histType.startswith('TH2'):
                    editFuncs.append(partial(self.editTH2,**editCfg))

            # Apply editings #
            for group in self.content[histName].keys():
                for systName,hist in self.content[histName][group].items():
                    if systName == 'nominal':
                        oldInt = hist.Integral()
                    for editFunc in editFuncs:
                        editFunc(hist)
                    if systName == 'nominal':
                        newInt = hist.Integral()
                logging.info(f'... group {group:20s} : Integral = {oldInt:9.3f} -> {newInt:9.3f} [{(oldInt-newInt)/(oldInt+1e-9)*100:5.2f}%]')

    def editTH1(self,h,val,x=None,i=None):
        assert (x is not None and i is None) or (x is None and i is not None)
        if x is not None:
            bins = [h.FindBin(x)]
        if i is not None:
            bins = [i]
        self.editHistByBinNumber(h,bins,val)

    def editTH2(self,h,val,x=None,y=None,i=None,j=None):
        assert ((x is not None or y is not None) and (i is None and j is None)) or \
               ((x is None and y is None) and (i is not None or j is not None))
        yAxis = h.GetYaxis()
        xAxis = h.GetXaxis()

        if x is not None or y is not None:
            if x is not None and y is None:
                initBin = h.FindBin(x,yAxis.GetBinCenter(1))
                bins = [initBin + k*(h.GetNbinsX()+2) for k in range(h.GetNbinsX()-1)]
            elif x is None and y is not None:
                initBin = h.FindBin(xAxis.GetBinCenter(1),y) 
                bins = [initBin + k for k in range(h.GetNbinsX())]
            else:
                bins = [h.FindBin(x,y)]
        if i is not None or j is not None:
            if i is not None and j is None:
                initBin = h.GetBin(i,1)
                bins = [initBin + k*h.GetNbinsX()+1 for k in range(h.GetNbinsX()-1)]
            elif i is None and j is not None:
                initBin = h.GetBin(1,j)
                bins = [initBin + k for k in range(h.GetNbinsX())]
            else:
                bins = [h.GetBin(i,j)]
        self.editHistByBinNumber(h,bins,val)

    @staticmethod
    def editHistByBinNumber(h,bins,val):
        assert isinstance(bins,list)
        for b in bins:
            h.SetBinContent(b,val)

    def applyRebinning(self):
        logging.info('Applying rebinning')
        self.checkForMissingNominals()

        # Apply rebinning schemes #
        for histName in self.content.keys():
            histType = self.content[histName][list(self.groups.keys())[0]]['nominal'].__class__.__name__ 
            logging.info(f'... {histName:20s} ({histType})')
            if histName in self.rebin.keys():
                rebinSchemes = self.rebin[histName]
                if not isinstance(rebinSchemes,list):
                    rebinSchemes = [rebinSchemes]
                for rebinScheme in rebinSchemes:
                    method = rebinScheme['method']
                    params = rebinScheme['params']
                    # 1D rebinnings #
                    if method == 'classic':
                        rebinFunc = self.rebinClassic
                    elif method == 'boundary':
                        rebinFunc = self.rebinBoundary
                    elif method == 'quantile':
                        rebinFunc = self.rebinInQuantile
                    elif method == 'threshold':
                        rebinFunc = self.rebinThreshold
                    # 2D rebinnings #
                    elif method == 'classic2d':
                        rebinFunc = self.rebinClassic2D
                    elif method == 'boundary2d':
                        rebinFunc = self.rebinBoundary2D
                    elif method == 'quantile2d':
                        rebinFunc = self.rebinInQuantile2D
                    elif method == 'threshold2d':
                        rebinFunc = self.rebinThreshold2D
                    # 2D Linearized #
                    elif method == 'linearize2d':
                        rebinFunc = self.rebinLinearize2D
                    elif method == 'rebinlinearize2d':
                        rebinFunc = self.rebinLinearizeSplit2D
                    # Error #
                    else:
                        raise RuntimeError("Could not understand rebinning method {} for histogram {}".format(method,histName))
                    logging.info(f'\t -> rebinning scheme : {method}')
                    rebinFunc(histName,params)
            else:
                logging.info('\t -> rebinning scheme not requested')

    def rebinGetHists(self,histName,groups_for_binning):
        if not set(groups_for_binning).issubset(set(self.groups.keys())):
            raise RuntimeError('Groups {'+','.join([g for g in groups_for_binning if g not in self.groups.keys()])+'} are not in group dict')
        hists_for_binning = [histDict['nominal'] for group,histDict in self.content[histName].items() if group in groups_for_binning if histDict['nominal'] is not None]
        return hists_for_binning

    def countPerHistName(self,histName):
        N = 0
        for group in self.content[histName].keys():
            N += len(self.content[histName][group])
        return N
        
    def rebinClassic(self,histName,params):
        """
            Using the classic ROOT rebinning (1D)
            histName: name of histogram in keys of self.content
            params : (list) [groups(int)]
        """
        assert isinstance(params,list)
        assert len(params) == 1
        assert isinstance(params[0],int)
        pbar = enlighten.Counter(total=self.countPerHistName(histName), desc='Progress', unit='histograms')
        for group,histDict in self.content[histName].items():
            for systName,hist in histDict.items():
                self.content[histName][group][systName].Rebin(params[0])
                pbar.update()

    def rebinClassic2D(self,histName,params):
        """
            Using the classic ROOT rebinning (2D)
            histName: name of histogram in keys of self.content
            params : (list) [groups(int),groups(int)]
        """
        assert isinstance(params,list)
        assert len(params) == 2
        assert isinstance(params[0],int)
        assert isinstance(params[1],int)
        pbar = enlighten.Counter(total=self.countPerHistName(histName), desc='Progress', unit='histograms')
        for group,histDict in self.content[histName].items():
            for systName,hist in histDict.items():
                self.content[histName][group][systName].Rebin(params[0],params[1])
                pbar.update()


    def rebinBoundary(self,histName,params):
        """
            Using the given hardcoded boundaries (1D)
            histName: name of histogram in keys of self.content
            params : (list) [boundaries(list)]
        """
        from Rebinning import Boundary
        assert isinstance(params,list)
        assert len(params) == 1
        assert isinstance(params[0],list)
        boundObj = Boundary(boundaries=params[0])
        pbar = enlighten.Counter(total=self.countPerHistName(histName), desc='Progress', unit='histograms')
        for group in self.content[histName].keys():
            for systName,hist in self.content[histName][group].items():
                self.content[histName][group][systName] = boundObj(hist) 
                pbar.update()

    def rebinBoundary2D(self,histName,params):
        """
            Using the given hardcoded boundaries (2D)
            histName: name of histogram in keys of self.content
            params : (list) [boundaries(list),boundaries(list)]
        """
        from Rebinning import Boundary2D
        assert isinstance(params,list)
        assert len(params) == 2
        assert isinstance(params[0],list)
        assert isinstance(params[1],list)
        boundObj = Boundary2D(params[0],params[1])
        pbar = enlighten.Counter(total=self.countPerHistName(histName), desc='Progress', unit='histograms')
        for group in self.content[histName].keys():
            for systName,hist in self.content[histName][group].items():
                self.content[histName][group][systName] = boundObj(hist) 
                pbar.update()


    def rebinInQuantile(self,histName,params):
        """
            Using the quantile rebinning (1D)
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
        hists_for_binning = self.rebinGetHists(histName,groups_for_binning)
        qObj = Quantile(hists_for_binning,quantiles)
        pbar = enlighten.Counter(total=self.countPerHistName(histName), desc='Progress', unit='histograms')
        for group in self.content[histName].keys():
            for systName,hist in self.content[histName][group].items():
                self.content[histName][group][systName] = qObj(hist) 
                pbar.update()

    def rebinInQuantile2D(self,histName,params):
        """
            Using the quantile rebinning (2D)
            histName: name of histogram in keys of self.content
            params : (list) [quantiles (list), quantiles (list), groups (list)]
        """
        from Rebinning import Quantile2D
        assert isinstance(params,list)
        assert len(params) == 3
        assert isinstance(params[0],list)
        assert isinstance(params[1],list)
        assert isinstance(params[2],list)
        qx = params[0]
        qy = params[1]
        groups_for_binning = params[2]
        hists_for_binning = self.rebinGetHists(histName,groups_for_binning)
        qObj = Quantile2D(hists_for_binning,qx,qy)
        pbar = enlighten.Counter(total=self.countPerHistName(histName), desc='Progress', unit='histograms')
        for group in self.content[histName].keys():
            for systName,hist in self.content[histName][group].items():
                self.content[histName][group][systName] = qObj(hist) 
                pbar.update()


    def rebinThreshold(self,histName,params):
        """
            Using the threshold rebinning (1D)
            histName: name of histogram in keys of self.content
            params : (list) [thresholds (list), extra contributions (list | str), number of bins (int), rsut (float), groups (list | str)]
        """
        from Rebinning import Threshold
        assert isinstance(params,list)
        assert len(params) == 5
        assert isinstance(params[0],list)
        assert isinstance(params[1],list) or isinstance(params[1],str)
        assert isinstance(params[2],int)
        assert isinstance(params[3],float)
        assert isinstance(params[4],list) or isinstance(params[1],str)
        thresholds          = params[0]
        extra               = params[1] # main processes to be kept above threshold
        nbins               = params[2] 
        rsut                = params[3] # relative stat. unc. threshold
        groups_for_binning  = params[4]
        hists_for_binning   = self.rebinGetHists(histName,groups_for_binning)
        hists_extra         = self.rebinGetHists(histName,extra)
        tObj = Threshold(hists_for_binning,thresholds,hists_extra,nbins,rsut)
        pbar = enlighten.Counter(total=self.countPerHistName(histName), desc='Progress', unit='histograms')
        for group in self.content[histName].keys():
            for systName,hist in self.content[histName][group].items():
                self.content[histName][group][systName] = tObj(hist) 
                pbar.update()
                
    def rebinThreshold2D(self,histName,params):
        """
            Using the threshold rebinning (2D)
            histName: name of histogram in keys of self.content
            params : (list) [thresholds (list), extra contributions (list), number of bins (int), rsut (float), groups (list | str)]
            params : (list) [thresholds x (list), thresholds y (list),extra contributions (list), number of bins (int) rsut (float), groups (list)]
        """
        from Rebinning import Threshold2D
        assert isinstance(params,list)
        assert len(params) == 4
        assert isinstance(params[0],list)
        assert isinstance(params[1],list)
        assert isinstance(params[2],list) or isinstance(params[2],str)
        assert isinstance(params[3],int) 
        assert isinstance(params[4],float)
        assert isinstance(params[5],list) or isinstance(params[5],str)
        threshx             = params[0]
        threshy             = params[1]
        extra               = params[2] # main processes to be kept above threshold
        nbins               = params[3] 
        rsut                = params[4] # relative stat. unc. threshold
        groups_for_binning  = params[5]
        hists_for_binning   = self.rebinGetHists(histName,groups_for_binning)
        hists_extra         = self.rebinGetHists(histName,extra)
        tObj = Threshold(hists_for_binning,threshx,threshy,hists_extra,nbins,rsut)
        pbar = enlighten.Counter(total=self.countPerHistName(histName), desc='Progress', unit='histograms')
        for group in self.content[histName].keys():
            for systName,hist in self.content[histName][group].items():
                self.content[histName][group][systName] = tObj(hist) 
                pbar.update()

    def rebinLinearize2D(self,histName,params):
        """
            Using the linearization of 2D hostogram
            histName: name of histogram in keys of self.content
            params : Major axis ('x' or 'y')
        """
        from Rebinning import Linearize2D
        assert isinstance(params,list)
        assert len(params) == 1
        assert params[0] in ['x','y']
        lObj = Linearize2D(major=params[0])
        pbar = enlighten.Counter(total=self.countPerHistName(histName), desc='Progress', unit='histograms')
        for group in self.content[histName].keys():
            for systName,hist in self.content[histName][group].items():
                self.content[histName][group][systName] = lObj(hist)
                pbar.update()
        if not hasattr(self,'plotLinearizeData'):
            self.plotLinearizeData = {histName:lObj.getPlotData()}
        else:
            self.plotLinearizeData[histName] = lObj.getPlotData()
         
    def rebinLinearizeSplit2D(self,histName,params):
        """
            Using the linearization of 2D hostogram
            histName: name of histogram in keys of self.content
            params : Major axis ('x' or 'y')
        """
        from Rebinning import LinearizeSplitRebin2D
        assert isinstance(params,list)
        assert len(params) == 6
        assert params[0] in ['x','y']       # Major binning
        assert isinstance(params[1],str)    # Name of class for major ax rebinning
        assert isinstance(params[2],list)   # Params of the major ax rebinning
        assert isinstance(params[3],str)    # Name of class for minor ax rebinning
        assert isinstance(params[4],list)   # Params of the minor ax rebinning
        assert isinstance(params[5],list)   # Groups to be used in the rebinning definition
        groups_for_binning = params[5]
        hists_for_binning  = self.rebinGetHists(histName,groups_for_binning)
        # Major #
        major_class = params[1]
        if major_class == "Boundary":
            major_params = params[2]
        elif major_class == "Quantile":
            major_params = params[2]
        elif major_class == "Threshold":
            major_params = params[2]
            major_params[2][1] = self.rebinGetHists(histName,params[2][1])
        else:
            major_class = None
            major_params = None
        # Minor #
        minor_class = params[3]
        if minor_class == "Boundary":
            minor_params = params[4]
        elif minor_class == "Quantile":
            minor_params = params[4]
        elif minor_class == "Threshold":
            minor_params = params[4]
            minor_params[1] = self.rebinGetHists(histName,params[4][1])
        else:
            minor_class = None
            minor_params = None
            
        lObj = LinearizeSplitRebin2D(major        = params[0],
                                     major_class  = major_class,
                                     major_params = major_params,
                                     minor_class  = minor_class,
                                     minor_params = minor_params,
                                     h            = hists_for_binning)
        pbar = enlighten.Counter(total=self.countPerHistName(histName), desc='Progress', unit='histograms')
        for group in self.content[histName].keys():
            for systName,hist in self.content[histName][group].items():
                self.content[histName][group][systName] = lObj(hist)
                pbar.update()

        # Save linearized info for plotIt #
        if not hasattr(self,'plotLinearizeData'):
            self.plotLinearizeData = {histName:lObj.getPlotData()}
        else:
            self.plotLinearizeData[histName] = lObj.getPlotData()
         


    def prepare_plotIt(self):
        logging.info("Preparing plotIt root files for categories")
        if self.regroup is not None:
            hist_conv = {}
            regroup_conv = {v:key for key,val in self.regroup.items() for v in val}
            for cat,hists in self.hist_conv.items():
                recat = regroup_conv[cat]
                if recat not in hist_conv.keys():
                    hist_conv[recat] = hists
                else:
                    hist_conv[recat].extend(hists)
        else:
            hist_conv = self.hist_conv

        categories = hist_conv.keys()
        for cat in hist_conv.keys():
            logging.info(f"... {cat}")

        # Initialize containers #
        content = {cat:{group:{} for group in self.groups.keys()} for cat in hist_conv.keys()}
        systematics = []

        # Make new directory with root files for plotIt #
        path_plotIt = os.path.join(self.datacardPath,'plotit')
        path_rootfiles = os.path.join(path_plotIt,'root')
        if not os.path.exists(path_rootfiles):
            os.makedirs(path_rootfiles)
            
        # Loop over root files to get content #
        datacardRoots = glob.glob(os.path.join(self.datacardPath,"*root"))
        for f in glob.glob(os.path.join(self.datacardPath,"*root")):
            category = os.path.basename(f).replace('.root','').replace(f'_{self.era}','')
            logging.debug(f"Looking at {f} -> category {category}")
            if category not in content.keys():
                continue
            F = ROOT.TFile(f)
            for key in F.GetListOfKeys():
                name = key.GetName()
                if '__' in name: # systematic histogram
                    group,systName = name.split('__')
                    systName = systName.replace('Up','up').replace('Down','down')
                    baseSyst = systName.replace('up','').replace('down','')
                else:   # nominal histogram
                    group = name
                    systName = 'nominal'
                    baseSyst = None
                content[category][group][systName] = copy.deepcopy(F.Get(name))
                if baseSyst is not None and baseSyst not in systematics:
                    systematics.append(baseSyst)
            F.Close()

        valid_content = True
        for cat in hist_conv.keys():
            for group in self.groups.keys():
                if len(content[cat][group]) == 0:
                    logging.error(f"Group {group} in category {cat} is empty")
                    valid_content = False
        if not valid_content:
            raise RuntimeError('Failed to run plotIt')

                    

        # Write to file #
        for group in self.groups.keys():
            F = ROOT.TFile(os.path.join(path_rootfiles,group+'.root'),'update')
            for histName in content.keys():
                for systName,hist in content[histName][group].items():
                    if hist is None:
                        continue
                    outName = histName if systName == "nominal" else f'{histName}__{systName}'
                    hist.Write(outName,ROOT.TObject.kOverwrite)
            F.Close()

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

        # Compute min-max values for plots #
        def getMinNonEmptyBins(h):
            hmin = math.inf
            for i in range(1,h.GetNbinsX()+1):
                if h.GetBinContent(i) > 0:
                    hmin = min(hmin,h.GetBinContent(i))
            return hmin
            
        histMax = {}
        histMin = {}
        for histName in content.keys():
            histMax[histName] = -math.inf
            histMin[histName] = 0.1
            hstack = ROOT.THStack(histName+"stack",histName+"stack")
            for group in content[histName].keys():
                if self.groups[group]['type'] == 'mc':
                    hstack.Add(content[histName][group]['nominal'])
                if self.groups[group]['type'] == 'signal':
                    if 'hide' in self.groups[group].keys() and self.groups[group]['hide']:
                        continue
                    histMin[histName] = min(histMin[histName],getMinNonEmptyBins(content[histName][group]['nominal']))
                    histMax[histName] = max(histMax[histName],content[histName][group]['nominal'].GetMaximum())
            histMax[histName] = max(hstack.GetMaximum(),histMax[histName])

        # Files informations #
        config['files'] = {}
        for group,gconfig in self.groups.items():
            if 'hide' in gconfig.keys() and gconfig['hide']:
                continue
            config['files'][group+'.root'] = {'cross-section'   : 1./lumi,
                                              'era'             : str(self.era),
                                              'generated-events': 1.,
                                              'type'            : gconfig['type']}
            if gconfig['type'] != 'signal':
                plotGroup = group if 'group' not in gconfig.keys() else gconfig['group']
                config['files'][group+'.root'].update({'group':plotGroup})
            else:
                config['files'][group+'.root'].update({k:v for k,v in gconfig.items() if k not in ['files','type']})

        # Groups informations #
        config['groups'] = {}
        for group,gconfig in self.groups.items():
            if gconfig['type'] != 'signal':
                plotGroup = group if 'group' not in gconfig.keys() else gconfig['group']
                if plotGroup not in config['groups'].keys():
                    config['groups'][plotGroup] = {k:v for k,v in gconfig.items() if k not in ['files','type']}


        # Plots informations #
        config['plots'] = {}
        for h1,h2 in hist_conv.items():
            # Check histogram dimension #
            hist = content[h1][list(self.groups.keys())[0]]['nominal']
            if isinstance(hist,ROOT.TH2F) or isinstance(hist,ROOT.TH2D):
                logging.error(f'Histogram {h1} is a TH2 and can therefore not be used in plotIt')
                continue
            # Get basic config in case not found in plots.yml
            baseCfg = {'log-y'              : 'both',
                       'x-axis'             : hist.GetXaxis().GetTitle(),
                       'x-axis-range'       : [hist.GetXaxis().GetBinLowEdge(1),hist.GetXaxis().GetBinUpEdge(hist.GetNbinsX())],
                       'y-axis'             : 'Events',
                       'y-axis-show-zero'   : True,
                       'ratio-y-axis'       : '#frac{Data}{MC}',
                       'show-ratio'         : True}
                       
            # Check if plots already in yaml #
            if isinstance(h2,list):
                if h2[0] in self.yaml_dict['plots'].keys():
                    config['plots'][h1] = self.yaml_dict['plots'][h2[0]]
                else:
                    config['plots'][h1] = baseCfg
            else:
                if h2 in self.yaml_dict['plots'].keys():
                    config['plots'][h1] = self.yaml_dict['plots'][h2]
                else:
                    config['plots'][h1] = baseCfg
            if 'VBF' in h1 or 'GGF' in h1:
                xmin, xmax = config['plots'][h1]['x-axis-range']
                config['plots'][h1]['blinded-range'] = [xmin,xmax] 
            else:
                config['plots'][h1].pop('blinded-range',None)
            if histMax[h1] != 0.:
                config['plots'][h1]['y-axis-range'] = [0.,histMax[h1]*1.5]
                config['plots'][h1]['log-y-axis-range'] = [histMin[h1],histMax[h1]*100]
                config['plots'][h1]['ratio-y-axis-range'] = [0.8,1.2]
            config['plots'][h1]['sort-by-yields'] = True
            config['plots'][h1]['show-overflow'] = True
            if 'labels' in config['plots'][h1].keys():
                del config['plots'][h1]['labels']

            if self.legend is not None:
                if 'position' in self.legend.keys():
                    assert len(self.legend['position']) == 4
                    config['plots'][h1]['legend-position'] = self.legend['position']
                if 'columns' in self.legend.keys():
                    config['plots'][h1]['legend-columns'] = self.legend['columns']

#            if hasattr(self,'plotLinearizeData'):
#                if h1 in self.plotLinearizeData.keys():
#                    margin_left = config['configuration']['margin-left']
#                    margin_right = 1-config['configuration']['margin-right']
#                    extraItems = self.plotLinearizeData[h1]
#                    for idx in range(len(extraItems['labels'])):
#                        extraItems['labels'][idx]['position'][0] = margin_left + extraItems['labels'][idx]['position'][0] * (margin_right-margin_left)
#                    for idx in range(len(extraItems['lines'])):
#                        extraItems['lines'][idx] = [[extraItems['lines'][idx],0.],[extraItems['lines'][idx],histMax[h1]*1.1]]
#                    config['plots'][h1].update(extraItems)

        if len(systematics) > 0:
            config['systematics'] = systematics
        return config

    @staticmethod
    def merge_plotIt(configs):
        mainConfig = {}
        for config in configs:
            for key in ['configuration','files','groups','plots']:
                if key not in mainConfig.keys():
                    mainConfig[key] = config[key]
                else:
                    mainConfig[key].update(**config[key])
            if 'systematics' in config:
                if 'systematics' not in mainConfig.keys():
                    mainConfig['systematics'] = []
                for syst in config['systematics']:
                    if syst not in mainConfig['systematics']:
                        mainConfig['systematics'].append(syst)
        return mainConfig

    def run_plotIt(self,config):
        # Write yaml file #
        path_plotIt = os.path.join(self.datacardPath,'plotit')
        path_yaml = os.path.join(path_plotIt,'plots.yml')
        with open(path_yaml,'w') as handle:
            yaml.dump(config,handle)
        logging.info("New yaml file for plotIt : %s"%path_yaml)

        # PlotIt command #
        path_pdf = os.path.join(path_plotIt,'plots_%s'%self.era)
        if not os.path.exists(path_pdf):
            os.makedirs(path_pdf)
        logging.info("plotIt command :")
        cmd = f"plotIt -i {path_plotIt} -o {path_pdf} -e {self.era} {path_yaml}"
        print (cmd)
        if self.which('plotIt') is not None:
            logging.info("Calling plotIt")
            exitCode,output = self.run_command(shlex.split(cmd),return_output=True)
            if exitCode != 0:
                if logging.root.level > 10:
                    for line in output:
                        logging.info(line.strip())
                logging.info('... failure (see log above)')
            else:
                logging.info('... success')
        else:
            logging.warning("plotIt not found")

    def run_combine(self,entries,additional_args={},debug=False):
        logging.info(f'Running combine on {self.outputDir}')
        if self.combineConfigs is None:
            logging.warning("No combine command in the yaml config file")
            return
        txtPaths = self.getTxtFilesPath()
        missingTxts = []
        for txtPath in txtPaths.values():
            if not os.path.exists(txtPath):
                missingTxts.append(txtPath)
        if len(missingTxts) > 0:
            logging.info('Missing following datacards :')
            for missingTxt in missingTxts:
                logging.info(f'... {missingTxt}')
            logging.info('Datacard not found, will produce it')
            self.run_production()
        
        if entries is None or (isinstance(entries,list) and len(entries)==0):
            entries = self.combineConfigs.keys()

        # Inference setup #
        setup_dir, setup_script = os.path.split(SETUP_PATH)
        script_dir = os.path.dirname(os.path.abspath(__file__))

        # Combine setup #
        combine_dir, env_script = os.path.split(COMBINE_ENV)
        combineModes = ['limits','gof','pulls_impact','prefit','postfit_b','postfit_s']
        slurmIdPerEntry = {}
        # Loop over all the entries inside the combine configuration #
        for entry in entries:
            if entry not in self.combineConfigs.keys():
                raise RuntimeError(f'Combine entry {entry} not in combine config in the yaml file')
            combineCfg = self.combineConfigs[entry]
            combineMode = combineCfg['mode']
            slurmIdPerEntry[entry] = []

            # Make subdirectory #
            logging.info(f'Running entry {entry} with mode {combineMode}')
            if combineMode not in combineModes:
                raise RuntimeError(f'Combine mode {combineMode} not understood')
            subdir = os.path.join(self.datacardPath,entry)
            if os.path.exists(subdir) and len(os.listdir(subdir)) != 0:
                logging.warning(f'Subdirectory {subdir} is not empty')
            if not os.path.exists(subdir):
                os.makedirs(subdir)

            if 'command' not in combineCfg.keys():
                logging.warning('No "command" entry for combine mode "{combineMode}"')

            # Select txt files #
            if not 'bins' in combineCfg.keys() \
                    and (('combine_bins' in combineCfg.keys() and combineCfg['combine_bins']) \
                    or   ('split_bins' in combineCfg.keys() and combineCfg['split_bins'])):
                raise RuntimeError(f'You did not specify the bins in the config entry {entry} but used `combine_ins` or `split_bins`')
            if 'bins' in combineCfg.keys() and list(txtPaths.keys()) == ['all']:
                raise RuntimeError(f'You asked for a single datacard text file, but filter on bins for the command mode {combineMode}')

            binsToUse = []
            if 'bins' in combineCfg.keys():
                if not ('combine_bins' in combineCfg.keys() and not combineCfg['combine_bins']): 
                    # Unless combine_bins == False : use all bins 
                    binsToUse.append(combineCfg['bins'])
                if ('split_bins' in combineCfg.keys() and combineCfg['split_bins']):
                    # If split_bins == True : do once per bin
                    for b in combineCfg['bins']:
                        binsToUse.append([b])
            else:
                binsToUse.append('all')

            # Loop over all the bin combinations (can be all of the one and/or one by one) #
            subdirBinPaths = []
            for binNames in binsToUse:
                # Make bin subdir #
                if len(binsToUse) == 1:
                    subdirBin = subdir
                    binSuffix = ''
                else:
                    if len(binNames) == 1:
                        binSuffix = binNames[0]
                    else:
                        binSuffix = 'combination'
                    subdirBin = os.path.join(subdir,binSuffix)
                    if not os.path.exists(subdirBin):
                        os.makedirs(subdirBin)
                    subdirBinPaths.append(subdirBin)
                
                # Check with argument (either debug or worker mode) #
                if len(binsToUse) > 1 and 'bin' in additional_args.keys():
                    if additional_args['bin'] != binSuffix:
                        continue

                # Check if root output is already in subdir #
                combinedTxtPath = os.path.join(subdirBin,'datacard.txt')
                workspacePath = os.path.join(subdirBin,'workspace.root')
                rootFiles = glob.glob(os.path.join(subdirBin,'*.root'))
                if len(rootFiles) == 0 or rootFiles == [workspacePath]:
                    logging.info(f'No root file found in {subdirBin}, will run the command for mode {combineMode}')
                    txtPathsForCombination = []
                    if binNames == 'all':
                        txtPathsForCombination = list(txtPaths.values())
                    else:
                        for binName in binNames:
                            if binName not in txtPaths.keys():
                                raise RuntimeError(f'{binName} not in hist content')
                            txtPathsForCombination.append(txtPaths[binName])
                    logging.info('Running on txt files :')
                    for txtPathComb in txtPathsForCombination:
                        logging.info(f'\t- {txtPathComb}')


                    # Combine txt files into one #
                    if not os.path.exists(combinedTxtPath):
                        logging.info('Producing combine datacard')
                        for txtPathComb in txtPathsForCombination:
                            logging.debug(f'... {txtPathComb}')
                        logging.debug(f'-> {combinedTxtPath}')
                        if os.path.exists(combinedTxtPath):
                            os.remove(combinedTxtPath)
                        self.combineCards(txtPathsForCombination,combinedTxtPath)
                        logging.info('... done')
                        
                    # Create workspace #
                    if not os.path.exists(workspacePath):
                        logging.info('Producing workspace')
                        workspaceCmd = f"cd {combine_dir}; "
                        workspaceCmd += f"env -i bash -c 'source {env_script} && cd {subdirBin} && text2workspace.py {combinedTxtPath} -o {workspacePath}'"
                        exitCode,output = self.run_command(workspaceCmd,return_output=True,shell=True)
                        if exitCode != 0:
                            if logging.root.level > 10:
                                for line in output:
                                    logging.info(line.strip())
                            raise RuntimeError('Could not produce the workspace, see log above')
                        logging.info('... done')

                    # List systematics for the pulls and impacts #
                    if combineMode == 'pulls_impact':
                        workspaceJson = workspacePath.replace('.root','.json')
                        if not os.path.exists(workspaceJson):
                            worspace_cmd = f"cd {setup_dir}; "
                            worspace_cmd += f"env -i bash -c 'source {setup_script} && cd {script_dir} && python helperWorkspace.py {workspacePath} {workspaceJson}'"
                            rc,output = self.run_command(worspace_cmd,shell=True,return_output=True)
                            if rc != 0: 
                                if logging.root.level > 10:
                                    for line in output:
                                        logging.info(line.strip())
                                raise RuntimeError(f'Could not convert {workspacePath} into {workspaceJson}, see log above')
                            else:
                                logging.info(f'Produced {workspaceJson}')
                        else:
                            logging.info(f'Loaded {workspaceJson}')

                        with open(workspaceJson,'r') as handle:
                            workspaceParams = json.load(handle)

                        # Get systematic names for the jobs #
                        if 'mc_stats' in combineCfg.keys() and combineCfg['mc_stats']:
                            systNames = [key for key in workspaceParams.keys() if key != 'r']
                            # Move mc stats to end of list #
                            systNames.sort(key=lambda x: 'prop_bin' in x) 
                        else:
                            systNames = [key for key in workspaceParams.keys() if key != 'r' and 'prop_bin' not in key]

                    # Submit mode #
                    if 'submit' in combineCfg.keys() and not self.worker:
                        params = combineCfg['submit']
                        args = {'combine':entry}
                        arrayIds = None
                        # Create slurm directories
                        slurmDir  = os.path.join(subdirBin,'batch')
                        logDir    = os.path.join(slurmDir,'logs')
                        outputDir = os.path.join(slurmDir,'output')
                        for directory in [slurmDir,logDir,outputDir]:
                            if not os.path.exists(directory):
                                os.makedirs(directory)

                        n_jobs = 1
                        if combineMode == 'gof':
                            toys = combineCfg['toys-per-job']
                            n_jobs = math.ceil(combineCfg['toys']/toys) + 1 # 1 : data, >1: toys (to match array ids)
                            def fileChecker(f,idx):
                                F = ROOT.TFile(f)
                                valid = True
                                if 'limit' not in [key.GetName() for key in F.GetListOfKeys()]:
                                    valid = False
                                    logging.debug(f'Tree limit not found in {f}')
                                else:
                                    tree = F.Get('limit')
                                    if idx == 1 and tree.GetEntries() != 1:
                                        logging.debug(f'Tree limit in {f} has {tree.GetEntries()} entries instead of 1')
                                        valid = False
                                    if idx > 1 and tree.GetEntries() != toys:
                                        logging.debug(f'Tree limit in {f} has {tree.GetEntries()} entries instead of {toys}')
                                        valid = False
                                F.Close()
                                return valid
                        elif combineMode == 'pulls_impact':
                            n_jobs = len(systNames) + 1 # 1 = initial fit, >1: all the systematics
                            def fileChecker(f,idx):
                                F = ROOT.TFile(f)
                                valid = True
                                if 'limit' not in [key.GetName() for key in F.GetListOfKeys()]:
                                    valid = False
                                    logging.debug(f'Tree limit not found in {f}')
                                else:
                                    tree = F.Get('limit')
                                    if tree.GetEntries() != 3:
                                        logging.debug(f'Tree limit in {f} has {tree.GetEntries()} entries instead of 3')
                                        valid = False
                                F.Close()
                                return valid

                        idxs = [_ for _ in range(1,n_jobs+1)]

                        subScript = os.path.join(slurmDir,'slurmSubmission.sh')

                        if not os.path.exists(subScript):
                            logging.info(f'{entry} : submitting {n_jobs} jobs')
                            jobArrayArgs = []
                            if n_jobs == 1:
                                subScript = self.writeSbatchCommand(slurmDir,params,args)
                            else:
                                for idx in idxs:
                                    subsubdir = os.path.join(outputDir,str(idx))
                                    if not os.path.exists(subsubdir):
                                        os.makedirs(subsubdir)
                                    jobArgs = {**args,'combine_args':f'idx={idx}'} 
                                    if binSuffix != '':
                                        jobArgs['combine_args'] += f' bin={binSuffix}'
                                    jobArrayArgs.append(jobArgs)
                                subScript = self.writeSbatchCommand(slurmDir,params,jobArrayArgs)
                        else:
                            arrayIds = []
                            logging.info(f'{entry}: found batch script, will look for unfinished jobs')
                            for idx in idxs:
                                subsubdir = os.path.join(outputDir,str(idx))
                                rootfile = glob.glob(os.path.join(subsubdir,'higgs*.root')) 
                                if len(rootfile) == 0:
                                    logging.debug(f'Root output not found in {subsubdir}')
                                    arrayIds.append(str(idx))
                                else:
                                    if not fileChecker(rootfile[0],idx):
                                        arrayIds.append(str(idx))
                        if arrayIds is None:
                            slurmCmd = f'sbatch {subScript}'
                        else:
                            if len(arrayIds) > 0:
                                logging.info('... will resubmit the following array ids : '+','.join(arrayIds))
                                if combineMode == 'pulls_impact':
                                    for arrayId in arrayIds:
                                        if int(arrayId) == 1: 
                                            logging.info(f'\t{arrayId:4s} -> initial fit')
                                        else:
                                            logging.info(f'\t{arrayId:4s} -> {systNames[int(arrayId)-2]}')
                                slurmCmd = f"sbatch --array={','.join(arrayIds)} {subScript}"
                            else:
                                logging.info('... all jobs have succeeded')
                                slurmCmd = ''
                        if slurmCmd != '' and not debug:
                            logging.debug(f'Submitting {subScript}')
                            rc,output = self.run_command(slurmCmd,return_output=True,shell=True)
                            slurm_id = None
                            for line in output:
                                if logging.root.level > 10:
                                    logging.info(line.strip())
                                numbers = re.findall(r'\d+',line.strip())
                                if len(numbers) == 1:
                                    slurm_id = numbers[0]
                                    if len(slurm_id) != 8:
                                        slurm_id = None
                            if slurm_id is None:
                                logging.error('Slurm job id could not be found')
                            else:
                                slurmIdPerEntry[entry].append(slurm_id)
                            if arrayIds is None or len(arrayIds) == len(idxs): 
                                # For single job we want to wait until something has finished before going to finalize
                                # For array job, as soon at least one job has run, produce finalize with intermediate results
                                continue

                    # Produce command #
                    combineCmd = combineCfg['command'].format(workspacePath)

                    if 'submit' in combineCfg.keys() and not self.worker: # pass to finalize part
                        combineCmd = ''
                    else:
                        if 'idx' in additional_args.keys():
                            subdirBin = os.path.join(subdirBin,'batch','output',additional_args['idx'])
                            if '--seed' not in combineCmd:
                                combineCmd += f" --seed {additional_args['idx']}"
                            if combineMode == 'gof':
                                if int(additional_args['idx']) > 1: # Toys case
                                    combineCmd += f" --toys {combineCfg['toys-per-job']}"
                            if combineMode == 'pulls_impact':
                                if int(additional_args['idx']) == 1: # initial fit 
                                    combineCmd += " --algo singles"
                                else: # One job per systematic
                                    systName = systNames[int(additional_args['idx'])-2] # Offset 0->1 (for array id) + #1 is for initial fit
                                    combineCmd += f" --algo impact -P {systName}  --floatOtherPOIs 1 --saveInactivePOI 1 "
                                if 'unblind' in combineCfg.keys() and combineCfg['unblind']:
                                    if '--toys' in combineCmd or '-t' in combineCmd:
                                        raise RuntimeError(f'You want to unblind entry {entry} but used `--toys` in the combine command, these are exclusive')
                                else:
                                    if not '--toys' in combineCmd and not '-t' in combineCmd:
                                        combineCmd += " --toys -1"

                    if combineCmd != '':
                        fullCombineCmd  = f"cd {combine_dir}; "
                        fullCombineCmd += f"env -i bash -c 'source {env_script} && cd {subdirBin} && {combineCmd}'"

                        logging.info('Extensive command is below')
                        logging.info(fullCombineCmd)

                        # Run command #
                        exitCode,output = self.run_command(fullCombineCmd,return_output=True,shell=True)
                        path_log = os.path.join(subdirBin,'log_combine.out')
                        with open(path_log,'w') as handle:
                            for line in output:
                                handle.write(line)

                        if exitCode != 0:
                            raise RuntimeError(f'Something went wrong in combine, see log in {path_log}') 
                        else:
                            logging.info('... done')
                else:
                    logging.warning(f"Already a root file in subdirectory {subdirBin}, will not run the combine command again (remove and run again if needed)")

                if self.worker:
                    continue

                # Per bin mode finalize #
                # Saving limit values in json #
                if combineMode == 'limits':
                    limits = {}
                    with open(os.path.join(subdirBin,'log_combine.out'),'r') as output:
                        for line in output:
                            l_str = None
                            r_str = None
                            if line.startswith('Observed'):
                                logging.info("\t"+line.strip())
                                r_str = line.split(':')[1]
                            elif line.startswith('Expected'):
                                logging.info("\t"+line.strip())
                                l_str,r_str = line.split(':')
                            else:
                                continue
                            if l_str is None: 
                                level = -1.
                            else:   
                                level = float(re.findall("\d+.\d+",l_str)[0])
                            limits[level] = float(re.findall("\d+.\d+",r_str)[0])
                            
                    path_limits = os.path.join(self.datacardPath,subdirBin,'limits.json')
                    with open(path_limits,'w') as handle:
                        json.dump(limits,handle,indent=4)
                    logging.info(f'Saved limits as {path_limits}') 

                # Producing GoF #
                if combineMode == 'gof':
                    if len(additional_args) > 0:
                        continue
                    # Hadding in case the file is not there
                    toys_file = os.path.join(subdirBin,'higgsCombine.GoodnessOfFit.toys.root')
                    data_file = os.path.join(subdirBin,'higgsCombine.GoodnessOfFit.data.root')
                    if not os.path.exists(toys_file):
                        hadd_cmd = ['hadd',toys_file]
                        for rootfile in glob.glob(os.path.join(subdirBin,'batch','output','*','higgs*root')):
                            if os.path.join(subdirBin,'batch','output','1','higgs') not in rootfile: # avoid taking data measurement in toys
                                hadd_cmd.append(rootfile)
                        rc = self.run_command(hadd_cmd)
                        if rc != 0:
                            raise RuntimeError("Hadd command `{' '.join(hadd_cmd)}` failed")
                        logging.debug(f'Created {toys_file}')
                    if not os.path.exists(data_file):
                        rootfile = glob.glob(os.path.join(subdirBin,'batch','output','1','*root'))
                        if len(rootfile) != 1:
                            raise RuntimeError(f'Cannot find file {rootfile.__repr__}()')
                        shutil.copy(rootfile[0],data_file)
                        logging.debug(f'Created {data_file}')

                    # Inspect root files to get data and toys test statistics #
                    def getLimitValues(f):
                        F = ROOT.TFile(f)
                        t = F.Get('limit')
                        limits = [event.limit for event in t]
                        F.Close()
                        return limits

                    limit_toys = getLimitValues(toys_file)
                    limit_data = getLimitValues(data_file)

                    # Find algorithm #
                    algorithm = None
                    for i,arg in enumerate(combineCfg['command'].split()):
                        if arg == '--algo': # Used --algo ...
                            algorithm = combineCmd.split()[i+1]
                        if '--algo' in arg and arg != '--algo': # Used --algo=...
                            algorithm = arg.replace('--algo','').replace('=','')
                    if algorithm is None:
                        raise RuntimeError('Could not understand algorithm in combine command, did you use `--algo=...` or `--algo ...` ?')

                    # Produce plots #
                    path_pdf = os.path.join(subdirBin,'gof.pdf')
                    content = {'paths'      : path_pdf,
                               'data'       : limit_data[0],
                               'toys'       : limit_toys,
                               'algorithm'  : algorithm,
                               'campaign'   : self.era}
                    if 'plotting' in combineCfg.keys():
                        content.update(combineCfg['plotting'])
                    path_json = os.path.join(subdirBin,'gof.json')
                    logging.info(f'Saved {os.path.basename(subdirBin)} data in {path_json}')
                    with open(path_json,'w') as handle:
                        json.dump(content,handle,indent=4)

                    gof_cmd = f"cd {setup_dir}; "
                    gof_cmd += f"env -i bash -c 'source {setup_script} && cd {script_dir} && python helperInference.py dhi.plots.gof.plot_gof_distribution {path_json}'"
                    rc,output = self.run_command(gof_cmd,shell=True, return_output=True)
                    if rc != 0: 
                        if logging.root.level > 10:
                            for line in output:
                                logging.info(line.strip())
                        logging.error('Failed to produce gof plot, see log above')
                    else:
                        logging.info(f'Produced {path_pdf}')

                # Pulls and impacts #
                if combineMode == 'pulls_impact':
                    data = {"params": []}
                    # Inspect root files to get poi values #
                    for idx in range(1,len(systNames)+2):
                        rootfiles = glob.glob(os.path.join(subdirBin,'batch','output',str(idx),'higgs*root'))
                        values = {'r':[-1,-1,-1]}
                        if idx != 1:
                            systName = systNames[idx-2]
                            values[systName] = [-1,-1,-1]
                        if len(rootfiles) != 0:
                            F = ROOT.TFile(rootfiles[0])
                            if 'limit' in [key.GetName() for key in F.GetListOfKeys()]:
                                tree = F.Get('limit')
                                values['r'] = [event.r for event in tree]
                                if idx != 1:
                                    values[systName] = [getattr(event,systName) for event in tree]
                                if len(values['r']) != 3:
                                    values['r'] = [-9999.,-9999.,-9999.]
                                    values[systName] = [-9999.,-9999.,-9999.]
                                    
                            F.Close()
                        if idx == 1:
                            data['POIs'] = [{"name": 'r', "fit": [values['r'][1],values['r'][0],values['r'][2]]}]
                        else:
                            d = {}
                            params = workspaceParams[systName]
                            d['name']     = systName
                            d["type"]     = params['type']
                            d["groups"]   = params['groups']
                            d["prefit"]   = params['prefit']
                            d["fit"]      = [values[systName][1],values[systName][0],values[systName][2]]
                            d["r"]        = [values['r'][1],values['r'][0],values['r'][2]]
                            d["impacts"]  = {
                                'r' : [
                                        d['r'][1] - d['r'][0],
                                        d['r'][2] - d['r'][1],
                                ]
                            }
                            d["impact_r"] = max(map(abs, d["impacts"]['r']))

                            data['params'].append(d)

                    path_pdf = os.path.join(subdirBin,'pulls_impact.pdf')
                    content = {
                        'paths'     : [path_pdf],
                        'data'      : data,
                        'poi'       : 'r',
                        'campaign'  : self.era,
                    }
                    if 'plotting' in combineCfg.keys():
                        content.update(combineCfg['plotting'])
                    path_json = os.path.join(subdirBin,'pulls_impact.json')
                    with open(path_json,'w') as handle:
                        json.dump(content,handle,indent=4)

                    pull_cmd = f"cd {setup_dir}; "
                    pull_cmd += f"env -i bash -c 'source {setup_script} && cd {script_dir} && python helperInference.py dhi.plots.pulls_impacts.plot_pulls_impacts {path_json}'"
                    rc,output = self.run_command(pull_cmd,shell=True, return_output=True)
                    if rc != 0: 
                        if logging.root.level > 10:
                            for line in output:
                                logging.info(line.strip())
                        logging.error('Failed to produce pulls and impact plot, see log above')
                    else:
                        logging.info(f'Produced {path_pdf}')

                # Producing prefit and postfit plots #
                if 'prefit' in combineMode or 'postfit' in combineMode:
                    config = combineCfg['plotting']
                    fitdiagFile = glob.glob(os.path.join(subdirBin,'fitDiagnostic*root'))
                    if len(fitdiagFile) == 0:
                        raise RuntimeError("Could not find any fitdiag file in subdir")
                    fitdiagFile = fitdiagFile[0]
                    config['main']['fitdiagnosis'] = fitdiagFile.replace('/auto','')
                    # Generate dat config file #
                    path_dict = os.path.join(subdirBin,'dict.dat')
                    self.makeDictForPostFit(path_dict,config)
                    # Run plotting #
#                    fitplot_cmd = f"cd {setup_dir}; "
#                    postfit_script = os.path.join(INFERENCE_PATH,'dhi','scripts','postfit_plots.py')
#                    fitplot_cmd += f"env -i bash -c 'source {setup_script} && cd {script_dir} && python {postfit_script} --plot_options_dict {path_dict} --output_folder {subdirBin} "
#                    if 'postfit' in combineMode:
#                        fitplot_cmd += " --doPostFit "
#                    if 'unblind' in combineCfg.keys() and combineCfg['unblind']:
#                        fitplot_cmd += " --unblind "
#                    fitplot_cmd += "'"
#
#                    rc,output = self.run_command(fitplot_cmd,shell=True, return_output=True)
#                    if rc != 0: 
#                        if logging.root.level > 10:
#                            for line in output:
#                                logging.info(line.strip())
#                        logging.error('Failed to produce fit plots, see log above')
#                    else:
#                        logging.info(f'Produced plots in {subdirBin}')

                    # Nuisances likelihoods #
                    if 'postfit' in combineMode:
                        fit_name = None
                        if combineMode == 'postfit_b':
                            fit_name = 'fit_b'
                        if combineMode == 'postfit_s':
                            fit_name = 'fit_s'
                        path_pdf = os.path.join(subdirBin,'nuisances.pdf')
                        resultFile = glob.glob(os.path.join(subdirBin,'higgs*root'))
                        if len(resultFile) == 0:
                            raise RuntimeError(f'Could not find result file in {subdirBin}')

                        content = {
                            'paths'                 : [path_pdf],
                            'poi'                   : 'r',
                            'fit_name'              : fit_name,
                            'workspace'             : resultFile[0],
                            'dataset'               : resultFile[0],
                            'fit_diagnostics_path'  : fitdiagFile,
                            'y_min'                 : 0.,
                            'y_max'                 : 5.,
                        }
                        path_json = os.path.join(subdirBin,'nuisances.json')
                        with open(path_json,'w') as handle:
                            json.dump(content,handle,indent=4)
                        getter = {
                            'workspace'     : {'type':'ROOT','name':'w'},
                            'dataset'       : {'type':'ROOT','name':'toys/toy_asimov'},
                        }
                        path_getter = os.path.join(subdirBin,'nuisances_getter.json')
                        with open(path_getter,'w') as handle:
                            json.dump(getter,handle,indent=4)

                        pull_cmd = f"cd {setup_dir}; "
                        pull_cmd += f"env -i bash -c 'source {setup_script} && cd {script_dir} && python helperInference.py dhi.plots.likelihoods.plot_nuisance_likelihood_scans {path_json} {path_getter}'"
                        rc,output = self.run_command(pull_cmd,shell=True, return_output=True)
                        if rc != 0: 
                            if logging.root.level > 10:
                                for line in output:
                                    logging.info(line.strip())
                            logging.error('Failed to produce nuisance plots, see log above')
                        else:
                            logging.info(f'Produced {path_pdf}')

#                        content_log = {
#                            'paths'                 : [os.path.join(subdirBin,'nuisances_log.pdf')],
#                            'poi'                   : 'r',
#                            'fit_name'              : fit_name,
#                            'workspace'             : resultFile[0],
#                            'fit_diagnostics_path'  : fitdiagFile,
#                            'y_log'                 : True,
#                        }

                    from fitdiag_plots import create_postfit_plots
                    from collections import OrderedDict
                    info_bin = eval(open(path_dict, "r").read())
                    folder = None
                    if combineMode == 'prefit':
                        folder = 'shapes_prefit'
                    elif combineMode == 'postfit_b':
                        folder = 'shapes_fit_b'
                    elif combineMode == 'postfit_s':
                        folder = 'shapes_fit_s'
                    else:
                        raise RuntimeError(f'fit mode {combineMode} not understood')

                    for key, binCfg in info_bin.items():
                        unblind = False
                        for plotCfg in config['plots']:
                            if plotCfg['name'] in key and 'unblind' in plotCfg.keys() and plotCfg['unblind']:
                                unblind = True
                        create_postfit_plots(path                    = subdirBin,
                                             fit_diagnostics_path    = fitdiagFile,
                                             normalize_X_original    = False,
                                             doPostFit               = 'postfit' in combineMode,
                                             divideByBinWidth        = False,
                                             folder                  = folder,
                                             bin                     = binCfg,
                                             binToRead               = key,
                                             unblind                 = unblind)


                        
                                        
            # Combined finalize mode for all bins (in case there was) #
            if len(subdirBinPaths) > 1 and not self.worker:
                # Goodness of fit : combine into one plot
                if combineMode == 'gof':
                    logging.info('GoF combination :')
                    path_pdf = os.path.join(subdir,'gof.pdf')
                    content = {'paths'      : path_pdf,
                               'n_bins'     : 32,
                               'data'       : [],
                               'algorithm'  : 'saturated',
                               'campaign'   : self.era}
                    for subdirBinPath in subdirBinPaths:
                        path_json = os.path.join(subdirBinPath,'gof.json')
                        if not os.path.exists(path_json):
                            logging.warning(f'Could not load {path_json}, will continue')
                            continue
                        logging.info(f'Loading {path_json}')
                        with open(path_json,'r') as handle:
                            subCont = json.load(handle)
                            content['data'].append({'data' : subCont['data'],
                                                    'toys' : subCont['toys'], 
                                                    'name' : os.path.basename(subdirBinPath)})
                    
                    path_json = os.path.join(subdir,'gof.json')
                    logging.info(f'Saved combination data in {path_json}')
                    with open(path_json,'w') as handle:
                        subCont = json.dump(content,handle,indent=4)

                    if len(content['data']) > 0:
                        gof_cmd = f"cd {setup_dir}; "
                        gof_cmd += f"env -i bash -c 'source {setup_script} && cd {script_dir} && python helperInference.py dhi.plots.gof.plot_gofs {path_json}'"
                        rc,output = self.run_command(gof_cmd,shell=True, return_output=True)
                        if rc != 0: 
                            if logging.root.level > 10:
                                for line in output:
                                    logging.info(line.strip())
                            logging.error('Failed to produce gof plot, see log above')
                        else:
                            logging.info(f'Produced {path_pdf}')

        # Slurm ID printout #
        if len([v for key,val in slurmIdPerEntry.items() for v in val]) > 0: # If any slurm Id has been registered
            logging.info('Following slurm job ids were submitted')
            for entry,slurm_ids in slurmIdPerEntry.items():
                logging.info(f'... Entry {entry} : '+' '.join(slurm_ids))

    def getTxtFilesPath(self):
        if self.textfiles is None:
            initTextFiles = {histName:os.path.join(self.datacardPath,f'{histName}_{self.era}.txt') 
                                for histName in self.hist_conv.keys()}
        elif isinstance(self.textfiles,str):
            if '{}' in self.textfiles:
                initTextFiles = {histName:os.path.join(self.datacardPath,
                                              self.textfiles.format(f'{histName}_{self.era}.txt'))
                                    for histName in self.hist_conv.keys()}
            else:
                initTextFiles = {'all',os.path.join(self.datacardPath,self.textfiles)}
        else:
            raise RuntimeError('Format of "textfiles" entry not understood')
        return initTextFiles

        
    def combineCards(self,initTextPaths,outputTextPath):
        for initTextPath in initTextPaths:
            if not os.path.exists(initTextPath):
                raise RuntimeError(f'File {initTextPath} not found')
        if len(initTextPaths) == 0:
            raise RuntimeError('No initial txt paths')
        # combine the datacards #
        combineCmd = "combineCards.py "
        for initTextPath in initTextPaths:
            if 'auto' in initTextPath:
                initTextPath = initTextPath.replace('/auto','')
            binName = os.path.basename(initTextPath).replace('.txt','')
            combineCmd += f" {binName}={initTextPath} "
        combine_dir, env_script = os.path.split(COMBINE_ENV)
        fullCombineCmd  = f"cd {combine_dir}; "
        fullCombineCmd += f"env -i bash -c 'source {env_script} && {combineCmd}'"
        rc, output = self.run_command(fullCombineCmd,return_output=True,shell=True)
        if rc != 0:
            raise RuntimeError(f'combineCards failed with {outputTextPath}')
        with open(outputTextPath,'w') as handle:
            for line in output:
                handle.write(line)

    def makeDictForPostFit(self,filePath,config):
        def colorTrad(s):
            if isinstance(s,int):
                return s
            if isinstance(s,str):
                return ROOT.TColor.GetColor(s)
            
        with open(filePath,'w') as handle:
            handle.write("{\n")
            for plotCfg in config['plots']:
                handle.write(f"\t'{plotCfg['name']}_{self.era}' : "+"{\n")
                for key,val in config['main'].items():
                    if isinstance(val,str):
                        handle.write(f"\t\t'{key}': '{val}',\n")
                    else:
                        handle.write(f"\t\t'{key}': {val},\n")
                handle.write(f"\t\t'header_legend' : '{plotCfg['header']}',\n")
                handle.write(f"\t\t'era' : {self.era},\n")
                handle.write("\t\t'align_cats': [")
                for b in plotCfg['bins']:
                    if b not in self.hist_conv.keys():
                        raise RuntimeError(f'Bin {b} not in defined histograms')
                    handle.write(f"'{b}_{self.era}',")
                handle.write("],\n")
                handle.write("\t\t'align_cats_labels'  :[ "+','.join([f"['{l}']" for l in plotCfg['labels']])+"],\n")
                handle.write("\t\t'align_cats_labelsX' :[ "+','.join([f"{l}" for l in plotCfg['labelpos']])+"],\n")
                handle.write("\t\t'procs_plot_options_bkg' : OrderedDict(\n\t\t\t[\n")
                for group,groupCfg in self.groups.items():
                    if groupCfg['type'] == 'mc':
                        handle.write(f"\t\t\t('{group}', "+"{")
                        color = colorTrad(groupCfg['fill-color'])
                        handle.write(f"'color' : {color}, ")
                        if 'fill-type' in groupCfg.keys():
                            handle.write(f"'fillStype' : {groupCfg['fill-type']}, ")
                        else:
                            handle.write("'fillStype' : 1001, ")
                        handle.write(f"'label': '{groupCfg['legend']}', ")
                        handle.write("'make border' : True}),\n")
                handle.write("\t\t\t]\n\t\t),\n")
                handle.write("\t\t'procs_plot_options_sig' : OrderedDict(\n\t\t\t[\n")
                for group,groupCfg in self.groups.items():
                    if groupCfg['type'] == 'signal':
                        if 'hide' in groupCfg.keys() and groupCfg['hide']:
                            continue
                        handle.write(f"\t\t\t('{group}', "+"{")
                        color = colorTrad(groupCfg['line-color'])
                        handle.write(f"'color' : {color}, ")
                        if 'fill-type' in groupCfg.keys():
                            handle.write(f"'fillStype' : {groupCfg['fill-type']}, ")
                        else:
                            handle.write("'fillStype' : 3351, ")
                        handle.write(f"'label': '{groupCfg['legend']}', ")
                        if 'scale' in groupCfg.keys():
                            handle.write(f"'scaleBy': {groupCfg['scale']}, ")
                        handle.write("'scaleBy' : 1.,")
                        handle.write("'make border' : True}),\n")
                handle.write("\t\t\t]\n\t\t),\n")
                handle.write("\t},")
            handle.write("}")

    def writeSbatchCommand(self,subdir,params={},args={}):
        # Generate file path #
        filepath = os.path.join(subdir,f'slurmSubmission.sh')
        scriptName = os.path.abspath(__file__)
        log_dir = os.path.join(subdir,'logs')
        if not os.path.exists(log_dir):
            os.makedirs(log_dir)

        # Check job type (array or not) #
        if isinstance(args,list): # Array submission
            for arg in args:
                if not isinstance(arg,dict):
                    raise RuntimeError('Argument needs to be a dict : '+arg.__repr__())
            N_arrays = len(args)
        elif isinstance(args,dict): # single job submission
            N_arrays = 0
        else:
            raise RuntimeError('Either a dict or a list needs to be used for args')


        with open(filepath,'w') as handle:
            # Write header #
            handle.write('#!/bin/bash\n\n')
            if N_arrays > 0:
                handle.write('#SBATCH --job-name=datacard_%a\n')
                handle.write(f'#SBATCH --output={os.path.join(log_dir,"log-%A_%a.out")}\n')
                handle.write(f'#SBATCH --array=1-{N_arrays}\n')
            else:
                handle.write('#SBATCH --job-name=datacard\n')
                handle.write(f'#SBATCH --output={os.path.join(log_dir,"log-%j.out")}\n')
            for key,val in params.items():
                handle.write(f'#SBATCH --{key}={val}\n')
            handle.write('\n')
            handle.write('\n')
                
            # command #
            if N_arrays > 0:
                handle.write('cmds=(\n')
                for arg in args:
                    handle.write(f'"python3 {scriptName} --yaml {self.configPath} --worker ')
                    if logging.root.level == 10:
                        handle.write(' -v ')
                    for key,val in arg.items():
                        handle.write(f'--{key} {val} ')
                    handle.write('"\n')
                handle.write(')\n')
                handle.write('${cmds[$SLURM_ARRAY_TASK_ID-1]}\n')
            else:
                handle.write(f'python3 {scriptName} --yaml {self.configPath} --worker ')
                if logging.root.level == 10:
                    handle.write(' -v ')
                for key,val in args.items():
                    handle.write(f'--{key} {val} ')
                handle.write("\n")

        logging.debug(f'Generated submission script {filepath}')
        return filepath

 
    @staticmethod
    def run_command(command,return_output=False,**kwargs):
        process = subprocess.Popen(command,universal_newlines=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,**kwargs)
        # Poll process for new output until finished #
        output = []
        while True:
            try:
                nextline = process.stdout.readline()
            except UnicodeDecodeError:
                continue
            if nextline == '' and process.poll() is not None:
                break
            logging.debug(nextline.strip())
#                sys.stdout.write(nextline)
#                sys.stdout.flush()
            if return_output:
                output.append(nextline)
        process.communicate()
        exitCode = process.returncode
        if return_output:
            return exitCode,output
        else:
            return exitCode


    @staticmethod
    def which(program):
        import os
        def is_exe(fpath):
            return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    
        fpath, fname = os.path.split(program)
        if fpath:
            if is_exe(program):
                return program
        else:
            for path in os.environ["PATH"].split(os.pathsep):
                exe_file = os.path.join(path, program)
                if is_exe(exe_file):
                    return exe_file
    
        return None



 

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Produce datacards')
    parser.add_argument('--yaml', action='store', required=True, type=str,
                        help='Yaml containing parameters')
    parser.add_argument('--pseudodata', action='store_true', required=False, default=False,
                        help='Whether to use pseudo data (data = sum of MC) [default = False]')
    parser.add_argument('--plotIt', action='store_true', required=False, default=False,
                        help='Browse datacard files and produce plots via plotIt \
                              (note : produced by default when datacards are produced) ') 
    parser.add_argument('--split', action='store', required=False, default=None, nargs='*',
                        help='If used without argument, will process all categories in serie, \
                              otherwise will process the categories given in argument')
    parser.add_argument('-j,','--jobs', action='store', required=False, default=None, type=int,
                        help='Number of parallel processes (only useful with `--split`)')
    parser.add_argument('--combine', action='store', required=False, default=None, nargs='*',
                        help='Run combine on the txt datacard only')
    parser.add_argument('--combine_args', action='store', required=False, default=[], nargs='*',
                        help='Additional args for the combine commands (used in jobs, keep for debug)')
    parser.add_argument('--custom', action='store', required=False, default=None, nargs='*',
                        help='Format the yaml file')
    parser.add_argument('-v','--verbose', action='store_true', required=False, default=False,
                        help='Verbose mode')
    parser.add_argument('--worker', action='store_true', required=False, default=False,
                        help='Force working locally')
    parser.add_argument('--debug', action='store_true', required=False, default=False,
                        help='Do not send jobs')
    args = parser.parse_args()

    logging.basicConfig(level   = logging.DEBUG,
                        format  = '%(asctime)s - %(levelname)s - %(message)s',
                        datefmt = '%m/%d/%Y %H:%M:%S')

    # Verbose level #
    if not args.verbose:
        logging.getLogger().setLevel(logging.INFO)

    # File checking #
    if args.yaml is None:
        raise RuntimeError("Must provide the YAML file")
    if not os.path.isfile(args.yaml):
        raise RuntimeError("YAML file {} is not a valid file".format(args.yaml))

    # Yaml parsing #
    if args.custom is not None:
        # Custom command line argument for parsing #
        formatting = {}
        for arg in args.custom:
            if '=' in arg:
                formatting[arg.split('=')[0]] = arg.split('=')[1]
            else:
                logging.warning(f'`--custom {arg}` will be ignored because no `=`')
        config = yaml.load({'filename':args.yaml,'formatting':formatting},
                      Loader=YMLIncludeLoader)
    else:
        # Classic parse #
        with open(args.yaml,'r') as handle:
            config = yaml.load(handle,Loader=YMLIncludeLoader)

    # Content checks #
    required_items = ['path','outputDir','yamlName','hist_conv','groups']
    if any(item not in config.keys() for item in required_items): 
        raise RuntimeError('Your config is missing the following items :'+ \
                ','.join([item for item in required_items if item not in config.keys()]))

    # Splitting #
    global_hist_conv = copy.deepcopy(config['hist_conv']) # Run all at once
    categoriesToRun = [config['hist_conv'].keys()]
        
    if args.split is not None:
        if len(args.split) == 0:  # run all sequantially
            categoriesToRun = [[cat] for cat in categoriesToRun[0]]
            if 'regroup' in config.keys():
                categoriesToRun = [val for val in config['regroup'].values()]
        else: # Run only requested
            if 'regroup' not in config.keys():
                if not all([item in config['hist_conv'].keys() for item in args.split]):
                    raise RuntimeError('Category(ies) requested in split not found : '+\
                            ','.join([item for item in args.split if item not in hist_conv.keys()]))
                categoriesToRun = [[cat] for cat in args.split]
            else:
                categoriesToRun = [[val for key in args.split for val in config['regroup'][key]]]

    # Instantiate #
    instances = []
    for icat,categories in enumerate(categoriesToRun):
        config['hist_conv'] = {cat:global_hist_conv[cat] for cat in categories}
        instances.append(DataCard(configPath      = os.path.abspath(args.yaml),
                                  worker          = args.worker,
                                  pseudodata      = args.pseudodata,
                                  logName         = f'log_{icat}.log',
                                  **config) )
    if args.combine is not None: 
        if args.split is not None:
            raise RuntimeError('Safer not to use `--split` with `--combine`')
        assert len(instances) == 1
        instance = instances[0]
        combine_args = {}
        for arg in args.combine_args:
            if '=' in arg:
                combine_args[arg.split('=')[0]] = arg.split('=')[1]
            else:
                logging.warning(f'`--combine_args {arg}` will be ignored because no `=`')
        instance.run_combine(args.combine,combine_args,debug=args.debug)
    else:
        if args.jobs is None:
            plotIt_configs = []
            for instance in instances:
                if not args.plotIt:
                    instance.run_production()
                plotIt_configs.append(instance.prepare_plotIt())
        else:
            def run(instance,methods):
                for method in methods:
                    assert hasattr(instance,method)
                    output = getattr(instance,method)()
                return output
                    
            methods = ['prepare_plotIt']
            if not args.plotIt:
                methods.insert(0,'run_production')
            print ('starting pool')
            with mp.Pool(processes=args.jobs) as pool:
                plotIt_configs = pool.starmap(run,[(instance,methods) for instance in instances])
            print ('finished pool')
                

        # Merge and run plotit (single thread anyway) #
        plotIt_config = DataCard.merge_plotIt(plotIt_configs)
        print ('merged')
        instances[0].run_plotIt(plotIt_config)
        print ('ran')

