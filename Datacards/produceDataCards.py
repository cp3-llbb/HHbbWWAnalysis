import os
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
import numpy as np
import math
from itertools import chain
from pprint import pprint
import numpy as np
import enlighten
import ROOT

from yamlLoader import YMLIncludeLoader

from IPython import embed
ROOT.gROOT.SetBatch(True) 
combine_env = os.path.join(os.path.dirname(os.path.abspath(__file__)),'HiggsAnalysis','CombinedLimit','env_standalone.sh')

class DataCard:
    def __init__(self,datacardName=None,configPath=None,path=None,yamlName=None,worker=None,groups=None,shapeSyst=None,normSyst=None,hist_conv=None,era=None,use_syst=False,root_subdir=None,DYnonclosure=None,pseudodata=False,rebin=None,textfiles=None,produce_plots=False,legend=None,combineConfigs=None,**kwargs):
        self.datacardName   = datacardName
        self.configPath     = configPath
        self.path           = path
        self.worker         = worker
        self.hist_conv      = hist_conv
        self.era            = str(era)
        self.use_syst       = use_syst
        self.root_subdir    = root_subdir
        self.pseudodata     = pseudodata
        self.DYnonclosure   = DYnonclosure
        self.rebin          = rebin
        self.textfiles      = textfiles
        self.produce_plots  = produce_plots
        self.legend         = legend
        self.combineConfigs    = combineConfigs

        self.groups = self.includeEntry(groups,'Groups')
        self.normSyst = self.includeEntry(normSyst,'Norm syst')
        self.shapeSyst = self.includeEntry(shapeSyst,'Shape syst')
               
        self.yaml_dict = self.loadYaml(self.path,yamlName)

        if self.pseudodata:
            logging.info('Will use pseudodata')
            self.groups = self.generatePseudoData(self.groups)
            if self.datacardName is not None:
                self.datacardName += "_pseudodata"
        self.content = {histName:{g:{} for g in self.groups.keys()} for histName in self.hist_conv.keys()}
        self.systPresent = {histName:{group:[] for group in self.groups.keys()} for histName in self.hist_conv.keys()}

        self.datacardPath = os.path.join(os.path.abspath(os.path.dirname(__file__)),'datacards',self.datacardName)


    def run_production(self):
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
        else: 
            raise RuntimeError(f"{name} format not understood")

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
             
    def loadYaml(self,paths,yamlNames):
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
        if self.datacardName is None:
            raise RuntimeError("Datacard name is not set")

        shapes = {histName:os.path.join(self.datacardPath,f"{histName}_{self.era}.root") for histName in self.content.keys()}

        # Save root file #
        if not os.path.exists(self.datacardPath):
            os.makedirs(self.datacardPath)
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
                        CMSNameUp   = f"{group}_{CMSName}Up"
                        CMSNameDown = f"{group}_{CMSName}Down"

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
                                raise RuntimeError("Could not find {} in systematic dict".format(baseSyst))
                            CMSName = self.shapeSyst[systName]
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
                            for line in systList:
                                writer.addFooter(line)
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
        #val     = 0.
        #valup   = 0.
        #valdown = 0.
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


    def applyDYNonClosure(self):
        logging.info('Applying DY non closure shift')
        from DYNonClosure import NonClosureSystematic
        dyObj = NonClosureSystematic.load_from_json(self.DYnonclosure['path'])
        for histName in self.content.keys():
            if histName not in self.DYnonclosure['hist'].keys():
                continue
            catName = self.DYnonclosure['hist'][histName]['entry']
            closureSystName = self.DYnonclosure['hist'][histName]['name']
            dyObj.SetCategory(catName)
            logging.info("... Histogram {:20s} -> applying non closure with key {:20s}".format(histName,catName))
            #dyObj.PlotFitToPdf('shapes_'+catName)
            group = self.DYnonclosure['group']

            # Add non closure systematics #
            if self.use_syst:
                ((h_shape1_up,h_shape1_down),(h_shape2_up,h_shape2_down)) = dyObj.GetNonClosureShape(self.content[histName][group]['nominal'])
                h_shape1_up.SetName(self.content[histName][group]['nominal'].GetName()+f'_{closureSystName}_shape1Up')
                h_shape1_down.SetName(self.content[histName][group]['nominal'].GetName()+f'_{closureSystName}_shape1Down')
                h_shape2_up.SetName(self.content[histName][group]['nominal'].GetName()+f'_{closureSystName}_shape2Up')
                h_shape2_down.SetName(self.content[histName][group]['nominal'].GetName()+f'_{closureSystName}_shape2Down')

                self.content[histName][group][f'{closureSystName}_shape1Up']   = h_shape1_up
                self.content[histName][group][f'{closureSystName}_shape1Down'] = h_shape1_down
                self.content[histName][group][f'{closureSystName}_shape2Up']   = h_shape2_up
                self.content[histName][group][f'{closureSystName}_shape2Down'] = h_shape2_down

            # Scale all DY histograms based on non closure #
            for systName in self.content[histName][group].keys():
                if closureSystName in systName:
                    continue # avoid applying twice
                dyObj(self.content[histName][group][systName])
                

    def applyRebinning(self):
        logging.info('Applying rebinning')
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

        
    def rebinClassic(self,histName,params):
        """
            Using the classic ROOT rebinning (1D)
            histName: name of histogram in keys of self.content
            params : (list) [groups(int)]
        """
        assert isinstance(params,list)
        assert len(params) == 1
        assert isinstance(params[0],int)
        for group,histDict in self.content[histName].items():
            for systName,hist in histDict.items():
                self.content[histName][group][systName].Rebin(params[0])

    def rebinClassic2D(self,histName,params):
        """
            Using the classic ROOT rebinning (2D)
            histName: name of histogram in keys of self.content
            params : (list) [groups(int),grousp(int)]
        """
        assert isinstance(params,list)
        assert len(params) == 2
        assert isinstance(params[0],int)
        assert isinstance(params[1],int)
        for group,histDict in self.content[histName].items():
            for systName,hist in histDict.items():
                self.content[histName][group][systName].Rebin(params[0],params[1])


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
        for group in self.content[histName].keys():
            for systName,hist in self.content[histName][group].items():
                self.content[histName][group][systName] = boundObj(hist) 

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
        for group in self.content[histName].keys():
            for systName,hist in self.content[histName][group].items():
                self.content[histName][group][systName] = boundObj(hist) 


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
        for group in self.content[histName].keys():
            for systName,hist in self.content[histName][group].items():
                self.content[histName][group][systName] = qObj(hist) 

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
        for group in self.content[histName].keys():
            for systName,hist in self.content[histName][group].items():
                self.content[histName][group][systName] = qObj(hist) 


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
        for group in self.content[histName].keys():
            for systName,hist in self.content[histName][group].items():
                self.content[histName][group][systName] = tObj(hist) 
                
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
        for group in self.content[histName].keys():
            for systName,hist in self.content[histName][group].items():
                self.content[histName][group][systName] = tObj(hist) 

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
        for group in self.content[histName].keys():
            for systName,hist in self.content[histName][group].items():
                self.content[histName][group][systName] = lObj(hist)
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
        for group in self.content[histName].keys():
            for systName,hist in self.content[histName][group].items():
                self.content[histName][group][systName] = lObj(hist)

        # Save linearized info for plotIt #
        if not hasattr(self,'plotLinearizeData'):
            self.plotLinearizeData = {histName:lObj.getPlotData()}
        else:
            self.plotLinearizeData[histName] = lObj.getPlotData()
         


    def preparePlotIt(self,suffix=''):
        logging.info("Preparing plotIt root files")

        # Make new directory with root files for plotIt #
        path_plotIt = os.path.join(self.datacardPath,'plotit')
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

        # Compute stack values #
        inMean = {'backgrounds':[]}
        histMax = {}
        for histName in self.content.keys():
            hstack = ROOT.THStack(histName+"stack",histName+"stack")
            for group in self.content[histName].keys():
                if self.groups[group]['type'] == 'mc':
                    hstack.Add(self.content[histName][group]['nominal'])
                elif self.groups[group]['type'] == 'signal':
                    if group not in inMean.keys():
                        inMean[group] = [self.content[histName][group]['nominal'].Integral()]
                    else:
                        inMean[group].append(self.content[histName][group]['nominal'].Integral())
            if hstack.GetNhists() > 0:
                inMean['backgrounds'].append(hstack.GetStack().Last().Integral())
                histMax[histName] = hstack.GetMaximum()
            else:
                inMean['backgrounds'] = [0.]
                histMax[histName] = 0.
        inMean = {k:sum(values)/len(values) for k,values in inMean.items()} 



        # Files informations #
        config['files'] = {}
        for group,gconfig in self.groups.items():
            if self.pseudodata and group == 'data_real':
                continue
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
                if inMean['backgrounds'] != 0.:
                    config['files'][group+'.root'].update({k:v for k,v in gconfig.items() if k not in ['files','type']})
                    config['files'][group+'.root'].update({'scale':0.01*inMean['backgrounds']/inMean[group]})

        # Groups informations #
        config['groups'] = {}
        for group,gconfig in self.groups.items():
            if self.pseudodata and group == 'data_real':
                continue
            if gconfig['type'] != 'signal':
                plotGroup = group if 'group' not in gconfig.keys() else gconfig['group']
                if plotGroup not in config['groups'].keys():
                    config['groups'][plotGroup] = {k:v for k,v in gconfig.items() if k not in ['files','type']}


        # Plots informations #
        config['plots'] = {}
        for h1,h2 in self.hist_conv.items():
            # Check histogram dimension #
            hist = self.content[h1][list(self.groups.keys())[0]]['nominal']
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
                config['plots'][h1]['blinded-range'] = [xmin,xmax] #[(xmax-xmin)/2,xmax]
            else:
                config['plots'][h1].pop('blinded-range',None)
            if histMax[h1] != 0.:
                config['plots'][h1]['y-axis-range'] = [0.,histMax[h1]*1.5]
                config['plots'][h1]['log-y-axis-range'] = [1.e-1,histMax[h1]*100]
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

            if hasattr(self,'plotLinearizeData'):
                if h1 in self.plotLinearizeData.keys():
                    margin_left = config['configuration']['margin-left']
                    margin_right = 1-config['configuration']['margin-right']
                    extraItems = self.plotLinearizeData[h1]
                    for idx in range(len(extraItems['labels'])):
                        extraItems['labels'][idx]['position'][0] = margin_left + extraItems['labels'][idx]['position'][0] * (margin_right-margin_left)
                    for idx in range(len(extraItems['lines'])):
                        extraItems['lines'][idx] = [[extraItems['lines'][idx],0.],[extraItems['lines'][idx],histMax[h1]*1.1]]
                    config['plots'][h1].update(extraItems)
                    #config['plots'][h1]['x-axis-hide-ticks'] = True
                    #config['plots'][h1]['blinded-range'] = [[0.5,1.],[1.5,2.]]

        if len(systematics) > 0:
            config['systematics'] = systematics

        # Write yaml file #
        path_yaml = os.path.join(path_plotIt,'plots.yml')
        with open(path_yaml,'w') as handle:
            yaml.dump(config,handle)
        logging.info("New yaml file for plotIt : %s"%path_yaml)

        # PlotIt command #
        path_pdf = os.path.join(path_plotIt,'plots_%s'%self.era)
        if not os.path.exists(path_pdf):
            os.makedirs(path_pdf)
        logging.info("plotIt command :")
        cmd = "plotIt -i {input} -o {output} -e {era} {yaml}".format(**{'input': path_plotIt,
                                                                        'output': path_pdf,   
                                                                        'era': self.era,
                                                                        'yaml': path_yaml})
        print (cmd)
        if self.which('plotIt') is not None:
            logging.info("Calling plotIt")
            rc = self.run_command(shlex.split(cmd))
            if rc == 0:
                logging.info('... success')
            else:
                logging.info('... failure')
        else:
            logging.warning("plotIt not found")

    def run_combine(self,entries):
        logging.info(f'Running combine on {self.datacardName}')
        if self.combineConfigs is None:
            logging.warning("No combine command in the yaml config file")
            return
        txtPaths = self.getTxtFilesPath()
        missingTxts = []
        for txtPath in txtPaths.values():
            if not os.path.exists(txtPath):
                missingTxts.append(txtPath)
        if len(missingTxts) > 0:
            logging.debug('Missing following datacards :')
            for missingTxt in missingTxts:
                logging.debug(f'... {missingTxt}')
            logging.info('Datacard not found, will produce it')
            self.run_production()
        
        if entries is None or (isinstance(entries,list) and len(entries)==0):
            entries = self.combineConfigs.keys()

        combine_dir, env_script = os.path.split(combine_env)
        combineModes = ['limits','prefit','postfit_b','postfit_s']
        for entry in entries:
            if entry not in self.combineConfigs.keys():
                raise RuntimeError(f'Combine entry {entry} not in combine config in the yaml file')
            combineCfg = self.combineConfigs[entry]
            combineMode = combineCfg['mode']

            # Make subdirectory #
            logging.info(f'Running mode {combineMode}')
            if combineMode not in combineModes:
                raise RuntimeError(f'Combine mode {combineMode} not understood')
            subdir = os.path.join(self.datacardPath,entry)
            if os.path.exists(subdir) and len(os.listdir(subdir)) != 0:
                logging.warning(f'Subdirectory {subdir} is not empty')
            if not os.path.exists(subdir):
                os.makedirs(subdir)



            # Check if root output is already in subdir #
            combinedTxtPath = os.path.join(subdir,'datacard.txt')
            workspacePath = os.path.join(subdir,'workspace.root')
            rootFiles = glob.glob(os.path.join(subdir,'*.root'))
            if len(rootFiles) == 0 or rootFiles == [workspacePath]:
                logging.info(f'No root file found in {subdir}, will run the command for mode {combineMode}')

                if 'command' not in combineCfg.keys():
                    logging.warning('No "command" entry for combine mode "{combineMode}"')

                # Select txt files #
                if 'bins' in combineCfg.keys():
                    if list(txtPaths.keys()) == ['all']:
                        raise RuntimeError(f'You asked for a single datacard text file, but filter on bins for the command mode {combineMode}')
                    txtPathsForCombination = []
                    for combBin in combineCfg['bins']:
                        if combBin not in txtPaths.keys():
                            raise RuntimeError(f'{combBin} not in hist content')
                        txtPathsForCombination.append(txtPaths[combBin])
                else:
                    txtPathsForCombination = list(txtPaths.values())

                # Submit mode #
                if 'submit' in combineCfg.keys() and not self.worker:
                    params = combineCfg['submit']
                    args = {'combine':entry}
                    subScript = self.writeSbatchCommand(subdir,params,args)
                    logging.info(f'Generated {subScript}')
                    rc,output = self.run_command(f'sbatch {subScript}',return_output=True,shell=True)
                    for line in output:
                        logging.info(line)
                    continue

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
                    workspaceCmd += f"env -i bash -c 'source {env_script} && cd {subdir} && text2workspace.py {combinedTxtPath} -o {workspacePath}'"
                    exitCode,output = self.run_command(workspaceCmd,return_output=True,shell=True)
                    if exitCode != 0:
                        for line in output:
                            logging.info(line)
                        raise RuntimeError('Could not produce the workspace, see log above')
                    logging.info('... done')
                
                # Produce command #
                combineCmd = combineCfg['command'].format(workspacePath)
#            if combineMode == 'fitdiag':
#                # Add specific string to output file to avoid overwriting
#                testName = 'Test'
#                seed = '123456'
#                mass = '120'
#                def findArg(cmd,arg):
#                    cmd_split = cmd.split(' ')
#                    return cmd_split[cmd_split.index(arg)+1]
#                if '-n ' in combineCmd:
#                    testName = findArg(combineCmd,'-n')
#                if '--name ' in combineCmd:
#                    testName = findArg(combineCmd,'--name')
#                if '-s ' in combineCmd:
#                    seed = findArg(combineCmd,'-s')
#                if '--seed ' in combineCmd:
#                    seed = findArg(combineCmd,'--seed')
#                if '-m ' in combineCmd:
#                    mass = findArg(combineCmd,'-m')
#                if '--mass ' in combineCmd:
#                    mass = findArg(combineCmd,'--mass')
#                if '--keyword-value' in combineCmd:
#                    raise RuntimeError('The combine argument "--keyword-value" is not supported in this script')
#
#                random_str = ''.join(random.choice(string.ascii_uppercase) for i in range(6))
#                combineCmd += f' --keyword-value WORD={random_str}' 
#                outputFile = f'higgsCombine{testName}.FitDiagnostics.mH{mass}.WORD{random_str}.{seed}.root'

                fullCombineCmd  = f"cd {combine_dir}; "
                fullCombineCmd += f"env -i bash -c 'source {env_script} && cd {subdir} && {combineCmd}'"

                logging.debug('Extensive command is below')
                logging.debug(fullCombineCmd)

                # Run command #
                exitCode,output = self.run_command(fullCombineCmd,return_output=True,shell=True)
                path_log = os.path.join(subdir,'log_combine.out')
                with open(path_log,'w') as handle:
                    for line in output:
                        handle.write(line)

                if exitCode != 0:
                    logging.critical(f'Something went wrong in combine, see log in {path_log}') 
                    continue
                else:
                    logging.info('... done')
            else:
                logging.warning(f"Already a root file in subdirectory {subdir}, will not run the combine command again (remove and run again if needed)")

            # Saving limit values in json #
            if combineMode == 'limits':
                limits = {}
                with open(os.path.join(subdir,'log_combine.out'),'r') as output:
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
                        
                path_limits = os.path.join(self.datacardPath,subdir,'limits.json')
                with open(path_limits,'w') as handle:
                    json.dump(limits,handle)
                logging.info(f'Saved limits as {path_limits}') 

            # Producing prefit and postfit plots #
            if 'prefit' in combineMode or 'postfit' in combineMode:
                config = combineCfg['plotting']
                fitdiagFile = glob.glob(os.path.join(subdir,'fitDiagnostic*root'))
                if len(fitdiagFile) == 0:
                    raise RuntimeError("Could not find any fitdiag file in subdir")
                fitdiagFile = fitdiagFile[0]
                # Generate dat config file #
                config['main'].update({'fitdiagnosis':fitdiagFile})
                path_dict = os.path.join(subdir,'dict.dat')
                self.makeDictForPostFit(path_dict,config)
                # Run plotting #
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
                    create_postfit_plots(path                    = subdir,
                                         fit_diagnostics_path    = fitdiagFile,
                                         normalize_X_original    = False,
                                         doPostFit               = 'postfit' in combineMode,
                                         divideByBinWidth        = False,
                                         folder                  = folder,
                                         bin                     = binCfg,
                                         binToRead               = key,
                                         unblind                 = unblind)
                                    

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
        combine_dir, env_script = os.path.split(combine_env)
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
        random_str = ''.join(random.choice(string.ascii_lowercase) for i in range(6))
        filepath = os.path.join(subdir,f'sbatch_{random_str}.sh')
        with open(filepath,'w') as handle:
            # Write header #
            handle.write('#!/bin/bash\n\n')
            handle.write('#SBATCH --job-name=datacard_prod\n')
            handle.write(f'#SBATCH --output={os.path.join(subdir,"log_"+random_str+"-%j.out")}\n')
            handle.write(f'#SBATCH --error={os.path.join(subdir,"log_"+random_str+"-%j.err")}\n')
            for key,val in params.items():
                handle.write(f'#SBATCH --{key}={val}\n\n')
            # Write command #
            handle.write(f'python3 produceDataCards.py --yaml {self.configPath} --worker ')
            for key,val in args.items():
                handle.write(f'--{key} {val} ')
            handle.write("\n")

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
                        help='Whether to use pseudo data (data = sum of MC)')
    parser.add_argument('--combine', action='store', required=False, default=None, nargs='*',
                        help='Run combine on the txt datacard only')
    parser.add_argument('-v','--verbose', action='store_true', required=False, default=False,
                        help='Verbose mode')
    parser.add_argument('--worker', action='store_true', required=False, default=False,
                        help='Force working locally without submitting jobs')
    args = parser.parse_args()

    logging.basicConfig(level=logging.DEBUG,format='%(asctime)s - %(levelname)s - %(message)s',datefmt='%m/%d/%Y %H:%M:%S')
    if not args.verbose:
        logging.getLogger().setLevel(logging.INFO)

    if args.yaml is None:
        raise RuntimeError("Must provide the YAML file")
    with open(args.yaml,'r') as handle:
        f = yaml.load(handle,Loader=YMLIncludeLoader)
    if not os.path.isfile(args.yaml):
        raise RuntimeError("YAML file {} is not a valid file".format(args.yaml))
    instance = DataCard(datacardName    = os.path.basename(args.yaml).replace('.yml',''),
                        configPath      = os.path.abspath(args.yaml),
                        worker          = args.worker,
                        pseudodata      = args.pseudodata,
                        **f)
    if args.combine is not None: 
        instance.run_combine(args.combine)
    else:
        instance.run_production()
        #instance.run_combine(args.combine)
