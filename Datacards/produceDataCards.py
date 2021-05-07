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
#ROOT.gROOT.SetBatch(True) # TODO: enable

class DataCard:
    def __init__(self,datacardName=None,path=None,yamlName=None,groups=None,shapeSyst=None,normSyst=None,hist_conv=None,era=None,use_syst=False,root_subdir=None,DYnonclosure=None,pseudodata=False,rebin=None,textfiles=None,produce_plots=False,legend=None,**kwargs):
        self.datacardName   = datacardName
        self.path           = path
        self.hist_conv      = hist_conv
        self.era            = str(era)
        self.use_syst       = use_syst
        self.shapeSyst      = shapeSyst
        self.normSyst       = normSyst
        self.root_subdir    = root_subdir
        self.pseudodata     = pseudodata
        self.DYnonclosure   = DYnonclosure
        self.rebin          = rebin
        self.textfiles      = textfiles
        self.produce_plots  = produce_plots
        self.legend         = legend

        if isinstance(groups,dict):
            self.groups = groups
        elif isinstance(groups,list) or isinstance(groups,tuple):
            assert all([isinstance(subgroups,dict) for subgroups in groups])
            self.groups = {key:val for subgroups in groups for key,val in subgroups.items()}
        else: 
            raise RuntimeError("Group format not understood")
        
        self.yaml_dict = self.loadYaml(self.path,yamlName)

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
        from txtwriter import Writer
        if self.datacardName is None:
            raise RuntimeError("Datacard name is not set")
        path_datacard = os.path.join(os.path.abspath(os.path.dirname(__file__)),'datacards',self.datacardName)

        shapes = {histName:os.path.join(path_datacard,f"{histName}_{self.era}.root") for histName in self.content.keys()}

        # Save root file #
        if not os.path.exists(path_datacard):
            os.makedirs(path_datacard)
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
                            raise RuntimeError("Could not find {} in systematic dict".format(systName))
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
            print (f"Saved file {shapes[histName]}")

        # Save txt file #
        if self.textfiles is None:
            self.textfiles = '{}.txt'
        if '{}' not in self.textfiles:
            print ('Will create a single datacard txt file for all the categories')
            writer = Writer([f'{histName}_{self.era}' for histName in self.content.keys()])
        else:
            print ('Will create one datacard txt file per categories')
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
                textPath = os.path.join(path_datacard,self.textfiles.format(f"{histName}_{self.era}"))
                writer.dump(textPath,os.path.basename(shapes[histName]))
                print (f"Saved file {textPath}")
        if '{}' not in self.textfiles:
            textPath = os.path.join(path_datacard,self.textfiles)
            writer.dump(textPath,[os.path.basename(shape) for shape in shapes.values()])
            print (f"Saved file {textPath}")

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


    def applyDYNonClosure(self):
        print ('Applying DY non closure shift')
        from DYNonClosure import NonClosureSystematic
        dyObj = NonClosureSystematic.load_from_json(self.DYnonclosure['path'])
        for histName in self.content.keys():
            if histName not in self.DYnonclosure['hist'].keys():
                continue
            catName = self.DYnonclosure['hist'][histName]['entry']
            closureSystName = self.DYnonclosure['hist'][histName]['name']
            dyObj.SetCategory(catName)
            print ("... Histogram {:20s} -> applying non closure with key {:20s}".format(histName,catName))
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
        print ('Applying rebinning')
        for histName in self.content.keys():
            if histName in self.rebin.keys():
                rebinSchemes = self.rebin[histName]
                if not isinstance(rebinSchemes,list):
                    rebinSchemes = [rebinSchemes]
                for rebinScheme in rebinSchemes:
                    method = rebinScheme['method']
                    params = rebinScheme['params']
                    print (f'... {histName:20s} -> rebinning scheme : {method}')
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
                    # Error #
                    else:
                        raise RuntimeError("Could not understand rebinning method {} for histogram {}".format(method,histName))
                    rebinFunc(histName,params)
            else:
                print ('... {:20s} -> rebinning scheme not requested'.format(histName))
        
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
        boundObj = Boundary(params[0])
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
        if not set(groups_for_binning).issubset(set(self.groups.keys())):
            raise RuntimeError('Groups {'+','.join([g for g in groups_for_binning if g not in self.groups.keys()])+'} for quantile are not in group dict')

        hists_for_binning = [histDict['nominal'] for group,histDict in self.content[histName].items() if group in groups_for_binning if histDict['nominal'] is not None]
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
        if not set(groups_for_binning).issubset(set(self.groups.keys())):
            raise RuntimeError('Groups {'+','.join([g for g in groups_for_binning if g not in self.groups.keys()])+'} for quantile are not in group dict')

        hists_for_binning = [histDict['nominal'] for group,histDict in self.content[histName].items() if group in groups_for_binning if histDict['nominal'] is not None]
        qObj = Quantile2D(hists_for_binning,qx,qy)
        for group in self.content[histName].keys():
            for systName,hist in self.content[histName][group].items():
                self.content[histName][group][systName] = qObj(hist) 


    def rebinThreshold(self,histName,params):
        """
            Using the threshold rebinning (1D)
            histName: name of histogram in keys of self.content
            params : (list) [thresholds (list), backgrounds (list), extra contributions (list), rsut (float)]
        """
        from Rebinning import Threshold
        assert isinstance(params,list)
        assert len(params) == 3
        assert isinstance(params[0],list)
        assert isinstance(params[1],list)
        assert isinstance(params[2],list)
        assert isinstance(params[3],float)
        thresholds  = params[0]
        backgrounds = params[1]
        extra       = params[2]
        rsut        = params[3] # relative stat. unc. threshold
        if not set(backgrounds).issubset(set(self.groups.keys())):
            raise RuntimeError('Groups {'+','.join([g for g in backgrounds if g not in self.groups.keys()])+'} for threshold are not in group dict')
        if not set(extra).issubset(set(self.groups.keys())):
            raise RuntimeError('Groups {'+','.join([g for g in extra if g not in self.groups.keys()])+'} for extra contributions are not in group dict')
        hists_for_binning = [histDict['nominal'] for group,histDict in self.content[histName].items() if group in backgrounds]
        hists_extra       = [histDict['nominal'] for group,histDict in self.content[histName].items() if group in extra]
        tObj = Threshold(hists_for_binning,thresholds,True,hists_extra,rsut)
        for group in self.content[histName].keys():
            for systName,hist in self.content[histName][group].items():
                self.content[histName][group][systName] = tObj(hist) 
                
    def rebinThreshold2D(self,histName,params):
        """
            Using the threshold rebinning (2D)
            histName: name of histogram in keys of self.content
            params : (list) [thresholds (list), thresholds (list), backgrounds (list), extra contributions (list), rsut (float)]
        """
        from Rebinning import Threshold2D
        assert isinstance(params,list)
        assert len(params) == 4
        assert isinstance(params[0],list)
        assert isinstance(params[1],list)
        assert isinstance(params[2],list)
        assert isinstance(params[3],list)
        assert isinstance(params[4],float)
        threshx     = params[0]
        threshy     = params[1]
        backgrounds = params[2]
        extra       = params[3]
        rsut        = params[4] # relative stat. unc. threshold
        if not set(backgrounds).issubset(set(self.groups.keys())):
            raise RuntimeError('Groups {'+','.join([g for g in backgrounds if g not in self.groups.keys()])+'} for threshold are not in group dict')
        if not set(extra).issubset(set(self.groups.keys())):
            raise RuntimeError('Groups {'+','.join([g for g in extra if g not in self.groups.keys()])+'} for extra contributions are not in group dict')
        hists_for_binning = [histDict['nominal'] for group,histDict in self.content[histName].items() if group in backgrounds]
        hists_extra       = [histDict['nominal'] for group,histDict in self.content[histName].items() if group in extra]
        tObj = Threshold(hists_for_binning,threshx,threshy,True,hists_extra,rsut)
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
         

    def preparePlotIt(self,suffix=''):
        print ("Preparing plotIt root files")
        # Make new directory with root files for plotIt #
        path_plotIt = os.path.join(os.path.abspath(os.path.dirname(__file__)),'datacards',self.datacardName,'plotit')
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
        config['groups'] = {}
        for group,gconfig in self.groups.items():
            if self.pseudodata and group == 'data_real':
                continue
            if gconfig['type'] != 'signal':
                plotGroup = group if 'group' not in gconfig.keys() else gconfig['group']
                if plotGroup not in config['groups'].keys():
                    config['groups'][plotGroup] = {k:v for k,v in gconfig.items() if k not in ['files','type']}


        config['plots'] = {}
        for h1,h2 in self.hist_conv.items():
            # Check histogram dimension #
            hist = self.content[h1][list(self.groups.keys())[0]]['nominal']
            if isinstance(hist,ROOT.TH2F) or isinstance(hist,ROOT.TH2D):
                print(f'Histogram {h1} is a TH2 and can therefore not be used in plotIt')
                continue
            print (hist.GetXaxis().GetNdivisions())
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
                config['plots'][h1]['blinded-range'] = [(xmax-xmin)/2,xmax]
            else:
                config['plots'][h1].pop('blinded-range',None)
            config['plots'][h1]['ratio-y-axis-range'] = [0.5,1.5]
            config['plots'][h1]['sort-by-yields'] = True
            config['plots'][h1]['show-overflow'] = True
            if 'labels' in config['plots'][h1].keys():
                del config['plots'][h1]['labels']

            if self.legend is not None and 'position' in self.legend.keys():
                assert len(self.legend['position']) == 4
                config['plots'][h1]['legend-position'] = self.legend['position']

            if hasattr(self,'plotLinearizeData'):
                if h1 in self.plotLinearizeData.keys():
                    margin_left = config['configuration']['margin-left']
                    margin_right = 1-config['configuration']['margin-right']
                    extraItems = self.plotLinearizeData[h1]
                    for idx in range(len(extraItems['labels'])):
                        extraItems['labels'][idx]['position'][0] = margin_left + extraItems['labels'][idx]['position'][0] * (margin_right-margin_left)
                    config['plots'][h1].update(extraItems)
                    config['plots'][h1]['x-axis-hide-ticks'] = True

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
        if self.which('plotIt') is not None:
            print ("Calling plotIt")
            rc = self.run_command(cmd.split(" "),print_output=True)
            if rc == 0:
                print ('... success')
            else:
                print ('... failure')
        else:
            print ("plotIt not found")

    @staticmethod
    def run_command(command,return_output=False,print_output=False):
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
        outputs = []
        while True:
            out = process.stdout.readline().strip().decode()
            if len(out) == 0 and process.poll() is not None:
                break
            if out and print_output:
                print (out)
            if return_output:
                outputs.append(out)
        rc = process.poll()
        if return_output:
            return rc,outputs
        else:
            return rc

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



 

class YMLIncludeLoader(yaml.SafeLoader): 
    """Custom yaml loading to support including config files. Use `!include (file)` to insert content of `file` at that position."""
    
    def __init__(self, stream):
        super(YMLIncludeLoader, self).__init__(stream)
        self._root = os.path.split(stream.name)[0]

    def include(self, node):
        filename = os.path.join(self._root, self.construct_scalar(node))
        with open(filename, 'r') as f:
            return yaml.load(f, YMLIncludeLoader)

YMLIncludeLoader.add_constructor('!include', YMLIncludeLoader.include)


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
        f = yaml.load(handle,Loader=YMLIncludeLoader)
    instance = DataCard(datacardName    = os.path.basename(args.yaml).replace('.yml',''),
                        pseudodata      = args.pseudodata,
                        **f)
