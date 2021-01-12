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
import numpy as np
import enlighten
import ROOT

ROOT.gROOT.SetBatch(True)

sys.path.append(os.path.abspath('../DYStudy'))
from CDFShift import CDFShift


class DataCard:
    def __init__(self,datacardName=None,path=None,yamlName=None,groups=None,hist_conv=None,era=None,use_syst=False,root_subdir=None,pseudodata=False,quantiles=None,DYEstimation=None,produce_plots=False,**kwargs):
        self.datacardName   = datacardName
        self.path           = path
        self.groups         = groups
        self.hist_conv      = hist_conv
        self.era            = str(era)
        self.use_syst       = use_syst
        self.root_subdir    = root_subdir
        self.pseudodata     = pseudodata
        self.quantiles      = quantiles
        self.DYEstimation   = DYEstimation
        self.produce_plots  = produce_plots
        
        self.yaml_dict = self.loadYaml(os.path.join(self.path,yamlName))

        if self.pseudodata:
            print ('Will use pseudodata')
            self.groups = self.generatePseudoData(self.groups)
            if self.datacardName is not None:
                self.datacardName += "_pseudodata"
        self.content = {k:{g:None for g in self.groups.keys()} for k in self.hist_conv.keys()}

        self.loopOverFiles()
#        if self.pseudodata:
#            self.roundFakeData()
        if self.DYEstimation is not None:
            self.produceDYShit()
        if self.quantiles is not None:
            self.rebinInQuantile()
            #self.rebinClassic()
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
        files = glob.glob(os.path.join(self.path,'results','*.root'))
        pbar = enlighten.Counter(total=len(files), desc='Progress', unit='files')
        for f in files:
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
            for group in groups:
                self.addSampleToGroup(copy.deepcopy(hist_dict),group)
                # deepcopy is needed so that if the histogram is in two groups
                # acting on one version will not change the other

        for histName in self.content.keys():
            for group,hist in self.content[histName].items():
                if self.pseudodata and group == 'data_real':
                    continue
                hist.SetName(hist.GetName()+'_'+group)

    def addSampleToGroup(self,hist_dict,group):
        for histname,hists in hist_dict.items():
            nominal = hists['nominal']
            if self.content[histname][group] is None:
                self.content[histname][group] = copy.deepcopy(nominal)
            else:
                self.content[histname][group].Add(nominal)
            if self.use_syst:
                for systName in hists.keys():
                    if systName == 'nominal':
                        continue
                    hist = hists[systName]
                    groupsyst = group + '_' + systName
                    if groupsyst not in self.content[histname].keys():
                        self.content[histname][groupsyst] = copy.deepcopy(hist)
                    else:
                        self.content[histname][groupsyst].Add(hist)


    def findGroup(self,sample):
        group_of_sample = []
        for group in self.groups.keys():
            files = self.groups[group]['files']
            if not isinstance(files,list):
                raise RuntimeError("Group %s does not consist in a list"%key)
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
            hist_dict[datacardname] = {}
            if not isinstance(histnames,list):
                histnames = [histnames]
            for histname in histnames:
                # Check #
                if not histname in list_histnames:
                    print ("Could not find hist %s in %s"%(histname,rootfile))
                    continue
                listsyst = [hn for hn in list_histnames if histname in hn and '__' in hn] if self.use_syst else []
                
                # Nominal histogram #
                h = self.getHistogram(f,histname,lumi,br,xsec,sumweight)
                if not 'nominal' in hist_dict[datacardname].keys():
                    hist_dict[datacardname]['nominal'] = copy.deepcopy(h)
                else:
                    hist_dict[datacardname]['nominal'].Add(h)
                # Systematic histograms #
                for syst in listsyst:
                    h = self.getHistogram(f,syst,lumi,br,xsec,sumweight)
                    systName = syst.split('__')[-1].replace('up','Up').replace('down','Down')
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

        return {'luminosity':lumi_dict,'samples':sample_dict,'plots':full_dict['plots']}

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
        for histName, gdict in self.content.items():
            filename = os.path.join(path_datacard,histName+'_'+self.era+'.root')

            f = ROOT.TFile(filename,'recreate')
            if self.root_subdir is not None:
                d = f.mkdir(self.root_subdir,self.root_subdir)
                d.cd()
            for group,hist in gdict.items():
                if self.pseudodata and group == 'data_real':
                    continue
                if hist is None:
                    continue
                for i in range(1,hist.GetNbinsX()+1):
                    #hist.SetBinError(i,0.) #TODO : check 
                    if hist.GetBinContent(i) < 0.:
                        hist.SetBinContent(i,0.) 
                hist.SetTitle(group)
                hist.SetName(group)
                hist.Write(group)
            f.Write()
            f.Close()
            print ("Saved file %s"%filename)

    def rebinInQuantile(self,qObjects=None):
        from QuantileBinning import Quantile
        for histName,gDict in self.content.items():
            nodes = [node for node in self.quantiles.keys() if node in histName]
            if len(nodes) != 1:
                raise RuntimeError("For histogram {} : found matches in quantile dict : "+','.join(nodes))
            else:
                node = nodes[0]
            groups_for_binning = self.quantiles[node]['groups']
            if not set(groups_for_binning).issubset(set(self.groups.keys())):
                raise RuntimeError('Groups {'+','.join([g for g in groups_for_binning if g not in self.groups.keys()])+'} for quantile are not in group dict')
            hists_for_binning = [hist for group,hist in gDict.items() if group in groups_for_binning if hist is not None]
            quantiles = np.array(self.quantiles[node]['quantile'])
            qObj = Quantile(hists_for_binning,quantiles)
            for group,hist in gDict.items():
                if hist is None:
                    continue
                self.content[histName][group] = qObj(hist)

    def rebinClassic(self):
        for histName,gDict in self.content.items():
            print (histName)
            for group in gDict.keys():
                if self.pseudodata and group=='data_real':
                    continue
                self.content[histName][group].Rebin(8)

    def produceDYShit(self):
        print ('Producing DY estimation')
        with open(self.DYEstimation['yaml'],'r') as handle:
            f = yaml.load(handle,Loader=yaml.FullLoader)

        if not self.pseudodata:
            with open('DYFits/ZPeak_2016.yml','r') as handle:
                norm_ZPeak = yaml.load(handle,Loader=yaml.FullLoader)
            with open('DYFits/ZVeto_2016.yml','r') as handle:
                norm_ZVeto = yaml.load(handle,Loader=yaml.FullLoader)

        f['quantiles'] = None
        f['produce_plots'] = False
        ZPeak_datacard = DataCard(**f,pseudodata=self.pseudodata)

        path_DY = os.path.join(os.path.abspath(os.path.dirname(__file__)),self.datacardName,'DYEstimation')
        if not os.path.exists(path_DY):
            os.makedirs(path_DY)

        for histName,gDict in self.content.items():
            print ("---------------------------")
            print (histName)
            if not histName in self.DYEstimation['shift'].keys():
                continue
            if histName not in self.DYEstimation['shift'].keys():
                continue
            proxyName = self.DYEstimation['shift'][histName]
            if proxyName not in ZPeak_datacard.content.keys():
                raise RuntimeError("Could not find histogram %s in Z Peak datacard"%proxyName)

            if not self.pseudodata:
                proxyPeak_norm = norm_ZPeak[proxyName]
                proxyVeto_norm = norm_ZVeto[proxyName]
                histPeak_norm  = norm_ZPeak[histName]
                #histVeto_norm  = norm_ZVeto[histName]

            histZVeto = None 
            histZPeak = None 
            proxyZVeto = None 
            proxyZPeak = None 
            if not set(self.DYEstimation['factors'].keys()).issubset(set(self.groups.keys())):
                raise RuntimeError('Group factors in DYestimation dict {'+','.join([g for g in self.DYEstimation['factors'].keys() if g not in self.groups.keys()])+'} are not in the groups dict')

            for group in gDict.keys():
                if group not in self.DYEstimation['factors'].keys():
                    continue
                factor = self.DYEstimation['factors'][group]
                hv = copy.deepcopy(self.content[histName][group])
                hz = copy.deepcopy(ZPeak_datacard.content[histName][group])
                pv = copy.deepcopy(self.content[proxyName][group])
                pz = copy.deepcopy(ZPeak_datacard.content[proxyName][group])
                
                if not self.pseudodata:
                    factor_hz = None
                    factor_pv = None
                    factor_pz = None
                    print (group,proxyName,histName)
                    #if group in histVeto_norm.keys():
                    #    factor_hv = histVeto_norm[group][1]/histVeto_norm[group][0]
                    if group in histPeak_norm.keys():
                        factor_hz = histPeak_norm[group][1]/histPeak_norm[group][0]
                    else:
                        factor_hz = 1.
                    if group in proxyVeto_norm.keys():
                        factor_pv = proxyVeto_norm[group][1]/proxyVeto_norm[group][0]
                    else:
                        factor_pv = 1. 
                    if group in proxyPeak_norm.keys():
                        factor_pz = proxyPeak_norm[group][1]/proxyPeak_norm[group][0]
                    else:
                        factor_pz = 1.
                    print (factor_hz,factor_pv,factor_pz)
                    hz.Scale(factor_hz)
                    pv.Scale(factor_pv)
                    pz.Scale(factor_pz)

                hv.Scale(factor)
                hz.Scale(factor)
                pv.Scale(factor)
                pz.Scale(factor)

                if histZVeto is None:
                    histZVeto = hv
                    histZPeak = hz
                    proxyZVeto = pv
                    proxyZPeak = pz
                else:
                    histZVeto.Add(hv)
                    histZPeak.Add(hz)
                    proxyZVeto.Add(pv)
                    proxyZPeak.Add(pz)

            print (histZPeak.Integral(),proxyZVeto.Integral(),proxyZPeak.Integral())
            if histZPeak.Integral() <= 0. or proxyZVeto.Integral() <= 0. or proxyZPeak.Integral() <= 0.:
                continue

            if self.pseudodata:
                norm = histZVeto.Integral() * histZPeak.Integral() / proxyZPeak.Integral()
            else:
                norm = proxyZVeto.Integral() * histZPeak.Integral() / proxyZPeak.Integral()


            shiftObj = CDFShift(proxyZPeak,histZPeak,proxyZVeto,norm=norm)
            sHist, sHistUp, sHistDown = shiftObj.returnHistAndSyst()

            shiftObj.rootsave['histZVeto'] = copy.deepcopy(histZVeto)
            shiftObj.rootsave['histZPeak'] = copy.deepcopy(histZPeak)
            shiftObj.rootsave['proxyZVeto'] = copy.deepcopy(proxyZVeto)
            shiftObj.rootsave['proxyZPeak'] = copy.deepcopy(proxyZPeak)
            shiftObj.saveToRoot(os.path.join(path_DY,histName+'.root'))
            #shiftObj.saveToPDF(histName+'.pdf')
            sHist.SetName(self.content[histName][self.DYEstimation['replace']].GetName())
            sHistUp.SetName(self.content[histName][self.DYEstimation['replace']].GetName()+'__DYsystUp')
            sHistDown.SetName(self.content[histName][self.DYEstimation['replace']].GetName()+'__DYsystDown')
            self.content[histName][self.DYEstimation['replace']] = sHist
            self.content[histName][self.DYEstimation['replace']+'__DYsystUp'] = sHistUp
            self.content[histName][self.DYEstimation['replace']+'__DYsystDown'] = sHistDown


            pdfname = os.path.join(path_DY,histName+'.pdf')
            C = ROOT.TCanvas("C","C",800,600)
            #C.Print(pdfname+'[')

            
            histZVeto.Draw("AP")
            histZPeak.Draw("same")
            proxyZVeto.Draw("same")
            proxyZPeak.Draw("same")

            #C.cd(0)
            C.Print(pdfname,'Title:Test')
#            C.Clear()
#            
#            legH = ROOT.TLegend(0.7,0.7,0.9,0.9)
#            legH.AddEntry(shiftObj.rootsave['h1'],'h1')
#            legH.AddEntry(shiftObj.rootsave['h2'],'h2')
#            legH.AddEntry(shiftObj.rootsave['h3'],'h3')
#            legH.AddEntry(shiftObj.rootsave['nh3'],'nh3')                                                                                                                                                      
#            shiftObj.rootsave['h1'].SetLineColor(601)
#            shiftObj.rootsave['h2'].SetLineColor(418)
#            shiftObj.rootsave['h3'].SetLineColor(634)
#            shiftObj.rootsave['nh3'].SetLineColor(619)
#            
#            shiftObj.rootsave['h1'].DrawClone()
#            shiftObj.rootsave['h2'].DrawClone()
#            shiftObj.rootsave['h3'].DrawClone()
#            shiftObj.rootsave['nh3'].DrawClone()
#            legH.Draw()
#            
#            C.Print(pdfname)
#            C.Clear()
#            
#            legC = ROOT.TLegend(0.7,0.7,0.9,0.9)
#            legC.AddEntry(shiftObj.rootsave['cdf1'],'cdf1')
#            legC.AddEntry(shiftObj.rootsave['cdf2'],'cdf2')
#            legC.AddEntry(shiftObj.rootsave['cdf3'],'cdf3')
#            legC.AddEntry(shiftObj.rootsave['ncdf3'],'ncdf3')
#            shiftObj.rootsave['cdf1'].SetLineColor(601)
#            shiftObj.rootsave['cdf2'].SetLineColor(418)
#            shiftObj.rootsave['cdf3'].SetLineColor(634)
#            shiftObj.rootsave['ncdf3'].SetLineColor(619)
#            
#            shiftObj.rootsave['cdf1'].Draw("H")
#            shiftObj.rootsave['cdf2'].Draw("same H")
#            shiftObj.rootsave['cdf3'].Draw("same H")
#            shiftObj.rootsave['ncdf3'].Draw("same H")
#            legC.Draw()
#            
#            C.Print(pdfname)
#            C.Clear()
#
#            legS = ROOT.TLegend(0.7,0.7,0.9,0.9)
#            legS.AddEntry(shiftObj.rootsave['shift'],'shift (nominal)')
#            legS.AddEntry(shiftObj.rootsave['shift_up'],'shift (up)')
#            legS.AddEntry(shiftObj.rootsave['shift'],'shift (down)')
#            shiftObj.rootsave['shift'].SetLineStyle(1)
#            shiftObj.rootsave['shift'].SetLineWidth(2)
#            shiftObj.rootsave['shift_up'].SetLineStyle(10)
#            shiftObj.rootsave['shift_down'].SetLineStyle(5)
#            legS.Draw()
#     
#            shiftObj.rootsave['shift'].Draw()
#            shiftObj.rootsave['shift_up'].Draw("same")
#            shiftObj.rootsave['shift_down'].Draw("same")
#     
#            C.Print(pdfname)
#            C.Clear()
#     
#            legE = ROOT.TLegend(0.7,0.7,0.9,0.9)
#            legE.AddEntry(shiftObj.rootsave['cdf_tot_err'],'Error (tot)')
#            legE.AddEntry(shiftObj.rootsave['cdf_bin_err'],'Error (bin)')
#            legE.AddEntry(shiftObj.rootsave['cdf_shift_err'],'Error (shift)')
#            shiftObj.rootsave['cdf_tot_err'].SetLineStyle(1)
#            shiftObj.rootsave['cdf_tot_err'].SetLineWidth(2)
#            shiftObj.rootsave['cdf_bin_err'].SetLineStyle(10)
#            shiftObj.rootsave['cdf_shift_err'].SetLineStyle(5)
#            legE.Draw()
#     
#            shiftObj.rootsave['cdf_tot_err'].Draw()
#            shiftObj.rootsave['cdf_bin_err'].Draw("same")
#            shiftObj.rootsave['cdf_shift_err'].Draw("same")
#     
#            C.Print(pdfname)
#            C.Clear()
#     
#            C.Print(pdfname+']')


 
    def preparePlotIt(self,suffix=''):
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
            for histName,gDict in self.content.items():
                for gr,hist in gDict.items():
                    if gr.split("__")[0] == group:
                        if len(gr.split("__")) == 1:
                            hist.Write(histName)
                        else:
                            syst = gr.split('__')[1]
                            systName = syst.replace('Up','').replace('Down','') 
                            if systName not in systematics:
                                systematics.append(systName)
                            hist.Write(histName+'__'+syst.replace('Up','up').replace('Down','down'))
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
            #config['files'][group+'.root'].update({k:v for k,v in gconfig.items() if k not in ['files','type']})
        config['groups'] = {}
        for group,gconfig in self.groups.items():
            if self.pseudodata and group == 'data_real':
                continue
            if gconfig['type'] != 'signal':
                config['groups'][group] = {k:v for k,v in gconfig.items() if k not in ['files','type']}
            if self.DYEstimation is not None and group == self.DYEstimation['replace']:
                if self.pseudodata:
                    config['groups'][group]['legend'] += ' (from pseudo-data)'
                else:
                    config['groups'][group]['legend'] += ' (from data)'


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
