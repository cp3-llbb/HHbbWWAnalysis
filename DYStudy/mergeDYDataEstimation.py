import os
import sys
import glob
import shutil
import yaml
import argparse
import subprocess
import ROOT
import time

class Merge:
    def __init__(self,path_base,path_estimation,btag,era,force=False):
        self.path_base          = path_base
        self.path_estimation    = path_estimation
        self.btag               = btag 
        self.era                = era
        self.force              = force

        assert self.btag == "1b" or self.btag == "2b"

        self.copyBase()
        self.renameHist()
        self.modifyYAML('compare')
        self.modifyYAML('replace')

    def copyBase(self):
        if self.path_base.endswith('/'):
            self.path_base = self.path_base[:-1]
        self.path_copy = self.path_base + '_DataDY' + self.btag 
        self.sample_dataDY = []
        self.sample_MCDY = []

        if os.path.exists(self.path_copy) and self.force:
            print ('Directory %s will be delete and recreated'%self.path_copy)
            shutil.rmtree(self.path_copy)
        if os.path.exists(self.path_copy):
            print ("Directory %s already exists"%self.path_copy)
        if not os.path.exists(self.path_copy):
            # Root files #
            os.makedirs(os.path.join(self.path_copy,'results'))
            print ('Created dir %s'%os.path.join(self.path_copy,'results'))
            # Copy from base #
            print ('Copying root files from base path %s'%self.path_base)
            for f in glob.glob(os.path.join(self.path_base,'results','*root')):
                try:
                    shutil.copyfile(f,os.path.join(self.path_copy,'results',os.path.basename(f)))
                    if 'DY' in os.path.basename(f) and not '__' in os.path.basename(f):
                        self.sample_MCDY.append(os.path.basename(f))
                except IOError:
                    print ("\t[WARNING] Could not copy file %s to directory %s"%(f,os.path.join(self.path_copy,'results')))
                print ('\t%s -> %s'%(f,os.path.join(self.path_copy,'results',os.path.basename(f))))
            # Hadd from DY estimation #
            print ('Hadding all data DY root files from estimation path %s'%self.path_estimation)
            list_files = glob.glob(os.path.join(self.path_estimation,'results','*root'))
            subprocess.call("hadd -f "+os.path.join(self.path_copy,'results','DataDY.root')+' '+' '.join([l for l in list_files if '__' not in l]),shell=True)
            self.sample_dataDY.append('DataDY.root')
            # Copy other useful stuff #
            try:
                shutil.copyfile(os.path.join(self.path_base,'plots.yml'),os.path.join(self.path_copy,'plots.yml'))
                self.path_yaml = os.path.join(self.path_copy,'plots.yml')
            except IOError:
                print ("\t[WARNING] Could not copy file %s to directory %s"%(os.path.join(self.path_base,'plots.yml'),os.path.join(self.path_copy,'plots.yml')))
            golden_json = glob.glob(os.path.join(self.path_base,'*txt'))
            if len(golden_json) != 1:
                raise RuntimeError('Found several txt files for golden json')
            else:
                golden_json = os.path.basename(golden_json[0])
            try:
                shutil.copyfile(os.path.join(self.path_base,golden_json),os.path.join(self.path_copy,golden_json))
            except IOError:
                print ("\t[WARNING] Could not copy file %s to directory %s"%(os.path.join(self.path_base,golden_json),os.path.join(self.path_copy,golden_json)))

        else:
            self.path_yaml = os.path.join(self.path_copy,'plots.yml')
            for f in glob.glob(os.path.join(self.path_copy,'results','*DY*.root')):
                if '__' in os.path.basename(f):
                    continue # Skeletons
                if os.path.basename(f).startswith('DataDY'):
                    self.sample_dataDY.append(os.path.basename(f))
                if os.path.basename(f).startswith('DY'):
                    self.sample_MCDY.append(os.path.basename(f))
        print ("DY MC to be replaced")
        for sample in self.sample_MCDY:
            print ('...',sample)
        print ("DY data estimation to used instead")
        for sample in self.sample_dataDY:
            print ('...',sample)

    @staticmethod
    def getListHistogram(root_path,class_name=None,exclude_str=None,include_systematics=False):
        root_file = ROOT.TFile(root_path)
        hist_names = []
        for key in root_file.GetListOfKeys():
            # Get name of object #
            kname = key.GetName()

            # If excluding string #
            if exclude_str and exclude_str in kname:
                continue

            # Check if object is a systematics hist #
            if not include_systematics and '__' in kname:
                continue

            # If only asking for a specific object class #
            if class_name:
                obj = root_file.Get(kname)
                if class_name not in obj.ClassName():
                    continue

            hist_names.append(kname)
        root_file.Close()
        return hist_names 


    def renameHist(self):
        # Find list of histograms in data corrected DY #
        hist_name = self.getListHistogram(os.path.join(self.path_copy,'results',self.sample_dataDY[0]),'TH1D','ElMu',True)

        # Loop over root files to rename #
        for sample in self.sample_dataDY:
            path = os.path.join(self.path_copy,'results',sample)
            for hname in hist_name:
                if self.btag == '1b':
                    if 'ResolvedNoBtag' in hname:
                        cmd = 'rootmv {0}:{1} {0}:{2}'.format(path,hname,hname.replace('ResolvedNoBtag','ResolvedOneBtag')) 
                        subprocess.call(cmd,shell=True)
                        print ("Renamed {} -> {}".format(hname,hname.replace('ResolvedNoBtag','ResolvedOneBtag')))
                        hname = hname.replace('ResolvedNoBtag','ResolvedOneBtag')
                    if '_leadjet_' in hname: # Only _leadjet_ -> _leadbjet_  (other _subleadjet can remain the same)
                        cmd = 'rootmv {0}:{1} {0}:{2}'.format(path,hname,hname.replace('leadjet','leadbjet')) 
                        subprocess.call(cmd,shell=True)
                        print ("Renamed {} -> {}".format(hname,hname.replace('leadjet','leadbjet')))
                        hname = hname.replace('leadjet','leadbjet')
                if self.btag == '2b':
                    if 'ResolvedNoBtag' in hname:
                        cmd = 'rootmv {0}:{1} {0}:{2}'.format(path,hname,hname.replace('ResolvedNoBtag','ResolvedTwoBtags')) 
                        subprocess.call(cmd,shell=True)
                        print ("Renamed {} -> {}".format(hname,hname.replace('ResolvedNoBtag','ResolvedTwoBtags')))
                        hname = hname.replace('ResolvedNoBtag','ResolvedTwoBtags')
                    if 'leadjet' in hname: # _leadjet_ -> _leadbjet_ + _subleadjet -> _subleadbjet_
                        cmd = 'rootmv {0}:{1} {0}:{2}'.format(path,hname,hname.replace('leadjet','leadbjet')) 
                        subprocess.call(cmd,shell=True)
                        print ("Renamed {} -> {}".format(hname,hname.replace('leadjet','leadbjet')))
                        hname = hname.replace('leadjet','leadbjet')

    def modifyYAML(self,mode):
        assert mode in ['compare','replace']
        with open(self.path_yaml,'r') as handle:
            config = yaml.load(handle,Loader=yaml.FullLoader)

        #----- Configuration -----#
        assert config["configuration"]["eras"][0] == self.era
        luminosity = float(config["configuration"]["luminosity"][self.era] )

        #----- Histograms -----#
        # Find all histograms in one of the data DY files #
        hist_names = self.getListHistogram(os.path.join(self.path_copy,'results',self.sample_dataDY[0]),'TH1D','ElMu',False)

        # Only keep in yaml dict the hist that are in data DY #
        list_plots = config['plots'].keys()
        hist_to_delete = [l for l in list_plots if l not in hist_names]
        for name in hist_to_delete:
            del config['plots'][name]


        #----- Files -----#

        # Samples already in the YAML #
        if mode == 'replace': # We replace DY MC by DY data estimation
            # Remove the DY MC files #
            for sample in self.sample_MCDY:
                del config['files'][sample]
        if mode == 'compare':
            # Remove all but the DY MC files #
            sample_to_delete = [s for s in config['files'].keys() if s not in self.sample_MCDY]
            for sample in sample_to_delete:
                del config['files'][sample]
            # Put Xsec and sum weight at 1 #
            for sample in self.sample_MCDY:
                config['files'][sample]['group'] = 'DY_MC'
                config['files'][sample]['stack-index'] = 0
            
        # Add the Data driven DY MC #
        for sample_dataDY in self.sample_dataDY:
            config['files'][sample_dataDY] = {
                                                'cross-section' : 1./luminosity,
                                                'era' : self.era,
                                                'generated-events' : 1.,
                                                'group' : 'DY' if mode == 'replace' else 'DY_data' ,
                                                'type' : 'mc',
                                                'stack-index' : 0 if mode == 'replace' else 1
                                             }
        #----- Group -----#
        if mode == 'compare':
            config['groups'] = {'DY_MC': {'line-color' : '#1a83a1',
                                          'line-width' : 2,
                                          'fill-color' : -1,
                                          'legend' : 'DY from MC',
                                          'order' : 1},
                                'DY_data' : {'fill-color' : '#0e0b6b',
                                             'line-width' : 2,
                                             'fill-color' : -1,
                                             'legend' : 'DY from data estimation',
                                             'order' : 2}}
    


        #---- Save new yaml file -----#
        if mode == 'replace':
            path_new_yaml = self.path_yaml.replace('.yml','_DYDataEstimation.yml')
        if mode == 'compare':
            path_new_yaml = self.path_yaml.replace('.yml','_compareDYEstimation.yml')
        with open(path_new_yaml,'w') as handle:
            yaml.dump(config,handle)

        print ('YAML files for plots modified as %s'%path_new_yaml)
        
        self.printCommand(path_new_yaml)
    
    def printCommand(self,path_yaml):
        outdir = path_yaml.replace('.yml','')+'_'+self.era
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        com = "plotIt -i {input} -o {output} -y -e {era} {yaml}".format(**{'input': self.path_copy,
                                                                           'output': outdir,
                                                                           'era': self.era,
                                                                           'yaml': path_yaml})

        print ("PlotIt command (must be used in a bamboo env) :")
        print (com)



parser = argparse.ArgumentParser(description='Merge usual bamboo output with the DY estimation from data')
parser.add_argument('--path_base', action='store', required=True, type=str,
                    help='Path containing the usual MC and data produced by bamboo')
parser.add_argument('--path_estimation', action='store', required=True, type=str,
                    help='Path containing DY estimation from data')
parser.add_argument('--era', action='store', required=True, type=str,
                    help='Era : 2016, 2017, 2018')
parser.add_argument('--b1', action='store_true', required=False, default=False,
                    help='1 btag region')
parser.add_argument('--b2', action='store_true', required=False, default=False,
                    help='2 btag region')
parser.add_argument('--force', action='store_true', required=False, default=False,
                    help='Force the recreation of the merged directory even if it already exists')

args = parser.parse_args()

if args.b1 and not args.b2:
    btag = '1b'    
if args.b2 and not args.b1:
    btag = '2b'

instance = Merge(args.path_base,
                 args.path_estimation,
                 btag,
                 args.era,args.force)

 
