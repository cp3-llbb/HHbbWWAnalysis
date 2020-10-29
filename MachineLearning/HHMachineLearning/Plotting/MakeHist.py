import os
import sys
import argparse
import glob
import logging
import copy
import yaml
import pprint
import traceback

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = 2000#[ROOT.kPrint, ROOT.kInfo]#, kWarning, kError, kBreak, kSysError, kFatal;

import Classes
from Classes import Plot_TH1, Plot_TH2, Plot_Ratio_TH1, Plot_Multi_TH1, Plot_ROC, LoopPlotOnCanvas, MakeROCPlot, ProcessYAML, Plot_Multi_ROC, MakeMultiROCPlot


def main():
    #############################################################################################
    # Options #
    #############################################################################################
    parser = argparse.ArgumentParser(description='From given set of root files, make different histograms in a root file')
    parser.add_argument('-m','--model', action='store', required=True, type=str, default='',
                  help='NN model to be used')
    parser.add_argument('-v','--verbose', action='store_true', required=False, default=False,
            help='Show DEGUG logging')
    opt = parser.parse_args() 

    # Logging #
    if opt.verbose:
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s - %(levelname)s - %(message)s') 
    else:
        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s - %(levelname)s - %(message)s')
    #############################################################################################
    # Select samples #
    #############################################################################################
    INPUT_DIR = os.path.join('/home/ucl/cp3/fbury/scratch/HHMachineLearning_output/',opt.model)
    #INPUT_DIR = os.path.asbpath(opt.model)
    logging.info('Taking inputs from %s'%INPUT_DIR)
    
    # For each directory, puth the path in dict the value is the list of histograms #
    class Plots():
        def __init__(self,name,override_params):
            self.name = name                            # Name of the subdir
            self.path = os.path.join(INPUT_DIR,name)    # Full path to files
            self.override_params = override_params      # Parameters to override in the configs
            self.list_histo = []                        # List of histograms to be filled

    # Select template #
    class Template:
        def __init__(self,tpl,class_name):
            self.tpl = tpl                          # Template ".yml.tpl" to be used
            self.class_name = class_name            # Name of one of the class to use (TH1,TH2,...)

    # Select template #
    class ROC:
        def __init__(self,tpl,class_name,def_name,plot_name):
            self.tpl = tpl                          # Template ".yml.tpl" to be used containing all the ROC curves parameters (they will all be plotted in the same figure)
            self.class_name = class_name            # Name of one of the class to use (ROC or multiROC)
            self.def_name = def_name                # Name of of the processing function 
            self.plot_name = plot_name              # Name of the plot -> .png
            self.list_instance = []
        def AddInstance(self,instance):
            self.list_instance.append(instance) # Contains the ROC configs listed in the tpl file
        def clearInstance(self):
            self.list_instance = []

    #///////////////      TO BE MODIFIED BY USER       ////////////////
    list_plots = [
                    Plots(name = 'test',
                          override_params = {}),
                                              
                 ]

    templates = [
                    ############       TH1       ##############
                    Template(tpl = 'TH1_event_weight.yml.tpl',
                             class_name = 'Plot_TH1'),
                    #########       TH1 Ratio       ###########
                    #Template(tpl = '.yml.tpl',
                    #         class_name = 'Plot_Ratio_TH1'),
                    #########       TH1 Multi       ###########
                    Template(tpl = 'TH1Multi_class.yml.tpl',
                             class_name = 'Plot_Multi_TH1'),
                    ############       TH2       ##############
                    #Template(tpl = '.yml.tpl',
                    #         class_name = 'Plot_TH2'),
               ] 

    rocs = [
                ########### Bivariate ROC ############
                #ROC(tpl = '.yml.tpl',
                #    class_name = 'Plot_ROC',
                #    def_name = 'MakeROCPlot',
                #    plot_name = ''),
                ########### Multivariate ROC ############
                ROC(tpl = 'ROCMulti_class.yml.tpl',
                    class_name = 'Plot_Multi_ROC',
                    def_name = 'MakeMultiROCPlot',
                    plot_name = 'Multiclass'),
           ]

    # Make the output dir #
    OUTPUT_PDF = os.path.join(os.getcwd(),'PDF',opt.model)
    if not os.path.exists(OUTPUT_PDF):
        os.makedirs(OUTPUT_PDF)

    #############################################################################################
    # Loop over the different paths #
    #############################################################################################
    # Loop over the files #
    for obj in list_plots:
        logging.info('Starting Plotting from %s subdir'%obj.name)
        files = sorted(glob.glob(obj.path+'/*.root'))

        if len(files) == 0:
            logging.error('Could not find %s at %s'%(obj.name,obj.path))

       # Instantiate all the ROCs #
        for roc in rocs:
            logging.debug('ROC template "%s" -> Class "%s" and process function %s '%(roc.tpl, roc.class_name, roc.def_name))
            roc.clearInstance()
            YAML = ProcessYAML(roc.tpl) # Contain the ProcessYAML objects
            YAML.Particularize(obj.name)
            for name,config in YAML.config.items():
                class_ = getattr(Classes, roc.class_name)
                logging.info('\tInitializing ROC %s'%(name))
                logging.debug('\t ... Containing')
                logging.debug(pprint.pformat(config))
                instance = class_(**config)
                roc.AddInstance(instance)

        # Loop over files #
        for f in files:
            fullname = os.path.basename(f).replace('.root','')
            logging.info('Processing weights from %s'%(fullname+'.root'))
            #if fullname.startswith('DY'):
            #    filename = 'Drell-Yan'
            #elif fullname.startswith('TT'):
            #    filename = 't#bar{t}'
            #elif fullname.startswith('HToZA'):
            #    filename = 'H#rightarrowZA'
            #else:
            filename = fullname
            ##############  ROC  section ################ 
            for roc in rocs:
                for inst_roc in roc.list_instance:
                    try:
                        valid = inst_roc.AddToROC(f)
                        if valid:
                            logging.info('\tAdded to ROC %s'%(inst_roc.title))
                    except Exception as e:
                        logging.warning('Could not add to ROC due to "%s"'%(e))
                        traceback.print_exc()
            
            ##############  HIST section ################ 
            # Loop over the templates #
            for template in templates: 
                logging.debug('Hist template "%s" -> Class "%s"'%(template.tpl, template.class_name))
                list_config = [] # Will contain the dictionaries of parameters
                YAML = ProcessYAML(template.tpl) # Contain the ProcessYAML objects
                params = {**{'filepath':f,'filename':filename},**obj.override_params}
                # Get the list of configs #
                YAML.Particularize(fullname)
                YAML.Override(params)
                # loop over the configs #
                for name,config in YAML.config.items():
                    try:
                        class_ = getattr(Classes, template.class_name)
                        logging.info('\tPlot %s'%(name))
                        instance = class_(**config)
                        instance.MakeHisto()
                        obj.list_histo.append(instance)
                    except Exception as e:
                        logging.warning('Could not plot %s due to "%s"'%(name,e))
                        traceback.print_exc()

        # Process ROCs #
        for roc in rocs:
            for inst_roc in roc.list_instance:
                try:
                    inst_roc.ProcessROC() 
                    logging.info('\tProcessed ROC %s'%(inst_roc.title))
                except Exception as e:
                    logging.warning('Could not process ROC due to "%s"'%(e))
                    traceback.print_exc()

        # Make ROCs graphs #
        for roc in rocs:
            try:
                def_ = getattr(Classes, roc.def_name)
                def_(roc.list_instance,name=os.path.join(OUTPUT_PDF,roc.plot_name+'_'+obj.name))
            except Exception as e:
                logging.warning('Could not plot ROC due to "%s"'%(e))
                traceback.print_exc()

    #############################################################################################
    # Save histograms #
    #############################################################################################
    for obj in list_plots:
        PDF_name = os.path.join(OUTPUT_PDF,obj.name) 
        if len(obj.list_histo) != 0:
            LoopPlotOnCanvas(PDF_name,obj.list_histo)
    logging.info("All Canvas have been printed")

if __name__ == "__main__":
    main()   
