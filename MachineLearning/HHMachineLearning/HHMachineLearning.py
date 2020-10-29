#!/usr/bin/env python

import re
import math
import glob
import csv
import yaml
import os
import sys
import pprint
import logging
import copy
import pickle
import json
#import psutil
from functools import reduce
import operator
import itertools
import matplotlib.pyplot as plt
if plt.rcParams['backend'] == 'TkAgg':
    raise ImportError("Change matplotlib backend to 'Agg' in ~/.config/matplotlib/matplotlibrc")

import argparse
import numpy as np
import pandas as pd

from sklearn import preprocessing
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder, OneHotEncoder

# Personal files #
    
def get_options():
    """
    Parse and return the arguments provided by the user.
    """
    parser = argparse.ArgumentParser(description='HHMachineLearning')

    # Scan, deploy and restore arguments #
    a = parser.add_argument_group('Scan, deploy and restore arguments')
    a.add_argument('--scan', action='store', required=False, type=str, default='',
        help='Name of the scan to be used')
    a.add_argument('--task', action='store', required=False, type=str, default='',
        help='Name of dict to be used for scan (Used by function itself when submitting jobs or DEBUG)')
    a.add_argument('--generator', action='store_true', required=False, default=False, 
        help='Whether to use a generator for the neural network')
    a.add_argument('--resume', action='store_true', required=False, default=False,
        help='Whether to resume the training of a given model (path in parameters.py)')

    # Splitting and submitting jobs arguments #
    b = parser.add_argument_group('Splitting and submitting jobs arguments')
    b.add_argument('--split', action='store', required=False, type=int, default=0,
        help='Number of parameter sets per jobs to be used for splitted training for slurm submission (if -1, will create a single subdict)')
    b.add_argument('--submit', action='store', required=False, default='', type=str,
        help='Whether to submit on slurm and name for the save (must have specified --split)')
    b.add_argument('--resubmit', action='store', required=False, default='', type=str,
        help='Whether to resubmit failed jobs given a specific path containing the jobs that succeded')
    b.add_argument('--debug', action='store_true', required=False, default=False,
        help='Debug mode of the slurm submission, does everything except submit the jobs')

    # Analyzing or producing outputs for given model (csv or zip file) #
    c = parser.add_argument_group('Analyzing or producing outputs for given model (csv or zip file)')
    c.add_argument('--report', action='store', required=False, type=str, default='',
        help='Name of the csv file for the reporting (without .csv)')
    c.add_argument('--model', action='store', required=False, nargs='+', type=str, default='',
        help='Loads the provided model name (without .zip and type, it will find them), list can be used in case of cross validation') 
    c.add_argument('--test', action='store_true', required=False, default=False,
        help='Applies the provided model (do not forget -o) on the test set and output the tree') 
    c.add_argument('--output', action='store', required=False, nargs='+', type=str, default=[], 
        help='Applies the provided model (do not forget -o) on the list of keys from sampleList.py (separated by spaces)') 

    # Concatenating csv files arguments #
    d = parser.add_argument_group('Concatenating csv files arguments')
    d.add_argument('--csv', action='store', required=False, type=str, default='',
        help='Whether to concatenate the csv files from different slurm jobs into a main one, \
              please provide the path to the csv files')

    # Additional arguments #
    f = parser.add_argument_group('Additional arguments')
    f.add_argument('-v','--verbose', action='store_true', required=False, default=False,
        help='Show DEGUG logging')
    f.add_argument('--nocache', action='store_true', required=False, default=False,
        help='Will not use the cache and will not save it')
    f.add_argument('--GPU', action='store_true', required=False, default=False,
        help='GPU requires to execute some commandes before')

    opt = parser.parse_args()

    if opt.split!=0 or opt.submit!='':
        if opt.scan!='' or opt.report!='':
            logging.critical('These parameters cannot be used together')  
            sys.exit(1)
    if opt.submit!='': # Need --output or --split arguments
        if opt.split==0 and len(opt.output)==0:
            logging.warning('In case of learning you forgot to specify --split')
            sys.exit(1)
    if opt.split!=0 and (opt.report!='' or opt.output!='' or opt.csv!='' or opt.scan!=''):
        logging.warning('Since you have specified a split, all the other arguments will be skipped')
    if opt.csv!='' and (opt.report!='' or opt.output!='' or opt.scan!=''):
        logging.warning('Since you have specified a csv concatenation, all the other arguments will be skipped')
    if opt.report!='' and (opt.output!='' or opt.scan!=''):
        logging.warning('Since you have specified a scan report, all the other arguments will be skipped')
    if (opt.test or len(opt.output)!=0) and opt.output == '': 
        logging.critical('You must specify the model with --output')
        sys.exit(1)
    if opt.generator:
        logging.info("Will use the generator")
    if opt.resume:
        logging.info("Will resume the training of the model")

    return opt

def main():
    #############################################################################################
    # Preparation #
    #############################################################################################
    # Get options from user #
    logging.basicConfig(level=logging.DEBUG,format='%(asctime)s - %(levelname)s - %(message)s',datefmt='%m/%d/%Y %H:%M:%S')
    opt = get_options()
    # Verbose logging #
    if not opt.verbose:
        logging.getLogger().setLevel(logging.INFO)


    # Private modules containing Pyroot #
    from NeuralNet import HyperModel
    from import_tree import LoopOverTrees
    from produce_output import ProduceOutput
    from make_scaler import MakeScaler
    from submit_on_slurm import submit_on_slurm
    from generate_mask import GenerateMask, GenerateSampleMasks, GenerateSliceIndices, GenerateSliceMask
    from split_training import DictSplit
    from concatenate_csv import ConcatenateCSV
    from threadGPU import utilizationGPU
    import parameters

    # Needed because PyROOT messes with argparse

    logging.info("="*98)
    logging.info("   _    _ _    _ __  __            _     _            _                           _             ")
    logging.info("  | |  | | |  | |  \/  |          | |   (_)          | |                         (_)            ")
    logging.info("  | |__| | |__| | \  / | __ _  ___| |__  _ _ __   ___| |     ___  __ _ _ __ _ __  _ _ __   __ _ ")
    logging.info("  |  __  |  __  | |\/| |/ _` |/ __| '_ \| | '_ \ / _ \ |    / _ \/ _` | '__| '_ \| | '_ \ / _` |")
    logging.info("  | |  | | |  | | |  | | (_| | (__| | | | | | | |  __/ |___|  __/ (_| | |  | | | | | | | | (_| |")
    logging.info("  |_|  |_|_|  |_|_|  |_|\__,_|\___|_| |_|_|_| |_|\___|______\___|\__,_|_|  |_| |_|_|_| |_|\__, |")
    logging.info("                                                                                           __/ |")
    logging.info("                                                                                          |___/ ")
    logging.info("="*98)

    # Make path model #
    path_model = os.path.join(parameters.main_path,'model')
    if not os.path.exists(path_model):
        os.mkdir(path_model)

    #############################################################################################
    # Splitting into sub-dicts and slurm submission #
    #############################################################################################
    if opt.submit != '':
        if opt.split != 0:
            DictSplit(opt.split,opt.submit,opt.resubmit)
            logging.info('Splitting jobs done')
        
        # Arguments to send #
        args = '' 
        if opt.generator:           args += ' --generator '
        if opt.generator:           args += ' --generator '
        if opt.GPU:                 args += ' --GPU '
        if opt.resume:              args += ' --resume '
        if opt.nocache:             args += ' --nocache'
        if opt.model!='':           args += ' --model '+opt.model+' '
        if len(opt.output)!=0:      args += ' --output '+ ' '.join(opt.output)+' '

        if opt.submit!='':
            logging.info('Submitting jobs with args "%s"'%args)
            if opt.resubmit:
                submit_on_slurm(name=opt.submit+'_resubmit',debug=opt.debug,args=args)
            else:
                submit_on_slurm(name=opt.submit,debug=opt.debug,args=args)
        sys.exit()

    #############################################################################################
    # CSV concatenation #
    #############################################################################################
    if opt.csv!='':
        logging.info('Concatenating csv files from : %s'%(opt.csv))
        dict_csv = ConcatenateCSV(opt.csv)

        sys.exit()

    #############################################################################################
    # Reporting given scan in csv file #
    #############################################################################################
    if opt.report != '':
        instance = HyperModel(opt.report)
        instance.HyperReport(parameters.eval_criterion)

        sys.exit()

    #############################################################################################
    # Output of given files from given model #
    #############################################################################################
    if opt.model != '' and len(opt.output) != 0:
        # Create directory #
        path_output = os.path.join(parameters.path_out,opt.model)
        if not os.path.exists(path_output):
            os.mkdir(path_output)

        # Instantiate #
        inst_out = ProduceOutput(model=os.path.join(parameters.path_model,opt.model),generator=opt.generator)
        # Loop over output keys #
        for key in opt.output:
            # Create subdir #
            path_output_sub = os.path.join(path_output,key+'_output')
            if not os.path.exists(path_output_sub):
                os.mkdir(path_output_sub)
            try:
                inst_out.OutputNewData(input_dir=samples_path,list_sample=samples_dict[key],path_output=path_output_sub)
            except Exception as e:
                logging.critical('Could not process key "%s" due to "%s"'%(key,e))
        sys.exit()
    #############################################################################################
    # Data Input and preprocessing #
    #############################################################################################
    # Memory Usage #
    #pid = psutil.Process(os.getpid())
    logging.info('Current pid : %d'%os.getpid())

    # Input path #
    logging.info('Starting tree importation')

    # Import variables from parameters.py
    variables = parameters.inputs+parameters.outputs+parameters.other_variables
    list_inputs  = parameters.inputs
    list_outputs = parameters.outputs

    # Load samples #
    with open (parameters.config,'r') as f:
        sampleConfig = yaml.load(f)

    if not opt.generator:
        if opt.nocache:
            logging.warning('No cache will be used nor saved')
        if os.path.exists(parameters.train_cache) and not opt.nocache:
            logging.info('Will load training data from cache')
            logging.info('... Training set : %s'%parameters.train_cache)
            train_all = pd.read_pickle(parameters.train_cache)
            if os.path.exists(parameters.test_cache) and not opt.nocache and not parameters.crossvalidation:
                logging.info('Will load testing data from cache')
                logging.info('... Testing  set : %s'%parameters.test_cache)
                test_all = pd.read_pickle(parameters.test_cache)
        else:
            # Import arrays #
            data_dict = {}
            xsec_dict = dict()
            event_weight_sum_dict = dict()
            for era in parameters.eras:
                with open(parameters.xsec_json.format(era=era),'r') as handle:
                    xsec_dict[era] = json.load(handle)
                with open(parameters.event_weight_sum_json.format(era=era),'r') as handle:
                    event_weight_sum_dict[era] = json.load(handle)

            for node in parameters.nodes:
                logging.info('Starting data importation for class %s'%node)
                strSelect = [f'{cat}_{channel}_{node}' for channel in parameters.channels for cat in parameters.categories]
                data_node = None

                for era in parameters.eras:
                    samples_dict = sampleConfig['sampleDict'][era]
                    if len(samples_dict.keys())==0:
                        logging.info('\tSample dict for era {} is empty'.format(era))
                        continue
                    list_sample = [sample for key in strSelect for sample in samples_dict[key]]
                    data_node_era = LoopOverTrees(input_dir                 = sampleConfig['sampleDir'],
                                                  variables                 = variables,
                                                  weight                    = parameters.weight,
                                                  list_sample               = list_sample,
                                                  cut                       = parameters.cut,
                                                  xsec_dict                 = xsec_dict,
                                                  event_weight_sum_dict     = event_weight_sum_dict,
                                                  lumi_dict                 = parameters.lumidict,
                                                  eras                      = era,
                                                  tree_name                 = parameters.tree_name,
                                                  additional_columns        = {'tag':node,'era':era},
                                                  stop                      = 300000) # TODO : remove 
                    data_node_era = data_node_era.sample(frac=1)[:300000] # TODO : remove 
                    if data_node is None:
                        data_node = data_node_era
                    else:
                        data_node = pd.concat([data_node,data_node_era],axis=0)
                    era_str = '{:5s} class in era {} : sample size = {:10d}'.format(node,era,data_node_era.shape[0])
                    if parameters.weight is not None:
                        era_str += ', weight sum = {:.3e} (with normalization = {:.3e})'.format(data_node_era[parameters.weight].sum(),data_node_era['event_weight'].sum())
                    logging.info(era_str)
                data_dict[node] = data_node
                all_eras_str = '{:5s} class for all eras : sample size = {:10d}'.format(node,data_node.shape[0])
                if parameters.weight is not None:
                    all_eras_str +=  ', weight sum = {:.3e} (with normalization = {:.3e})'.format(data_node[parameters.weight].sum(),data_node['event_weight'].sum())
                logging.info(all_eras_str)
            list_inputs  = [var.replace('$','') for var in parameters.inputs]
            list_outputs = [var.replace('$','') for var in parameters.outputs]

            # Weight equalization #
            N = sum([data.shape[0] for data in data_dict.values()])/len(data_dict)
            for node, data in data_dict.items():
                if parameters.weight is not None:
                    sum_weights = data['event_weight'].sum()
                    logging.info('Sum of weight for %s samples : %.2e'%(node,sum_weights))
                    data['learning_weights'] = data['event_weight']/sum_weights*N
                    logging.info('\t -> After equalization : %0.2e (factor %0.2e)'%(data['learning_weights'].sum(),1.e5/sum_weights))
                else:
                    data['learning_weights'] = pd.Series([1.]*data.shape[0],index=data.index)

            # Data splitting #
            train_dict = {}
            test_dict = {}
            for node,data in data_dict.items():
                if parameters.crossvalidation: # Cross-validation
                    if parameters.splitbranch not in data.columns:
                        raise RuntimeError('Asked for cross validation mask but cannot find the slicing array')
                    try:
                        data['mask'] = (data[parameters.splitbranch] % parameters.N_slices).to_numpy()
                        # Will contain numbers : 0,1,2,...N_slices-1
                    except ValueError:
                        logging.critical("Problem with the masking")
                        raise ValueError
                else: # Classic separation
                    mask = GenerateMask(data.shape[0],parameters.suffix+'_'+node)
                    try:
                        train_dict[node] = data[mask==True]
                        test_dict[node]  = data[mask==False]
                    except ValueError:
                        logging.critical("Problem with the mask you imported, has the data changed since it was generated ?")
                        raise ValueError
             
            if parameters.crossvalidation:
                train_all = pd.concat(data_dict.values(),copy=True).reset_index(drop=True)
                test_all = pd.DataFrame(columns=train_all.columns) # Empty to not break rest of script
            else:
                train_all = pd.concat(train_dict.values(),copy=True).reset_index(drop=True)
                test_all  = pd.concat(test_dict.values(),copy=True).reset_index(drop=True)
            del data_dict 
            if not parameters.crossvalidation:
                del train_dict, test_dict
            #logging.info('Current memory usage : %0.3f GB'%(pid.memory_info().rss/(1024**3)))

            # Randomize order, we don't want only one type per batch #
            random_train = np.arange(0,train_all.shape[0]) # needed to randomize x,y and w in same fashion
            np.random.shuffle(random_train) # Not needed for testing
            train_all = train_all.iloc[random_train]
              
            # Add target #
            label_encoder = LabelEncoder()
            onehot_encoder = OneHotEncoder(sparse=False)
            label_encoder.fit(parameters.nodes)
            # From strings to labels #
            train_integers = label_encoder.transform(train_all['tag']).reshape(-1, 1)
            if not parameters.crossvalidation:
                test_integers = label_encoder.transform(test_all['tag']).reshape(-1, 1)
            # From labels to strings #
            onehotobj = onehot_encoder.fit(np.arange(len(list_outputs)).reshape(-1, 1))
            train_onehot = onehotobj.transform(train_integers)
            if not parameters.crossvalidation:
                test_onehot = onehotobj.transform(test_integers)
            # From arrays to pd DF #
            train_cat = pd.DataFrame(train_onehot,columns=label_encoder.classes_,index=train_all.index)
            if not parameters.crossvalidation:
                test_cat = pd.DataFrame(test_onehot,columns=label_encoder.classes_,index=test_all.index)
            # Add to full #
            train_all = pd.concat([train_all,train_cat],axis=1)
            if not parameters.crossvalidation:
                test_all = pd.concat([test_all,test_cat],axis=1)

            # Preprocessing #
            # The purpose is to create a scaler object and save it
            # The preprocessing will be implemented in the network with a custom layer
            if opt.scan!='': # If we don't scan we don't need to scale the data
                MakeScaler(train_all,list_inputs) 

          # Caching #
            if not opt.nocache:
                train_all.to_pickle(parameters.train_cache)
                logging.info('Data saved to cache')
                logging.info('... Training set : %s'%parameters.train_cache)
                if not parameters.crossvalidation:
                    test_all.to_pickle(parameters.test_cache)
                    logging.info('... Testing  set : %s'%parameters.test_cache)
         
        logging.info("Sample size seen by network : %d"%train_all.shape[0])
        #logging.info('Current memory usage : %0.3f GB'%(pid.memory_info().rss/(1024**3)))
        if parameters.crossvalidation: 
            N = train_all.shape[0]
            logging.info('Cross-validation has been requested on set of %d events'%N)
            for i in range(parameters.N_models):
                slices_apply , slices_eval, slices_train = GenerateSliceIndices(i)
                logging.info('... Model %d :'%i)
                for slicename, slices in zip (['Applied','Evaluated','Trained'],[slices_apply , slices_eval, slices_train]):
                    selector = np.full((train_all.shape[0]), False, dtype=bool)
                    selector = GenerateSliceMask(slices,train_all['mask'])
                    n = train_all[selector].shape[0]
                    logging.info('     %10s on %10d [%3.2f%%] events'%(slicename,n,n*100/N)+' (With mask indices : ['+','.join([str(s) for s in slices])+'])')
        else:
            logging.info("Sample size for the output  : %d"%test_all.shape[0])
    else:
        logging.info("You asked for generator so no data input has been done")
        list_samples = [os.path.join(sampleConfig['sampleDir'],sample) for era in parameters.eras for samples in sampleConfig['sampleDict'][era].values() for sample in samples ]
        # Produce mask if not cross val #
        if not parameters.crossvalidation:
            logging.info("Will generate masks for each sample")
            GenerateSampleMasks(list_samples,parameters.suffix)
        # Produce scaler #
        MakeScaler(list_inputs  = list_inputs,
                   generator    = True,
                   batch        = parameters.output_batch_size,
                   list_samples = list_samples) 
            
        train_all = None
        test_all = None
    list_inputs  = [var.replace('$','') for var in parameters.inputs]
    list_outputs = [var.replace('$','') for var in parameters.outputs]

    #############################################################################################
    # DNN #
    #############################################################################################
    # Start the GPU monitoring thread #
    if opt.GPU:
        thread = utilizationGPU(print_time = 900,
                                print_current = False,
                                time_step=0.01)
        thread.start()

    if opt.scan != '':
        instance = HyperModel(opt.scan,list_inputs,list_outputs)
        if parameters.crossvalidation:
            for i in range(parameters.N_models):
                logging.info("*"*80)
                logging.info("Starting training of model %d"%i)
                instance.HyperScan(data=train_all,
                                   task=opt.task,
                                   generator=opt.generator,
                                   resume=opt.resume,
                                   model_idx=i)
                instance.HyperDeploy(best='eval_error')
        else:
            instance.HyperScan(data=train_all,
                               task=opt.task,
                               generator=opt.generator,
                               resume=opt.resume)
            instance.HyperDeploy(best='eval_error')

    if opt.GPU:
        # Closing monitor thread #
        thread.stopLoop()
        thread.join()
        
    if len(opt.model) != 0: 
        # Make path #
        output_name = "test" 
        model_name = opt.model[0]
        if parameters.crossvalidation:
            model_name = model_name[:-1]
        path_output = os.path.join(parameters.path_out,model_name,output_name)
        if not os.path.exists(path_output):
            os.makedirs(path_output)

        # Instance of output class #
        inst_out = ProduceOutput(model=[os.path.join(parameters.main_path,'model',model) for model in opt.model],
                                 generator=opt.generator,
                                 list_inputs=list_inputs)

        # Use it on test samples #
        if opt.test:
            logging.info('  Processing test output sample  '.center(80,'*'))
            if parameters.crossvalidation: # in cross validation the testing set in inside the training DF
                inst_out.OutputFromTraining(data=train_all,path_output=path_output)
            else:
                inst_out.OutputFromTraining(data=test_all,path_output=path_output)
            logging.info('')
             
   
if __name__ == "__main__":
    main()
