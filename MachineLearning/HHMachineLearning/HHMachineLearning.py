#!/usr/bin/env python

import re
import math
import glob
import csv
import os
import sys
import pprint
import logging
import copy
import pickle
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
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import OneHotEncoder

# Personal files #
    
def get_options():
    """
    Parse and return the arguments provided by the user.
    """
    parser = argparse.ArgumentParser(description='HHMachineLearning')

    # Scan, deploy and restore arguments #
    a = parser.add_argument_group('Scan, deploy and restore arguments')
    a.add_argument('-s','--scan', action='store', required=False, type=str, default='',
        help='Name of the scan to be used (modify scan parameters in NeuralNet.py)')
    a.add_argument('-task','--task', action='store', required=False, type=str, default='',
        help='Name of dict to be used for scan (Used by function itself when submitting jobs or DEBUG)')
    a.add_argument('--generator', action='store_true', required=False, default=False, 
        help='Wether to use a generator for the neural network')
    a.add_argument('--resume', action='store_true', required=False, default=False,
        help='Wether to resume the training of a given model (path in parameters.py)')

    # Splitting and submitting jobs arguments #
    b = parser.add_argument_group('Splitting and submitting jobs arguments')
    b.add_argument('-split','--split', action='store', required=False, type=int, default=0,
        help='Number of parameter sets per jobs to be used for splitted training for slurm submission (if -1, will create a single subdict)')
    b.add_argument('-submit','--submit', action='store', required=False, default='', type=str,
        help='Wether to submit on slurm and name for the save (must have specified --split)')
    b.add_argument('-resubmit','--resubmit', action='store', required=False, default='', type=str,
        help='Wether to resubmit failed jobs given a specific path containing the jobs that succeded')
    b.add_argument('-debug','--debug', action='store_true', required=False, default=False,
        help='Debug mode of the slurm submission, does everything except submit the jobs')

    # Analyzing or producing outputs for given model (csv or zip file) #
    c = parser.add_argument_group('Analyzing or producing outputs for given model (csv or zip file)')
    c.add_argument('-r','--report', action='store', required=False, type=str, default='',
        help='Name of the csv file for the reporting (without .csv)')
    c.add_argument('-m','--model', action='store', required=False, type=str, default='',                                                                                                          
        help='Loads the provided model name (without .zip and type, it will find them)') 
    c.add_argument('--test', action='store_true', required=False, default=False,
        help='Applies the provided model (do not forget -o) on the test set and output the tree') 
    c.add_argument('-o','--output', action='store', required=False, nargs='+', type=str, default=[], 
        help='Applies the provided model (do not forget -o) on the list of keys from sampleList.py (separated by spaces)') 

    # Concatenating csv files arguments #
    d = parser.add_argument_group('Concatenating csv files arguments')
    d.add_argument('-csv','--csv', action='store', required=False, type=str, default='',
        help='Wether to concatenate the csv files from different slurm jobs into a main one, \
              please provide the path to the csv files')

    # Concatenating csv files arguments #
    e = parser.add_argument_group('Physics arguments')
    e.add_argument('--resolved1b', action='store_true', required=False, default=False,
       help='Resolved topology with 1 bjet')
    e.add_argument('--resolved2b', action='store_true', required=False, default=False,
       help='Resolved topology with 2 bjets')
    e.add_argument('--boosted', action='store_true', required=False, default=False,
       help='Boosted topology')


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
    from generate_mask import GenerateMask
    from split_training import DictSplit
    from concatenate_csv import ConcatenateCSV
    from sampleList import samples_dict_2016, samples_dict_2017, samples_dict_2018, samples_path
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
        args = ' ' # Do not forget the spaces after each arg!
        if opt.resolved1b:          args += ' --resolved1b '
        if opt.resolved2b:          args += ' --resolved2b '
        if opt.boosted:             args += ' --boosted '
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

    lumidict = {'2016':35922,'2017':41529.152060112,'2018':59740.565201546}

    if opt.nocache:
        logging.warning('No cache will be used not saved')
    if os.path.exists(parameters.train_cache) and os.path.exists(parameters.test_cache) and not opt.nocache:
        logging.info('Will load data from cache')
        train_all = pd.read_pickle(parameters.train_cache)
        test_all = pd.read_pickle(parameters.test_cache)
    else:
        # Import arrays #
        nodes = ['TT','DY','ST','VVV','TTVX','H','Rare','HH']
        channels = ['ElEl','MuMu','ElMu']
        data_dict = {}
        for node in nodes:
            list_sample = []
            if opt.resolved1b:
                strSelect = ['resolved_1b_{}_{}'.format(channel,node) for channel in channels]
            if opt.resolved2b:
                strSelect = ['resolved_2b_{}_{}'.format(channel,node) for channel in channels]
            if opt.boosted:
                strSelect = ['boosted_{}_{}'.format(channel,node) for channel in channels]

            data_node = None

            for era,samples_dict in zip(['2016','2017','2018'],[samples_dict_2016,samples_dict_2017,samples_dict_2018]):
                if len(samples_dict.keys())==0:
                    logging.info('Sample dict for era {} is empty'.format(era))
                    continue
                if node != 'HH':
                    xsec_json = parameters.xsec_json.format(era=era)
                    event_weight_sum_json = parameters.event_weight_sum_json.format(era=era)
                else:
                    xsec_json = None
                    event_weight_sum_json = None
                list_sample = [sample for key in strSelect for sample in samples_dict[key]]
                data_node_era = LoopOverTrees(input_dir                 = samples_path,
                                              variables                 = variables,
                                              weight                    = parameters.weights,
                                              list_sample               = list_sample,
                                              cut                       = parameters.cut,
                                              xsec_json                 = xsec_json,
                                              event_weight_sum_json     = event_weight_sum_json,
                                              luminosity                = lumidict[era],
                                              additional_columns        = {'tag':node,'era':era})
                if data_node is None:
                    data_node = data_node_era
                else:
                    data_node = pd.concat([data_node,data_node_era],axis=0)
                logging.info('{} Sample size for era {} : {}'.format(node,era,data_node_era.shape[0]))
            data_dict[node] = data_node
            logging.info('{} Sample size for all eras : {}'.format(node,data_node.shape[0]))

        # Weight equalization #
        if parameters.weights is not None:
            for node, data in data_dict.items():
                sum_weights = data[parameters.weights].sum()
                logging.info('Sum of weight for %s samples : %.2e'%(node,sum_weights))
                data['learning_weights'] = data[parameters.weights]/sum_weights*1e5
                logging.info('\t -> After equalization : %0.2e'%data['learning_weights'].sum())

        # Data splitting #
        train_dict = {}
        test_dict = {}
        for node,data in data_dict.items():
            mask = GenerateMask(data.shape[0],parameters.suffix+'_'+node)
            try:
                train_dict[node] = data[mask==True]
                test_dict[node]  = data[mask==False]
            except ValueError:
                logging.critical("Problem with the mask you imported, has the data changed since it was generated ?")
                raise ValueError
            
        del data_dict
         
        train_all = pd.concat(train_dict.values(),copy=True).reset_index(drop=True)
        test_all  = pd.concat(test_dict.values(),copy=True).reset_index(drop=True)
        del train_dict, test_dict
        #logging.info('Current memory usage : %0.3f GB'%(pid.memory_info().rss/(1024**3)))

        # Randomize order, we don't want only one type per batch #
        random_train = np.arange(0,train_all.shape[0]) # needed to randomize x,y and w in same fashion
        np.random.shuffle(random_train) # Not needed for testing
        train_all = train_all.iloc[random_train]
          
        # Add target #
        label_encoder = LabelEncoder()
        onehot_encoder = OneHotEncoder(sparse=False)
        label_encoder.fit(train_all['tag'])
        # From strings to labels #
        train_integers = label_encoder.transform(train_all['tag']).reshape(-1, 1)
        test_integers = label_encoder.transform(test_all['tag']).reshape(-1, 1)
        # From labels to strings #
        train_onehot = onehot_encoder.fit_transform(train_integers)
        test_onehot = onehot_encoder.fit_transform(test_integers)
        # From arrays to pd DF #
        train_cat = pd.DataFrame(train_onehot,columns=label_encoder.classes_,index=train_all.index)
        test_cat = pd.DataFrame(test_onehot,columns=label_encoder.classes_,index=test_all.index)
        # Add to full #
        train_all = pd.concat([train_all,train_cat],axis=1)
        test_all = pd.concat([test_all,test_cat],axis=1)

        # Preprocessing #
        # The purpose is to create a scaler object and save it
        # The preprocessing will be implemented in the network with a custom layer
        if opt.scan!='': # If we don't scan we don't need to scale the data
            MakeScaler(train_all,list_inputs) 

        if not opt.nocache:
            train_all.to_pickle(parameters.train_cache)
            test_all.to_pickle(parameters.test_cache)
            logging.info('Data saved to cache')
     
    list_inputs  = [var.replace('$','') for var in parameters.inputs]
    list_outputs = [var.replace('$','') for var in parameters.outputs]
    logging.info("Sample size seen by network : %d"%train_all.shape[0])
    logging.info("Sample size for the output  : %d"%test_all.shape[0])
    #logging.info('Current memory usage : %0.3f GB'%(pid.memory_info().rss/(1024**3)))

    #############################################################################################
    # DNN #
    #############################################################################################
    if opt.GPU:
        # Start the GPU monitoring thread #
        thread = utilizationGPU(print_time = 900,
                                print_current = False,
                                time_step=0.01)
        thread.start()

    if opt.scan != '':
        instance = HyperModel(opt.scan)
        instance.HyperScan(data=train_all,
                           list_inputs=list_inputs,
                           list_outputs=list_outputs,
                           task=opt.task,
                           generator=opt.generator,
                           resume=opt.resume)
        instance.HyperDeploy(best='eval_error')

    if opt.GPU:
        # Closing monitor thread #
        thread.stopLoop()
        thread.join()
        
    if opt.model!='': 
        # Make path #
        output_name = "test" 
        path_output = os.path.join(parameters.path_out,opt.model,output_name)
        if not os.path.exists(path_output):
            os.makedirs(path_output)

        # Instance of output class #
        inst_out = ProduceOutput(model=os.path.join(parameters.main_path,'model',opt.model),
                                 generator=opt.generator,
                                 list_inputs=list_inputs)

        # Use it on test samples #
        if opt.test:
            logging.info('  Processing test output sample  '.center(80,'*'))
            inst_out.OutputFromTraining(data=test_all,path_output=path_output)
            logging.info('')
             
   
if __name__ == "__main__":
    main()
