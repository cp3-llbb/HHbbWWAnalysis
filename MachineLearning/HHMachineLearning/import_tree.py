import glob
import os
import sys
import logging
import json
import re
import collections
import copy

import array
import numpy as np
import pandas as pd

import parameters
from root_numpy import tree2array, rec2array
from ROOT import TChain, TFile, TTree


###############################################################################
# Tree2Pandas#
###############################################################################

def Tree2Pandas(input_file, variables, weight=None, cut=None, xsec=None, event_weight_sum=None, luminosity=None, tree_name='tree',start=None,stop=None):
    """
    Convert a ROOT TTree to a pandas DF
    """
    variables = [var for var in variables if not var.startswith("$")]
    variables = copy.copy(variables) # Otherwise will add the weight and have a duplicate branch

    # Check for repetitions in variables -> makes root_numpy crash #
    repeated_var = [item for item, count in collections.Counter(variables).items() if count > 1]
    if len(repeated_var) != 0:
        logging.critical('There are repeated variables')
        for var in repeated_var:
            logging.critical('... %s'%var)
        raise RuntimeError("Repeated arguments for importing data")

    # Get root tree, check if exists first #
    if not os.path.exists(input_file):
        logging.warning("File %s does not exist"%input_file)
        return None
    file_handle = TFile.Open(input_file)
    if not file_handle.GetListOfKeys().Contains(tree_name):
        #logging.warning("Could not find tree %s in %s"%(tree_name,input_file))
        logging.debug("Could not find tree %s in %s"%(tree_name,input_file))
        return None
    tree = file_handle.Get(tree_name)
    N = tree.GetEntries()
    logging.debug('\tNumber of events : %d'%N)

    # Read the tree and convert it to a numpy structured array
    if weight is not None:
        variables += [weight]
    data = tree2array(tree, branches=variables, selection=cut, start=start, stop=stop)

    # Convert to pandas dataframe #
    df = pd.DataFrame(data)

     # Reweighting #
    relative_weight = 1
    if weight is not None and xsec is not None and event_weight_sum is not None:
        if luminosity is None:
            luminosity = 1
        relative_weight = xsec * luminosity / event_weight_sum
        logging.debug('\t\tReweighting requested')
        logging.debug('\t\t\tCross section : %0.5f'%xsec)
        logging.debug('\t\t\tEvent weight sum : %0.2f'%event_weight_sum)
        logging.debug('\t\t\tLuminosity : %0.2f'%luminosity)
        logging.debug('\t\tRelative weight %0.3e'%relative_weight)
        df['cross_section'] = np.ones(df.shape[0])*xsec
        df['luminosity'] = np.ones(df.shape[0])*luminosity
        df['event_weight_sum'] = np.ones(df.shape[0])*event_weight_sum
   
    if weight is not None:
        df['event_weight'] = df[weight]*relative_weight
    else:
        df['event_weight'] = np.ones(df.shape[0])

    # Only part of tree #
    if start is not None or stop is not None:
        ni = start if start is not None else 0
        nf = stop if stop is not None else N
        logging.debug("Reading from {} to {} in input tree (over {} entries)".format(ni,nf,N))
    file_handle.Close()

    return df

###############################################################################
# LoopOverTrees #
###############################################################################

def LoopOverTrees(input_dir, variables, weight=None, additional_columns={}, cut=None, xsec_dict=None, event_weight_sum_dict=None, lumi_dict=None, eras=None, tree_name='tree', list_sample=None, start=None, stop=None):
    """
    Loop over ROOT trees inside input_dir and process them using Tree2Pandas.
    """
    # Check if directory #
    if not os.path.exists(input_dir):
        raise RuntimeError("%s does not exist"%input_dir)

    logging.debug("Accessing directory : "+input_dir)

    # Xsec #
    if xsec_dict is not None and not isinstance(xsec_dict,dict):
        raise NotImplementedError('Cannot handle xsec not being a dict')
    # Event weight sum #
    if event_weight_sum_dict is not None and not isinstance(event_weight_sum_dict,dict):
        raise NotImplementedError('Cannot handle event weight sum not being a dict')
    if eras is None and (xsec_dict is None or event_weight_sum_dict is None):
        raise RuntimeError('If you plan to use xsec and even weight sum you need to provide the eras (either one value or a list with one element per sample)')
        

    # Wether to use a given sample list or loop over files inside a dir #
    if list_sample is None:
        list_sample = glob.glob(os.path.join(input_dir,"*.root"))
    else:
        list_sample = [os.path.join(input_dir,s) for s in list_sample]

    # Start and stop #
    if isinstance(start,list) and isinstance(stop,list):
        if len(start) != len(list_sample):
            raise RuntimeError("Start events list does not match the list samples")
        if len(stop) != len(list_sample):
            raise RuntimeError("Stop events list does not match the list samples")

    # Loop over the files #
    first_file = True
    all_df = pd.DataFrame() 
    for i,sample in enumerate(list_sample):
        logging.debug("\tAccessing file : %s"%sample)
        sample_name = os.path.basename(sample)

        # Get era if given #
        if isinstance(eras,str):
            era = eras
        elif isinstance(eras,list):
            if len(eras) != len(list_sample):
                raise RuntimeError("List of eras has length %d and sample list has length %d"(len(eras),len(list_sample)))
            era = eras[i]

        # Cross section #
        xsec = None
        if xsec_dict is not None and sample_name in xsec_dict[era].keys():
            xsec = xsec_dict[era][sample_name]

        # Event weight sum #
        event_weight_sum = None
        if event_weight_sum_dict is not None and sample_name in event_weight_sum_dict[era].keys():
            event_weight_sum = event_weight_sum_dict[era][sample_name]

        # Luminority #
        luminosity = None
        if lumi_dict is not None:
            luminosity = lumi_dict[era]

        # Start #
        ni = None
        if start is not None:
            if isinstance(start,list):
                ni = start[i]
            else:
                ni = start

        # Stop #
        nf = None
        if stop is not None:
            if isinstance(stop,list):
                nf = stop[i]
            else:
                nf = stop

        # Get the data as pandas df #
        df = Tree2Pandas(input_file                 = sample,
                         variables                  = variables,
                         weight                     = weight,
                         cut                        = cut,
                         xsec                       = xsec,
                         event_weight_sum           = event_weight_sum,
                         luminosity                 = luminosity,
                         tree_name                  = tree_name,
                         start                      = ni,
                         stop                       = nf) 
        if df is None:
            continue

        # Register sample name #
        df['sample'] = pd.Series([sample_name.replace('.root','')]*df.shape[0])
        
        # Register additional columns #
        if len(additional_columns.keys()) != 0:
            for key,val in additional_columns.items():
                if isinstance(val,list):
                    if len(val) != len(list_sample):
                        raise RuntimeError('Value list %s you want to add has len %d while there are %d samples'%(key,len(val),len(list_sample)))
                    df[key] = pd.Series([val[i]]*df.shape[0])
                else:
                    df[key] = pd.Series([val]*df.shape[0])

        # Concatenate into full df #
        if first_file:
            all_df = df
            first_file = False
        else:
            all_df = pd.concat([all_df,df])
        all_df = all_df.reset_index(drop=True) # Otherwise there will be an index repetition for each file
    return all_df
