import os
import sys
import yaml
import math
import copy
import argparse
import random
import logging
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import tensorflow as tf
from lbn import LBNLayer
from sklearn.metrics import f1_score
from sklearn.preprocessing import LabelEncoder, OneHotEncoder

sys.path.append('..')
import parameters
import Operations
from import_tree import LoopOverTrees
from input_plots import InputPlots

class Analyze:
    def __init__(self,suffix,N=None):
        self.inputs      = parameters.inputs
        self.LBN_inputs  = parameters.LBN_inputs
        self.outputs     = parameters.outputs
        self.mask_onehot = parameters.mask_op
        self.config      = parameters.config
        self.path_out    = os.path.join(os.path.abspath("."),suffix)
        self.N           = N
        self.suffix      = suffix
        self.cacheName   = 'cache_df.pkl'

        self.custom_objects =  {'LBNLayer':LBNLayer} 
        self.custom_objects.update({name:getattr(Operations,name) for name in dir(Operations) if name.startswith('op')})

        if not os.path.exists(self.path_out):
            os.makedirs(self.path_out)

        self.getDF()

    def getDF(self):
        variables = list(set(self.inputs+self.outputs+self.LBN_inputs))
        with open (self.config,'r') as f:
            sampleConfig = yaml.load(f)

        if not os.path.exists(self.cacheName):
            if self.N is None:
                raise RuntimeError("You need to specify a number of events")
            # Import arrays #
            data_dict = {}
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
                                                  list_sample               = list_sample,
                                                  cut                       = parameters.cut,
                                                  eras                      = era,
                                                  tree_name                 = parameters.tree_name,
                                                  additional_columns        = {'tag':node,'era':era},
                                                  stop                      = self.N) 
                    if data_node is None:
                        data_node = data_node_era
                    else:
                        data_node = pd.concat([data_node,data_node_era],axis=0)
                data_dict[node] = data_node
            self.data = pd.concat(data_dict.values(),copy=True).reset_index(drop=True)
            self.saveDF()
        else:
            self.loadDF()
            self.N = self.data.shape[0]

        self.inputs  = [inp.replace('$','') for inp in parameters.inputs]
        self.outputs = [out.replace('$','') for out in parameters.outputs]

        # Make one hot encoding of the output #

        # Add target #
        label_encoder = LabelEncoder()
        onehot_encoder = OneHotEncoder(sparse=False)
        label_encoder.fit(self.outputs)
        # From strings to labels #
        integers = label_encoder.transform(self.data['tag']).reshape(-1, 1)
        # From labels to strings #
        onehot_encoder.fit(np.arange(len(self.outputs)).reshape(-1, 1))
        onehot = onehot_encoder.transform(integers)
        # From arrays to pd DF #
        onehot_df = pd.DataFrame(onehot,columns=label_encoder.classes_,index=self.data.index)
        # Add to full #
        self.data = pd.concat([self.data,onehot_df],axis=1)
            
        logging.info("... Number of events : %d"%self.data.shape[0])

    def saveDF(self):
        self.data.to_pickle(self.cacheName)
        logging.info("Saved data to cache")

    def loadDF(self):
        self.data = pd.read_pickle(self.cacheName)
        logging.info("Loaded data from cache")

    def makeCorrelationMatrix(self):
        logging.info("Starting the correlation matrix plot")
        # Compute the correlation matrix
        inputs = [self.inputs[i] for i in range(len(self.inputs)) if not self.mask_onehot[i]]
        corr = self.data[inputs].astype(np.float32).corr()
        corr = corr.fillna(0.)

        # Generate a mask for the upper triangle
        mask = np.zeros_like(corr, dtype=np.bool)
        mask[np.triu_indices_from(mask)] = True

        # Set up the matplotlib figure
        f, ax = plt.subplots(figsize=(11, 9))
        f.subplots_adjust(left=0.1, right=1., top=0.95, bottom=0.15)

        # Generate a custom diverging colormap
        cmap = sns.diverging_palette(250, 15, s=75, l=40, n=10, center='light', as_cmap=True)

        # Draw the heatmap with the mask and correct aspect ratio
        g = sns.heatmap(corr, mask=mask, cmap=cmap,  center=0, xticklabels=True, yticklabels=True,
                    square=True, linewidths=.5)

        g.set_yticklabels(g.get_ymajorticklabels(), size = int(50/math.sqrt(corr.shape[1])))
        g.set_xticklabels(g.get_xmajorticklabels(), size = int(50/math.sqrt(corr.shape[0])))

        name = os.path.join(self.path_out,'correlation_matrix.png')
        plt.savefig(name)
        logging.info("... saved as %s"%name)
#        import IPython
#        IPython.embed()

    def inputPlots(self):
        InputPlots(self.data,self.inputs,self.path_out)

    def loadModel(self,basemodel):
        # Need to find basemodel.json and basemodel.h5 #
        self.basemodel = basemodel
        if os.path.exists(basemodel+'.zip'):
            logging.info("Loading model from zip file %s.zip"%basemodel)
            self.model = talos.Restore(basemodel+'.zip',custom_objects=self.custom_objects).model 
        elif os.path.exists(basemodel+'.json') and os.path.exists(basemodel+'.h5'):
            with open(basemodel+'.json','r') as f:
                loaded_model_json = f.read()
            self.model = tf.keras.models.model_from_json(loaded_model_json,custom_objects=self.custom_objects)
            self.model.load_weights(basemodel+".h5")
        else:
            raise RuntimeError('Could not find the suitable format for %s'%basemodel)

    def PerformutationFeatureImportance(self):
        data = self.data.sample(n=100000,axis=0)
        inputsLL = np.hsplit(data[self.inputs].astype(np.float32).values,len(self.inputs))
        has_LBN = any([l.__class__.__name__ == 'LBNLayer' for l in self.model.layers])
        if has_LBN:
            inputsLBN = data[self.LBN_inputs].astype(np.float32).values.reshape(-1,4,len(self.LBN_inputs)//4)
            inputs = (inputsLL,inputsLBN)
        else:
            inputs = inputsLL
        logging.info("Producing output for true score")
        outputs = self.model.predict(inputs,batch_size=10000,verbose=0)
        targets = data[self.outputs].values
        true_F1_score = f1_score(targets.argmax(1),outputs.argmax(1),average='macro')

        permutations = 10
        f1_scores = np.zeros(len(self.inputs))
        f1_scores_err = np.zeros(len(self.inputs))
        logging.info("Producing output for permutation scores")
        for idxPerm,inputName in enumerate(self.inputs):
            logging.info("Looking at input %d/%d"%(idxPerm,len(self.inputs)))
            inputs_perm = copy.deepcopy(inputsLL)
            perm_f1_scores = []
            for perm in range(permutations):
                logging.info("... Permutation %d/%d"%(perm,permutations))
                np.random.shuffle(inputs_perm[idxPerm])
                if has_LBN:
                    outputs_perm = self.model.predict((inputs_perm,inputsLBN),batch_size=10000,verbose=0)
                else:
                    outputs_perm = self.model.predict(inputs_perm,batch_size=10000,verbose=0)
                perm_f1_scores.append(f1_score(targets.argmax(1),outputs_perm.argmax(1),average='macro'))
            perm_f1_scores = np.array(perm_f1_scores)
            f1_scores[idxPerm] = abs(perm_f1_scores.mean()-true_F1_score)
            f1_scores_err[idxPerm] = perm_f1_scores.std()
        
        idxSort = np.flip(np.argsort(f1_scores))
        inputNames = np.array(self.inputs,dtype=object)[idxSort]
        f1_scores = f1_scores[idxSort]
        f1_scores_err = f1_scores_err[idxSort]

        fig, ax = plt.subplots(figsize=(8,15))
        plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.1)
        y_pos = np.arange(len(self.inputs))
        ax.barh(y_pos,f1_scores,xerr=f1_scores_err,align='center')
        ax.set_yticks(y_pos)
        ax.set_yticklabels(inputNames,size=int(80/math.sqrt(len(self.inputs))))
        ax.invert_yaxis()  # labels read top-to-bottom
        ax.set_xlabel('Importance',fontsize=18)
        name = os.path.join(self.path_out,'feature_importance.png')
        fig.savefig(name)
        logging.info("... saved as %s"%name)


if __name__ == "__main__":
    # Argparse #
    parser = argparse.ArgumentParser(description='Analyzer : variables and/or network')

    # Main args #
    parser.add_argument('-N','--number', action='store', required=False, type=int, default=0,
        help='Number of events max per file')
    parser.add_argument('--suffix', action='store', required=True, type=str, default = '',
        help='Suffix for the png plots (default is "")')
    parser.add_argument('--model', action='store', required=False, type=str, default = None,
        help='DNN model')
    parser.add_argument('--correlation', action='store_true', required=False, default=False,
        help='Produce correlation matrix plot')
    parser.add_argument('--input_plots', action='store_true', required=False, default=False,
        help='Produce input feature plots')
    parser.add_argument('-v','--verbose', action='store_true', required=False, default=False,
        help='Verbose mode')

    args = parser.parse_args()

    # Verbose #
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG,format='%(asctime)s - %(levelname)s - %(message)s',datefmt='%m/%d/%Y %H:%M:%S')
    else:
        logging.basicConfig(level=logging.INFO,format='%(asctime)s - %(levelname)s - %(message)s',datefmt='%m/%d/%Y %H:%M:%S')

    # Instance #
    instance = Analyze(suffix=args.suffix,N=args.number)
    if args.correlation:
        instance.makeCorrelationMatrix()
    if args.input_plots:
        instance.inputPlots()
    if args.model is not None:
        instance.loadModel(args.model)
        instance.PerformutationFeatureImportance()
        
