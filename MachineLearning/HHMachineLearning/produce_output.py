import os
import copy
import sys
import logging
import numpy as np
import pandas as pd
from root_numpy import array2root

from NeuralNet import HyperModel
from import_tree import Tree2Pandas
from generate_mask import GenerateSliceIndices, GenerateSliceMask
from data_generator import DataGenerator
import parameters


class ProduceOutput:
    def __init__(self,model,generator=False,list_inputs=None):
        self.model = model              # Model of the NeuralNet
        self.list_inputs = list_inputs
        self.generator = generator
        if self.list_inputs is None:
            self.list_inputs = copy.deepcopy(parameters.inputs) 

    def OutputFromTraining(self,data,path_output,output_name=None):
        """
            Get the output of the model from the test set
            This is data separated from the training
            If output_name is specified, the whole data will be written in 'output_name'.root
                if not, the samples in the dataframe are used to split into different files with names 'sample'.root
        """
        if not self.generator:
            inputs = data[self.list_inputs].values
            if len(self.model) == 1: # classic training
                instance = HyperModel(self.model[0])
                output = instance.HyperRestore(inputs)
            else:   # cross validation
                output=None
                for model_idx,model in enumerate(self.model):
                    instance = HyperModel(model)
                    apply_idx, _, _ = GenerateSliceIndices(model_idx)
                    apply_mask = GenerateSliceMask(apply_idx,data['mask']) 
                    model_out = instance.HyperRestore(inputs[apply_mask],generator=self.generator,model_idx=model_idx,return_inputs=True)
                    if output is None:
                        output = model_out
                    else:
                        output = np.concatenate((output,model_out),axis=0)
            assert output.shape[0] == data.shape[0]
            output_df = pd.DataFrame(output,columns=[('output_%s'%o).replace('$','') for o in parameters.outputs],index=pd.RangeIndex(start=0,stop=output.shape[0]))
            full_df = pd.concat([data,output_df],axis=1)
            self.SaveToRoot(full_df,path_output,output_name)
        else:
            if len(self.model) == 1: # classic training
                output_generator = DataGenerator(path = parameters.config,
                                                 inputs = parameters.inputs,
                                                 outputs = parameters.outputs,
                                                 cut = parameters.cut,
                                                 weight  = parameters.weight,
                                                 batch_size = parameters.output_batch_size,
                                                 state_set = 'output')
                instance = HyperModel(self.model[0])
                for i in range(len(output_generator)):
                    inputs,targets,weight,other = output_generator.__getitem__(i,True)
                    output = instance.HyperRestore(inputs)
                    output_df = pd.DataFrame(output,columns=[('output_%s'%o).replace('$','') for o in parameters.outputs],index=pd.RangeIndex(start=0,stop=output.shape[0]))
                    full_df = pd.concat([inputs,targets,other,weight,output_df],axis=1)
                    self.SaveToRoot(full_df,path_output,output_name,out_idx='_slice%d'%i)
            else:   # cross validation
                output=None
                for model_idx,model in enumerate(self.model):
                    instance = HyperModel(model)
                    output_generator = DataGenerator(path = parameters.config,
                                                     inputs = parameters.inputs,
                                                     outputs = parameters.outputs,
                                                     cut = parameters.cut,
                                                     weight  = parameters.weight,
                                                     batch_size = parameters.output_batch_size,
                                                     state_set = 'output',
                                                     model_idx = model_idx)
                    for i in range(len(output_generator)):
                        inputs,targets,weight,other = output_generator.__getitem__(i,True)
                        output = instance.HyperRestore(inputs)
                        output_df = pd.DataFrame(output,columns=[('output_%s'%o).replace('$','') for o in parameters.outputs],index=pd.RangeIndex(start=0,stop=output.shape[0]))
                        full_df = pd.concat([inputs,targets,other,weight,output_df],axis=1)
                        self.SaveToRoot(full_df,path_output,output_name,out_idx='_model%d_slice%d'%(model_idx,i))

 
    def SaveToRoot(self,df,path_output,output_name=None,out_idx=''):
        # Get the unique samples as a list #
        if output_name is None:
            sample_list = list(df[parameters.split_name].unique())

            # Loop over samples #
            for sample in sample_list:
                sample_df = df.loc[df[parameters.split_name]==sample] # We select the rows corresponding to this sample

                # Remove tag and sample name (info in target as bool) #
                sample_df = sample_df.drop('tag',axis=1)
                sample_df = sample_df.drop('sample',axis=1)

                # From df to numpy array with dtype #
                sample_output = sample_df.to_records(index=False,column_dtypes='float64')
                sample_output.dtype.names = parameters.make_dtype(sample_output.dtype.names)# because ( ) and . are an issue for root_numpy
                sample_output_name = os.path.join(path_output,sample+out_idx+'.root')

                # Save as root file #
                array2root(sample_output,sample_output_name,mode='recreate')
                logging.info('Output saved as : '+sample_output_name)
        else:
            # From df to numpy array with dtype #
            full_output = df.to_records(index=False,column_dtypes='float64')
            full_output.dtype.names = parameters.make_dtype(full_output.dtype.names)# because ( ) and . are an issue for root_numpy
            full_output_name = os.path.join(path_output,output_name)
            array2root(full_output,full_output_name,mode='recreate')
            logging.info('Output saved as : '+full_output_name)

         
    def OutputNewData(self,input_dir,list_sample,path_output,variables=None):
        """
            Given a model, produce the output 
            The Network has never seen this data !
        """
        # Loop over datasets #
        logging.info('Input directory : %s'%input_dir)
        for f in list_sample: 
            name = os.path.basename(f)
            full_path = os.path.join(input_dir,f)
            logging.info('Looking at %s'%f)

            # Get the data #
            if variables is None:
                var = parameters.inputs+parameters.outputs+parameters.other_variables
            else:
                var = copy.deepcopy(variables) # Avoid bug where variables is changed at each new file

            if self.generator:
                data = None
            else:
                data = Tree2Pandas(input_file=full_path,
                                   variables=var,
                                   weight=parameters.weights,
                                   cut = parameters.cut,
                                   reweight_to_cross_section=False)
                    
                if data.shape[0]==0:
                    logging.info('\tEmpty tree')
                    continue # Avoids empty trees

            self.OutputFromTraining(data=data,path_output=path_output,output_name=name)
