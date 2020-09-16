import numpy as np
import os
import logging

import parameters

def GenerateMask(N,name):
    path_mask = os.path.join(parameters.main_path,'mask_'+name)
    if not os.path.exists(path_mask+'.npy'):                     
        mask = np.full((N,), False, dtype=bool)     
        size = parameters.training_ratio+parameters.evaluation_ratio
        mask[:int(size*N)] = True                         
        # False => Evaluation set, True => Training set           
        np.random.shuffle(mask)                                     
        # Save #
        np.save(path_mask,mask)                                 
        logging.info('Mask not found at '+path_mask+' -> Has been generated')
    else:                                                        
        mask = np.load(path_mask+'.npy')     
        logging.info('Mask found at '+path_mask+'.npy')

    return mask

def GenerateSliceIndices(model_idx):
    Nm = parameters.N_models
    assert model_idx < Nm
    Nt = parameters.N_train
    Ne = parameters.N_eval
    Na = parameters.N_apply
    Ns = parameters.N_slices
    model_idx *= Ns/Nm
    apply_idx = [int(model_idx+i) for i in range(0,Na)]
    eval_idx = [int ((apply_idx[-1]+1+i)%Ns) for i in range(0,Ne)]
    train_idx = [i for i in range(0,Ns) if i not in apply_idx and not i in eval_idx]
    return apply_idx,eval_idx,train_idx

def GenerateSliceMask(slices,mask):
    selector = np.full((mask.shape[0]), False, dtype=bool)
    for s in slices:
        selector = selector | (mask == s)
    return selector

