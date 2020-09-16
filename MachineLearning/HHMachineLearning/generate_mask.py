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
