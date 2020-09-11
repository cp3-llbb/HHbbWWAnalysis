import numpy as np
import os
import logging

import parameters                                                                                                                                                                                      

def GenerateMask(N,name):
    path_mask = os.path.join(parameters.main_path,'mask_'+name)
    if not os.path.exists(path_mask+'.npy'):                     
        mask = np.full((N,), False, dtype=bool)     
        size = parameters.training_ratio+parameters.validation_ratio
        mask[:int(size*N)] = True                         
        np.random.shuffle(mask)                                     
        np.save(path_mask,mask)                                 
        # False => Evaluation set, True => Training set           
        logging.info('Mask not found at '+path_mask+' -> Has been generated')
    else:                                                        
        mask = np.load(path_mask+'.npy')     
        logging.info('Mask found at '+path_mask+'.npy')

    return mask
