# Defines the parameters that users might need to change
# Must be included manually in each script

# /!\ Two types of files need to be filled with users parameters 
#           - parameters.py (mort important one)
#           - sampleList.py (on what samples to run)
#           (optionnaly NeuralNet.py for early_stopping etc)
import multiprocessing
import os
from tensorflow.keras.losses import binary_crossentropy, mean_squared_error, logcosh, categorical_crossentropy
from tensorflow.keras.optimizers import RMSprop, Adam, Nadam, SGD            
from tensorflow.keras.activations import relu, elu, selu, softmax, tanh, sigmoid
from tensorflow.keras.regularizers import l1,l2 

from grouped_entropy import GroupedXEnt

##################################  Path variables ####################################

main_path = os.path.abspath(os.path.dirname(__file__))
path_out = '/nfs/scratch/fynu/gsaha/HHMachineLearning_output_res/'
path_model = os.path.join(main_path,'model')

##############################  Datasets proportion   #################################
crossvalidation = True

# Classic training #
# -> For crossvalidation == False
training_ratio = 0.7    # Training set sent to keras
evaluation_ratio = 0.1  # Evaluation set sent to autom8
output_ratio = 0.2      # Output set for plotting later
assert training_ratio + evaluation_ratio + output_ratio == 1 
# Will only be taken into account for the masks generation, ignored after

# Cross-validation #
# -> For crossvalidation == True
N_models = 5  # Number of models to train
N_train  = 3  # Number of slices on which to train
N_eval   = 1  # Number of slices on which to evaluate the model
N_apply  = 1  # Number of slices on which the model will be applied for uses in analysis
N_slices = N_train+N_eval+N_apply
splitbranch = 'event' # Will split the training based on "event % N_slices" 
# /!\ Don't forget to add "splitbranch" in other_variables

if N_slices % N_models != 0: # will not work otherwise
    raise RuntimeError("N_slices [{}] % N_models [{}] should be == 0".format(N_slices,N_models))
if N_apply != N_slices/N_models: # Otherwise the same slice can be applied on several networks
    raise RuntimeError("You have asked {} models that should be applied {} times, the number of slices should be {} but is {}+{}+{}".format(N_models,N_apply,N_models*N_apply,N_train,N_eval,N_apply))

############################### Slurm parameters ######################################
partition = 'gpu'  # Def, cp3, cp3-gpu or gpu
QOS = 'normal' # cp3, normal, cp3-gpu
time = '0-04:00:00' # days-hh:mm:ss
mem = '90000' # ram in MB
tasks = 1 # Number of threads(as a string) (not parallel training for classic mode)
cpus = 8
gpus = 1
additional_options = ""
workers = 20
split_per_model = True # in case of cross validation, to send one job per model or not

##################################  Naming ######################################
# Physics Config #
config = os.path.join(os.path.abspath(os.path.dirname(__file__)),'sampleListDL_res.yml')
lumidict = {2016:35922,2017:41529.152060112,2018:59740.565201546}
eras = [2016,2017,2018] # To enable or disable eras, add or remove from this list

categories = ["Resolved1Btag","Resolved2Btag","Boosted1Btag"]
channels = ['ElEl','MuMu','ElMu']

# Better put them in alphabetical order
nodes = ['DY','GGF','H','Rare','ST','TT','TTVX','VVV']

weight_groups = [
                  (1.0, ('DY')),
                  (1.0, ('GGF')),
                  (1.0, ('H')),
                  (1.0, ('Rare')),
                  (1.0, ('ST')),
                  (1.0, ('TT')),
                  (1.0, ('TTVX')),
                  (1.0, ('VVV')),
                ]


quantile = 0.95 # repeat weights with too high learning weights

# Input plots options #
#node_colors = {
#            'Ewk'  : '#610596',
#            'GGF'   : '#288a24',
#            'H'     : '#06b894',
#            'Top'    : '#cc7a16',
#            'WJets' : '#d95564',
#             }

# Tree name #
tree_name = 'Events'

# scaler and mask names #
suffix = 'resonant' 
# scaler_name -> 'scaler_{suffix}.pkl'  If does not exist will be created 
# mask_name -> 'mask_{suffix}_{sample}.npy'  If does not exist will be created 
scaler_name = 'scaler_'+suffix+'_'.join([str(era) for era in eras])+'.pkl'
scaler_path = os.path.join(main_path,scaler_name)

# Data cache #
train_cache = os.path.join(path_out,'train_cache_'+suffix+'_'+'_'.join([str(era) for era in eras])+'.pkl')
test_cache = os.path.join(path_out,'test_cache_'+suffix+'_'+'_'.join([str(era) for era in eras])+'.pkl')

# Meta config info #
xsec_json = os.path.join(main_path,'xsec_{era}.json')
event_weight_sum_json = os.path.join(main_path,'event_weight_sum_{era}.json')

# Training resume #
resume_model = ''

# Output #
output_batch_size = 100000
split_name = 'tag' # 'sample' or 'tag' : criterion for output file splitting

##############################  Evaluation criterion   ################################

eval_criterion = "eval_error" # either val_loss or eval_error or val_acc

##############################  Model callbacks ################################
# Early stopping to stop learning after some criterion 
early_stopping_params = {'monitor'   : 'val_loss',  # Value to monitor
                         'min_delta' : 0.0001,          # Minimum delta to declare an improvement
                         'patience'  : 50,          # How much time to wait for an improvement
                         'verbose'   : 1,           # Verbosity level
                         'restore_best_weights': False,
                         'mode'      : 'min'}       # Mode : 'auto', 'min', 'max'

# Reduce LR on plateau : if no improvement for some time, will reduce lr by a certain factor
reduceLR_params = {'monitor'    : 'val_loss',   # Value to monitor
                   'factor'     : 0.5,          # Multiplicative factor by which to multiply LR
                   'min_lr'     : 1e-5,         # Minimum value for LR
                   'min_delta'  : 0.0001,       # Minimum delta to declare an improvement
                   'patience'   : 20,            # How much time to wait for an improvement
                   'cooldown'   : 0,            # How many epochs before starting again to monitor
                   'verbose'    : 1,            # Verbosity level
                   'mode'      : 'min'}         # Mode : 'auto', 'min', 'max'


    
#################################  Scan dictionary   ##################################
# /!\ Lists must always contain something (even if 0), otherwise 0 hyperparameters #
# Classification #
#p = { 
#    'lr' : [0.01,0.001], 
#    'first_neuron' : [64,128,256,512,1024],
#    'activation' : [relu],
#    'dropout' : [0.,0.1,0.2],
#    'hidden_layers' : [2,3,4,5,6], # does not take into account the first layer
#    'output_activation' : [softmax],
#    'l2' : [0,0.01,0.1],
#    'optimizer' : [Adam],  
#    'epochs' : [50],   
#    'batch_size' : [10000], 
#    'n_particles' : [16],
#    'loss_function' : [grouped_loss] , #  [categorical_crossentropy]
#}
p = { 
    'lr' : [0.01], 
    'first_neuron' : [512],
    'activation' : [relu],
    'dropout' : [0.01],
    'hidden_layers' : [4], # does not take into account the first layer
    'output_activation' : [softmax],
    'l2' : [1e-3],
    'optimizer' : [Adam],  
    'epochs' : [500],   
    'batch_size' : [50000], 
    'n_particles' : [10],
    'loss_function' : [categorical_crossentropy],
}


repetition = 1 # How many times each hyperparameter has to be used 

###################################  Variables   ######################################

cut = 'total_weight > 0 && MC_weight <10000 && total_weight < 10e10'
#cut = 'MC_weight > 0'

weight = 'total_weight'
#weight = None

# Input branches (combinations possible just as in ROOT #
#/!\ onehot variables need to be at the beginning of the list (checked later)
inputs = [
            ###########################
            ####----- Resolved ----####
            ###########################
            # Onehot #
            '$era@op_era',
            'l1_charge@op_charge',
            'l2_charge@op_charge',
            'l1_pdgId@op_pdgid',
            'l2_pdgId@op_pdgid',
            # Parameter #
            '$param',
            # LL variables #
            'l1_E',          
            'l1_Px',    
            'l1_Py',    
            'l1_Pz',    
            'l2_E',     
            'l2_Px',    
            'l2_Py',    
            'l2_Pz',    
            'j1_E',
            'j1_Px',
            'j1_Py', 
            'j1_Pz', 
            #'j1_btag',
            'j2_E',   
            'j2_Px',  
            'j2_Py',  
            'j2_Pz',  
            #'j2_btag',
            'j3_E',   
            'j3_Px',  
            'j3_Py',  
            'j3_Pz',  
            #'j3_btag',
            'j4_E',   
            'j4_Px',  
            'j4_Py',  
            'j4_Pz',  
            #'j4_btag',
            'fatjet_E',       
            'fatjet_Px',      
            'fatjet_Py',      
            'fatjet_Pz',      
            'fatjet_tau1',    
            'fatjet_tau2',    
            'fatjet_tau3',    
            'fatjet_tau4',    
            'fatjet_softdrop',
            'met_E', 
            'met_Px',
            'met_Py',
            'met_Pz',
            # HL variables #
            'm_bb_bregcorr',
            'ht', 
            'min_dr_jets_lep1',       
            'min_dr_jets_lep2',       
            'min_deta_jets_lep1',       
            'min_deta_jets_lep2',       
            'm_ll',                   
            'dr_ll',                  
            'deta_ll',                  
            'min_dr_jet',             
            'min_dphi_jet',           
            'min_deta_jet',
            'm_hh_simplemet_bregcorr',
            'met_ld',                 
            'dr_bb',                  
            'min_dr_leps_b1',         
            'min_dr_leps_b2',         
            'lep1_conept',            
            'lep2_conept',            
            'mww_simplemet',          
            #'vbf_tag',                
            'boosted_tag',            
            'dphi_met_dilep',         
            'dphi_met_dibjet',        
            'dr_dilep_dijet',         
            'dr_dilep_dibjet',        
            'm_T',
            'dr_lep1_dibjet',
            'dr_lep2_dibjet',
            'cosThetaS_Hbb',
    ]

operations = [inp.split('@')[1] if '@' in inp else None  for inp  in  inputs]
check_op = [(o is not None)*1 for o in operations]
if check_op != sorted(check_op,reverse=True):
    raise RuntimeError('Onehot inputs need to be at the beginning of the inputs list')
mask_op = [len(inp.split('@'))==2 for inp  in  inputs]
inputs = [inp.split('@')[0] for inp  in  inputs]

#### Resolved ####
LBN_inputs = [
              'l1_E','l1_Px','l1_Py','l1_Pz',
              'l2_E','l2_Px','l2_Py','l2_Pz',
              'j1_E','j1_Px','j1_Py','j1_Pz',
              'j2_E','j2_Px','j2_Py','j2_Pz',
              'fatjet_E','fatjet_Px','fatjet_Py','fatjet_Pz',
             ]

#### Boosted ####
#LBN_inputs = [
#              'lep_E','lep_Px','lep_Py','lep_Pz',
#              'fatbj_E','fatbj_Px','fatbj_Py','fatbj_Pz',
#              'wj1_E','wj1_Px','wj1_Py','wj1_Pz',
#              'wj2_E','wj2_Px','wj2_Py','wj2_Pz',
#             ]
#


#### /!\ format = [E,Px,Py,Pz] per particle


assert len(LBN_inputs)%4 == 0


# Output branches #

outputs = [
            '$DY',
            '$GGF',
            '$H',
            '$Rare',
            '$ST',
            '$TT',
            '$TTVX',
            '$VVV',
          ] 

# Other variables you might want to save in the tree #
other_variables = [
                    'event',
                    'ls',
                    'run',
                    'MC_weight',
                 ]


################################  dtype operation ##############################

def make_dtype(list_names): # root_numpy does not like . and ()
    list_dtype = [(name.replace('.','_').replace('(','').replace(')','').replace('-','_minus_').replace('*','_times_')) for name in list_names]
    return list_dtype
        



                                



