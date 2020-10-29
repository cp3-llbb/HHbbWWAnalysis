# Defines the parameters that users might need to change
# Must be included manually in each script

# /!\ Two types of files need to be filled with users parameters 
#           - parameters.py (mort important one)
#           - sampleList.py (on what samples to run)
#           (optionnaly NeuralNet.py for early_stopping etc)
import multiprocessing
import os
from keras.losses import binary_crossentropy, mean_squared_error, logcosh, cosine_proximity, categorical_crossentropy
from keras.optimizers import RMSprop, Adam, Nadam, SGD            
from keras.activations import relu, elu, selu, softmax, tanh, sigmoid
from keras.regularizers import l1,l2 

##################################  Path variables ####################################

main_path = os.path.abspath(os.path.dirname(__file__))
path_out = '/home/ucl/cp3/fbury/scratch/HHMachineLearning_output/'
path_model = os.path.join(main_path,'model')

##############################  Datasets proportion   #################################
crossvalidation = False

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
partition = 'Def'  # Def, cp3 or cp3-gpu
QOS = 'normal' # cp3 or normal
time = '0-08:00:00' # days-hh:mm:ss
mem = '9000' # ram in MB
tasks = '1' # Number of threads(as a string) (not parallel training for classic mode)
workers = 20

##################################  Naming ######################################
# Physics Config #
config = os.path.join(os.path.abspath(os.path.dirname(__file__)),'sampleListSL.yml')
lumidict = {'2016':35922,'2017':41529.152060112,'2018':59740.565201546}
eras = ['2016']
#eras = ['2016','2017','2018'] # To enable or disable eras, add or remove from this list

categories = ['resolved1b3j','resolved2b2j']

# Better put them in alphabetical order
nodes = ['DY', 'H', 'HH', 'Rare', 'ST', 'TT', 'TTVX', 'VVV','WJets']
channels = ['El','Mu']

# Tree name #
tree_name = 'Events'

# scaler and mask names #
suffix = 'resolved' 
# scaler_name -> 'scaler_{suffix}.pkl'  If does not exist will be created 
# mask_name -> 'mask_{suffix}_{sample}.npy'  If does not exist will be created 

# Data cache #
train_cache = os.path.join(path_out,'train_cache.pkl' )
test_cache = os.path.join(path_out,'test_cache.pkl' )

# Meta config info #
xsec_json = os.path.join(main_path,'background_{era}_xsec.json')
event_weight_sum_json = os.path.join(main_path,'background_{era}_event_weight_sum.json')

# Training resume #
resume_model = ''

# Output #
output_batch_size = 1000
split_name = 'tag' # 'sample' or 'tag' : criterion for output file splitting

##############################  Evaluation criterion   ################################

eval_criterion = "eval_error" # either val_loss or eval_error or val_acc

##############################  Model callbacks ################################
# Early stopping to stop learning after some criterion 
early_stopping_params = {'monitor'   : 'val_loss',  # Value to monitor
                         'min_delta' : 0.001,          # Minimum delta to declare an improvement
                         'patience'  : 20,          # How much time to wait for an improvement
                         'verbose'   : 1,           # Verbosity level
                         'restore_best_weights':True,
                         'mode'      : 'min'}       # Mode : 'auto', 'min', 'max'

# Reduce LR on plateau : if no improvement for some time, will reduce lr by a certain factor
reduceLR_params = {'monitor'    : 'val_loss',   # Value to monitor
                   'factor'     : 0.1,          # Multiplicative factor by which to multiply LR
                   'min_lr'     : 1e-6,         # Minimum value for LR
                   'patience'   : 5,           # How much time to wait for an improvement
                   'cooldown'   : 5,            # How many epochs before starting again to monitor
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
#    'epochs' : [100],   
#    'batch_size' : [10000], 
#    'loss_function' : [categorical_crossentropy] 
#}
p = { 
    'lr' : [0.01], 
    'first_neuron' : [1024],
    'activation' : [relu],
    'dropout' : [0.1],
    'hidden_layers' : [6], # does not take into account the first layer
    'output_activation' : [softmax],
    'l2' : [0.001],
    'optimizer' : [Adam],  
    'epochs' : [100],   
    'batch_size' : [10000], 
    'loss_function' : [categorical_crossentropy] 
}


repetition = 1 # How many times each hyperparameter has to be used 

###################################  Variables   ######################################

cut = 'MC_weight > 0'

weight = 'total_weight'
#weight = None

# Input branches (combinations possible just as in ROOT #
inputs = [
    '$era@onehot_era',
    'lep_E',
    'lep_Px',
    'lep_Py',
    'lep_Pz',
    'lep_charge@onehot_charge',
    'lep_pid@onehot_pdgid',
    'j1_E',
    'j1_Px',
    'j1_Py',
    'j1_Pz',
    'j1_bTagDeepFlavB',
    'j2_E',
    'j2_Px',
    'j2_Py',
    'j2_Pz',
    'j2_bTagDeepFlavB',
    'j3_E',
    'j3_Px',
    'j3_Py',
    'j3_Pz',
    #'j3_bTagDeepFlavB',
    'j4_E',
    'j4_Px',
    'j4_Py',
    'j4_Pz',
    #'j4_bTagDeepFlavB',
    #'nAk4Jets',
    #'nAk4BJets',
    #'METpt',
    #'METphi',
    #'lep_Px',
    #'lep_Py',
    #'lep_Pz',
    #'lep_E',
    #'lep_pt',
    #'lep_eta',
    #'lepmet_DPhi',
    #'lepmet_pt',
    #'lep_MT',
    #'MET_LD',
    #'j1_pt',
    #'j2_pt',
    #'j3_pt',
    #'j4_pt',
    #'j1j2_DR',
    #'j2j3_DR',
    #'j1j4_DR',
    #'j3j4_DR',
    #'j1j2_M',
    #'j3j4_M',
    #'j1MetDPhi',
    #'j3MetDPhi',
    #'j1LepDR',
    #'j3LepDR',
    #'minJetDR',
    #'minDR_lep_allJets',
    #'w1w2_MT',
    #'HT2_lepJetMet',
    #'HT2R_lepJetMet',
    #'mT_top_3particle',
    #'HWW_Mass',
    #'HWW_Simple_Mass',
    #'HWW_dR',
    #'cosThetaS_Hbb',
    #'cosThetaS_Wjj_simple',
    #'cosThetaS_WW_simple_met',
    #'cosThetaS_HH_simple_met',
         ]
onehots = [inp.split('@')[1] if '@' in inp else 'onehot_unit'  for inp  in  inputs]
mask_onehot = [len(inp.split('@'))==2 for inp  in  inputs]
inputs = [inp.split('@')[0] for inp  in  inputs]

# Output branches #
outputs = [
            '$DY',
            '$H',
            '$HH',
            '$Rare',
            '$ST',
            '$TT',
            '$TTVX',
            '$VVV',
            '$WJets',
          ] 

# Other variables you might want to save in the tree #
other_variables = [
                    'event',
                    'ls',
                    'run',
                 ]


################################  dtype operation ##############################

def make_dtype(list_names): # root_numpy does not like . and ()
    list_dtype = [(name.replace('.','_').replace('(','').replace(')','').replace('-','_minus_').replace('*','_times_')) for name in list_names]
    return list_dtype
        



                                



