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
path_out = '/nfs/scratch/fynu/fbury/HHMachineLearning_output/'
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
partition = 'gpu'  # Def, cp3, cp3-gpu or gpu
QOS = 'normal' # cp3, normal, cp3-gpu
time = '0-08:00:00' # days-hh:mm:ss
mem = '80000' # ram in MB
tasks = 1 # Number of threads(as a string) (not parallel training for classic mode)
cpus = 1
gpus = 1
additional_options = ""
workers = 20
split_per_model = True # in case of cross validation, to send one job per model or not

##################################  Naming ######################################
# Physics Config #
config = os.path.join(os.path.abspath(os.path.dirname(__file__)),'sampleListSL_SM.yml')
lumidict = {2016:35922,2017:41529.152060112,2018:59740.565201546}
#eras = [2016,2017,2018] # To enable or disable eras, add or remove from this list
eras = [2018]

categories = ['resolved2b2Wj','resolved2b1Wj','resolved2b0Wj','resolved1b2Wj','resolved1b1Wj','resolved1b0Wj','resolved0b']
#categories = ['boosted2b2Wj','boosted2b1Wj', 'boosted2b0Wj']
channels = ['El','Mu']

# Better put them in alphabetical order
nodes = ['Ewk','GGF','H','Top','VBF','WJets']

weight_groups = [
                  (1.0/10, ('GGF')),
                  (1.0/10, ('VBF')),
                  (1.0, ('H')),
                  (1.0, ('Ewk')),
                  (1.0, ('Top')),
                  (1.0, ('WJets')),
                ]


quantile = 0.95 # repeat weights with too high learning weights

#    def loss(self): # -> DL
#        return GroupedXEnt(
#            group_ids=[
#                (1.0, [0, 1]),  # Signal
#                (1.0, [2, 3, 4, 5, 6, 7, 8]),  # Bkg
#                (1.0, [0]),  # GGF
#                (1.0, [1]),  # VBF
#                (1.0, [3]),  # DY
#                (1.0, [5]),  # ttbar
#                (1.0, [2, 4, 6, 7, 8]),  # Other
#            ]
#        )#

# Input plots options #
node_colors = {
            'Ewk'  : '#610596',
            'GGF'   : '#288a24',
            'H'     : '#06b894',
            'Top'    : '#cc7a16',
            'VBF'   : '#8f0a1e',
            'WJets' : '#d95564',
             }

# Tree name #
tree_name = 'Events'

# scaler and mask names #
#suffix = 'boosted' 
suffix = 'resolved' 
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
early_stopping_params = {'monitor'   : 'val_categorical_accuracy',  # Value to monitor
                         'min_delta' : 0.0001,          # Minimum delta to declare an improvement
                         'patience'  : 50,          # How much time to wait for an improvement
                         'verbose'   : 1,           # Verbosity level
                         'restore_best_weights': False,
                         'mode'      : 'max'}       # Mode : 'auto', 'min', 'max'

# Reduce LR on plateau : if no improvement for some time, will reduce lr by a certain factor
reduceLR_params = {'monitor'    : 'val_categorical_accuracy',   # Value to monitor
                   'factor'     : 0.5,          # Multiplicative factor by which to multiply LR
                   'min_lr'     : 1e-5,         # Minimum value for LR
                   'min_delta'  : 0.001,       # Minimum delta to declare an improvement
                   'patience'   : 20,            # How much time to wait for an improvement
                   'cooldown'   : 10,            # How many epochs before starting again to monitor
                   'verbose'    : 1,            # Verbosity level
                   'mode'      : 'max'}         # Mode : 'auto', 'min', 'max'


    
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
    'first_neuron' : [256],
    'activation' : [relu],
    'dropout' : [0.001],
    'hidden_layers' : [4], # does not take into account the first layer
    'output_activation' : [softmax],
    'l2' : [1e-6],
    'optimizer' : [Adam],  
    'epochs' : [500],   
    'batch_size' : [50000], 
    'n_particles' : [10],
    'loss_function' : [categorical_crossentropy],
}


repetition = 1 # How many times each hyperparameter has to be used 

###################################  Variables   ######################################

cut = 'total_weight > 0 && MC_weight <10000'
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
            #'$era@op_era',
            'lep_pdgId@op_pdgid',
            'lep_charge@op_charge',
            'JPAcat@op_resolved_jpacat',
            # JPA values #
            'L2_2b2Wj',
            'L2_2b1Wj',
            'L2_2b0Wj',
            'L2_1b2Wj',
            'L2_1b1Wj',
            'L2_1b0Wj',
            'L2_0b',
            # LL variables #
#            'METpt',               
            'METpx',               # discard               
            'METpy',               # discard
            'METenergy',           # discard
            'lep_Px',              # discard
            'lep_Py',              # discard
            'lep_Pz',              # discard
            'lep_E',               # discard     
#           'lep_pt',
#           'lep_eta',             # discard
            'bj1_Px',              # discard
            'bj1_Py',              # discard
            'bj1_Pz',              # discard
            'bj1_E',               # discard
#           'bj1_pt',
#            'bj1_eta',             # discard
           'bj1_bTagDeepFlavB',
            'bj2_Px',              # discard
            'bj2_Py',              # discard
            'bj2_Pz',              # discard
            'bj2_E',               # discard
#           'bj2_pt',
#            'bj2_eta',             # discard
           'bj2_bTagDeepFlavB',
            'wj1_Px',              # discard
            'wj1_Py',              # discard
            'wj1_Pz',              # discard
            'wj1_E',               # discard
#           'wj1_pt',
#            'wj1_eta',             # discard
           'wj1_bTagDeepFlavB',
            'wj2_Px',              # discard
            'wj2_Py',              # discard
            'wj2_Pz',              # discard
            'wj2_E',               # discard
#            'wj2_pt',
#            'wj2_eta',             # discard
           'wj2_bTagDeepFlavB ',
#            'nAk4BJets',           
#            'nAk8BJets',           # discard
            'VBF_tag',
           # HL variables #
            'lepmet_DPhi',
            'lepmet_pt',
            'lep_MT',
            'MET_LD',
            'hT',
            'bj1LepDR',
            'bj2LepDR',
            'bj1MetDPhi',
            'bj2MetDPhi',
#            'bj1LepDPhi',
#            'bj2LepDPhi',
            'minDR_lep_allJets',
            'bj1bj2_pt',
            'wj1wj2_pt',
#            'bj1bj2_M',
            'cosThetaS_Hbb',
            'mT_top_3particle',
            'wj1LepDR',
            'wj2LepDR',
#            'wj1LepDPhi',
#            'wj2LepDPhi',
            'wj1MetDPhi',
            'wj2MetDPhi',
            'wj1wj2_M',
#            'w1w2_MT',
#            'HWW_Mass',
            'HWW_Simple_Mass',
            'HWW_dR',
            'cosThetaS_Wjj_simple',
            'cosThetaS_WW_simple_met ',
            'cosThetaS_HH_simple_met',
            'angleBetWWPlane',
            'angleBetHWPlane',
            'bj1bj2_DR',
            'bj1wj1_DR',
            'bj1wj2_DR',
            'bj2wj1_DR',
            'wj1wj2_DR',
#            'bj1bj2_DPhi',
#            'bj2wj1_DPhi',
#            'wj1wj2_DPhi',
#            'bj1wj2_DPhi',
#            'bj1wj1_DPhi',
            'VBFj1pt',
            'VBFj2pt',
#            'VBFj1eta',
#            'VBFj2eta',
            'VBFj1j2dEta',
            'VBFj1j2dPhi',
            'VBFj1j2invM',
            'zeppenfeldVar',
            'minJetDR',
            'minLepJetDR',
            'HbbLepDR',
#            'HbbLepDPhi',
#            'HbbMetDPhi',
            'HWWLepDR',
#            'HWWLepDPhi',
            'HWWMetDPhi',
#            'HT2_lepJetMet',
#            'HT2R_lepJetMet',

#            ############################
#            ####-----  Boosted  ----####
#            ############################
#            # Onehot #
#            '$era@op_era',
#            'lep_pdgId@op_pdgid',
#            'lep_charge@op_charge',
#            'JPAcat@op_boosted_jpacat',
#            # JPA values #
#            'L2_Hbb2Wj',
#            'L2_Hbb1Wj',
#            'L2_Hbb0Wj',
#            # LL variables #
#            'METpt',               
#            'METpx',               # discard               
#            'METpy',               # discard
#            'lep_Px',              # discard
#            'lep_Py',              # discard
#            'lep_Pz',              # discard
#            'lep_E',               # discard    
#            'fatbj_Px',
#            'fatbj_Py',
#            'fatbj_Pz',
#            'fatbj_E',
#            'fatbj_softdropMass',
#            'fatbj_tau1',
#            'fatbj_tau2',
#            'fatbj_tau3',
#            'fatbj_tau4',
#            #'fatbj_btagDeepB',
#            'wj1_Px',
#            'wj1_Py',
#            'wj1_Pz',
#            'wj1_E',
#            'wj1_bTagDeepFlavB',
#            'wj2_Px',
#            'wj2_Py',
#            'wj2_Pz',
#            'wj2_E',
#            'wj2_bTagDeepFlavB',
#            'VBF_tag',
#            'VBFj1pt',
#            'VBFj2pt',
#            'VBFj1j2dEta',
#            'VBFj1j2invM',
#            # HL variables #
#            'hT',
#            'lepmet_DPhi',
#            'lepmet_pt',
#            'lep_MT',
#            'MET_LD',
#            'fatbj_lepDR',
#            'fatbj_Wj1DR',
#            'fatbj_wj2DR',
#            #'fatbjSub1_lepDR',
#            #'fatbjSub2_lepDR',
#            #'fatbjSub1_Wj1DR',
#            #'fatbjSub2_Wj1DR',
#            #'fatbjSub1_Wj2DR',
#            #'fatbjSub2_Wj2DR',
#            'wj1_lepDR',
#            'wj2_lepDR',
#            'wj1wj2_pt',
#            'wj1wj2DR',
#            'wj1wj2invM',
#            'mT_top_3particle',
#            'WWplaneAngle_withMET', 
#            'HWplaneAngle',
#            'HWW_Simple_Mass',
#            'HWW_dR',
#            'cosThetaS_Hbb',
#            'cosThetaS_Wjj_simple',
#            'cosThetaS_WW_simple_met',
#            'cosThetaS_HH_simple_met',
#            #'minSubJetLepDR',
#            'MT_W1W2',
#            'zeppenfeldVar',
#            #'HT2_lepJetMet',
#            #'HT2R_lepJetMet',

    ]

operations = [inp.split('@')[1] if '@' in inp else None  for inp  in  inputs]
check_op = [(o is not None)*1 for o in operations]
if check_op != sorted(check_op,reverse=True):
    raise RuntimeError('Onehot inputs need to be at the beginning of the inputs list')
mask_op = [len(inp.split('@'))==2 for inp  in  inputs]
inputs = [inp.split('@')[0] for inp  in  inputs]

#### Resolved ####
LBN_inputs = [
              'lep_E','lep_Px','lep_Py','lep_Pz',
              'bj1_E','bj1_Px','bj1_Py','bj1_Pz',
              'bj2_E','bj2_Px','bj2_Py','bj2_Pz',
              'wj1_E','wj1_Px','wj1_Py','wj1_Pz',
              'wj2_E','wj2_Px','wj2_Py','wj2_Pz',
             ]

#### Boosted ####
#LBN_inputs = [
#              'lep_E','lep_Px','lep_Py','lep_Pz',
#              'fatbj_E','fatbj_Px','fatbj_Py','fatbj_Pz',
#              'wj1_E','wj1_Px','wj1_Py','wj1_Pz',
#              'wj2_E','wj2_Px','wj2_Py','wj2_Pz',
#             ]
#
#

#### /!\ format = [E,Px,Py,Pz] per particle


assert len(LBN_inputs)%4 == 0


# Output branches #

outputs = [
            '$Ewk',
            '$GGF',
            '$H',
            '$Top',
            '$VBF',
            '$WJets',
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
        



                                



