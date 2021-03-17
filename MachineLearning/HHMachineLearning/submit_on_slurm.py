#! /bin/env python

import copy
import os
import datetime
import sys
import glob
import logging

# Slurm configuration
from CP3SlurmUtils.Configuration import Configuration
from CP3SlurmUtils.SubmitWorker import SubmitWorker
from CP3SlurmUtils.Exceptions import CP3SlurmUtilsException

# Personal files #
import parameters

def submit_on_slurm(name,args,debug=False):
    # Check arguments #
    GPU = args.find("--GPU") != -1
    output = args.find("--output") != -1

    config = Configuration()
    config.sbatch_partition = parameters.partition
    config.sbatch_qos = parameters.QOS
    config.sbatch_chdir = parameters.main_path
    config.sbatch_time = parameters.time
    config.sbatch_additionalOptions = [parameters.additional_options]
    config.sbatch_memPerCPU = parameters.mem
    if parameters.partition == 'cp3-gpu':
        config.sbatch_additionalOptions += ['--export=NONE']
        #if parameters.cpus > 1:
        #    config.sbatch_additionalOptions += ["--cpus-per-gpu={}".format(parameters.cpus)]
        #config.sbatch_additionalOptions += ['--mem-per-gpu={}'.format(parameters.mem)]
    elif parameters.partition == 'gpu':
        config.sbatch_additionalOptions += ['--gres=gpu:TeslaV100:{}'.format(parameters.gpus),'--export=NONE']
        #config.sbatch_additionalOptions += ['--mem-per-gpu={}'.format(parameters.mem)]
        #if parameters.cpus > 1:
        #    config.sbatch_additionalOptions += ["--cpus-per-gpu={}".format(parameters.cpus)]
        #if parameters.cpus > 1:
        #    config.sbatch_additionalOptions += ["--cpus-per-gpu={}".format(parameters.cpus)]
        #    config.sbatch_additionalOptions += ["--cpus-per-task={}".format(parameters.cpus)]
            #config.sbatch_additionalOptions += ["-c {}".format(parameters.cpus)]
    else:
        if parameters.tasks > 1:
            config.sbatch_additionalOptions += ["-n={}".format(parameters.tasks)]
        if parameters.cpus > 1:
            config.sbatch_additionalOptions += ["--cpus-per-task={}".format(parameters.cpus)]
        
    config.inputSandboxContent = []
    config.useJobArray = True
    config.inputParamsNames = []
    config.inputParams = []
    if output:
        config.inputParamsNames += ["--verbose"]
        config.inputParams += [[""]]
    if not output:
        config.inputParamsNames += ['scan','task']
        if parameters.crossvalidation and parameters.split_per_model:
            config.inputParamsNames += ['modelId']

    config.payload = ""

    if parameters.partition == 'cp3-gpu':
        config.payload += "export PYTHONPATH=/python3/lib/python3.6/site-packages/:$PYTHONPATH\n" # GPU tf
        config.payload += "export PYTHONPATH=/root6/lib:$PYTHONPATH\n" # ROOT
        config.payload += "module load cp3\n" # needed on gpu to load slurm_utils
        config.payload += "module load python/python36_sl7_gcc73\n" 
        config.payload += "module load slurm/slurm_utils\n"
    if parameters.partition == 'gpu':
        config.payload += "module load releases/2019b_test \n"
        config.payload += "module load cp3\n" # needed on gpu to load slurm_utils
        config.payload += "module load root/6.12.04-sl7_gcc73 \n"
        config.payload += "module load root_numpy \n"
        config.payload += "module load TensorFlow \n"
        config.payload += "module load slurm/slurm_utils\n"
        
    config.payload += "python3 {script} "
    if not output:
        config.payload += "--scan ${{scan}} --task ${{task}} "
    if parameters.crossvalidation and parameters.split_per_model:
        config.payload += "--modelId ${{modelId}}"
    config.payload += args

    timestamp = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    out_dir = parameters.main_path

    slurm_config = copy.deepcopy(config)
    slurm_working_dir = os.path.join(out_dir,'slurm',name+'_'+timestamp)

    slurm_config.batchScriptsDir = os.path.join(slurm_working_dir, 'scripts')
    slurm_config.inputSandboxDir = slurm_config.batchScriptsDir
    slurm_config.stageoutDir = os.path.join(slurm_working_dir, 'output')
    slurm_config.stageoutLogsDir = os.path.join(slurm_working_dir, 'logs')
    slurm_config.stageoutFiles = ["*.csv","*.zip","*.png"]

    slurm_config.payload = config.payload.format(script=os.path.join(out_dir,"HHMachineLearning.py"))

    if not output:
        for f in glob.glob(os.path.join(parameters.main_path,'split',name,'*.pkl')):
            task = os.path.basename(f)
            if parameters.crossvalidation and parameters.split_per_model:
                for N in range(parameters.N_models):
                    slurm_config.inputParams.append([name,task,N])
            else:
                slurm_config.inputParams.append([name,task])

    # Submit job!

    logging.info("Submitting job...")
    if not debug:
        submitWorker = SubmitWorker(slurm_config, submit=True, yes=True, debug=False, quiet=False)
        submitWorker()
        logging.info("Done")
    else:
        logging.info("Number of jobs : %d"%len(slurm_config.inputParams))
        logging.info(slurm_config.payload)
        logging.info(slurm_config.inputParamsNames)
        for inputParam in slurm_config.inputParams:
            logging.info(inputParam)
        logging.info('... don\'t worry, jobs not sent')

