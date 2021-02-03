#!/bin/bash
module purge
export PYTHONPATH=/python3/lib/python3.6/site-packages/:$PYTHONPATH
export PYTHONPATH=/root6/lib:$PYTHONPATH
module load python/python36_sl7_gcc73
module load slurm/slurm_utils
