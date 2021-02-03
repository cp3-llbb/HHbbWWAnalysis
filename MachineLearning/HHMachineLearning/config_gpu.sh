#!/bin/bash
module purge
module load releases/2019b_test
module load root/6.12.04-sl7_gcc73
module load root_numpy
module load TensorFlow
module load slurm/slurm_utils

