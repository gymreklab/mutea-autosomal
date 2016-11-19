#!/bin/bash

# Directories
SCRIPTSDIR=/home/mag50/workspace/cteam/
LOGDIR=/groups/reich/melissa/mutation_rate_project/analysis_round1/simulations_joint/log
OUTDIR=/groups/reich/melissa/mutation_rate_project/analysis_round1/simulations_joint/
DATADIR=/groups/reich/melissa/mutation_rate_project/analysis_round1/simulations_joint/data

# Params
NUMSIM=100
SAMPLES=600 # SGDP sample size. will get half this
NEFF=10000

# Preprocessed filenames
PREFIX=simjoint
NUMPROC=5

# Params
LENC=0.3
STEPPARAM=0.9
MUTMU="-6 -5 -4 -3 -2"
LENCOEFF="0 0.01 0.02 0.03"

PAIRS="--diploid"
STUTTER="0.0,0.0,1.0"