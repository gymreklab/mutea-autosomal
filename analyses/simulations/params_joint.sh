#!/bin/bash

# Directories
SCRIPTSDIR=.
MUTEADIR=/home/mag50/workspace/mutea-autosomal/
LOGDIR=/groups/reich/melissa/mutation_rate_project/analysis_round1/simulations/log
OUTDIR=/groups/reich/melissa/mutation_rate_project/analysis_round1/simulations
DATADIR=/groups/reich/melissa/mutation_rate_project/analysis_round1/simulations/simdata

# Params
NUMSIM=1000
SAMPLES=600 # SGDP sample size. will get half this
NUMBS=0
NUMPROC=5
NEFF=100000

# Preprocessed filenames
PREFIX=treejoint
OUTFILE=${OUTDIR}/${PREFIX}_asdhet.bed
VCFFILE=${OUTDIR}/${PREFIX}_asdhet.vcf
TRUTHFILE=${OUTDIR}/${PREFIX}_truth.bed
LOCFILE=${OUTDIR}/${PREFIX}_loci.bed
STRSDFILE=${OUTDIR}/${PREFIX}_strsd.bed
SUMMFILE=${OUTDIR}/${PREFIX}_summ.bed

STEPPARAMS="0.95"
MUTS="0.00001 0.000001 0.0000001 0.00000001"
LENCS="0.3"

PAIRS="--diploid"
STUTTER="0.0,0.0,1.0"
COVERAGE=10
VCFARGS=""
#MAXLOCI="10 100 500 1000"
MAXLOCI="1 10 20 30 40 50 60 70 80 90 100 500 750 1000"