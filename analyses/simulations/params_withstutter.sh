#!/bin/bash

source params_base.sh

# Directories
SCRIPTSDIR=.
MUTEADIR=${HOME}/workspace/mutea-autosomal/
LOGDIR=${BASEDIR}/simulations/log
OUTDIR=${BASEDIR}/simulations/
DATADIR=${BASEDIR}/simulations/simdata

# Params
NUMSIM=10
SAMPLES=600 # SGDP sample size. will get half this
NUMBS=0
NUMPROC=5
NEFF=100000

# Preprocessed filenames
PREFIX=treestutter
VCFFILE=${OUTDIR}/${PREFIX}_asdhet.vcf
OUTFILE=${OUTDIR}/${PREFIX}_asdhet.bed
TRUTHFILE=${OUTDIR}/${PREFIX}_truth.bed
LOCFILE=${OUTDIR}/${PREFIX}_loci.bed
STRSDFILE=${OUTDIR}/${PREFIX}_strsd.bed
SUMMFILE=${OUTDIR}/${PREFIX}_summ.bed

STEPPARAMS="0.7 0.95"
MUTS="0.01 0.001 0.0001 0.00001 0.000001 0.0000001 0.00000001"
LENCS="0.1 0.3 0.5 0.7"

PAIRS="--diploid"
STUTTER="0.1,0.05,0.9"
COVERAGE=5
VCFARGS="--reestimate-gt"
