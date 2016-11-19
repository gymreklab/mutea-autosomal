#!/bin/bash

MUTEADIR=${HOME}/workspace/mutea-anutosomal/
BASEDIR=/groups/reich/melissa/mutation_rate_project/analysis_round1/autosomal_estimates
DATADIR=/groups/reich/melissa/mutation_rate_project/sgdp_asdt_vcf
LOBREF=/home/mag50/hs37d5_v3.0.3/lobstr_v3.0.2_hg19_ref_nochr.bed

NUMLINES=5000 # Break to batches of this many
LOCDIR=${BASEDIR}/getloci/batches/
OUTDIR=${BASEDIR}/batches_ml/
LOGDIR=${BASEDIR}/log_ml/
TMPLOC=/n/scratch2/gymrek


