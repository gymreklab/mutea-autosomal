#!/bin/bash

ACCOUNT=ddp268
MUTEADIR=${HOME}/workspace/mutea-autosomal/

# On orchestra:
#BASEDIR=/groups/reich/melissa/mutation_rate_project/analysis_round1/autosomal_estimates
#DATADIR=/groups/reich/melissa/mutation_rate_project/sgdp_asdt_vcf
#LOBREF=/home/mag50/hs37d5_v3.0.3/lobstr_v3.0.2_hg19_ref_nochr.bed
#TMPLOC=/n/scratch2/gymrek

# On comet:
BASEDIR=/oasis/projects/nsf/ddp268/mgymrek/mutea-autosomal/autosomal_estimates
DATADIR=/oasis/projects/nsf/ddp268/mgymrek/mutea-autosomal/sgdp_asdt_vcf
LOBREF=/oasis/projects/nsf/ddp268/mgymrek/dbase/lobstr/lobstr_v3.0.2_hg19_ref_nochr.bed
TMPLOC=/oasis/scratch/comet/mgymrek/temp_project

NUMLINES=5000 # Break to batches of this many
LOCDIR=${BASEDIR}/getloci/batches/

OUTDIR=${BASEDIR}/batches_ml/
LOGDIR=${BASEDIR}/log_ml/

OUTDIRSMM=${BASEDIR}/batches_smm/
LOGDIRSMM=${BASEDIR}/log_smm/

