#!/bin/bash

source params_base.sh

#ASDHET=/groups/reich/melissa/mutation_rate_project/sgdp_asdt_vcf/
ASDHET=${BASEDIR}/sgdp_asdt_vcf/
OUTPREFIX=${BASEDIR}/autosomal_validation/codis_estimates_ml
LOCI=codis_loci.bed
NUMPROC=1
TRUTHFILE2=codis_truth_np.bed
SCALE=0.5