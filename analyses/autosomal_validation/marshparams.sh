#!/bin/bash

source params_base.sh

#ASDHET=/groups/reich/melissa/mutation_rate_project/sgdp_asdt_vcf/
ASDHET=${BASEDIR}/sgdp_asdt_vcf/
OUTPREFIX=${BASEDIR}/autosomal_validation/sunetal_estimates_ml
LOCI=sun_etal_msat_data_hg19_lobSTRcoords.bed
NUMPROC=5
TRUTHFILE2=marsh_truth_np.bed
SCALE=0.5