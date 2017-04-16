#!/bin/bash

PARAMFILE=$1
source $PARAMFILE

./filter_and_scale.py \
    --ests ${BASEDIR}/autosomal_estimates_smm.bed.gz \
    --vcf ${DATADIR}/ \
    | bgzip -c > \
    ${BASEDIR}/autosomal_estimates_smm_scaled.bed.gz
tabix -p bed ${BASEDIR}/autosomal_estimates_smm_scaled.bed.gz
zcat ${BASEDIR}/autosomal_estimates_smm_scaled.bed.gz | grep PASS | bgzip -c > ${BASEDIR}/autosomal_estimates_smm_filtered.bed.gz
tabix -p vcf ${BASEDIR}/autosomal_estimates_smm_filtered.bed.gz
