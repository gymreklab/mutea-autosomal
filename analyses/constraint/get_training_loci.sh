#!/bin/bash

PARAMFILE=$1
source ${PARAMFILE}

zcat ${BASEDIR}/constraint/autosomal_perlocus_observed.bed.gz | \
    head -n 1 > ${BASEDIR}/constraint/autosomal_perlocus_train_intergenic.bed
zcat ${BASEDIR}/constraint/autosomal_perlocus_observed.bed.gz | grep -v start | \
    intersectBed -a stdin -b ${INTERGENIC} >> \
    ${BASEDIR}/constraint/autosomal_perlocus_train_intergenic.bed
gzip ${BASEDIR}/constraint/autosomal_perlocus_train_intergenic.bed
