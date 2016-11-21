#!/bin/bash

PARAMFILE=$1
source $PARAMFILE

tmpdir=$(mktemp -d -p ${TMPLOC})

zcat ${BASEDIR}/constraint/lobSTR_ref_GRCh37_properties_filtered.tab | grep -v chrom | \
    intersectBed -a ${BASEDIR}/autosomal_estimates/perlocus/autosomal_estimates_ml_filtered.bed.gz -b stdin -wa -wb | \
    cut -f 10-12 --complement > ${tmpdir}/tmp
echo "chrom,start,end,ml_mu,ml_beta,ml_p,ml_mu_stderr,numsamples,strfilter,motif,length,uninterrupted_length,recomb,gc,entropy,reptiming,featurefilter" | \
    sed 's/,/\t/g' > ${BASEDIR}/constraint/autosomal_perlocus_observed.bed
cat ${tmpdir}/tmp >> ${BASEDIR}/constraint/autosomal_perlocus_observed.bed
gzip ${BASEDIR}/constraint/autosomal_perlocus_observed.bed
