#!/bin/bash

PARAMFILE=$1
source $PARAMFILE

${MUTEADIR}/mutea-auto/main_autosomal.py \
    --asdhet ${VCF} --vcf \
    --loci ${LOCI} \
    --numproc 5 \
    --max_mu 0.01 \
    --stderrs fisher \
    --min_samples 50 \
    --max_samples 10000 \
    --eststutter ${OUTPREFIX}_stutter_h1.tab \
    --use-sample-pairs ${SAMPLEPAIRS} \
    --out ${OUTPREFIX}_h1.tab \
    --output-likelihood

${MUTEADIR}/mutea-auto/main_autosomal.py \
    --asdhet ${VCF} --vcf \
    --loci ${LOCI} \
    --numproc 5 \
    --max_mu 0.01 \
    --max_beta 0.1 \
    --stderrs fisher \
    --min_samples 50 \
    --max_samples 10000 \
    --eststutter ${OUTPREFIX}_stutter_h0.tab \
    --use-sample-pairs ${SAMPLEPAIRS} \
    --out ${OUTPREFIX}_h0.tab \
    --output-likelihood
