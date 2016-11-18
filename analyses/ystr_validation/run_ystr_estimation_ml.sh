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
    --eststutter ${OUTPREFIX}_stutter.tab \
    --use-sample-pairs ${SAMPLEPAIRS} \
    --out ${OUTPREFIX}_ml.tab