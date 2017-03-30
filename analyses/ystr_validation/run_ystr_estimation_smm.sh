#!/bin/bash

PARAMFILE=$1
source $PARAMFILE

${MUTEADIR}/mutea-auto/main_autosomal.py \
    --asdhet ${VCF} --vcf \
    --loci ${LOCI} \
    --numproc 5 \
    --min_samples 50 \
    --eststutter ${OUTPREFIX}_stutter.tab \
    --use-sample-pairs ${SAMPLEPAIRS} \
    --smm \
    --out ${OUTPREFIX}_smm.tab