#!/bin/bash

PARAMFILE=$1
source $PARAMFILE

${MUTEADIR}/main_autosomal.py \
    --asdhet ${OUTDIR}/${PREFIX}_asdhet.vcf.gz --vcf \
    --loci ${OUTDIR}/${PREFIX}_loci.bed \
    --min_samples 50 \
    --stderrs fisher \
    --numproc ${NUMPROC} \
    --out ${PREFIX}_maxlik.tab

