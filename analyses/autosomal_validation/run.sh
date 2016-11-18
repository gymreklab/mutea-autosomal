#!/bin/bash

PARAMFILE=$1
source $PARAMFILE

MINSAMPLES=50

OUTFILE=${OUTPREFIX}.bed
${MUTEADIR}/mutea-auto/main_autosomal.py \
    --asdhet ${ASDHET} --vcf \
    --asdhetdir \
    --loci ${LOCI} \
    --min_samples ${MINSAMPLES} \
    --numproc ${NUMPROC} \
    --stderrs fisher \
    --eststutter ${OUTPREFIX}_maxlik_stutterparams.tab \
    --out ${OUTFILE}