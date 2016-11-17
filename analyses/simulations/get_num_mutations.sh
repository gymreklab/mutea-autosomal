#!/bin/bash

PARAMFILE=$1
source $PARAMFILE

${SCRIPTSDIR}/MutationRateEstimatorTMRCA.py \
    --asdhet ${OUTDIR}/${PREFIX}_asdhet.bed.gz \
    --loci ${OUTDIR}/${PREFIX}_loci.bed \
    --mutinfo \
    --out ${OUTDIR}/${PREFIX}_mutinfo.tab