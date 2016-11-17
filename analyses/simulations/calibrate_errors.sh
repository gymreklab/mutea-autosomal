#!/bin/bash

PARAMFILE=$1
source $PARAMFILE

# Make truth file
cat ${OUTDIR}/${PREFIX}_truth.bed | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $4 "\t" $4}' > ${OUTDIR}/${PREFIX}_truth2.bed
python2.7 ~/workspace/MUTEA/calibrate_errors.py \
    --ests ${OUTDIR}/${PREFIX}_maxlik_calcstutter.tab \
    --truth ${OUTDIR}/${PREFIX}_truth2.bed > ${OUTDIR}/${PREFIX}_calibrate_errors.tab
