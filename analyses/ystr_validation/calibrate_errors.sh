#!/bin/bash

PARAMFILE=$1
source $PARAMFILE

python2.7 ${MUTEADIR}/mutea-auto/calibrate_errors.py \
    --ests ${OUTPREFIX}_ml.tab \
    --truth ${TRUTHFILE} > ${OUTPREFIX}_calibrate_errors.tab

python2.7 ${MUTEADIR}/mutea-auto/calibrate_errors.py \
    --ests ${OUTPREFIX}_ml.tab \
    --scale ${SCALE} \
    --truthnp ${TRUTHFILE2} > ${OUTPREFIX}_calibrate_errors_np.tab