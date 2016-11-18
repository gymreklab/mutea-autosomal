#!/bin/bash

PARAMFILE=$1
source $PARAMFILE

python2.7 ~/workspace/MUTEA/calibrate_errors.py \
    --ests ${OUTPREFIX}.bed \
    --scale ${SCALE} \
    --truthnp ${TRUTHFILE2} > ${OUTPREFIX}_calibrate_errors_np.tab
