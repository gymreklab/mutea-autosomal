#!/bin/bash

PARAMFILE=$1
source $PARAMFILE

for STEPPARAM in ${STEPPARAMS}
do
    strsd=$(${SCRIPTSDIR}/get_strsd.py ${STEPPARAM})
    for STR_MUT_RATE in ${MUTS}
    do
	for LENC in ${LENCS}
	do
	    for ML in ${MAXLOCI}
	    do
		logmu=$(echo $STR_MUT_RATE | awk '{print log($1)/log(10)}')
		echo $logmu $LENC $strsd $ML $(cat ${OUTDIR}/joint/${PREFIX}_joint_${STEPPARAM}_${STR_MUT_RATE}_${LENC}_${ML}.tab | cut -f 2-)
	    done
	done
    done
done