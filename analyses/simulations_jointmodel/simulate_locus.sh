#!/bin/bash

logmu=$1
lencoeff=$2
paramfile=$3

source $paramfile


STRSD=$(${SCRIPTSDIR}/simulations/get_strsd.py ${STEPPARAM})

for sim in $(seq 1 $NUMSIM)
do
    adjmu=$(echo $logmu | awk -v"lencoeff=$lencoeff" -v"sim=$sim" -v "numsim=$NUMSIM" \
	'{print 10**($1+lencoeff*(sim-numsim/2))}')
    echo "[run_joint_simulation.sh] Simulate mutation rate... $adjmu"
    python2.7 ${SCRIPTSDIR}/simulations/SimulateSTRMutationTree.py \
	--Neff ${NEFF} \
	--samples ${SAMPLES} \
	--strmutrate ${adjmu} \
	--numsim 1 \
	--lenconstraint ${LENC} \
	--strsd ${STRSD} \
	--strstepparam ${STEPPARAM} \
	${PAIRS} \
	--out ${DATADIR}/${PREFIX}_logmu${logmu}_lencoeff${lencoeff}_${sim} \
	--method discreteMultiOU \
	--stutter ${STUTTER}
done
