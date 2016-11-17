#!/bin/bash

PARAMFILE=$4
source $PARAMFILE

STR_MUT_RATE=$1
LENC=$2
STEPPARAM=$3
STRSD=$(${SCRIPTSDIR}/get_strsd.py ${STEPPARAM})

# Run simulation
echo "[simulate_locus.sh] Simulate mutations on trees..."
python2.7 ${SCRIPTSDIR}/SimulateSTRMutationTree.py \
    --Neff ${NEFF} \
    --samples ${SAMPLES} \
    --strmutrate ${STR_MUT_RATE} \
    --numsim ${NUMSIM} \
    --lenconstraint ${LENC} \
    --strsd ${STRSD} \
    --strstepparam ${STEPPARAM} \
    ${PAIRS} \
    --out ${DATADIR}/${PREFIX}_sim_mutrate${STR_MUT_RATE}_stepparam${STEPPARAM}_lenc${LENC} \
    --method discreteMultiOU