#!/bin/bash

PARAMFILE=$1
source $PARAMFILE

# Step 1: simulate all loci
step1jobs=""
counter=0
for STEPPARAM in ${STEPPARAMS}
do
    for STR_MUT_RATE in ${MUTS}
    do
	for LENC in ${LENCS}
	do
	    jobname="${PREFIX}${counter}"
	    bsub -q short -W 4:00 -J ${jobname} \
		-oo ${LOGDIR}/${PREFIX}_${STR_MUT_RATE}_${STEPPARAM}_lenc${LENC}.log \
		-eo ${LOGDIR}/${PREFIX}_${STR_MUT_RATE}_${STEPPARAM}_lenc${LENC}.err \
		./simulate_locus.sh \
		${STR_MUT_RATE} ${LENC} ${STEPPARAM} ${PARAMFILE}
	    step1jobs="${step1jobs} ended(${jobname})"
	    counter=$((counter+1))
	done
    done
done

# Step 2: Get estimation files
step1jobs=$(echo ${step1jobs} | sed 's/ /\&\&/g')
bsub -q short -W 1:00 -J ${PREFIX}.est -w "${step1jobs}" \
    -eo ${LOGDIR}/${PREFIX}.est.err \
    -oo ${LOGDIR}/${PREFIX}.est.out \
    ./get_locus_files.sh ${PARAMFILE}

# Step 3: Run estimation - maxlik - nostutter
bsub -q short -W 12:00 -J ${PREFIX}.ml -w "ended(${PREFIX}.est)" \
    -eo ${LOGDIR}/${PREFIX}.ml.err \
    -oo ${LOGDIR}/${PREFIX}.ml.out \
    -n ${NUMPROC} \
    ./run_maxlik.sh ${PARAMFILE}

# Step 3: Run estimation - maxlik - nostutter
bsub -q short -W 12:00 -J ${PREFIX}.ml -w "ended(${PREFIX}.est)" \
    -eo ${LOGDIR}/${PREFIX}.ml.stutter.err \
    -oo ${LOGDIR}/${PREFIX}.ml.stutter.out \
    -n ${NUMPROC} \
    ./run_maxlik_stutter.sh ${PARAMFILE}
