#!/bin/bash

PARAMFILE=$1
source $PARAMFILE

# Step 1: simulate all loci
step1jobs="afterok"
counter=0
for STEPPARAM in ${STEPPARAMS}
do
    for STR_MUT_RATE in ${MUTS}
    do
	for LENC in ${LENCS}
	do
	    jobname="${PREFIX}${counter}"
#	    bsub -q short -W 4:00 -J ${jobname} \
#		-oo ${LOGDIR}/${PREFIX}_${STR_MUT_RATE}_${STEPPARAM}_lenc${LENC}.log \
#		-eo ${LOGDIR}/${PREFIX}_${STR_MUT_RATE}_${STEPPARAM}_lenc${LENC}.err \
	    jid=$(sbatch --job-name=${jobname} --time=60 --partition=compute --account=${ACCOUNT} --get-user-env \
		--error=${LOGDIR}/${PREFIX}_${STR_MUT_RATE}_${STEPPARAM}_lenc${LENC}.err \
		--output=${LOGDIR}/${PREFIX}_${STR_MUT_RATE}_${STEPPARAM}_lenc${LENC}.log \
		./simulate_locus.sh \
		${STR_MUT_RATE} ${LENC} ${STEPPARAM} ${PARAMFILE} | awk '{print $NF}')
#	    step1jobs="${step1jobs} ended(${jobname})"
	    step1jobs="${step1jobs}:${jid}"
	    counter=$((counter+1))
	done
    done
done

echo "Step 1 jobs: ${step1jobs}"

# Step 2: Get estimation files
#step1jobs=$(echo ${step1jobs} | sed 's/ /\&\&/g')
#bsub -q short -W 1:00 -J ${PREFIX}.est -w "${step1jobs}" \
#    -eo ${LOGDIR}/${PREFIX}.est.err \
#    -oo ${LOGDIR}/${PREFIX}.est.out \
jidest=$(sbatch --job-name=${PREFIX}.est --time=60 --partition=compute \
    --account=${ACCOUNT} --get-user-env --dependency=${step1jobs} \
    --error=${LOGDIR}/${PREFIX}.est.err \
    --output=${LOGDIR}/${PREFIX}.est.out \
    ./get_locus_files.sh ${PARAMFILE} | awk '{print $NF}')

# Step 3: Run estimation - maxlik - nostutter
#bsub -q short -W 12:00 -J ${PREFIX}.ml -w "ended(${PREFIX}.est)" \
#    -eo ${LOGDIR}/${PREFIX}.ml.err \
#    -oo ${LOGDIR}/${PREFIX}.ml.out \
#    -n ${NUMPROC} \
jidml=$(sbatch --job-name=${PREFIX}.ml --time=200 --partition=compute \
    --account=${ACCOUNT} --get-user-env --dependency=afterok:${jidest} \
    --mincpus=${NUMPROC} \
    --error=${LOGDIR}/${PREFIX}.ml.err \
    --output=${LOGDIR}/${PREFIX}.est.out \
    ./run_maxlik.sh ${PARAMFILE} | awk '{print $NF}')

# Step 4: Run estimation - maxlik - stutter
#bsub -q short -W 12:00 -J ${PREFIX}.ml -w "ended(${PREFIX}.est)" \
#    -eo ${LOGDIR}/${PREFIX}.ml.stutter.err \
#    -oo ${LOGDIR}/${PREFIX}.ml.stutter.out \
#    -n ${NUMPROC} \
jidmls=$(sbatch --job-name=${PREFIX}.ml.stutter  --time 200 --partition=compute \
    --account=${ACCOUNT} --get-user-env --dependency=afterok:${jidest} \
    --mincpus=${NUMPROC} \
    --error=${LOGDIR}/${PREFIX}.ml.stutter.err \
    --output=${LOGDIR}/${PREFIX}.ml.stutter.out \
    ./run_maxlik_stutter.sh ${PARAMFILE} | awk '{print $NF}')
