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
	    # Get list of loci
	    locfile="${OUTDIR}/loci/loci_${STEPPARAM}_${STR_MUT_RATE}_${LENC}.bed"
	    cat ${PREFIX}_truth.bed | intersectBed -a stdin -b ${PREFIX}_strsd.bed -wa -wb | \
		awk -v"mu=${STR_MUT_RATE}" -v "lenc=${LENC}" -v "strsd=${strsd}" \
		'($4==mu && $5==lenc && $9==strsd)' | cut -f 1-3 > ${locfile}
	    for ML in ${MAXLOCI}
	    do
                # Run joint estimation
#		bsub -q short -W 4:00 \
#		    -eo ${LOGDIR}/joint_${STEPPARAM}_${STR_MUT_RATE}_${LENC}_${ML}.err \
#		    -oo ${LOGDIR}/joint_${STEPPARAM}_${STR_MUT_RATE}_${LENC}_${ML}.out \
		sbatch --partition=compute --account=${ACCOUNT} --get-user-env --time=60 \
		    --error=${LOGDIR}/joint_${STEPPARAM}_${STR_MUT_RATE}_${LENC}_${ML}.err \
		    --output=${LOGDIR}/joint_${STEPPARAM}_${STR_MUT_RATE}_${LENC}_${ML}.out \
		    /home/mag50/workspace/MUTEA/main_autosomal.py \
		    --asdhet ${OUTDIR}/${PREFIX}_asdhet.vcf.gz --vcf \
		    --locus_priors ${OUTDIR}/${PREFIX}_maxlik.tab \
		    --loci ${locfile} \
		    --joint --use_locus_means \
		    --maxloci ${ML} \
		    --out ${OUTDIR}/joint/${PREFIX}_joint_${STEPPARAM}_${STR_MUT_RATE}_${LENC}_${ML}.tab
	    done
	done
    done
done