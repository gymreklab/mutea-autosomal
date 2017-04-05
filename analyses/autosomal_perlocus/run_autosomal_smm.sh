#!/bin/bash

PARAMFILE=$1
source ${PARAMFILE}

#    bsub -q long -W 500:00 \
#	-oo ${LOGDIRSMM}/${chrom}.${batch}.log \
#	-eo ${LOGDIRSMM}/${chrom}.${batch}.err \
#	-J ${locfile}ml \
# 	--eststutter ${OUTDIRSMM}/${chrom}.${batch}_stutter.tab \

#    sbatch --job-name=${locfile}smm --time=100 --partition=compute --mem=2000 \
#	--account=${ACCOUNT} --get-user-env \
#	--error=${LOGDIRSMM}/${chrom}.${batch}.err \
#	--output=${LOGDIRSMM}/${chrom}.${batch}.log \

for locfile in $(ls ${LOCDIR}) 
do
    chrom=$(echo $locfile | cut -f 1 -d'.')
    batch=$(echo $locfile | cut -f 2 -d'.')
    cmd="${MUTEADIR}/mutea-auto/main_autosomal.py \
	--asdhet ${DATADIR}/sgdp_asdt_chr${chrom}.vcf.gz --vcf \
	--loci ${LOCDIR}/${locfile} \
	--min_samples 50 \
	--smm \
	--usestutter ${BASEDIR}/perlocus/autosomal_stutter_ml.bed \
	--out ${OUTDIRSMM}/${chrom}.${batch}_combined_smm.estimates.tab 2> ${LOGDIRSMM}/${chrom}.${batch}.err"
    echo $cmd
done | xargs -I% -n1 -P5 sh -c "%"
