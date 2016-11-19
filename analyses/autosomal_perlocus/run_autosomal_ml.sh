#!/bin/bash

PARAMFILE=$1
source ${PARAMFILE}

for locfile in $(ls ${LOCDIR})
do
    chrom=$(echo $locfile | cut -f 1 -d'.')
    batch=$(echo $locfile | cut -f 2 -d'.')
    bsub -q long -W 500:00 \
	-oo ${LOGDIR}/${chrom}.${batch}.log \
	-eo ${LOGDIR}/${chrom}.${batch}.err \
	-J ${locfile}ml \
	${MUTEADIR}/mutea-auto/main_autosomal.py \
	--asdhet ${DATADIR}/sgdp_asdt_chr${chrom}.vcf.gz --vcf \
	--eststutter ${OUTDIR}/${chrom}.${batch}_stutter.tab \
	--loci ${LOCDIR}/${locfile} \
	--min_samples 50 \
	--out ${OUTDIR}/${chrom}.${batch}_combined.estimates.tab
done