#!/bin/bash

PARAMFILE=$1
source ${PARAMFILE}

mkdir -p ${LOCDIR}

tmpdir=$(mktemp -d -p ${TMPLOC})

for chrom in $(seq 1 22)
do
    echo $chrom
    cat ${LOBREF} | awk -v"chrom=$chrom" '($1==chrom)' | cut -f 1-3 > ${tmpdir}/loc_${chrom}.bed
    split -l ${NUMLINES} -d ${tmpdir}/loc_${chrom}.bed ${LOCDIR}/${chrom}"."
done