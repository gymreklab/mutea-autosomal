#!/bin/bash

PARAMFILE=$1
source $PARAMFILE

set -e

die()
{
    BASE=$(basename -- "$0")
    echo "$BASE: error: $@">&2
    exit 1
}

rm -f ${OUTFILE}
rm -f ${VCFFILE}
rm -f ${TRUTHFILE}
rm -f ${LOCFILE}
rm -f ${STRSDFILE}
rm -f ${OUTFILE}.gz
rm -f ${OUTFILE}.gz.tbi
rm -f ${COVARFILE}.gz
rm -f ${COVARFILE}.gz.tbi
rm -f ${SUMMFILE}

# Get VCF header
cat vcf_header_template.vcf > ${VCFFILE}
echo "#CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT"  $(seq 1 $(echo "$SAMPLES/2" | bc -l)) | sed 's/ /\t/g' | sed 's/,/\t/g' >> ${VCFFILE}

echo "chrom,start,end,treenum,mutrate,strsd,lenc,num_mutations,total_branch_length,faction_boundary,eff_mut_rate" | sed 's/,/\t/g' > ${SUMMFILE}
locstart=1
for STEPPARAM in ${STEPPARAMS}
do
    for STR_MUT_RATE in ${MUTS}
    do
	for LENC in ${LENCS}
	do
	    echo ${STEPPARAM} ${STR_MUT_RATE} ${LENC}
	    # Get ASDHET - deprecated
	    cat ${DATADIR}/${PREFIX}_sim_mutrate${STR_MUT_RATE}_stepparam${STEPPARAM}_lenc${LENC}.asd_vs_tmrca.tab | grep -v tmrca | \
		awk -v"start=$locstart" -v"numcalls=$numcalls" \
		'{print "Z" "\t" start*1000+$1 "\t" start*1000+1+$1 "\t" $5  "\t" "s"$2"_"$3 "\t" $10 "\t" $11 "\t" "1" "\t" $4}' >> ${OUTFILE}
	    # Get VCF
	    python2.7 ${SCRIPTSDIR}/GenerateVCF.py \
		--asdhet ${DATADIR}/${PREFIX}_sim_mutrate${STR_MUT_RATE}_stepparam${STEPPARAM}_lenc${LENC}.asd_vs_tmrca.tab \
		--stutter ${STUTTER} \
		--locstart ${locstart} \
		--coverage ${COVERAGE} ${VCFARGS} | grep -v "\#" >> ${VCFFILE}
	    # Get truth
	    cat ${DATADIR}/${PREFIX}_sim_mutrate${STR_MUT_RATE}_stepparam${STEPPARAM}_lenc${LENC}.asd_vs_tmrca.tab | grep -v tmrca | \
		cut -f 1 | sort | uniq | \
		awk -v"start=${locstart}" -v"mu=${STR_MUT_RATE}" -v"beta=${LENC}" \
		'{print "Z" "\t" start*1000+$1 "\t" start*1000+1+$1 "\t" mu "\t" beta}' >> ${TRUTHFILE}
            # Get locfile
	    cat ${DATADIR}/${PREFIX}_sim_mutrate${STR_MUT_RATE}_stepparam${STEPPARAM}_lenc${LENC}.asd_vs_tmrca.tab | grep -v tmrca | \
		cut -f 1 | sort | uniq | \
		awk -v"start=${locstart}" \
		'{print "Z" "\t" start*1000+$1 "\t" start*1000+1+$1}' >> ${LOCFILE} 
	    # Get strsdfile
	    strsd=$(${SCRIPTSDIR}/get_strsd.py ${STEPPARAM})
	    cat ${DATADIR}/${PREFIX}_sim_mutrate${STR_MUT_RATE}_stepparam${STEPPARAM}_lenc${LENC}.asd_vs_tmrca.tab | grep -v tmrca | \
		cut -f 1 | sort | uniq | \
		awk -v"start=${locstart}" -v"strsd=$strsd" \
		'{print "Z" "\t" start*1000+$1 "\t" start*1000+1+$1 "\t" strsd}' >> ${STRSDFILE}
	    # Get summary file
	    cat ${DATADIR}/${PREFIX}_sim_mutrate${STR_MUT_RATE}_stepparam${STEPPARAM}_lenc${LENC}.summary.tab | grep -v treenum | \
		awk -v"start=${locstart}" '{print "Z" "\t" start*1000+$1 "\t" start*1000+1+$1 "\t" $0}' >> ${SUMMFILE}
	    locstart=$((locstart+1))
	done
    done
done

# Sort and index bed file and covar file (already sorted)
cat ${OUTFILE} | bgzip -c > ${OUTFILE}.gz
tabix -p bed ${OUTFILE}.gz
cat ${VCFFILE} | vcf-sort | bgzip -c >${VCFFILE}.gz
tabix -p vcf ${VCFFILE}.gz
