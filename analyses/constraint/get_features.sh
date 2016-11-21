#!/bin/bash

PARAMFILE=$1
source ${PARAMFILE}

##################################
# Get features for each STR
##################################

tmpdir=$(mktemp -d -p ${TMPLOC})

# Motif, length, Uninterrupted track length
cat ${LOBREF} | ./canon 15 > ${tmpdir}/lobstr_ref_canon.bed
./get_motif_length_features.py ${tmpdir}/lobstr_ref_canon.bed ${REFDB} > ${tmpdir}/lobstr_ref_features.bed

# Recombination, GC, entropy
cat ${PROPFILE} | grep -v entropy | cut -d',' -f 1-3,5,8,9 | sed 's/,/\t/g' | uniq > ${tmpdir}/lobstr_ref_localfeatures.bed

# Replication timing
cat ${TORFILE} | sed 's/chr//' | intersectBed -a ${tmpdir}/lobstr_ref_canon.bed -b stdin -wa -wb | \
    cut -f 1-3,20 | datamash -g1,2,3 mean 4 > ${tmpdir}/lobstr_ref_reptiming.bed

# Combine
echo "chrom,start,end,motif,length,uninterrupted_length,recomb,gc,entropy,reptiming" | \
    sed 's/,/\t/g' > ${BASEDIR}/lobSTR_ref_GRCh37_properties.tab
intersectBed -a ${tmpdir}/lobstr_ref_features.bed -b ${tmpdir}/lobstr_ref_localfeatures.bed -wa -wb -f 1 | \
    cut -f 7-9 --complement | \
    intersectBed -a stdin -b ${tmpdir}/lobstr_ref_reptiming.bed -wa -wb -f 1| \
    cut -f 10-12 --complement >> ${BASEDIR}/lobSTR_ref_GRCh37_properties.tab

# Filter
./filter_features.py ${BASEDIR}/lobSTR_ref_GRCh37_properties.tab > ${BASEDIR}/lobSTR_ref_GRCh37_properties_filtered.tab
