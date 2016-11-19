#!/bin/bash

set -e

LOGMU=$1
lencoeff=$2
PARAMFILE=$3

source $PARAMFILE

# Output filenames
simprefix=${PREFIX}.${LOGMU}.${lencoeff}
OUTFILE=${simprefix}.asdhet.bed
LOCFILE=${simprefix}.loci.bed
FEATUREFILE=${simprefix}.features.bed

rm -f $LOCFILE
rm -f $OUTFILE
rm -f $FEATUREFILE

for sim in $(seq 1 $NUMSIM)
do
    echo $sim
    cat ${DATADIR}/${PREFIX}_logmu${LOGMU}_lencoeff${lencoeff}_${sim}.asd_vs_tmrca.tab | grep -v tmrca | \
	awk -v"sim=${sim}" '{print "Z" "\t" sim "\t" sim+1 "\t" $5  "\t" "s"$2"_"$3 "\t" $10 "\t" $11 "\t" "1" "\t" $4}' >> ${OUTFILE}
    echo "$sim" | awk '{print "Z" "\t" $1 "\t" $1+1}' >> ${LOCFILE}
    echo "$sim" | awk '{print "Z" "\t" $1 "\t" $1+1 "\t" $1}' >> ${FEATUREFILE}
done

# Sort and index bed file and covar file (already sorted)
cat ${OUTFILE} | bgzip -c > ${OUTFILE}.gz
tabix -p bed ${OUTFILE}.gz
