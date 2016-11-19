#!/bin/bash

PARAMFILE=$1
source ${PARAMFILE}

# Estimates
echo "gathering estimates..."
cat ${LOGDIR}/*.err | grep PROGRESS | awk '(NF==9)' | sed 's/PROGRESS: //' | grep -v "PROGRESS: PROGRESS:" | \
    bedtools sort -i stdin | uniq | bgzip -c > ${BASEDIR}/autosomal_estimates_ml.bed.gz
tabix -p bed ${BASEDIR}/autosomal_estimates_ml.bed.gz

# Stutter
echo "gathering stutter..."
cat ${LOGDIR}/*.err | \
    grep -v RuntimeWarning | grep -v fisher | grep -v "Unable" | grep -v "logprobs" | grep -v weights | grep -v "numpy\.max" | \
    grep "Estimated parameters:" -A 1 | \
    awk 'NR%3{printf "%s ",$0;next;}1' | grep PROGRESS | grep -v Estimating | \
    sed 's/Estimated parameters: //' | \
    sed 's/PROGRESS: //' | \
    sed 's/P_GEOM=//' | sed 's/P_DOWN=//' | sed 's/P_UP=//' | sed 's/,//g' | \
    awk '{print $4 "\t" $5 "\t" $6 "\t" $3 "\t" $2 "\t" $1}' | \
    awk '(NF==6)' | \
    bedtools sort -i stdin | uniq | bgzip -c > ${BASEDIR}/autosomal_stutter_ml.bed.gz
tabix -p bed ${BASEDIR}/autosomal_stutter_ml.bed.gz