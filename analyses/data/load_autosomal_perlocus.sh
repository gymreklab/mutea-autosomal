#!/bin/bash

source params.sh

#scp ${HOST}:${HOSTBASEDIR}/autosomal_estimates/perlocus/autosomal_estimates_ml.bed.gz ${DATADIR}/autosomal_estimates/perlocus/
#scp ${HOST}:${HOSTBASEDIR}/autosomal_estimates/perlocus/autosomal_stutter_ml.bed.gz ${DATADIR}/autosomal_estimates/perlocus/
#scp ${HOST}:${HOSTBASEDIR}/autosomal_estimates/perlocus/autosomal_estimates_ml_filtered.bed.gz ${DATADIR}/autosomal_estimates/perlocus/

scp ${HOST}:${HOSTBASEDIR}/autosomal_estimates/*v2*.gz ${DATADIR}/autosomal_estimates/perlocus/

