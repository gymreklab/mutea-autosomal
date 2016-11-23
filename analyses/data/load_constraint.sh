#!/bin/bash

source params.sh

scp ${HOST}:${HOSTBASEDIR}/constraint/lobSTR_ref_GRCh37_properties_filtered.tab.gz ${DATADIR}/constraint/
scp ${HOST}:${HOSTBASEDIR}/constraint/autosomal_perlocus_estimates.bed ${DATADIR}/constraint/

