#!/bin/bash

source params.sh

scp ${HOST}:${HOSTBASEDIR}/autosomal_validation/codis_estimates_ml* ${DATADIR}/autosomal_validation/
scp ${HOST}:${HOSTBASEDIR}/autosomal_validation/sunetal_estimates_ml*  ${DATADIR}/autosomal_validation/

