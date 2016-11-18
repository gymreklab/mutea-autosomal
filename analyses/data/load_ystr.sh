#!/bin/bash

source params.sh

scp ${HOST}:${HOSTBASEDIR}/ystr_validation/ystrs_sgdp_ml.tab ${DATADIR}/ystr_validation
scp ${HOST}:${HOSTBASEDIR}/ystr_validation/ystrs_1kg_ml.tab ${DATADIR}/ystr_validation

