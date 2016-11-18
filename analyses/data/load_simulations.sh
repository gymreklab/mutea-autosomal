#!/bin/bash

source params.sh

scp ${HOST}:${HOSTBASEDIR}/simulations/tree_maxlik.tab ${DATADIR}/simulations/
scp ${HOST}:${HOSTBASEDIR}/simulations/tree_truth.bed ${DATADIR}/simulations/
scp ${HOST}:${HOSTBASEDIR}/simulations/tree_strsd.bed ${DATADIR}/simulations/
scp ${HOST}:${HOSTBASEDIR}/simulations/joint_estimates.tab ${DATADIR}/simulations/
