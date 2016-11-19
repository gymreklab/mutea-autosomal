#!/bin/bash

source params.sh

scp ${HOST}:${HOSTBASEDIR}/simulations/tree_maxlik.tab ${DATADIR}/simulations/
scp ${HOST}:${HOSTBASEDIR}/simulations/tree_truth.bed ${DATADIR}/simulations/
scp ${HOST}:${HOSTBASEDIR}/simulations/tree_strsd.bed ${DATADIR}/simulations/
scp ${HOST}:${HOSTBASEDIR}/simulations/tree_stutter.tab ${DATADIR}/simulations/
scp ${HOST}:${HOSTBASEDIR}/simulations/tree_calibrate_errors.tab ${DATADIR}/simulations/
scp ${HOST}:${HOSTBASEDIR}/simulations/treestutter_calibrate_errors.tab ${DATADIR}/simulations/
scp ${HOST}:${HOSTBASEDIR}/simulations/joint_estimates.tab ${DATADIR}/simulations/

scp ${HOST}:${HOSTBASEDIR}/simulations/treestutter_maxlik.tab ${DATADIR}/simulations/
scp ${HOST}:${HOSTBASEDIR}/simulations/treestutter_truth.bed ${DATADIR}/simulations/
scp ${HOST}:${HOSTBASEDIR}/simulations/treestutter_strsd.bed ${DATADIR}/simulations/
scp ${HOST}:${HOSTBASEDIR}/simulations/treestutter_maxlik_calcstutter.tab ${DATADIR}/simulations/
scp ${HOST}:${HOSTBASEDIR}/simulations/treestutter_stutter.tab ${DATADIR}/simulations/
