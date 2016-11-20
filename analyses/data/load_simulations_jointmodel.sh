#!/bin/bash

source params.sh

scp ${HOST}:${HOSTBASEDIR}/simulations_joint/simulation_joint_results.txt ${DATADIR}/simulations_joint/
