#!/bin/bash
source params_base.sh

############# Per locus ############
# Run simulations for each condition and get estimates
./run_simulate_locus.sh params.sh
./run_simulate_locus.sh params_withstutter.sh

# Calibrate standard errors
./calibrate_errors.sh params.sh
./calibrate_errors.sh params_withstutter.sh

############# Joint ############
./run_simulate_locus.sh params_joint.sh
./run_joint.sh params_joint.sh
./combine_results.sh params_joint.sh | awk -F" " '(NF==8)' | sed 's/ loci//' > ${BASEDIR}/simulations/joint_estimates.tab
