#!/bin/bash

############# Per locus ############
# Run simulations for each condition and get estimates
./run_simulate_locus.sh params.sh
./run_simulate_locus.sh params_withstutter.sh

# Get number of mutations
./get_num_mutations.sh params.sh
cat tree_mutinfo.tab | grep -v chrom | intersectBed -a stdin -b tree_truth.bed -wa -wb | awk '{print $5+$6 "\t" $14}' | sort -g -k 2 | datamash -g 2 median 1

# Calibrate standard errors
./calibrate_errors.sh params.sh
./calibrate_errors.sh params_withstutter.sh

############# Joint ############
./run_simulate_locus.sh params_joint.sh
./run_joint.sh params_joint.sh
./combine_results.sh params_joint.sh | awk -F" " '(NF==8)' | sed 's/ loci//' > joint_estimates.tab
