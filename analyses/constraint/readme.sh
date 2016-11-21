#!/bin/bash

##############################################################
# Get genome-wide file with of features + observed rates
./get_features.sh params.sh
./combine_features_estimates.sh params.sh
##############################################################

##############################################################
# Get set of loci to train on
./get_training_loci.sh params.sh
##############################################################

##############################################################
# Build model - TODO
runipy BuildConstraintModel.ipynb
##############################################################

##############################################################
# Simulations to get null - TODO
./simulate_constraint_nulls.sh params.sh
##############################################################
