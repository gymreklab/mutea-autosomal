#!/bin/bash

##############################################################
# Get genome-wide file with of features + observed rates
./get_features.sh params.sh

./combine_features_estimates.sh params.sh
##############################################################

# Get set of loci to train on
head -n 1 autosomal_perlocus_observed.bed > autosomal_perlocus_train_intergenic.bed
cat autosomal_perlocus_observed.bed | grep -v start | \
    intersectBed -a stdin -b ../../../autosomal_byclass/explore_classes/lobSTR_ref_GRCh37_intergenic.bed >> \
    autosomal_perlocus_train_intergenic.bed

gzip autosomal_perlocus_observed.bed
gzip autosomal_perlocus_train_intergenic.bed
# Train model on intergenic - see randomforest

# Get Z-scores per gene
#./calc_zscores.py autosomal_perlocus_observed.bed randomforest/autosomal_model_rf > autosomal_perlocus_zscores.bed
    