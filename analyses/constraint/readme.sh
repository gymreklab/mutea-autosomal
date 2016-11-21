#!/bin/bash

##############################################################
# Get genome-wide file with of features + observed rates
./get_features.sh params.sh

cat features/lobSTR_ref_GRCh37_properties_filtered.tab | grep -v chrom | \
    intersectBed -a ../perlocus/autosomal_estimates_ml_filtered.bed.gz -b stdin -wa -wb | \
    cut -f 10-12 --complement > tmp
echo "chrom,start,end,ml_mu,ml_beta,ml_p,ml_mu_stderr,numsamples,strfilter,motif,length,uninterrupted_length,recomb,gc,entropy,reptiming,featurefilter" | \
    sed 's/,/\t/g' > autosomal_perlocus_observed.bed
cat tmp >> autosomal_perlocus_observed.bed
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
    