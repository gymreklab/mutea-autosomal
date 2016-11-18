#!/bin/bash
source params_base.sh

MUTEADIR=${HOME}/workspace/mutea-autosomal/
#VCF=/groups/reich/melissa/mutation_rate_project/validation/fromThomas/final_callset.sorted.vcf.gz
VCF=/oasis/projects/nsf/ddp268/mgymrek/mutea-autosomal/ystr_validation/data/final_callset.sorted.vcf.gz
SAMPLEPAIRS=${BASEDIR}/1kg_sample_pairs.tab
OUTPREFIX=ystrs_1kg
LOCI=ystr_coords_forest_chr.bed
TRUTHFILE=willems_etal_truth_1kg.tab
TRUTHFILE2=ystrs_ballantyne_truth_np_chr.bed
SCALE=0.95
