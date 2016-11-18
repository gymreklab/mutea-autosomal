#!/bin/bash
source params_base.sh

MUTEADIR=${HOME}/workspace/mutea-autosomal/
#VCF=/groups/reich/melissa/mutation_rate_project/str_genotypes/round2/sgdp_memstrs_Y.filt.sorted.vcf.gz
VCF=/oasis/projects/nsf/ddp268/mgymrek/mutea-autosomal/ystr_validation/data/sgdp_memstrs_Y.filt.sorted.vcf.gz
SAMPLEPAIRS=sgdp_sample_pairs.tab
OUTPREFIX=ystrs_sgdp
LOCI=ystr_coords_forest.bed
TRUTHFILE=willems_etal_truth_sgdp.tab
TRUTHFILE2=ystrs_ballantyne_truth_np.bed
SCALE=1.60
