#!/bin/bash

BASEDIR=/n/scratch2/gymrek
DATADIR=/groups/reich/melissa/mutation_rate_project/sgdp_asdt_vcf
NUMSIM=10
BATCHSIZE=1
PREDFILE=test.bed
MUTEADIR=~/workspace/mutea-autosomal

# Run MUTEA on the locus
#${MUTEADIR}/mutea-auto/main_autosomal.py \
#    --asdhet ${DATADIR}/sgdp_asdt_chr2.vcf.gz --vcf \
#    --eststutter test.stutter.tab \
#    --loci test_loci.bed \
#    --out test.estimates.tab

# Simulate new haplotype values based on predicted, automatically batch
./simulate_constraint_nulls.py \
    --asdhet ${DATADIR}/sgdp_asdt_chr2.vcf.gz \
    --pred ${PREDFILE} \
    --numsim ${NUMSIM} \
    --scale ${SCALE} \
    --out test.vcf
bgzip -f test.vcf
tabix -p vcf test.vcf.gz

# Run MUTEA on simulated data
${MUTEADIR}/mutea-auto/main_autosomal.py \
    --asdhet test.vcf.gz --vcf \
    --loci test.vcf.info \
    --out test_sim.estimates.tab

# Aggregate, scale, get pvalue and updated zscore
./aggregate_sims.py \
    test.vcf.info \
    test_sim.estimates.tab \
    test.estimates.tab \
    ${SCALE}
