#!/bin/bash

BASEDIR=/n/scratch2/gymrek
DATADIR=/oasis/projects/nsf/ddp268/mgymrek/mutea-autosomal/sgdp_asdt_vcf
NUMSIM=10
BATCHSIZE=1
PREDFILE=

# Batch loci - TODO

# Simulate new haplotype values based on predicted, automatically batch - TODO
mkdir -p ${BASEDIR}/constraint/batches
./simulate_constraint_nulls.py \
    --asdhet ${DATADIR}/ \
    --pred ${PREDFILE} \
    --numsim ${NUMSIM} \
    --outdir ${BASEDIR}/constraint/batches



