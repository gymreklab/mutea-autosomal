#!/usr/bin/env python
"""
Draw samples from predicted mutation rates
Simulate haplotypes, output vcfpy file
"""

import argparse
import glob
import pandas as pd
import vcf

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--asdhet", help="Directory containing asdhet files", required=True, type=str)
    parser.add_argument("--pred", help="File with predicted + observed values for each locus", required=True, type=str)
    parser.add_argument("--numsim", help="Number of simulations for each locus", required=True, type=int)
    parser.add_argument("--batchsize", help="Number of loci to include in each batch output", required=True, type=int)
    parser.add_argument("--outdir", help="Output directory", required=True, type=str)
    args = parser.parse_args()

    # Get VCF data
    asdhet = glob.glob(args.asdhet + "/*.vcf.gz")

    # Load predictions
    pred = pd.read_csv(args.pred, sep="\t")

    # For each predicted locus, do simulations
    batchnum = 0
    numinbatch = 0
    for i in range(pred.shape[0]):
        RunLocusSimulations(pred.chrom.values[i], pred.start.values[i], pred.end.values[i], \
                                asdhet, args.numsim, args.outdir)

if __name__ == "__main__":
    main()
