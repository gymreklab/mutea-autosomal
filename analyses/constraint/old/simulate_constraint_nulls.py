#!/usr/bin/env python
"""
Draw samples from predicted mutation rates
Simulate haplotypes, output vcfpy file
"""

import argparse
import pandas as pd
import numpy as np
import sys
import vcf

sys.path.append("../simulations/")
import STRMutationGenerator

def ERROR(msg):
    sys.stderr.write("[simulate_constraint_nulls.py] " + msg.strip() + "\n")

def RunLocusSimulations(chrom, start, end, vcf_reader, numsim, pred_data, vcf_writer, info_writer, scale=1):
    ERROR("%s:%s-%s"%(chrom, start, end))
    # Fetch VCF record
    try:
        record = vcf_reader.fetch(chrom, start)
    except: return
    ident = "."
    # Get VCF info common to all records
    REF = str(record.REF)
    FORMAT = "GT:GB:TMRCA"
    for i in range(numsim):
        INFO = "MOTIF=%s;END=%s"%(record.INFO["MOTIF"], i+(record.INFO["END"]))
        alt_alleles = []
        # Get mutation generator
#        beta = pred_data["ml_beta"]
#        p = pred_data["ml_p"]
        # Just set beta and p to reasonable values
        beta = 0.3
        p = 0.8 
        pred_mu = pred_data["pred_mu"]
        pred_mu_se = pred_data["pred_mu_se"]
        mu_draw = np.random.normal(loc=pred_mu, scale=pred_mu_se) + np.log10(scale) # Scale back
        smg = STRMutationGenerator.STRMutationGenerator(None, str_mut_rate=10**mu_draw, \
                                                            str_step_param=p, method="discreteMultiOU", \
                                                            length_constraint=beta)
        # Get sample data
        SAMPLES = []
        for sample in record:
            tmrca = sample["TMRCA"]
            if tmrca == -1: # Will ignore these
                gt = [0, 0]
                gb = [0, 0]
            else:
                a1, a2 = SimulateMutations(tmrca, smg)
                gb = [a1, a2]
                gt = []
                if a1 != 0:
                    if a1 not in alt_alleles: alt_alleles.append(a1)
                    gt.append(alt_alleles.index(a1))
                else: gt.append(0)
                if a2 != 0:
                    if a2 not in alt_alleles: alt_alleles.append(a2)
                    gt.append(alt_alleles.index(a2))
                else: gt.append(0)
            sdata = "%s:%s:%s"%("/".join(map(str, gt)), "/".join(map(str, gb)), tmrca)
            SAMPLES.append(sdata)
        ALT = GetAltAlleles(alt_alleles, REF, record.INFO["MOTIF"])
        vcf_writer.write("\t".join(map(str, [chrom, i+start, ident, REF, ",".join(ALT), 0, "PASS", INFO, FORMAT]+SAMPLES))+"\n")
        info_writer.write("\t".join(map(str, [chrom, start+i, end+i, i, start, end]))+"\n")

def SimulateMutations(tmrca, smg):
    a1 = smg.GetSTRForT(tmrca)
    a2 = smg.GetSTRForT(tmrca)
    return a1, a2

# Just has to be correct length...
def GetAltAlleles(lengths, ref, motif):
    alleles = []
    for l in lengths:
        if l > 0: newlen = ref + motif*(l/len(motif))
        else: newlen = ref[0:-1*l]
        alleles.append(newlen)
    return alleles

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--asdhet", help="Asdhet file", required=True, type=str)
    parser.add_argument("--pred", help="File with predicted + observed values for each locus", required=True, type=str)
    parser.add_argument("--numsim", help="Number of simulations for each locus", required=True, type=int)
    parser.add_argument("--out", help="Output file", required=True, type=str)
    parser.add_argument("--scale", help="Scale mutation rate", required=True, type=float)
    args = parser.parse_args()

    # Load predictions
    pred = pd.read_csv(args.pred, sep="\t")

    # Get reader and writer
    vcf_reader = vcf.Reader(open(args.asdhet, "rb"))
    vcf_writer = vcf.Writer(open(args.out, "w"), vcf_reader)
    vcf_writer.close()
    vcf_writer = open(args.out, "a")
    info_writer = open(args.out + ".info", "w")

    # For each predicted locus, do simulations
    for i in range(pred.shape[0]):
        RunLocusSimulations(pred.chrom.values[i], pred.start.values[i], pred.end.values[i], \
                                vcf_reader, args.numsim, pred.iloc[i,:], vcf_writer, info_writer, scale=args.scale)
    vcf_writer.close()
    info_writer.close()

if __name__ == "__main__":
    main()
