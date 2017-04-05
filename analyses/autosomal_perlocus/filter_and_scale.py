#!/usr/bin/env python

LOBREF="/storage/resources/datasets/lobstr/hg19/lobstr_v3.0.2_hg19_ref_nochr.bed"
import argparse
import os
import numpy as np
import pandas as pd
import vcf
import sys

SCALE=0.472423187416 # Scale mutation rate
GAMMA=1.2 # Scale standard error

def LoadMLData(datafile):
    data = pd.read_csv(datafile, sep="\t", names=["chrom","start","end","est_logmu_ml","est_beta_ml","est_pgeom_ml","est_logmu_stderr","numsamples_ml"])
    data["est_beta_eff_ml"] = data.apply(lambda x: x.est_beta_ml/((2-x.est_pgeom_ml)/x.est_pgeom_ml**2), 1)
    data["start"] = data["start"].apply(int)
    data["end"] = data["end"].apply(int)
    data["chrom"] = data["chrom"].apply(int)
    return data[["chrom","start","end","est_logmu_ml","est_logmu_stderr", "est_beta_eff_ml","est_beta_ml", "est_pgeom_ml", "numsamples_ml"]]

def LoadVCFReaders(vcfdir):
    vcffiles = [item for item in os.listdir(vcfdir) if item.endswith(".vcf.gz")]
    readers = []
    for vf in vcffiles:
        reader = vcf.Reader(open(os.path.join(vcfdir, vf), "rb"))
        readers.append(reader)
    return readers

def GetVCFFilter(chrom, start, end, readers):
    for reader in readers:
        try:
            record = reader.fetch(str(int(chrom)), start, end)
            for r in record:
                if len(r.FILTER) > 1:
                    return ":".join(r.FILTER)
                else: return "PASS"
        except ValueError: continue
    return "None"

def main():
    parser = argparse.ArgumentParser("script to filter and scale per locus estimates")
    parser.add_argument("--ests", help="Input autosomal estimates", required=True, type=str)
    parser.add_argument("--ref", help="lobSTR strinfo", required=False, default=LOBREF)
    parser.add_argument("--vcf", help="Directory of VCF files. If present, annotate output with VCF filter.", type=str, required=False)
    args = parser.parse_args()

    # Load autosomal estimates
    data = LoadMLData(args.ests)
    sys.stderr.write("before filter %s\n"%data.shape[0])

    # Load VCF readers per chromosome
    if args.vcf:
        vcfreaders = LoadVCFReaders(args.vcf)
        data["vcffilter"] = data.apply(lambda x: GetVCFFilter(x.chrom, x.start, x.end, vcfreaders), 1)

    # Filter to contain only perfect repeats - REMOVED later use uninterrupted tract length
    data = data #pd.merge(data, lobperfect, on=["chrom","start","end"])
    sys.stderr.write("after filter perfect %s\n"%data.shape[0])

    # Set stderr nan to -1 so we know who they are
    data.ix[np.isnan(data["est_logmu_stderr"]), "est_logmu_stderr"] = -1 # set to -1, will filter later
    data.ix[np.isinf(data["est_logmu_stderr"]), "est_logmu_stderr"] = -1 # set to -1, will filter later
    sys.stderr.write("after filter stderr %s\n"%data.shape[0])

    # Scale mutation rates using CODIS scaling factor
    data["est_logmu_ml"] = data["est_logmu_ml"] - np.log10(SCALE)

    # Scale stderrs
    data["est_logmu_stderr"] = data.apply(lambda x: x.est_logmu_stderr*abs(x.est_logmu_ml)*GAMMA, 1)

    # Output results
    data[["chrom","start","end","est_logmu_ml","est_beta_ml","est_pgeom_ml","est_logmu_stderr","numsamples_ml","vcffilter"]].to_csv(sys.stdout, sep="\t", index=False, header=False)

if __name__ == "__main__":
    main()
