#!/usr/bin/env python

"""
Estimate autosomal mutation rates using MUTEA
"""

import argparse
import glob
import joblib
import numpy as np
import os
import scipy.optimize
import scipy.stats
import sys

import locus
import jointlocus

# For profiling
from pycallgraph import PyCallGraph
from pycallgraph.output import GraphvizOutput

def LoadJointPriors(priorfile):
    priors = map(lambda x: float(x.strip()), open(priorfile, "r").readlines())
    return priors

def LoadLocusFeatures(locusfeatures, locfile, use_features=None):
    if locusfeatures is None: return None
    # Filter to contain loci in locfile
    loci = []
    with open(locfile, "r") as f:
        for line in f:
            chrom, start, end = line.strip().split()[0:3]
            chrom = str(chrom)
            start = int(start)
            end = int(end)
            loci.append((chrom, start, end))
    features = {}
    numfeatures = 0
    with open(locusfeatures, "r") as f:
        for line in f:
            items = line.strip().split()
            chrom = items[0]
            start = int(items[1])
            end = int(items[2])
            fvals = map(float, items[3:])
            if use_features is not None:
                fvals = [fvals[i] for i in map(int, use_features.split(","))]
            numfeatures = len(fvals)
            if (chrom, start, end) in loci:
                features[(chrom, start, end)] = fvals
    sys.stderr.write("Loaded features for %s loci...\n"%len(features.keys()))
    if numfeatures == 0: return features
    # Center each feature
    feature_means = []
    for i in range(numfeatures):
        feature_means.append(np.mean([features[key][i] for key in features]))
    for key in features:
        features[key] = [features[key][i]-feature_means[i] for i in range(numfeatures)]
    sys.stdout.write("# Feature means: %s\n"%",".join(map(str, feature_means)))
    sys.stdout.flush()
    return features

def LoadPriors(locuspriors, use_locus_means, locfile):
    if locuspriors is None: return None
    # Filter to contain loci in locfile
    loci = []
    with open(locfile, "r") as f:
        for line in f:
            chrom, start, end = line.strip().split()[0:3]
            chrom = str(chrom)
            start = int(start)
            end = int(end)
            loci.append((chrom, start, end))
    priors = {}
    with open(locuspriors, "r") as f:
        for line in f:
            chrom, start, end, logmu, beta, pgeom = line.strip().split()[0:6]
            start = int(start)
            end = int(end)
            logmu = float(logmu)
            beta = float(beta)
            pgeom = float(pgeom)
            if (chrom, start, end) in loci:
                priors[(chrom, start, end)] = (logmu, beta, pgeom)
    sys.stderr.write("Loaded priors for %s loci...\n"%len(priors.keys()))
    if use_locus_means:
        mean_beta = np.mean([value[1] for value in priors.values()])
        mean_p = np.mean([value[2] for value in priors.values()])
        for locus in priors: priors[locus] = (priors[locus][0], mean_beta, mean_p)
    return priors

def LoadStutter(sfile):
    if sfile is None: return None
    st = {}
    with open(sfile, "r") as f:
        for line in f:
            if len(line.strip().split()) != 6: continue
            chrom, start, end, up, down, pgeom  = line.strip().split()
            start = int(start)
            end = int(end)
            up = float(up)
            down = float(down)
            pgeom = float(pgeom)
            st[(chrom, start, end)] = (up, down, pgeom)
    return st

def LoadLoci(locfile, datafiles, minsamples, maxsamples, \
                 locuspriors, locusfeatures, use_features, \
                 stderrs, jackknife_blocksize, isvcf, eststutter, usestutter, \
                 use_sample_pairs, use_locus_means, usesmm, \
                 debug=False):
    sys.stderr.write("Loading priors\n")
    priors = LoadPriors(locuspriors, use_locus_means, locfile)
    sys.stderr.write("Loading features\n")
    features = LoadLocusFeatures(locusfeatures, locfile, use_features=use_features)
    sys.stderr.write("Loading stutter\n")
    stutter = LoadStutter(usestutter)
    loci = []
    with open(locfile, "r") as f:
        for line in f:
            chrom, start, end = line.strip().split()[0:3]
            chrom = str(chrom)
            start = int(start)
            end = int(end)
            loc = locus.Locus(chrom, start, end, datafiles, minsamples, maxsamples, \
                                  stderrs, isvcf, eststutter, _debug=debug)
            loc.jkblocksize = jackknife_blocksize
            loc.usesmm = usesmm
            key = (chrom, start, end)
            if priors is not None and key not in priors: continue
            if features is not None and key not in features: continue
            if stutter is not None and key not in stutter: continue
            if priors is not None and key in priors:
                loc.LoadPriors(*priors[key])
            if features is not None and key in features:
                loc.LoadFeatures(features[key])
                sys.stderr.write("%s,%s\n"%(str(key), str(features[key])))
            if use_sample_pairs is not None:
                loc.use_sample_pairs = True
                loc.pair_file = use_sample_pairs
            if stutter is not None and key in stutter:
                loc.LoadStutter(*stutter[key])
            loci.append(loc)
    return loci

def ERROR(msg):
    sys.stderr.write(msg.strip() + "\n")
    sys.exit(1)
    
def MSG(msg):
    sys.stderr.write(msg.strip() + "\n")

def RunLocus(locus, args=None):
    if args.only_stutter:
        locus.LoadData()
    elif locus.usesmm:
        locus.SMM(debug=args.debug)
    else:
        locus.MaximizeLikelihood(mu_bounds=(args.min_mu, args.max_mu), \
                                     beta_bounds=(args.min_beta, args.max_beta), \
                                     pgeom_bounds=(args.min_pgeom, args.max_pgeom), \
                                     lencoeff_bounds=(args.min_lencoeff, args.max_lencoeff), \
                                     debug=args.debug)
    # Print intermediate results to stderr so we can track
    outline = locus.GetResultsString(include_likelihood=args.output_likelihood)
    # Print stutter results
    if locus.eststutter is not None:
        stutterline = locus.GetStutterString()
    else: stutterline = None
    sys.stderr.write("PROGRESS: " + outline)
    return outline, stutterline

def main():
    ############################
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--debug", help="Print helpful debug messages", action="store_true")
    parser.add_argument("--profile", help="Profile code performance", action="store_true")
    parser.add_argument("--maxloci", help="Only process this many loci (for debugging)", default=np.inf, type=int)
    parser.add_argument("--numproc", help="Number of processors", default=1, type=int)

    # I/O data options
    parser.add_argument("--asdhet", help="ASD-Het file. Must be indexed bed file. See help for columns.", type=str, required=True) 
    parser.add_argument("--asdhetdir", help="Is asdhet a directory", action="store_true")
    parser.add_argument("--vcf", help="Input is a VCF file (not asdhet)", action="store_true")
    parser.add_argument("--loci", help="Bed file with loci to process. First three columns are chrom, start, end", type=str, required=True)
    parser.add_argument("--out", help="Output file (default stdout)", type=str, required=False)
    parser.add_argument("--use-sample-pairs", help="Run estimtes on pairs of haploids (e.g. Y-STRs). Input sample1,sample2,tmrca", type=str, required=False)

    # Stutter estimation options
    parser.add_argument("--eststutter", help="Estimate stutter model from reads", type=str)
    parser.add_argument("--usestutter", help="Use previously estimated stutter model", type=str)
    parser.add_argument("--only-stutter", help="Only estimate stutter noise.", action="store_true")

    # Filtering options
    parser.add_argument("--min_samples", help="Don't process if less than this many samples", type=int, default=50)
    parser.add_argument("--max_samples", help="Only process this many samples (for debugging)", type=int, default=1000000)

    # Per-locus Estimation options
    parser.add_argument("--min_mu", required=False, type=float, default=0.00000001, help="Lower optimization boundary for mu.")
    parser.add_argument("--max_mu", required=False, type=float, default=0.05, help="Upper optimization boundary for mu.")
    parser.add_argument("--min_pgeom", required=False, type=float, default=0.7, help="Lower optimization boundary for pgeom.")
    parser.add_argument("--max_pgeom", required=False, type=float, default=1.0, help="Upper optimization boundary for pgeom.")
    parser.add_argument("--min_beta", required=False, type=float, default=0.0, help="Lower optimization boundary for beta.")
    parser.add_argument("--max_beta", required=False, type=float, default=0.9, help="Upper optimization boundary for beta.")
    parser.add_argument("--stderrs", required=False, type=str, default="fisher", help="Method to calc stderrs. Options: fisher, jackknife, both.")
    parser.add_argument("--jackknife_blocksize", required=False, type=int, default=10, help="Jackknife block size.")
    parser.add_argument("--min_lencoeff", required=False, type=float, default=0.0, help="Lower optimization boundary for lencoeff.")
    parser.add_argument("--max_lencoeff", required=False, type=float, default=0.0, help="Upper optimization boundary for lencoeff.")
    parser.add_argument("--smm", help="Assume an SMM model, perform simple linear fit", action="store_true")

    # Joint estimation options
    parser.add_argument("--joint", help="Estimate parameters jointly across loci", action="store_true")
    parser.add_argument("--locus_priors", help="Per-locus results to use as priors. Default: draw from uniform", type=str)
    parser.add_argument("--locus_features", help="Tab file of chrom, start, end, feature1, feature2, ...", type=str)
    parser.add_argument("--joint_priors", help="File with priors. One line per parameter.", type=str)
    parser.add_argument("--use_features", help="List of feature numbers to use. (default all)", type=str)
    parser.add_argument("--use_locus_means", help="Use mean of per-locus beta and p, rather than individual locus estimates", action="store_true")

    # Other options
    parser.add_argument("--output-likelihood", help="Add a column to the output with the model likelihood", action="store_true")
    args = parser.parse_args()
    ############################
    # Check options
    if not os.path.exists(args.asdhet): ERROR("%s does not exist"%args.asdhet)
    if not os.path.exists(args.loci): ERROR("%s does not exist"%args.loci)
    if args.numproc < 1: ERROR("%s must be >=1"%args.numproc)
    if args.locus_priors is not None and not os.path.exists(args.locus_priors): ERROR("%s does not exist"%args.locus_priors)
    if args.locus_features is not None and not os.path.exists(args.locus_features): ERROR("%s does not exist"%args.locus_features)
    ############################

    # Get list of asdhet files
    asdhet = []
    if args.asdhetdir:
        if args.vcf:
            asdhet = glob.glob(args.asdhet + "/*.vcf.gz")
        else: asdhet = glob.glob(args.asdhet + "/*.bed.gz")
    else: asdhet = [args.asdhet]

    # Load loci to process
    loci = LoadLoci(args.loci, asdhet, args.min_samples, args.max_samples, \
                        args.locus_priors, args.locus_features, args.use_features, \
                        args.stderrs, args.jackknife_blocksize, args.vcf, args.eststutter, \
                        args.usestutter, \
                        args.use_sample_pairs, args.use_locus_means, args.smm, \
                        debug=args.debug)
    if len(loci) > args.maxloci: loci = loci[:args.maxloci]

    # Get output
    if args.out is None: output = sys.stdout
    else: output = open(args.out, "w")
    if args.eststutter is not None: stutter = open(args.eststutter, "w")

    # Load priors
    if args.joint and args.joint_priors is not None:
        jpriors = LoadJointPriors(args.joint_priors)
    else: jpriors = None

    ####### Profiling #####
    if args.profile:
        with PyCallGraph(output=GraphvizOutput()):
            for locus in loci: RunLocus(locus, args=args)
        sys.exit(1)
    #######################

    # Run estimation
    if args.joint:
        jlocus = jointlocus.JointLocus(loci, _ires=None, _numproc=args.numproc)
        jlocus.SetPriors(jpriors)
        MSG("[main_autosomal.py] Processing joint locus with %s loci... (%s max)"%(len(jlocus.loci), args.maxloci))
        jlocus.MaximizeLikelihood(mu_bounds=(args.min_mu, args.max_mu), \
                                      sd_bounds=(None, None), \
                                      beta_bounds=(args.min_beta, args.max_beta), \
                                      pgeom_bounds=(args.min_pgeom, args.max_pgeom), \
                                      lencoeff_bounds=(args.min_lencoeff, args.max_lencoeff), \
                                      drawmu=False, \
                                      debug=args.debug)
        jlocus.PrintResults(output)
    else:
        outlines = joblib.Parallel(n_jobs=args.numproc, verbose=50)(joblib.delayed(RunLocus)(locus, args=args) for locus in loci)
        for res in outlines:
            l, s = res
            if not args.only_stutter:
                output.write(l)
                output.flush()
            if s is not None and args.eststutter is not None:
                stutter.write(s)
                stutter.flush()

if __name__ == "__main__":
    main()


