#!/usr/bin/env python
"""
Estimate per locus mutation rates

asdhet has columns:
 chrom
 start
 end
 tmrca (posterior estimate)
 sample
 a1
 a2
 period
 asd (not corrected for period)
This file must be a .bed.gz file that has been indexed by tabix
"""

import argparse
import os
import glob
import pandas as pd
import numpy as np
import sys
import tabix
import ModelEstimatorTMRCA
import multiprocessing

debug = False
CURVEFIT="curvefit"
MAXLIK="maxlik"
BNLR="bnlr"
ESTIMATORS = [CURVEFIT, MAXLIK, BNLR]
LINEERROR = "error"

def ERROR(msg):
    sys.stderr.write(msg.strip() + "\n")
    sys.exit(1)
    
def MSG(msg):
    sys.stderr.write(msg.strip() + "\n")
 
def GetSinglePairSamples(samples):
    """
    Prune samples for pairwise chrY analysis so that each 
    individual is present in only a single pair
    """
    # example sample pair (unfortunate use of "_" to split...): LP6005442-DNA_E07_LP6005443-DNA_H05
    pairs_to_use = set()
    used_samples = set()
    for s in samples:
        s1 = "_".join(s.split("_")[0:2])
        s2 = "_".join(s.split("_")[2:])
        if s1 not in used_samples and s2 not in used_samples:
            pairs_to_use.add(s)
            used_samples.add(s1)
            used_samples.add(s2)
    return pairs_to_use

def PrepareData(i, loci, truth, \
                    prior_mu, prior_beta, priors,\
                    strsd, strsds, \
                    asdhet, unit, min_samples, \
                    singlepair, maxT, covar_matrix, \
                    samples_to_use):
    """
    Pull out data needed for estimation
    (chrom, start, end, prior_mu, prior_beta, strsd, asds, tmrcas)
    returns dictionary with all of these things
    """
    locdata = {}
    if i == -1:
        locdata["chrom"] = "ALL"
        locdata["start"] = -1
        locdata["end"] = -1
    else:
        locdata["chrom"] = loci.chrom.values[i]
        locdata["start"] = loci.start.values[i]
        locdata["end"] = loci.end.values[i]
    locdata["true_mu"] = None
    locdata["true_beta"] = None
    if truth is not None:
        x = truth[(truth["chrom"]==locdata["chrom"]) & \
                      (truth["start"]==locdata["start"]) & \
                      (truth["end"]==locdata["end"])]
        if x.shape[0] == 1:
            locdata["true_mu"] = x.true_mu.values[0]
            locdata["true_beta"] = x.true_beta.values[0]
    locdata["prior_mu"] = prior_mu
    locdata["prior_beta"] = prior_beta
    if priors is not None:
        x = priors[(priors["chrom"]==locdata["chrom"]) & \
                       (priors["start"]==locdata["start"]) & \
                       (priors["end"]==locdata["end"])]
        if x.shape[0] == 1:
            locdata["prior_mu"] = x.prior_mu.values[0]
            locdata["prior_beta"] = x.prior_beta.values[0]
    locdata["strsd"] = strsd
    if strsds is not None:
        x = strsds[(strsds["chrom"]==locdata["chrom"]) & \
                       (strsds["start"]==locdata["start"]) & \
                       (strsds["end"]==locdata["end"])]
        if x.shape[0] == 1:
            locdata["strsd"] = x["strsd"].values[0]
    # Get ASD and tmrca
    locdata["tmrcas"] = []
    locdata["asds"] = []
    locdata["genotypes"] = []
    locdata["samples"] = []
    if i == -1:
        records = []
        for j in range(loci.shape[0]):
            chrom = loci.chrom.values[j]
            start = loci.start.values[j]
            end = loci.end.values[j]
            for ah in asdhet:
                try:
                    r = list(ah.query(str(chrom), int(start), int(end)))
                    records.extend(r)
                except: pass
    else:
        records = []
        for ah in asdhet:
            try:
                r = list(ah.query(str(locdata["chrom"]), int(locdata["start"]), int(locdata["end"])))
                records.extend(r)
            except: pass
    if singlepair:
        samples_to_keep = GetSinglePairSamples([r[4] for r in records])
    for record in records:
        if singlepair and record[4] not in samples_to_keep: continue # for singlepair
        if len(samples_to_use) > 0 and record[4] not in samples_to_use: continue # for overall samples list provided
        tmrca_r = int(float(record[3]))
        if tmrca_r > maxT: continue
        a1_r = float(record[5])
        a2_r = float(record[6])
        period_r = int(record[7])
        s_r = record[4]
        if unit and (a1_r%period_r != 0 or a2_r%period_r != 0):
            continue
        locdata["tmrcas"].append(tmrca_r)
        locdata["asds"].append(((a1_r/period_r)-(a2_r/period_r))**2)
        locdata["genotypes"].append((a1_r, a2_r))
        locdata["period"] = period_r
        locdata["samples"].append(s_r)
    # Get covariance matrix
    if str(covar_matrix) == "identity":
        locdata["covar"] = "identity"
    else: 
        covar = np.identity(len(locdata["samples"]))
        for cv in covar_matrix:
            try:
                c = cv.query(str(locdata["chrom"]), int(locdata["start"]), int(locdata["end"]))
                for record in c:
                    ind1 = locdata["samples"].index(record[3])
                    ind2 = locdata["samples"].index(record[4])
                    val = record[5]
                    covar[ind1, ind2] = val
                    covar[ind2, ind1] = val
            except tabix.TabixError: pass
        locdata["covar"] = covar
    if len(locdata["asds"]) < min_samples:
        MSG("Skipping %s:%s-%s. Only %s samples"%(locdata["chrom"], locdata["start"], locdata["end"], len(locdata["asds"])))
        return None
    else:
        MSG("%s:%s-%s, %s samples"%(locdata["chrom"], locdata["start"], locdata["end"], len(locdata["asds"])))
    return locdata

def RunLocusMutinfo(locdata):
    """
    Return:
    output_line
    """
    chrom = locdata["chrom"]
    start = locdata["start"]
    end = locdata["end"]
    period = locdata["period"]
    numsamples = len(locdata["asds"])
    genotypes = locdata["genotypes"]
    tmrcas = locdata["tmrcas"]
    allele_to_count = {}
    for gt in genotypes:
        allele_to_count[gt[0]] = allele_to_count.get(gt[0], 0) + 1
        allele_to_count[gt[1]] = allele_to_count.get(gt[1], 0) + 1
    alleles = [(allele_to_count[a], a) for a in allele_to_count.keys()]
    major_allele = alleles[0]
    num_alt_alleles = 0
    num_recurrent_alleles = 0
    recurrent_alleles = []
    single_alleles = []
    other_alleles = []
    major_allele = alleles[0][1]
    for i in range(1, len(alleles)):
        allele = alleles[i][1]
        tmrcas_hom = [tmrcas[j] for j in range(len(tmrcas)) if allele in genotypes[j] and genotypes[j][0]==genotypes[j][1]]
        tmrcas_het = [tmrcas[j] for j in range(len(tmrcas)) if allele in genotypes[j] and genotypes[j][0]!=genotypes[j][1]]
        label=""
        if len(tmrcas_hom) == 0 or len(tmrcas_het) == 0 or allele%period != 0:
            other_alleles.append(allele)
            label="other"
            if allele%period == 0: num_alt_alleles += 1
        elif max(tmrcas_hom) < min(tmrcas_het):
            single_alleles.append(allele)
            label="single"
            num_alt_alleles += 1
        else:
            num_recurrent_alleles += 1
            recurrent_alleles.append(allele)
            num_alt_alleles += 1
            label="recurrent"
        # LOG
        sys.stderr.write("%s:%s-%s allele=%s, tmrcas_hom=%s, tmrcas_het=%s, label=%s\n"%(chrom, start, end, allele, str(tmrcas_hom), str(tmrcas_het), label))
    return "\t".join(map(str, [chrom, int(start), int(end), numsamples, num_recurrent_alleles, num_alt_alleles, "'"+",".join(map(str, recurrent_alleles))+"'", "'"+",".join(map(str, single_alleles))+"'", "'" + ",".join(map(str, other_alleles)) +"'", major_allele])) + "\n"

def RunLocus(locdata):
    """
    Return:
    output_line,
    [bootstrap_output_lines]
    """
    if locdata["mutinfo"]: return RunLocusMutinfo(locdata), []
    locus = str(locdata["chrom"])+":"+str(locdata["start"])
    if debug: sys.stderr.write("[RunLocus] Running locus %s\n"%locus)
    # Make estimator
    if locdata["estimator_type"] == CURVEFIT:
        estimator = ModelEstimatorTMRCA.LineFitEstimator()
        estimator.SetNumBS(locdata["numBS"])
        estimator.SetNorm(locdata["norm"])
        estimator.SetEstEff(locdata["output_eff"])
        estimator.SetCovar(locdata["covar"])
    elif locdata["estimator_type"] == MAXLIK:
        estimator = ModelEstimatorTMRCA.MLGridEstimator()
    elif locdata["estimator_type"] == BNLR:
        estimator = ModelEstimatorTMRCA.BNLREstimator()
        estimator.options["progress_bar"] = False
        estimator.options["distribution"] = locdata["distribution"]
    else:
        ERROR("%s not implemented"%locdata["estimator_type"])

    # Run estimator
    estimator.SetPriors(mu=locdata["prior_mu"], beta=locdata["prior_beta"])
    estimator.SetParams(strsd=locdata["strsd"])
    estimator.LoadData(asds=locdata["asds"], tmrcas=locdata["tmrcas"])
    if debug: sys.stderr.write("[RunLocus] Estimator %s\n"%locus)
    estimator.Predict()
    # Get output line
    output_line = "\t".join(map(str, [locdata["chrom"], locdata["start"], locdata["end"], len(locdata["asds"]), \
                                          locdata["prior_mu"], locdata["prior_beta"], \
                                          locdata["true_mu"], locdata["true_beta"], estimator.params["strsd"], \
                                          estimator.ests["mu"], estimator.ests_stderrs["mu"], \
                                          estimator.ests["beta"], estimator.ests_stderrs["beta"]]))
    if locdata["estimator_type"] == CURVEFIT and locdata["output_kl"]:
        output_line = output_line + "\t" + "\t".join(map(str, [estimator.ests["k"], estimator.ests_stderrs["k"], \
                                                                   estimator.ests["L"], estimator.ests_stderrs["L"]]))
    if locdata["estimator_type"] == BNLR and locdata["output_ci"]:
        output_line = output_line + "\t" + "\t".join(map(str, [estimator.ests_ci["mu"][0], estimator.ests_ci["mu"][1], \
                                                                   estimator.ests_ci["beta"][0], estimator.ests_ci["beta"][1]]))
    bootstrap_lines = ModelEstimatorTMRCA.GetBootstrapLines(locdata["chrom"], locdata["start"], locdata["end"], \
                                                                estimator.bootstrap)
    estimator.Clear()
    if debug: sys.stderr.write("[RunLocus] Returning line %s\n"%(locus))
    return (output_line + "\n", bootstrap_lines)

class Writer(multiprocessing.Process):
    def __init__(self, results_queue, outfile, bsfile):
        multiprocessing.Process.__init__(self)
        self.results_queue = results_queue
        self.outfile = outfile
        self.bsfile = bsfile
        self.kill_received = False

    def run(self):
        while not self.kill_received:
            try:
                res = self.results_queue.get()
                if res is None:
                    self.results_queue.task_done()
                    break
                next_line, bs_lines = res
                if next_line != LINEERROR:
                    self.outfile.write(next_line)
                    self.outfile.flush()
                if self.bsfile is not None:
                    for line in bs_lines:
                        self.bsfile.write(line)
                self.results_queue.task_done()
            except (KeyboardInterrupt, IOError):
                break

class Consumer(multiprocessing.Process):
    def __init__(self, task_queue, results_queue):
        multiprocessing.Process.__init__(self)
        self.task_queue = task_queue
        self.results_queue = results_queue
        self.kill_received = False

    def run(self):
        while True:
            if debug: sys.stderr.write("[Consumer] waiting for task\n")
            next_locdata = self.task_queue.get()
            if next_locdata is None:
                self.task_queue.task_done()
                break
            try:
                output_line, bs_lines = RunLocus(next_locdata)
            except:
                output_line = LINEERROR
                bs_lines = []
                sys.stderr.write("[Consumer]: Error processing %s:%s\n"%(next_locdata["chrom"], next_locdata["start"]))
            if debug: sys.stderr.write("[Consumer] Got output line\n")
            self.task_queue.task_done()
            self.results_queue.put([output_line, bs_lines])
            if debug: sys.stderr.write("[Consumer] Looking for next task\n")
        return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--debug", help="Print helpful debug messages", action="store_true")
    parser.add_argument("--asdhet", help="ASD-Het file. Must be indexed bed file. See help for columns.", type=str, required=True)
    parser.add_argument("--asdhetdir", help="Is asdhet a directory", action="store_true")
    parser.add_argument("--loci", help="Bed file with loci to process. First three columns are chrom, start, end", type=str, required=True)
    parser.add_argument("--truth", help="Bed file with ground truth data. Columns are chrom, start, end, true_mu, true_beta", type=str, required=False)
    parser.add_argument("--unit", help="Filter calls that are not multiples of the repeat unit", action="store_true")
    parser.add_argument("--estimator", help="Estimator type. Options: %s"%",".join(ESTIMATORS), type=str)
    parser.add_argument("--min-samples", help="Minimum number of samples with calls to use a locus", type=int, default=50)
    parser.add_argument("--prior-mu", help="Prior for mutation rate", type=float, default=0.001)
    parser.add_argument("--prior-beta", help="Prior for beta", type=float, default=0.3)
    parser.add_argument("--priors", help="File with priors for each locus. Format: chrom, start, end, mu_prior, beta_prior", type=str, required=False)
    parser.add_argument("--strsd", help="Step size stddev", type=float, default=1.0)
    parser.add_argument("--strsds", help="File with strsd for each locus. Format: chrom, start, end, strsd", type=str, required=False)
    parser.add_argument("--out", help="Output file (default stdout)", type=str, required=False)
    parser.add_argument("--numBS", help="Number of bootstrap samples to estimate stderrs (curvefit)", type=int, default=100)
    parser.add_argument("--estimate-step", help="Estimate strsd", action="store_true")
    parser.add_argument("--numproc", help="Number of processors to use", type=int, default=1)
    parser.add_argument("--output-kl", help="Output k and L parameters", action="store_true")
    parser.add_argument("--output-eff", help="Output sigma^2*mu and sigma^2/beta rather than mu and beta", action="store_true")
    parser.add_argument("--output-ci", help="Output 95% confidence intervals", action="store_true")
    parser.add_argument("--aggregate", help="Aggregate all loci to one estimate", action="store_true")
    parser.add_argument("--numjobs", help="Only try this many jobs (for debugging)", type=int, default=-1)
    parser.add_argument("--distribution", help="Distribution to use for PyMC (bnlr only)", type=str, default="normal")
    parser.add_argument("--singlepair", help="For chrY analysis, allow each individual to be part of only one pair", action="store_true")
    parser.add_argument("--norm", help="Norm to use for NLS (curvefit only). Default: L2", type=str, default="L2")
    parser.add_argument("--output-bs", help="Output bootstrap samples to this file", type=str)
    parser.add_argument("--maxT", help="Filter points with TMRCA greater than this", type=int, default=100000000)
    parser.add_argument("--mutinfo", help="Output mutation info rather than estimating mutation rates", action="store_true")
    parser.add_argument("--covar-matrix", help="Weight regression. Default: identity. Otherwise file with covariance matrix info.", type=str, default="identity")
    parser.add_argument("--samples", help="Only run on list of samples in this file", type=str)
    args = parser.parse_args()

    if args.debug: debug=True
    # Check inputs
    if args.estimator not in ESTIMATORS and not args.mutinfo:
        ERROR("--estimator must be one of %s"%",".join(ESTIMATORS))
    if not os.path.exists(args.asdhet) and "*" not in args.asdhet:
        ERROR("%s does not exist"%args.asdhet)
    if not os.path.exists(args.loci):
        ERROR("%s does not exist"%args.loci)
    if args.truth and not os.path.exists(args.truth):
        ERROR("%s does not exist"%args.truth)
    if args.output_kl and args.estimator != CURVEFIT:
        ERROR("Can only output kl for %s"%CURVEFIT)
    if args.output_ci and args.estimator != BNLR:
        ERROR("Confidence intervals only implemented for %s"%BNLR)
    if args.norm not in ["L1","L2"]:
        ERROR("Invalid norm type %s"%args.norm)

    # Load loci
    loci = pd.read_csv(args.loci, usecols=range(3), names=["chrom","start","end"], sep="\t")
    MSG("Loading %s loci"%loci.shape[0])

    # Load truth if exists
    truth = None
    if args.truth:
        truth = pd.read_csv(args.truth, sep="\t", names=["chrom","start","end","true_mu","true_beta"])

    # Load priors if exists
    priors = None
    if args.priors:
        priors = pd.read_csv(args.priors, sep="\t", names=["chrom","start","end","prior_mu","prior_beta"])

    # Load strsds if exists:
    strsds = None
    if args.strsds:
        strsds = pd.read_csv(args.strsds, sep="\t", names=["chrom","start","end","strsd"])
    
    # Load list of samples
    if args.samples:
        samples_to_use = [line.strip() for line in open(args.samples, "r").readlines()]
        if len(samples_to_use) < args.min_samples:
            ERROR("Samples list is already less than min number of samples required (%s vs. %s)\n"%(len(samples_to_use), args.min_samples))
    else: samples_to_use = []

    # Get reader for asdhet file(s)
    asdhet = []
    if "*" not in args.asdhet:
        if args.asdhetdir:
            asdhetfiles =  glob.glob(args.asdhet + "/*.bed.gz")
            for f in asdhetfiles:
                asdhet.append(tabix.open(f))
        else:
            asdhet = [tabix.open(args.asdhet)]
    else:
        asdhetfiles = glob.glob(args.asdhet)
        for f in asdhetfiles:
            asdhet.append(tabix.open(f))
    if len(asdhet) == 0:
        ERROR("No asdhet files found")
    else:
        MSG("Using %s asdhet files"%len(asdhet))

    # Get reader for covar files(s)
    if str(args.covar_matrix) != "identity":
        covar = []
        if "*" not in args.covar_matrix:
            covar = [tabix.open(args.covar_matrix)]
        else:
            covarfiles = glob.glob(args.covar_matrix)
            for f in covarfiles:
                covar.append(tabix.open(f))
    else: covar = "identity"
    

    if not args.mutinfo:
        if args.output_eff: effs = "_eff"
        else: effs = ""
        header = ["chrom","start","end","numsamples",\
                      "prior_mu","prior_beta",\
                      "true_mu","true_beta","strsd",\
                      "est_mu%s"%effs,"est_mu%s_stderr"%effs,\
                      "est_beta%s"%effs,"est_beta%s_stderr"%effs]
        if args.output_kl:
            header.extend(["est_k","est_k_stderr","est_L","est_L_stderr"])
        if args.output_ci:
            header.extend(["est_mu_low", "est_mu_high", "est_beta_low", "est_beta_high"])
    else:
        header = ["chrom","start","end","numsamples","num_recurrent_alleles","num_alt_alleles","recurrent_alleles","single_alleles","other_alleles","major_allele"]
    if args.out is not None: f = open(args.out, "w")
    else: f = sys.stdout
    f.write("\t".join(header)+"\n")
    if args.output_bs is not None:
        g = open(args.output_bs, "w")
        bsheader = ["chrom","start","end","mu","beta"]
        g.write("\t".join(bsheader)+"\n")
    else: g = None
    if args.aggregate:
        locdata = PrepareData(-1, loci, truth, args.prior_mu, args.prior_beta, priors, \
                                   args.strsd, strsds, asdhet, args.unit, args.min_samples, args.singlepair, args.maxT, covar, samples_to_use)
        if locdata is not None:
            locdata["estimator_type"] = args.estimator
            locdata["numBS"] = args.numBS
            locdata["output_kl"] = args.output_kl
            locdata["output_ci"] = args.output_ci
            locdata["output_eff"] = args.output_eff
            locdata["distribution"] = args.distribution
            locdata["norm"] = args.norm
            # run and write output
            output_line, bs_lines = RunLocus(locdata)
            f.write(output_line)
            f.flush()
            if args.output_bs is not None:
                for line in bs_lines:
                    g.write(line)
        f.close()
        if args.output_bs: g.close()
        sys.exit(0)
    # Get queues ready
    tasks = multiprocessing.JoinableQueue()
    results = multiprocessing.JoinableQueue()
    # Get consumers ready
    consumers = [Consumer(tasks, results) for i in xrange(args.numproc)]
    for w in consumers: w.start()
    writer = Writer(results, f, g)
    writer.start()
    # Set up job for each locus
    sys.stderr.write("Preparing jobs\n")
    if args.numjobs != -1:
        loci = loci[0:min(len(loci)-1, args.numjobs)]
    numjobs = 0
    for i in range(len(loci)):
        try:
            #if i%1000 == 0: 
            sys.stderr.write("Preparing job %s out of %s\n"%(i, len(loci)))
            locdata = PrepareData(i, loci, truth, args.prior_mu, args.prior_beta, priors, \
                                      args.strsd, strsds, asdhet, args.unit, args.min_samples, args.singlepair, args.maxT, covar, samples_to_use)
            if locdata is not None:
                locdata["estimator_type"] = args.estimator
                locdata["numBS"] = args.numBS
                locdata["output_kl"] = args.output_kl
                locdata["output_ci"] = args.output_ci
                locdata["output_eff"] = args.output_eff
                locdata["distribution"] = args.distribution
                locdata["norm"] = args.norm
                locdata["mutinfo"] = args.mutinfo
                numjobs = numjobs + 1
                tasks.put(locdata)
        except KeyboardInterrupt:
            break
    if debug: sys.stderr.write("Done putting tasks in task queue\n")
    # Add a poison pill for each consumer
    for i in xrange(len(consumers)):
        tasks.put(None)
    # Wait for tasks to finish
    tasks.join()
    results.put(None)
    results.join()
    f.close()
