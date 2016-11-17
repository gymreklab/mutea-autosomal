import argparse
import numpy as np
import random
import sys
import scipy.stats

def LOG(msg):
    sys.stderr.write(msg.strip() + "\n")

def GetStutter(allele, stutter):
    upprob, downprob, geomp = map(float, stutter.split(","))
    # Draw step size from geomtric distribution
    step_size = np.random.geometric(p=geomp)
    r = random.random()
    # < u: up, < u+d: down, else no mutation
    if r < upprob:
        return allele + step_size
    elif r < upprob + downprob:
        return allele - step_size
    else: return allele

def GetReads(gt1, gt2, stutter, coverage):
    reads = []
    for i in range(coverage):
        if random.random() > 0.5:
            reads.append(GetStutter(gt1, stutter))
        else: reads.append(GetStutter(gt2, stutter))
    return reads

def GetAllReadsString(allreads):
    d = {}
    for a in set(allreads):
        d[a] = allreads.count(a)
    return ";".join(["%s|%s"%(a, d[a]) for a in d])

p_up = 0.01
p_down = 0.01
p_geom = 0.99
stutter_probs = [scipy.stats.geom(p_geom).pmf(i) for i in range(1,100)]
def GetTransitionProb(allele, read):
    if allele == read: return 1-p_up-p_down
    stepprob = stutter_probs[abs(allele-read)-1]
    if read > allele: return p_up*stepprob
    else: return p_down*stepprob

def GetLogLik(a1, a2, reads):
    loglik = 0
    for read in reads:
        prob1 = GetTransitionProb(a1, read)
        prob2 = GetTransitionProb(a2, read)
        prob = 0.5*prob1+0.5*prob2
        loglik += np.log10(prob)
    return loglik

def ReEstimateGT(reads, alleles):
    best_gt = (None, None)
    max_lik = -1*np.inf
    for i in xrange(len(alleles)):
        for j in xrange(len(alleles)):
            loglik = GetLogLik(alleles[i], alleles[j], reads)
            if loglik > max_lik:
                best_gt = (alleles[i], alleles[j])
                max_lik = loglik
    return best_gt

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--asdhet", help="ASD het tab file from simulation", type=str, required=True)
    parser.add_argument("--stutter", help="Stutter model (<upprob>,<downprob>,<geomp>)", type=str, default="0,0,1")
    parser.add_argument("--reestimate-gt", help="Restimate genotypes using PCR free stutter model", action="store_true")
    parser.add_argument("--coverage", help="Number of reads to simulate", type=int, default=10)
    parser.add_argument("--chrom", help="Use this chromosome", type=str, default="Z")
    parser.add_argument("--locstart", help="Use start coordinate locstart*1000+treenum", type=int, default=0)
    parser.add_argument("--haploid", help="For debugging. Use TMRCA=0 and only take one allele", action="store_true")
    parser.add_argument("--max-samples", help="Only use this many samples", type=int, default=1000)
    args = parser.parse_args()

    # Load asd vs tmrca info. Just need list of (g1, g2, tmrca) for each treenum
    sampledata  = {}
    allelesizes = {}
    with open(args.asdhet, "r") as f:
        for line in f:
            if "treenum" in line: continue # header
            items = line.strip().split()
            treenum = int(items[0])
            gt1 = int(items[9])
            gt2 = int(items[10])
            tmrca = float(items[4])
            if treenum not in sampledata:
                sampledata[treenum] = []
                allelesizes[treenum] = set()
            if args.haploid:
                sampledata[treenum].append((gt1, gt1, -1))
                sampledata[treenum].append((gt2, gt2, -1))
            else: sampledata[treenum].append((gt1, gt2, tmrca))
            allelesizes[treenum].add(gt1)
            allelesizes[treenum].add(gt2)

    # Write one VCF entry per locus
    for tree in sampledata:
        if len(sampledata[tree])>args.max_samples: sampledata[tree] = sampledata[tree][0:args.max_samples]
        # Get locus info
        min_allele = min(allelesizes[tree])
        reflen = -1*min_allele+1
        # Get VCF fields
        chrom = args.chrom
        ident = "."
        REF="A"*reflen
        ALT=[]
        start = args.locstart*1000+tree
        end = start+len(REF)
        for a in allelesizes[tree]:
            if a != 0:
                ALT.append("A"*(reflen+a))
        ALT.sort()
        alt_lengths = [len(item)-len(REF) for item in ALT]
        INFO = "MOTIF=A;END=%s"%end
        FORMAT = "GT:ALLREADS:TGB:GB:TMRCA"
        # Get samples
        SAMPLES = []
        for sample in sampledata[tree]:
            if sample[0] == 0: gt1 = 0
            else: gt1 = alt_lengths.index(sample[0])+1
            if sample[1] == 0: gt2 = 0
            else: gt2 = alt_lengths.index(sample[1])+1
            allreads = GetReads(sample[0], sample[1], args.stutter, args.coverage)
            allreads_string = GetAllReadsString(allreads)
            if args.reestimate_gt:
                err1, err2 = ReEstimateGT(allreads, [0]+alt_lengths)
                if err1 is None or err2 is None:
                    sys.stderr.write("Error getting best GT %s, %s %s,%s\n"%(allreads, [0]+alt_lengths, sample[0], sample[1]))
                    sys.exit(1)
                if err1 == 0: estgt1 = 0
                else:
                    if err1 not in alt_lengths:
                        ALT.append("A"*(reflen+err1))
                        alt_lengths.append(err1)
                    estgt1 = alt_lengths.index(err1)+1
                if err2 == 0: estgt2 = 0
                else:
                    if err2 not in alt_lengths:
                        ALT.append("A"*(reflen+err2))
                        alt_lengths.append(err2)
                    estgt2 = alt_lengths.index(err2)+1
            else:
                err1, err2 = sample[0], sample[1]
                estgt1, estgt2 = gt1, gt2
            sdata = "%s/%s:%s:%s/%s:%s/%s:%s"%(estgt1, estgt2, allreads_string, sample[0], sample[1], err1, err2, sample[2])
            SAMPLES.append(sdata)
        sys.stdout.write("\t".join(map(str, [chrom, start, ident, REF, ",".join(ALT), 0, "PASS", INFO, FORMAT]+SAMPLES))+"\n")
            
    
