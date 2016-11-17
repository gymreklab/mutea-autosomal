import argparse
from Bio.PopGen.SimCoal.Template import generate_simcoal_from_template
from Bio.PopGen.SimCoal.Controller import FastSimCoalController
import Bio.Phylo
from cStringIO import StringIO
import dendropy
import os
import random
from STRMutationGenerator import *
import numpy as np
import sys
import tempfile
import vcf

def LOG(msg):
    sys.stderr.write(msg.strip() + "\n")

##### General parameters #####
OUTPREFIX = ""
DEBUG = False
DIPLOID = False
INDEPENDENT = False

##### STR paramters #####
STR_MUT_RATE = 0.001
STR_STEP_PARAM = 0
STR_ALLELE_SD = 1
LENGTH_CONSTRAINT = 0
LENCOEFF = 0
ROOTNODE = 0

##### Coalescent parameters #####
NEFF = 20000 # this is really 2*Neff
SAMPLES = 100
NUMSIM = 1
METHOD = "SMML" # SMML, OU

def GetTMRCA(tree, leaf1, leaf2):
    """
    Get TMRCA
    """
    return tree.distance(leaf1, leaf2)/2

def GetNewickStrings(treefile):
    """
    Get list of strings for coalescent trees in newick format
    """
    newickstrings = []
    treelines = [line for line in open(treefile, "r").readlines() if "NumGen_tree" in line]
    for t in treelines:
        newickstrings.append(t.split("=")[1].strip())
    return newickstrings

def GetSeqData(arpfile):
    """
    Return dictionary of sample -> SNPs, STR
    """
    seqinfo = {}
    f = open(arpfile, "r")
    line = f.readline()
    while "SampleData" not in line and line != "":
        line = f.readline()
    line = f.readline()
    if line.strip() == "": return {}
    while line.strip() != "":
        sample, dummy, gt = line.strip().split("\t")
        snps, strgt = map(str.strip, gt.split())
        sample = sample.split("_")[1] + ".1"
        seqinfo[sample] = {"SNP": snps, "STR": int(strgt)}
        line = f.readline()
    f.close()
    return seqinfo

def GetHet(snp1_fsc, snp2_fsc):
    """
    Return % of bases differing
    """
    return len([i for i in range(len(snp1_fsc)) if snp1_fsc[i] != snp2_fsc[i]])*1.0/len(snp1_fsc)

def GetBranches(tree, leaf1, leaf2):
    """
    Return list of branches on path between leaf1 and leaf2
    """
    path1 = tree.get_path(leaf1)
    path2 = tree.get_path(leaf2)
    unique = []
    for c in path1+path2:
        if (path1+path2).count(c) == 1:
            unique.append(str(hash(c))+":"+str(c.branch_length))
    return unique

def GenerateErrors(allele, model):
    """
    model is upprob, downprob, geomp
    """
    upprob, downprob, geomp = model
    # Draw step size from geometric distribution
    step_size = np.random.geometric(p=geomp)
    r = random.random()
    # < u: up, < u+d: down, else no mutation
    if r < upprob:
        return allele + step_size
    elif r < upprob + downprob:
        return allele - step_size
    else: return allele

def ApplyStutter(nodeToSTR, model):
    nodeToSTR_stutter = {}
    for key in nodeToSTR:
        nodeToSTR_stutter[key] = GenerateErrors(nodeToSTR[key], model)
    return nodeToSTR_stutter

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", help="Prefix for output files", type=str, required=True)
    parser.add_argument("--Neff", help="Effective population size", type=int)
    parser.add_argument("--root", help="Root node value", default=0, type=int)
    parser.add_argument("--samples", help="Number of samples", type=int)
    parser.add_argument("--strmutrate", help="STR mutation rate", type=float)
    parser.add_argument("--strstepparam", help="STR step parameter", type=float)
    parser.add_argument("--strsd", help="Standard deviation of step size distribution", type=float)
    parser.add_argument("--lenconstraint", help="Length constraint parameter", type=float)
    parser.add_argument("--lencoeff", help="Mutation rate length coefficient", type=float, default=0.0)
    parser.add_argument("--numsim", help="Number of trees to build", type=int)
    parser.add_argument("--diploid", help="Don't take all pairs. Generate diploid individuals", action="store_true")
    parser.add_argument("--method", help="How to generate mutations. Default: SMML", type=str)
    parser.add_argument("--independent", help="Treat each pair as independent (just use tree to get TMRCA)", action="store_true")
    parser.add_argument("--notree", help="Don't get tmrca from tree (only for --independent)", action="store_true")
    parser.add_argument("--use-same-tree", help="Use same tree for all simulations. Output one line per pair", action="store_true")
    parser.add_argument("--stutter", help="Stutter model (<upprob>,<downprob>,<geomp>)", type=str, default="0,0,1")
    parser.add_argument("--debug", help="Debug", action="store_true")

    args = parser.parse_args()
    if args.out: OUTPREFIX = args.out
    if args.Neff: NEFF = args.Neff
    if args.samples: SAMPLES = args.samples
    if args.strmutrate: STR_MUT_RATE = args.strmutrate
    if args.strstepparam: STR_STEP_PARAM = args.strstepparam
    if args.lencoeff: LENCOEFF = args.lencoeff
    if args.root is not None: ROOTNODE = args.root
    if args.strsd: STR_ALLELE_SD = args.strsd
    if args.lenconstraint: LENGTH_CONSTRAINT = args.lenconstraint
    if args.numsim is not None: NUMSIM = args.numsim
    if args.method is not None: METHOD = args.method
    INDEPENDENT = args.independent
    NOTREE = args.notree
    DIPLOID = args.diploid
    DEBUG = args.debug
    STUTTERMODEL = map(float, args.stutter.split(","))

    # Write log file
    f = open(OUTPREFIX + ".log", "w")
    f.write("NEFF: %s\n"%NEFF)
    f.write("SAMPLES: %s\n"%SAMPLES)
    f.write("STRMUTRATE: %s\n"%STR_MUT_RATE)
    f.write("STRSTEPPARAM: %s\n"%STR_STEP_PARAM)
    f.write("STRSD: %s\n"%STR_ALLELE_SD)
    f.write("LENCONSTRAINT: %s\n"%LENGTH_CONSTRAINT)
    f.write("LENCOEFF: %s\n"%LENCOEFF)
    f.write("ROOT: %s\n"%ROOTNODE)
    f.write("NUMSIM: %s\n"%NUMSIM)
    f.write("DIPLOID: %s\n"%DIPLOID)
    f.write("METHOD: %s\n"%METHOD)
    f.write("IND: %s\n"%INDEPENDENT)
    f.write("ONETREE: %s\n"%args.use_same_tree)
    f.close()

    # Generate coalescent tree
    tmpdir = tempfile.mkdtemp(prefix='ousim.')
    LOG("Generate coalescent trees")
    generate_simcoal_from_template("simple", [(1,[("SNP", [100, 0, 0.01]),
                                                  ("MICROSAT", [1, 0, STR_MUT_RATE, STR_STEP_PARAM, 0])])],
                                   [("sample_size", [SAMPLES]),("pop_size", [NEFF])], out_dir=tmpdir)
    filename = os.path.join(tmpdir,"simple_%s_%s.par"%(NEFF, SAMPLES))
    fsc = FastSimCoalController(bin_name="fastsimcoal")
    os.chdir(tmpdir)
    if args.use_same_tree:
        fsc.run_fastsimcoal(filename, 1, opts={"tree": True})
    else: fsc.run_fastsimcoal(filename, NUMSIM, opts={"tree": True})

    # Read in coalescent trees
    treefile = os.path.join(tmpdir,"simple_%s_%s/simple_%s_%s_1_true_trees.trees"%(NEFF, SAMPLES, NEFF, SAMPLES))
    newickstrings = GetNewickStrings(treefile)
    phylotrees = [Bio.Phylo.read(StringIO(ns), "newick") for ns in newickstrings]

    # Read in SNP/STR data
    seqdata = []
    if args.use_same_tree: numtrees = 1
    else: numtrees = NUMSIM
    for i in range(numtrees):
        arpfile = os.path.join(tmpdir, "simple_%s_%s/simple_%s_%s_1_%s.arp"%(NEFF, SAMPLES, NEFF, SAMPLES, i+1))
        seqdata.append(GetSeqData(arpfile))
    
    if args.use_same_tree:
        f = open(OUTPREFIX + ".asd_vs_tmrca.tab", "w")
        header = ["sample1", "sample2", "tmrca", "branches", "snps"] + map(lambda x: "asd_%s"%x, range(NUMSIM))
        f.write("\t".join(header)+"\n")
        LOG("Evolve STRs")
        pt = phylotrees[0]
        treeseqdata = seqdata[0]
        strgen = STRMutationGenerator(pt, str_mut_rate=STR_MUT_RATE, str_step_param=STR_STEP_PARAM, \
                                          method=METHOD, str_allele_sd=STR_ALLELE_SD, \
                                          length_constraint=LENGTH_CONSTRAINT, lencoeff=LENCOEFF, root=ROOTNODE)
        nodeToSTR = []
        for i in range(NUMSIM):
            nodeToSTR.append(strgen.EvolveSTR())
        # Process each pair
        leaves = list(pt.find_elements(terminal=True))
#        leaf_ind = range(len(leaves))
#        random.shuffle(leaf_ind)
#        leaf_indices1 = [leaf_ind[i] for i in range(len(leaf_ind)) if leaf_ind[i]%2 == 0]
#        leaf_indices2 = [leaf_ind[i] for i in range(len(leaf_ind)) if leaf_ind[i]%2 == 1]
#        leaf_pairs = zip(leaf_indices1, leaf_indices2)
        leaf_pairs = []
        for i in range(len(leaves)):
            for j in range(len(leaves)):
                if i != j:
                    leaf_pairs.append((i,j))
        # Randomly choose leaves
        leaf_pairs = [random.choice(leaf_pairs) for i in range(SAMPLES)]
        for i,j in leaf_pairs:
            leaves[i].confidence = leaves[i] # added 04/27/15
            leaves[j].confidence = leaves[j]
            tmrca = GetTMRCA(pt, leaves[i], leaves[j])
            branches = GetBranches(pt, leaves[i], leaves[j])
            asds = []
            for ns in range(NUMSIM):
                str1 = nodeToSTR[ns][leaves[i]]
                str2 = nodeToSTR[ns][leaves[j]]
                asds.append((str2-str1)**2)
            snps1 = treeseqdata[str(leaves[i].confidence)]["SNP"]
            snps2 = treeseqdata[str(leaves[j].confidence)]["SNP"]
            # Need indi, indj, tij, bij, asdij(1...numsim), snpsij
            snps = "".join(map(lambda x: str(sum(x)), zip(map(int, list(snps1)), map(int, list(snps2)))))
            f.write("\t".join(map(str, [i, j, tmrca, ",".join(map(str, branches)), snps, "\t".join(map(str, asds))]))+"\n")
        f.close()
    else:
        # File to output summary info
        sumf = open(OUTPREFIX + ".summary.tab", "w")
        sumf.write("\t".join(["treenum","mutrate","sd","lenc","num_mutations","total_branch_length","fraction_boundaries","eff_mut_rate"])+"\n")

        # Get STR mutations
        LOG("Evolve STRs")
        nodeToSTR = {}
        strgens = [STRMutationGenerator(pt, str_mut_rate=STR_MUT_RATE, str_step_param=STR_STEP_PARAM, lencoeff=LENCOEFF, root=ROOTNODE, method=METHOD, \
                                            str_allele_sd=STR_ALLELE_SD, length_constraint=LENGTH_CONSTRAINT) for pt in phylotrees]
        eff_mut_rates = []
        num_mutations = []
        mutations = []
        total_branch_lengths = []
        fraction_boundaries = []
        for sg in strgens:
            nodeToSTR[sg.tree] = ApplyStutter(sg.EvolveSTR(), STUTTERMODEL)
            emr = sg.effective_mut_rate
            eff_mut_rates.append(emr)
            num_mutations.append(sg.totalMutations)
            total_branch_lengths.append(sg.totalTime)
            fraction_boundaries.append(sg.fractionBoundaries)
            mutations.append(sg.allMutations)
        LOG("Mean effective mutation rate: %s:"%np.mean(eff_mut_rates))

        # Get for each pair, the TMRCA and ASD
        LOG("Get TMRCA/ASD")
        f = open(OUTPREFIX + ".asd_vs_tmrca.tab", "w")
        f.write("\t".join(["treenum","sample1","sample2","asd","tmrca","fastsimcoal.asd","het","snps1","snps2","str1","str2"])+"\n")
        treenum = 0
        for p in range(len(phylotrees)):
            pt = phylotrees[p]
            treeseqdata = seqdata[p]
            leaves =  list(pt.find_elements(terminal=True))
            if DIPLOID:
                leaf_ind = range(len(leaves))
                random.shuffle(leaf_ind)
                leaf_indices1 = [leaf_ind[i] for i in range(len(leaf_ind)) if leaf_ind[i]%2 == 0]
                leaf_indices2 = [leaf_ind[i] for i in range(len(leaf_ind)) if leaf_ind[i]%2 == 1]
                leaf_pairs = zip(leaf_indices1, leaf_indices2)
            else:
                leaf_pairs = []
                for i in range(len(leaves)):
                    for j in range(len(leaves)): leaf_pairs.append((i,j))
            for i,j in leaf_pairs:
                leaves[i].confidence = leaves[i] # added 04/27/15
                leaves[j].confidence = leaves[j]
                tmrca = GetTMRCA(pt, leaves[i], leaves[j])
                if INDEPENDENT:
                    if NOTREE: # ignore tmrca from tree. obviously faster to not deal with trees, but easier to code here
                        tmrca = random.randint(0, 100000)
                    str1 = strgens[p].GetSTRForT(tmrca)
                    str2 = strgens[p].GetSTRForT(tmrca)
                else:
                    str1 = nodeToSTR[pt][leaves[i]]
                    str2 = nodeToSTR[pt][leaves[j]]
                asd = (str2-str1)**2
                str1_fsc = treeseqdata[str(leaves[i].confidence)]["STR"]
                str2_fsc = treeseqdata[str(leaves[j].confidence)]["STR"]
                fastsimcoal_asd = (str1_fsc-str2_fsc)**2
                snp1_fsc = treeseqdata[str(leaves[i].confidence)]["SNP"]
                snp2_fsc = treeseqdata[str(leaves[j].confidence)]["SNP"]
                het = GetHet(snp1_fsc, snp2_fsc)
                f.write("\t".join(map(str, [treenum, leaves[i].confidence, leaves[j].confidence, asd, tmrca, fastsimcoal_asd, het, snp1_fsc, snp2_fsc, str1, str2])) + "\n")
            sumf.write("\t".join(map(str, [treenum, STR_MUT_RATE, STR_ALLELE_SD, LENGTH_CONSTRAINT, \
                                               num_mutations[p], total_branch_lengths[p], fraction_boundaries[p], \
                                               eff_mut_rates[p]]))+"\n")
            treenum = treenum + 1
        f.close()
        sumf.close()

"""
    # Write out alleles at each locus
    f = open(OUTPREFIX + ".alleles.tab", "w")
    f.write("\t".join(["locus","STR"])+"\n")
    for p in range(len(phylotrees)):
        pt = phylotrees[p]
        leaves = list(pt.find_elements(terminal=True))
        for i in range(len(leaves)):
            strgt = nodeToSTR[pt][leaves[i]]
            f.write("\t".join(map(str, [p, strgt]))+"\n")
    f.close()
"""
