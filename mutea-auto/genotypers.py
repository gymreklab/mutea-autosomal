import math
import numpy

from scipy.misc  import logsumexp
from scipy.stats import geom
import sys
import geom_stutter_em

def ERROR(msg):
    sys.stderr.write(msg.strip() + "\n")
    sys.exit(1)
    
def MSG(msg):
    sys.stderr.write(msg.strip() + "\n")

# Classes that genotype an STR given a set of reads

'''
Correctly genotypes all samples and assigns the correct allele a posterior probability of 1.0
'''
class ExactGenotyper:
    def __init__(self):
        return

    def train(self, sample_read_counts, min_allele, max_allele):
        return True

    def get_genotype_posteriors(self, read_count_dict, true_genotype, stutter_model=None):
        return {true_genotype:1.0}

    def __str__(self):
        return "ExactGenotyper"

'''
Only returns posteriors if all of the reads support the same allele, in which case that
allele has a posterior of 1.0
'''
class AllGenotyper:
    def __init__(self):
        return

    def train(self, sample_read_counts, min_allele, max_allele):
        return True

    def get_genotype_posteriors(self, read_count_dict, true_genotype=None, stutter_model=None):
        if len(read_count_dict) == 0:
            exit("Cannot generate posteriors if no read lengths are provided")
        elif len(read_count_dict) == 1:
            return {read_count_dict.keys()[0]:1.0}
        else:
            return {}

    def __str__(self):
        return "AllGenotyper"

''' 
Computes the posteriors based on the fraction of reads supporting each allele
'''
class FractionGenotyper:
    def __init__(self):
        return

    def train(self, sample_read_counts, min_allele, max_allele):
        return True

    def get_genotype_posteriors(self, read_count_dict, true_genotype=None, stutter_model=None):
        if len(read_count_dict) == 0:
            exit("Cannot generate posteriors if no read lengths are provided")
        total_reads = sum(read_count_dict.values())
        posteriors = {}
        for allele,count in read_count_dict.items():
            posteriors[allele] = 1.0*count/total_reads
        return posteriors

    def __str__(self):
        return "FractionGenotyper"

'''
Computes the posteriors based on the exact stutter model and a uniform prior
'''
class ExactStutterGenotyper:
    def __init__(self):
        return

    def train(self, sample_read_counts, min_allele, max_allele):
        self.min_n = min_allele
        self.max_n = max_allele
        return True

    def get_genotype_posteriors(self, read_count_dict, true_genotype, stutter_model):
        log_prior  = math.log(1.0/(self.max_n-self.min_n+1))
        log_priors = dict(zip(range(self.min_n,self.max_n+1), (self.max_n-self.min_n+1)*[log_prior]))
        return stutter_model.get_genotype_posteriors(log_priors, read_count_dict)
        
    def __str__(self):
        return "ExactStutterGenotyper"

'''
Computes the posteriors based on an EM-estimated stutter model and a uniform prior
'''
class EstStutterGenotyper:
    def __init__(self, _diploid=False):
        self.diploid = _diploid
        return

    def construct_matrix(self, down, up, p_geom, min_allele, max_allele):
        self.log_down     = numpy.log(down)
        self.log_eq       = numpy.log(1.0-down-up)
        self.log_up       = numpy.log(up)
        self.p_geom       = p_geom
        self.min_allele   = min_allele
        self.max_allele   = max_allele
        self.nalleles     = self.max_allele - self.min_allele + 1
        self.stutter_dist = geom(self.p_geom)

        # Construct matrix where each row contains the stutter transition probabilites for a particular allele
        for j in xrange(self.nalleles):
            allele_probs = numpy.hstack(([self.log_down + self.stutter_dist.logpmf(j-x) for x in range(0, j)],
                                         [self.log_eq],
                                         [self.log_up   + self.stutter_dist.logpmf(x-j) for x in range(j+1, self.nalleles)]))
            if j == 0:
                step_probs = allele_probs
            else:
                step_probs = numpy.vstack((step_probs, allele_probs))
        if self.nalleles == 1:
            step_probs = numpy.expand_dims(step_probs,axis=0)
        self.step_probs = step_probs

    def create(self, down, up, p_geom, min_allele, max_allele):
        self.construct_matrix(down, up, p_geom, min_allele, max_allele)
    
    def train(self, sample_read_counts, min_allele, max_allele, debug=False):
        MSG("Estimating stutter model...")
        status, eff_coverage, down, up, p_geom, new_LL = geom_stutter_em.run_EM(sample_read_counts, diploid=self.diploid)
        if status == geom_stutter_em.CONVERGED:
            self.status = "CONVERGED"
        elif status == geom_stutter_em.ITERATION_LIMIT:
            self.status = "ITERATION_LIMIT"
        elif status == geom_stutter_em.COVERAGE_LIMIT:
            self.status = "COVERAGE_LIMIT"
        else:
            exit("ERROR: Invalid EM status: %d"%(status))
        self.eff_coverage = eff_coverage
        if status != geom_stutter_em.CONVERGED:
            self.LL        = "N/A" 
            self.est_up    = "N/A"
            self.est_down  = "N/A"
            self.est_pgeom = "N/A"
            return False
        else:
            self.LL        = new_LL
            self.est_up    = up
            self.est_down  = down
            self.est_pgeom = p_geom
        
        MSG("\tEstimated parameters: P_GEOM=%f, P_DOWN=%f, P_UP=%f"%(p_geom, down, up))
        self.construct_matrix(down, up, p_geom, min_allele, max_allele)
        return True

    def get_genotype_posteriors(self, read_count_dict, true_genotype=None, stutter_model=None, errprob=0):
        read_counts_array = numpy.zeros(self.nalleles)
        for length,count in read_count_dict.items():
            if length < self.min_allele or length > self.max_allele:
                exit("ERROR: Allele length is outside of stutter model's range")
            read_counts_array[length-self.min_allele] = count
        if self.diploid:
            num_gts = self.nalleles**2
            allele_sizes = range(self.min_allele, self.min_allele+self.nalleles)
            LLs = numpy.zeros(num_gts)
            gtind = 0
            keys = []
            for a1 in xrange(self.nalleles):
                for a2 in xrange(self.nalleles):
                    if a1 == a2: LLs[gtind] = numpy.log(1)
                    else: LLs[gtind] = numpy.log(0.5)
                    step_probs1 = self.step_probs[a1,:]
                    step_probs2 = self.step_probs[a2,:]
                    step_probs = numpy.logaddexp(step_probs1, step_probs2)+numpy.log(0.5)+numpy.log(1-errprob)
                    step_probs = numpy.logaddexp(numpy.log(errprob), step_probs)
                    LLs[gtind] += numpy.sum(read_counts_array*step_probs)
                    keys.append((allele_sizes[a1], allele_sizes[a2]))
                    gtind += 1
            posteriors = numpy.exp(LLs-logsumexp(LLs))
            # Combine (a, b) and (b, a)
            posts = dict(zip(keys, posteriors))
            posts_new = {}
            for key in posts:
                a, b = key
                if a < b: posts_new[(a, b)] = posts_new.get((a, b), 0) + posts[key]
                else: posts_new[(b, a)] = posts_new.get((b, a), 0) + posts[key]
            return posts_new
        else:
            LLs        = numpy.sum(self.step_probs*read_counts_array, axis=1)
            posteriors = numpy.exp(LLs-logsumexp(LLs))
            return dict(zip(range(self.min_allele, self.max_allele+1), posteriors))

    def __str__(self):
        return "EstStutterGenotyper"

'''
Does not compute posteriors. Instead, simple class that returns a dictionary
of read counts associated with each sample. Used when simulating PCR stutter
'''
class ReadCounter:
    def __init__(self):
        return

    def train(self, sample_read_counts, min_allele, max_allele):
        return True

    def get_genotype_posteriors(self, read_count_dict, true_genotype=None, stutter_model=None):
        return read_count_dict

    def __str__(self):
        return "ReadCounter"

def get_genotyper(name):
    if name == "EXACT":
        return ExactGenotyper()
    elif name == "ALL":
        return AllGenotyper()
    elif name == "FRACTION":
        return FractionGenotyper()
    elif name == "EXACT_STUTTER":
        return ExactStutterGenotyper()
    elif name == "EST_STUTTER":
        return EstStutterGenotyper()
    else:
        exit("ERROR: Invalid genotyper requested: %s"%(genotyper_name))
