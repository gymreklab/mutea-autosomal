import Bio.Phylo
import math
import numpy as np
import random
import sys

def weighted_choice(choices):
    """
    http://stackoverflow.com/questions/3679694/a-weighted-version-of-random-choiec
    """
    total = sum(w for c, w in choices)
    r = random.uniform(0, total)
    upto = 0
    for c, w in choices:
        if upto + w > r:
            return c
        upto += w
    assert False, "Shouldn't get here"

class STRMutationGenerator(object):
    def __init__(self, tree, str_mut_rate=0.001, str_step_param=0.9, method="SMML",
                 str_allele_sd=1, length_constraint=0, lencoeff=0, root=0):
        self.tree = tree
        self.str_mut_rate = str_mut_rate
        self.str_step_param = str_step_param
        self.method = method
        self.str_allele_sd = str_allele_sd
        self.length_constraint = length_constraint
        self.lencoeff = lencoeff
        self.root = root
        self.nodeToSTR = {}
        if tree is not None: self.parent = self.GetParents(tree)
        self.effective_mut_rate = None
        self.fractionBoundaries = None
        self.allMutations = None

    def GetParents(self, tree):
        """
        Dictionary of clade->parent
        """
        parents = {}
        clades = tree.find_elements(target=Bio.Phylo.Newick.Clade, order="preorder")
        for clade in clades:
            for subc in clade:
                parents[subc] = clade
        return parents

    def GetSTRMutation(self, starting_allele):
        """
        Get STR mutation event
        Return new_allele, boundary(T/F)
        """
        boundary = False
        if self.method == "SMML":
            mutlen = np.random.geometric(p=self.str_step_param)
            if random.random() < 0.5: mutlen = mutlen*-1
            scale = 1 + self.length_constraint*1.0/self.str_allele_sd
            new_allele = (starting_allele+mutlen)*scale
        elif self.method == "OU":
            # draw mut len from distribution with mean -lenc*starting_allele, variance in change str_allele_sd**2
            mutlen_exp = starting_allele*-1*self.length_constraint
            mutlen_sd = self.str_allele_sd
            mutlen = np.random.normal(loc=mutlen_exp, scale=mutlen_sd)
            new_allele = starting_allele + mutlen
        elif self.method == "discreteOU": # D.W.C.Miao equation 4
            # Get probability of downprob and upprob
            u = (self.str_allele_sd**2  - self.length_constraint*starting_allele)*0.5
            d = (self.str_allele_sd**2  + self.length_constraint*starting_allele)*0.5
            # Boundary states
            if u > 1 or d < 0:
                u = 1
                d = 0
            if u < 0 or d > 1:
                u = 0
                d = 1
            if u == 0 or u == 1:
                boundary = True
            # Draw step size (-1, 0, or 1)
            mutlen = weighted_choice([(-1, d), (0, 1-(d+u)), (1, u)])
            new_allele = starting_allele + mutlen
        elif self.method == "discreteMultiOU": # my extension of D.W.C. Miao for multiple step sizes
            # get up and down probs
            denom = 2*(2-self.str_step_param)
            u = (self.str_allele_sd**2*self.str_step_param**2 - (2-self.str_step_param)*self.length_constraint*self.str_step_param*starting_allele)*1.0/denom
            d = (self.str_allele_sd**2*self.str_step_param**2 + (2-self.str_step_param)*self.length_constraint*self.str_step_param*starting_allele)*1.0/denom
            # Boundary states
            if u > 1 or d < 0:
                u = 1
                d = 0
            if u < 0 or d > 1:
                u = 0
                d = 1
            if u == 0 or u == 1:
                boundary = True
            # Get step magnitude from geometric distribution
            step = np.random.geometric(p=self.str_step_param)
            mutlen = weighted_choice([(-1*step, d), (0, 1-(d+u)), (step, u)])
            new_allele = starting_allele + mutlen
        else: new_allele = starting_allele # no method specified, so no mutation
        return new_allele, boundary

    def GetSTRForT(self, tmrca, start_allele = 0):
        """
        Evolve STR on two branches with given tmrca. Get ASD
        """
        current_allele = start_allele
        tau = tmrca
        T = np.random.exponential(scale=1/self.str_mut_rate, size=1)[0]
        while tau > 0:
            if T < tau:
                before = current_allele
                current_allele, boundary = self.GetSTRMutation(current_allele)
                tau = tau - T
            else:
                break
        return current_allele

    def ProcessEdge(self, clade):
        """
        Evolve STRs along an edge
        Return:
        nummutations, allmutations, numboundary
        """
        numMutations = 0
        allmutations = []
        numboundary = 0
        current_allele = self.nodeToSTR.get(self.parent.get(clade, None), self.root) # root
        current_mu = 10**(math.log10(self.str_mut_rate)+self.lencoeff*(current_allele))
        tau = clade.branch_length
        while tau > 0:
            # Draw time from exponential distribution
            T = np.random.exponential(scale=1/current_mu, size=1)[0]
            if T < tau:
                before = current_allele
                current_allele, b = self.GetSTRMutation(current_allele)
                current_mu = 10**(math.log10(self.str_mut_rate)+self.lencoeff*(current_allele))
                if current_allele != before:
                    numMutations += 1
                    allmutations.append("%s->%s"%(before, current_allele))
                    numboundary += int(b)
                tau = tau - T
                if self.lencoeff != 0: # Redraw T to reflect updated mutation rate
                    T = np.random.exponential(scale=1/current_mu)
            else: break
        self.nodeToSTR[clade] = current_allele
        return numMutations, allmutations, numboundary

    def EvolveSTR(self):
        """
        Evolve STR on the tree. return NodeToSTR
        """
        totalTime = 0 # keep track of number of total mutations
        totalMutations = 0
        numboundary = 0
        allmutations = []
        self.nodeToSTR = {}
        clades = self.tree.find_elements(target=Bio.Phylo.Newick.Clade, order="preorder")
        for clade in clades:
            if clade.branch_length is not None:
                tm, m, boundary = self.ProcessEdge(clade)
                totalMutations += tm
                allmutations += m
                numboundary += boundary
                totalTime += clade.branch_length
        self.effective_mut_rate = totalMutations*1.0/totalTime
        self.totalMutations = totalMutations
        self.totalTime = totalTime
        self.allMutations = allmutations
        if totalMutations == 0:
            self.fractionBoundaries = 0
        else: self.fractionBoundaries = numboundary*1.0/totalMutations
        return self.nodeToSTR
