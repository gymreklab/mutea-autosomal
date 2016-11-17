"""
Store info for joint estimation across loci
"""

import joblib
import numpy as np
import random
import scipy.integrate
import sys
import time
from numpy.linalg import inv

def ERROR(msg):
    sys.stderr.write(msg.strip() + "\n")
    sys.exit(1)
    
def MSG(msg):
    sys.stderr.write(msg.strip() + "\n")

def GetLocusNegLL(locus, drawmu, params, feature_param_index, numfeatures, mu, sd, lencoeff, \
                      mu_bounds, beta_bounds, pgeom_bounds, lencoeff_bounds, mean_length=0, ires=100, debug=False):
    # Adjust mu for features
    adj_mu = mu + sum([locus.features[j]*params[j+feature_param_index] for j in range(numfeatures)])
    if drawmu:
        adj_sd = sd + sum([locus.features[j]*params[j+feature_param_index+numfeatures] for j in range(numfeatures)])
        # Draw samples for mu
        musamples = np.random.normal(loc=adj_mu, scale=adj_sd, size=ires)
    else:
        musamples = [adj_mu]
    # Approximate integration of P(D|mu)P(u)du
    nll_values = [GetLocusNegLogLikelihood(locus, mval, lencoeff, mean_length, mu_bounds, \
                                               beta_bounds, pgeom_bounds, lencoeff_bounds, debug=debug) for mval in musamples]
    if np.inf in nll_values:
        raise LLException("Negative infinite neg ll")
    nll = GetProbAvg(nll_values)
    return nll

def GetLocusNegLogLikelihood(locus, mu, lencoeff, mean_length, mu_bounds, \
                                 beta_bounds, pgeom_bounds, lencoeff_bounds, debug=False):
    beta = np.random.uniform(*beta_bounds)
    pgeom = np.random.uniform(*pgeom_bounds)
    if locus.prior_beta is not None: beta = locus.prior_beta
    if locus.prior_pgeom is not None: pgeom = locus.prior_pgeom
    val = locus.NegativeLogLikelihood(mu, beta, pgeom, lencoeff, range(len(locus.data)), \
                                          mu_bounds=mu_bounds, beta_bounds=beta_bounds, pgeom_bounds=pgeom_bounds, \
                                          lencoeff_bounds=lencoeff_bounds, \
                                          mut_model=None, allele_range=None, optimizer=None, mean_length=mean_length, joint=True)
    return val

def GetProbAvg(nll_values):
    logsump = reduce(lambda x, y: np.logaddexp(x, y), map(lambda x: -1*x, nll_values))
    return -1*(logsump - np.log(len(nll_values)))

class LLException(Exception):
    def __init__(self,*args,**kwargs):
        Exception.__init__(self,*args,**kwargs)

class JointLocus:
    def __init__(self, _locilist, _ires=10, _numproc=1):
        self.loci = _locilist
        self.best_res = None
        self.numiter = 2
        self.method = "Nelder-Mead"
        self.max_cycle_per_iter = 250
        self.ires = _ires
        self.stderrs = []
        self.numproc = _numproc
        self.priors = None
        self.mean_length = self.GetMeanLength()

    def GetMeanLength(self):
        lengths = []
        for locus in self.loci:
            l = (locus.end - locus.start + 1)
            lengths.append(l)
        return np.mean(lengths)

    def SetPriors(self, _priors):
        self.priors = _priors

    def callback_function(self, val):
        sys.stderr.write("Current parameters: %s\n"%(str(val)))

    def NegativeLogLikelihood(self, params, numfeatures, drawmu, mu_bounds, sd_bounds, \
                                  beta_bounds, pgeom_bounds, lencoeff_bounds, debug=False):
        # Order of params: [mu0, [lencoeff], mu_coeff1, mu_coeff2, mu_coeff3...] if drawmu=False
        # Order of params: [mu0, sd0, mu_coeff1, coeff2, ... sd_coeff1, sd_coeff2,...] if drawmu=True
        mu = params[0]
        sd = None
        feature_param_index = 1
        if lencoeff_bounds[0] != lencoeff_bounds[1]:
            lencoeff = params[1]
            feature_param_index += 1
        else: lencoeff = 0
        if drawmu:
            sd = params[feature_param_index]
            feature_param_index += 1
        # Check bounds
        if mu < np.log10(mu_bounds[0]) or mu > np.log10(mu_bounds[1]): return np.inf
        if drawmu:
            if sd < sd_bounds[0] or sd > sd_bounds[1]: return np.inf
        # Loop over each locus
        if debug: MSG("Starting likelihood calculation... %s loci"%len(self.loci))
        start_time = time.time()
        try:
            locnlogliks = joblib.Parallel(n_jobs=self.numproc)(joblib.delayed(GetLocusNegLL)(self.loci[i], drawmu, params, \
                                                                                                 feature_param_index, numfeatures, mu, sd, lencoeff, \
                                                                                                 mu_bounds, beta_bounds, pgeom_bounds, lencoeff_bounds, \
                                                                                                 mean_length=self.mean_length, ires=self.ires, debug=debug) \
                                                                   for i in range(len(self.loci)))
            result = sum(locnlogliks)
        except LLException:
            if debug: MSG("Exception in likelihood calculation, returning np.inf")
            result = np.inf
#        except:
#            ERROR("Unexpected exception in NegativeLogLikelihood")
        end_time = time.time()
        if debug: MSG("Done with likelihood calculation... %s Params: %s, value %s"%((end_time-start_time), params, result))
        return result

    def LoadData(self):
        toremove = []
        for locus in self.loci:
            locus.LoadData()
            if len(locus.data) <= locus.minsamples:
                toremove.append(locus)
        for locus in toremove: self.loci.remove(locus)

    def MaximizeLikelihood(self, mu_bounds=None, sd_bounds=None, beta_bounds=None, \
                               pgeom_bounds=None, lencoeff_bounds=None, drawmu=False, debug=False):
        if len(self.loci) == 0: return None
        MSG("Loading data...")
        # Load data for each locus
        self.LoadData()

        MSG("Setting up likelihood... %s loci"%len(self.loci))
        # How many params? mu + numfeatures
        numfeatures = len(self.loci[0].features)

        # Likelihood function
        fn = (lambda x: self.NegativeLogLikelihood(x, numfeatures, drawmu, mu_bounds, sd_bounds, \
                                                       beta_bounds, pgeom_bounds, lencoeff_bounds, debug=debug))

        MSG("Optimize likelihood...")
        # Optimize likelihood
        if debug: callback = self.callback_function
        else: callback = None
        best_res = None
        if self.priors is None:
            for i in xrange(self.numiter):
                if debug: MSG("########### Likelihood iteration %s ########"%(i))
                while True:
                    x0 = [random.uniform(np.log10(mu_bounds[0]), np.log10(mu_bounds[1]))]
                    if lencoeff_bounds[0] != lencoeff_bounds[1]:
                        x0.append(random.uniform(*lencoeff_bounds))
                    if drawmu: x0.append(random.uniform(*sd_bounds))
                    x0.extend([0 for i in range(numfeatures)]) # start coeff features at 0
                    if drawmu: x0.extend([0 for i in range(numfeatures)]) 
                    if not np.isnan(fn(x0)): break
                res = scipy.optimize.minimize(fn, x0, callback=callback, method=self.method, \
                                                  options={'maxiter': self.max_cycle_per_iter, 'xtol': 0.001, 'ftol':0.001})
                if best_res is None or (res.success and res.fun < best_res.fun):
                    best_res = res
        else:
            x0 = self.priors
            explen = numfeatures*(1+int(drawmu))+1
            if lencoeff_bounds[0] != lencoeff_bounds[1]: explen += 1
            if len(x0) != explen: ERROR("Priors is wrong length (%s vs. %s)"%(len(x0), explen))
            res = scipy.optimize.minimize(fn, x0, callback=callback, method=self.method, \
                                               options={'maxiter': self.max_cycle_per_iter, 'xtol': 0.001, 'ftol':0.001})
            best_res = res                                   
        self.best_res = best_res

        # Calculate stderr
        self.CalculateStdErrors(drawmu=drawmu, mu_bounds=mu_bounds, \
                                    beta_bounds=beta_bounds, pgeom_bounds=pgeom_bounds, \
                                    lencoeff_bounds=lencoeff_bounds, \
                                    sd_bounds=sd_bounds)

    def PartialDerivative(self, func, var=0, n=1, point=[], dx=1e-4):
        args = point[:]
        def wraps(x):
            args[var] = x
            return func(args)
        return scipy.misc.derivative(wraps, point[var], n=n, dx=dx)

    def GetLogLikelihoodSecondDeriv(self, dim1, dim2, numfeatures, drawmu, \
                                        mu_bounds=None, beta_bounds=None, pgeom_bounds=None, lencoeff_bounds=None, sd_bounds=None, dx=1e-4):
        point = list(self.best_res.x)
        deriv1_fnc = (lambda y: self.PartialDerivative(lambda x: -1*self.NegativeLogLikelihood(x, numfeatures, drawmu, \
                                                                                                   mu_bounds, \
                                                                                                   sd_bounds, \
                                                                                                   beta_bounds, \
                                                                                                   pgeom_bounds, \
                                                                                                   lencoeff_bounds), \
                                                           var=dim1, n=1, point=y, dx=dx))
        deriv2 = self.PartialDerivative(deriv1_fnc, var=dim2, n=1, point=point, dx=dx)
        return deriv2

    def GetFisherInfo(self, drawmu=False, mu_bounds=None, beta_bounds=None, pgeom_bounds=None, sd_bounds=None, lencoeff_bounds=None):
        if self.best_res is None: return
        numfeatures = len(self.loci[0].features)
        numparams = numfeatures*(1+int(drawmu))+1
        if lencoeff_bounds[0] != lencoeff_bounds[1]: numparams += 1
        fisher_info = np.zeros((numparams, numparams))
        dfs = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1]
        for i in range(numparams):
            for j in range(numparams):
                found_deriv = False
                for df in dfs:
                    deriv = -1*self.GetLogLikelihoodSecondDeriv(i, j, numfeatures, drawmu, \
                                                                    mu_bounds=mu_bounds, \
                                                                    beta_bounds=beta_bounds, \
                                                                    pgeom_bounds=pgeom_bounds, \
                                                                    lencoeff_bounds=lencoeff_bounds, \
                                                                    sd_bounds=sd_bounds, dx=df)
                    if not np.isnan(deriv):
                        found_deriv = True
                        fisher_info[i, j] = deriv
                        break
                if not found_deriv: fisher_info[i, j] = np.nan
        return fisher_info

    def CalculateStdErrors(self, drawmu=False, mu_bounds=None, beta_bounds=None, pgeom_bounds=None, sd_bounds=None, lencoeff_bounds=None):
        fisher_info = self.GetFisherInfo(drawmu=drawmu, mu_bounds=mu_bounds, beta_bounds=beta_bounds, pgeom_bounds=pgeom_bounds, sd_bounds=sd_bounds, lencoeff_bounds=lencoeff_bounds)
        try:
            self.stderrs = list(np.sqrt(np.diagonal(inv(fisher_info))))
        except np.linalg.linalg.LinAlgError:
            MSG("Could not get stderrs from fisher_info matrix \n%s"%fisher_info)
            self.stderrs = [np.nan for i in range(fisher_info.shape[0])]

    def PrintResults(self, out):
        if self.best_res is None: return
        out.write("\t".join(map(str, ["JOINT"]+list(self.best_res.x) + self.stderrs + ["%s loci"%len(self.loci)]))+"\n")
        out.flush()
