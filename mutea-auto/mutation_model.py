import numpy
from scipy.stats import norm, geom
from numpy.linalg import matrix_power
from numpy import linalg as LA

class OUGeomSTRMutationModel:
    def __init__(self, p_geom, mu, beta, allele_range, max_step=None, len_coeff=0.0, a0=0):
        if p_geom <= 0.0 or p_geom > 1:
            exit("Invalid geometric distribution probability = %f"%p_geom)

        if mu > 0.1 or mu <= 0:
            exit("Invalid mutation rate mu = %f"%mu)

        if beta < 0 or beta > 1.0:
            exit("Invalid length constraint beta = %f"%beta)

        if max_step is None:
            prob       = p_geom
            tot_rem    = 1.0
            max_step   = 0
            while tot_rem > 1e-6:
                tot_rem  -= prob
                prob     *= (1-p_geom)
                max_step += 1
            self.max_step = max_step
        else:
            self.max_step = max_step
        
        self.allele_range   = allele_range
        self.mu             = mu
        self.beta           = beta
        self.p_geom         = p_geom
        self.len_coeff      = len_coeff
        self.a0             = a0
        self.init_transition_matrix()

        # Memoized matrix info
        self.prev_forward_tmrca  = None
        self.forward_matrix      = None
        
        # Based on properties of exact geometric distribution
        self.step_size_variance = (2-p_geom)/(p_geom**2)

    def init_transition_matrix(self):
        max_step     = self.max_step
        max_n        = self.allele_range + max_step
        min_n        = -self.allele_range - max_step
        N            = 2*self.allele_range + 2*max_step + 1
        trans_matrix = numpy.zeros((N, N))
        steps        = numpy.arange(1, 2*max_n+1, 1)
        step_probs   = geom.pmf(steps, self.p_geom)

        # Fill in transition matrix                                                                                                                                             
        for i in xrange(min_n+1, max_n):
            up_prob      = min(1, max(0, 0.5*(1-self.beta*self.p_geom*i)))
            down_prob    = 1.0-up_prob
            lrem         = sum(step_probs[i-min_n-1:])
            rrem         = sum(step_probs[max_n-i-1:])
            mutrate      = 10**(numpy.log10(self.mu)+self.len_coeff*(i-self.a0))
            trans_matrix[:,i-min_n] = numpy.hstack((numpy.array([mutrate*down_prob*lrem]),
                                                    mutrate*down_prob*numpy.array(step_probs[:i-min_n-1][::-1]),
                                                    numpy.array([1-mutrate]),
                                                    mutrate*up_prob*numpy.array(step_probs[:max_n-i-1]),
                                                    numpy.array([mutrate*up_prob*rrem])))

        # Add boundaries to prevent probability leakage
        trans_matrix[:,0]     = 0
        trans_matrix[0,0]     = 1
        trans_matrix[:,N-1]   = 0
        trans_matrix[N-1,N-1] = 1
                        
        # Convert to matrix, or numpy functions won't do appropriate thing, and save for later use
        self.trans_matrix = numpy.matrix(trans_matrix)

        # Save various matrix-related variables
        self.min_n = min_n
        self.max_n = max_n
        self.N     = N

    def compute_forward_matrix(self, tmrca, allow_bigger_err=True):
        if tmrca != self.prev_forward_tmrca:
            self.prev_forward_tmrca = tmrca
            self.forward_matrix     = self.trans_matrix**tmrca
            numpy.clip(self.forward_matrix, 1e-10, 1.0, out=self.forward_matrix)

    def get_forward_matrix(self):
        return self.forward_matrix

    def get_forward_str_probs(self, start_allele, tmrca):
        self.compute_forward_matrix(tmrca)
        vec = numpy.zeros((self.N,1))                                   
        vec[start_allele-self.min_n] = 1           
        return numpy.array(self.forward_matrix.dot(vec).transpose())[0]

    def plot_transition_probabilities(self, x_vals, output_pdf, window_width=5):
        fig, axes = plt.subplots(1, len(x_vals), sharex=True, sharey=True)
        if len(x_vals) == 1:
            axes = [axes]
        for i in xrange(len(x_vals)):
            trans_probs = numpy.array(self.trans_matrix[:,x_vals[i]-self.min_n].transpose())[0]
            trans_probs = trans_probs[max(0, x_vals[i]-window_width-self.min_n): min(len(trans_probs), x_vals[i]+window_width+1-self.min_n)]
            x = range(max(0, x_vals[i]-window_width-self.min_n), min(self.max_n-self.min_n+2, x_vals[i]+window_width+1-self.min_n))
            x = numpy.array(x) + self.min_n - x_vals[i]
            axes[i].bar(x, trans_probs, align="center")
            axes[i].set_title("Allele=%d"%(x_vals[i]))
        map(lambda x: x.xaxis.set_ticks_position('bottom'), axes)
        map(lambda x: x.yaxis.set_ticks_position('left'), axes)
        map(lambda x: x.set_ylim((0.0, 0.8)), axes)
        axes[0].set_ylabel("Transition probabilities")
        axes[len(x_vals)/2].set_xlabel("STR Step Size")
        output_pdf.savefig(fig)
        plt.close(fig)
