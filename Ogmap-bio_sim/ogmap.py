#@+leo-ver=5-thin
#@+node:gte.20161001183816.3: * @file ogmap.py
#@+at
# Probability field:
#    A tool for assigning, to any element in some domain, a probability distribution.
#    Estimated by local influence methods from a survey: 
#    random samples from a `real' distribution at elements of the domain.
# 
# Node class:
#    Contains, most importantly, an estimate of a (cumulative) probability distribution:
#      which means functions for computing the cumulative probability of a random sample,
#      and the quantile of a probability
#    Also specifies where the node is
#    Possibly had a random sample from a real underlying cdf
# 
# Relevance class:
#    The influence that one element of the domain has for inferences about the cdf at another element
#    Computed from a measure of distance between Nodes, 
#    and a measure of how relevance decreases with increasing distance.
# 
# ProbField class:
#    The method for estimating a probability distribution at a Node.
#@@c
#@@language python
#@+<<imports>>
#@+node:gte.20161003082105.2935: ** <<imports>>
from __future__ import division
from numpy import *
import numpy
from numpy.random import shuffle, random_sample
import pandas as pd
import os
import random
from bisect import bisect_left as bisect
from scipy.stats import kstest, norm, spearmanr
from scipy.special import gamma as gammafn
from scipy.interpolate import interp1d # as TP
#if not '/home/gte/.ipython' in sys.path:  sys.path.append('/home/gte/.ipython')
set_printoptions(precision=3, suppress=True) # for printing np arrays
#@-<<imports>>
#@+<<gadgets>>
#@+node:gte.20161011080630.1380: ** <<gadgets>>
#@+at
# n unequal ordered points make n+1 intervals (two semi-infinite).
# bisect(_left) works nicely with this: it returns values from 0 to n inclusive,
# such that bisect(List, List[i])==i
# cum0sum() is made to work nicely as well by prepending a 0 to cumsum
# 
# ordi() makes points in the middle of n intervals between 0 and 1.
# 
#@@c
cum0sum = lambda vec: hstack([[0], cumsum(vec)])
ordi = lambda n: (.5+arange(n))/n  # linspace(0.5/n, 1-0.5/n, n)
quantiles=lambda n: (1+arange(n)) / (1.+n)  # linspace(1/(1.+n), 1-1/(1.+n), n)
Glogspace = lambda a,b,n: exp(linspace(log(a),log(b),n))
#@-<<gadgets>>
#@+others
#@+node:gte.20161001184024.1: ** Node class
#@+at
# A Node is basically a container for a cumulative probability distribution, i.e. for the functions that return a cumulative probability given a statistic vlue, and a quantile given a cumulatie probability.  It also has a guide attribute, used for ocmputing the cdf functions.
# 
# Strictly speaking, this is a class for a step-function cdf, with steps of size "steps" at locations "stepsat", with Prob() and Quantile() methods appropriate thereto.  Ogive mapping in general allows for other kernels than Heaviside along the probability axis, raising the possibility of other subclasses of a generic cdfclass.
# 
# Also stepsat is a class attribute rather than instance attribute; one could contemplate a superclass with the methods and named subclasses, created as needed in specific applications, with different stepsat vectors (e.g. surveys and subsurveys).  Not implemented here.
# 
#@@c
class Node(object):
   def __init__(self, guide):
       if isinstance(guide, (float,int)): guide=[guide] 
       self.guide = array(guide)

   @classmethod
   def set_stepsat(cls, stepsat):
       cls.stepsat = stepsat

   @classmethod
   def conform(cls, PF): # a ProbField
       cls.set_stepsat(PF.observed)
    
   def set_step_sizes(self, raw_steps):
      ogive = cum0sum(raw_steps)
      norm = ogive[-1]
      if norm==0.: norm=1.  # no info gives degenerate cdf
      self.steps = raw_steps/norm ; self.ogive = ogive/norm
      self.information = norm
      # raw_steps after redundancy adjust: maybe a bad idea?
   
   def little_set_step_sizes(self, steps, ogive):
       self.steps = steps ; self.ogive = ogive  # stepsizes instead of steps?

   def Prob(self, value):   
      ''' The probability that a random variable will be less than that value '''
      return self.ogive[bisect(self.stepsat, value)]

   def Quantile(self, prob):
       ''' The value at which that much cumulative probability is exceeded '''
       return self.stepsat[bisect(self.ogive[1:-1], prob)]  # why the restricted range?

   def meva(self, retval=False):
       self.mean =me= inner(self.steps, self.stepsat)
       msq = inner(self.steps, self.stepsat*self.stepsat)
       self.var = (msq-me*me) #/ (1-1./self.divers)  for bias correction
       self.divers = 1. / inner(self.steps, self.steps)
       if retval: return self.mean, self.var, self.divers

   def medex(self):  # index of the median
       return bisect(self.ogive[1:-1], 0.5)
#@+node:gte.20170801081339.1609: *3* Survey Node
class SurveyNode(Node):
  def __init__(self, guide, observed, rank):
    if isinstance(guide, (float,int)): guide=[guide] 
    self.guide = array(guide)
    self.observed=observed  ;  self.rank=rank

  def ownprobs(self):
    self.fullOwnProb = self.ogive[self.rank] + self.steps[self.rank]/2.
    self.jackOwnProb = self.ogive[self.rank] / (1.-self.steps[self.rank])   

  def diagnostics(self):
     #  full probability, jackknifed probability, diversity, information, covariate
     # takes the last guide for 1-d covariate ease of making nearby pairs
     self.meva()
     SO=self.ogive ; SS=self.steps
     return [SO[self.rank] + SS[self.rank]/2., SO[self.rank] / (1.-SS[self.rank]), \
               self.divers, self.information, self.guide[-1]  ]
#@+node:gte.20161122130825.1489: ** survey array from DataFrame
def SurveyArrayFromDataFrame(data_frame, observed_name, guide_names):
   columns = [observed_name] + guide_names
   return data_frame[columns].values

#@+at
# Just use 
# SurveyArrayFromDataFrame(pd.read_csv('file_name', ...)
# 
# def SurveyArrayFromCSV(csv_file, observed_name, guide_names, \
#                 delim_whitespace=False):
#    data_frame = pd.read_csv(csv_file, \
#       delim_whitespace=delim_whitespace)
#    columns = [observed_name] + guide_names
#    return data_frame[columns].values
# 
# How to sort a DataFrame on one column, kind of like an array
#     seed(423) ; L=range(len(survey)) ; 
#     shuffle(L)  # same shuffle every time
#     survey=survey.iloc[L]   # works on a DataFrame
#     self.survey_array = survey_array[observed.argsort()]
# 
#@+node:gte.20161001185821.1: ** Relevance class
#@+at
# This is the distances-first, then single shape model.
# There is potentially a shaped_distances, then add relevances model.
# 
#@@c
class Relevance(object):
   def __init__(self, kernel, distance, bandwidths):
      self.kernel=kernel ; self.distance=distance
      self.bandwidths=bandwidths
      self.numerical_init(bandwidths)

   def __call__(self, node_t, node_s):
      return self.kernel( self.distance(node_t.guide, node_s.guide) )

   def numerical_init(self, bandwidths):
       bandwidths = self.kernel.get_shape(bandwidths)
       self.distance.get_scale(bandwidths)
       #self.redundancy_and_cdfs_of_survey_nodes()

#@+<<kernel catalogue>>
#@+node:gte.20161003141412.1849: *3* <<kernel catalogue>>
#@+at
# The Kernel class corresponds to the usual 1-d kernel in kernel regression: a decreasing function of a scalar distance. By convention it is scaled so that Kernel(0)=1 and Kernel(1)=0.5.  
#@@c
class Kernel(object):
    def get_shape(self, bandwidth):
        ''' Take the first bandwidth parameter for the kernel shape and return the rest for distance scales. '''
        self.shape = bandwidth[0]
        return bandwidth[1:]

# The defining feature of the Cauchy distribution is its tails: so long that is has infinite variance.  subcauchy() isn't that extreme: it has thinner tails because its shape exponent is usually greater than 2; but it is still sort of like Cauchy both in its functional form and in that its tails are thicker than exponential or compact support distributions.
class subcauchy(Kernel):
    def __call__(self, distance):
        return 1./(1.+distance**self.shape)

class subcau_exp(Kernel):
    def __call__(self, distance):
        if distance<1.: return 1./(1.+distance**self.shape)
        return max(1.e-6, 0.5*exp(self.shape/2.*(1-distance)) )

# logistic function that eventually decays exponentially.  The shape parameter here is the offset from 0; an offset of 100 is approx like a subcauchy shape of 3.
class logist(Kernel):
    def __call__(self, distance):
      so = self.shape
      return (so + 1)/(so + (so+2)**distance)

class debug_kernel(Kernel):
    ''' It is a subclass of the kernel class only for the sake of intuition; everything is overridden.  There is no shape parameter; so, a dummy call to get_shape(). '''
    def get_shape(self, bandwidth):
        return bandwidth
    def __call__(self, distance):
        if distance<0.5: return 1.
        if distance<1.5: return 0.5
        if distance<2.5: return 0.25
        return 0.

class square_kernel(Kernel):
    def get_shape(self, bandwidth):
        return bandwidth
    def __call__(self, distance):
        return 1. if distance<1. else 0.

def display_kernel(kernel,shape,until):
    ks=kernel()
    ks.shape=shape
    distance = linspace(0,until,101)
    plot(distance, [sqrt(ks(d)) for d in distance])

#@-<<kernel catalogue>>
#@+<<distance catalogue>>
#@+node:gte.20161003141412.1851: *3* <<distance catalogue>>
#@+at
# Distance is a *scalar* function of two array guide variables.
# scaled so that relevance is 0.5 at the desired separation.
# 
# Possibilities:
#     native or ordinal
#     L1 or L2 combination
#     aspect ratio vector
#     native = linear, log, ...
# 
#@@c
class Distance(object):
    pass

class linear_distance(Distance):
    # x and y are arrays of length 1
    def __call__(self, x, y):
        return abs(x-y)[0]/self.scale
    def get_scale(self, bandwidth):
        self.scale = bandwidth[0]

class unscaled_linear_distance(linear_distance):
    def get_scale(self, bandwidth):
        self.scale = 1

class L1_distance(Distance):
    #  better to have separate scales for each direction
    def __call__(self, x, y):
        return sum([abs(x[i]-y[i]) for i in range(len(x))])/self.scale
    def get_scale(self, bandwidth):
        self.scale = bandwidth[0]

class zero_distance(Distance):
    def __call__(self, x, y):
        return 0.
    def get_scale(self, bandwidth):
        pass    

class euclid_distance(Distance):
    # isotropic: one scale for all directions
    def __call__(self, x, y):
        return sqrt(sum((x-y)**2)) / self.scale
    def get_scale(self, bandwidth):
        self.scale = bandwidth[0]    

#@+at
# This is just a form of L1.  Nothing makes y a rank variable.
# # column_stack([*,*])
# class euc_ord_distance(Distance):
#     def get_scale(self, bandwidth):
#         self.euc_scale, self.ord_scale, self.euc_weight = bandwidth
#     def __call__(self, x, y):
#         return self.euc_weight*abs(x[0]-y[0])/self.euc_scale \
#              + (1-self.euc_weight)*abs(x[1]-y[1])/self.ord_scale
#@@c
class log_distance(Distance):
    def get_scale(self, bandwidth):
        self.scale = bandwidth[0]
    def __call__(self, x, y):
        return abs(log(x)-log(y))[0]/self.scale

class ordinal_1d_distance(Distance):
    def __init__(self, PF, guide_col=-1):
        guide=PF.survey_array[:,guidecol]
        asguide = argsort(guide)
    def get_scale(self, bandwidth):
        self.scale = bandwidth[0]
    def getindex(self, x):
        self.inx = bisect(guide[asguide], x) - 1
        if x != guide[asguide][self.inx]: self.inx-= 0.5
        # the ordinal value is stored with the distance function, not the point
    def __call__(self, x, y):
        ' x is real, y is observed, ranky is known: interpolate to get rankx'
        ' x, y observed: go to rankx, ranky instead'
        pass
    
#@-<<distance catalogue>>

#@+node:gte.20161001185842.1: ** ProbField class
#@+at
# What needs to go here?
# 
# Broadly, ProbField uses a survey and a relevance function to compute the cdf at a node.
# 
# __init__() stores the survey and relevance function
# 
# Redundancy (clustering) of survey nodes: this is a property of the whole field, relevance plus survey; 
# it is used in computing influence from (pointwise) relevance.
# 
# __call__(self, node)  constructs the cdf (i.e. the steps and ogive) for the node.
# 
# Diagnostics of the whole field:  Is it a field the survey could plausibly be independent random samples from?
# Is it a field that uses the data widely (high diversity)?
# Adjust bandwidths to get good diagnostics.
# 
# Resample (resurvey) the field :  for bootstrap confidence bounds on statistics.
# Do the bootstrap methods use the resampled values or indices?
# 
# survey_array has observed first then guides
# 
#@@c
class ProbField(object):
   def __init__(self, survey_array, relevance_function, \
                adjust_for_redundancy=True):
      # shuffled_and_sorted=False (might want an option to skipp the next few lines that do this)
      # randomize
      sa=survey_array
      L=range(len(sa))
      shuffle(L) ; sa=sa[L]
      # sort by observed
      self.survey_array = sa = sa[sa[:,0].argsort()]
      self.survey = [SurveyNode(s[1:], s[0], inx) for inx, s in enumerate(sa)]
      self.samples = self.observed = sa[:,0]
      self.length=len(self.observed) ; self.forlen = range(self.length)
      self.put_cdf_step_locations_in_Node_class()
      self.REL = relevance_function
      self.adjust_for_redundancy=adjust_for_redundancy
      self.redundancy_and_cdfs_of_survey_nodes()

   def numerical_init(self, bandwidths):
      self.REL.numerical_init(bandwidths) 
      self.redundancy_and_cdfs_of_survey_nodes()

   def redundancy_and_cdfs_of_survey_nodes(self):  
      ISR = intra_survey_relevance = eye(self.length)
      for i,s in enumerate(self.survey):
        for k,t in enumerate(self.survey[i:]):
          ISR[i,k+i] = ISR[k+i,i] = self.REL(s,t)
      self.redundancy = sum(ISR, axis=0) \
         if self.adjust_for_redundancy else ones(self.length)
      for i,s in enumerate(self.survey):
         s.redundancy = self.redundancy[i]
      for s in self.survey: self(s)

   def put_cdf_step_locations_in_Node_class(self):
      Node.set_stepsat(self.observed)

   def influence(self, survey_node, target_node):
      return self.REL(survey_node, target_node) / survey_node.redundancy

   def __call__(self, node):
      raw_weights = [ self.influence(s_node, node) for s_node in self.survey ]
      node.set_step_sizes(array(raw_weights))
   
   def call__(self, node):
      raw_weights = [ self.influence(s_node, node) for s_node in self.survey ]
      ogive = cum0sum(raw_weights)
      norm = ogive[-1]
      if norm==0.: norm=1.  # no info gives degenerate cdf
      node.set_step_sizes(array(raw_weights)/norm, ogive/norm)

   #@+others
   #@+node:gte.20161005065757.1764: *3* diagnostics
   def generic_cost(self, bandwidth):
       self.numerical_init(bandwidth)
       self.make_diagnostics()
       return self.cost_kernel()

   def chatty_generic_cost(self, bandwidth):
       self.numerical_init(bandwidth)
       self.make_diagnostics()
       cost = self.cost_kernel()
       print array(bandwidth), '%.2f' % cost
       return cost

   def set_cost_kernel(self, function):
       self.cost_kernel = function
       # choose among worse_cost, crude_worse_cost, full_cost, jack_cost, ?? silent_cost_function ??
       
   def make_diagnostics(self):
       # related to how efficient and consistent the use of survey data is
       self.diagnostic_array  = qa = array([node.diagnostics() for node in self.survey])
       self.guide_sorted_diagnostic_array = array(qa[qa[:,-1].argsort()])
       # 1-d 

   def full_cost(self):
       fc = self.diagnostic_array[:,0]
       fullunif = kstest(abs(fc-0.5), lambda x : 2*x)
       fullindep = kstest([abs(fc[p]-fc[q]) for (p,q) in self.raw_pairs ], lambda x : x*(2-x))
       P = fullunif[1] * fullindep[1]  ; full = P * (1-log(P))
       if full<=0: return 1000.
       return -log(full)

   def jack_cost(self):
       jc = self.diagnostic_array[:,1]
       jackunif = kstest(abs(jc-0.5), lambda x : 2*x)
       jackindep = kstest([abs(jc[p]-jc[q]) for (p,q) in self.raw_pairs ], lambda x : x*(2-x))
       P = jackunif[1] *jackindep[1]  ;jack = P * (1-log(P))
       if jack<=0: return 1000.
       return -log(jack)

   def worse_cost(self):
       fc = self.diagnostic_array[:,0]
       fullunif = kstest(abs(fc-0.5), lambda x : 2*x)
       print fullunif
       fullindep = kstest([abs(fc[p]-fc[q]) for (p,q) in self.raw_pairs ], lambda x : x*(2-x))
       P = fullunif[1] * fullindep[1]  ; full = P * (1-log(P))
       if full<=0: return 1000.  
       jc = self.diagnostic_array[:,1]
       jackunif = kstest(abs(jc-0.5), lambda x : 2*x)
       jackindep = kstest([abs(jc[p]-jc[q]) for (p,q) in self.raw_pairs ], lambda x : x*(2-x))
       P = jackunif[1] *jackindep[1]  ;jack = P * (1-log(P))
       if jack<=0: return 1000.
       try: 
          co=self.cost_option 
          if co=='jack': return -log(jack)
          if co=='full': return -log(full)
       except:
           pass
       return -log(min(jack, full))

   def cost_components(self):
       fc = self.diagnostic_array[:,0]
       jc = self.diagnostic_array[:,1]
       fullunif = kstest(abs(fc-0.5), lambda x : 2*x)
       fullindep = kstest([abs(fc[p]-fc[q]) for (p,q) in self.raw_pairs ], lambda x : x*(2-x))
       jackunif = kstest(abs(jc-0.5), lambda x : 2*x)
       jackindep = kstest([abs(jc[p]-jc[q]) for (p,q) in self.raw_pairs ], lambda x : x*(2-x))
       return (fullunif[0], fullindep[0], jackunif[0],jackindep[0])
       
   def crude_worse_cost(self):
       fc = self.diagnostic_array[:,0]
       jc = self.diagnostic_array[:,1]
       fullunif = kstest(abs(fc-0.5), lambda x : 2*x)
       fullindep = kstest([abs(fc[p]-fc[q]) for (p,q) in self.raw_pairs ], lambda x : x*(2-x))
       #P = fullunif[1] * fullindep[1]  ; full = P * (1-log(P))
       full = 2*fullunif[0] + fullindep[0]
       if full<=0: return 1000.    
       jackunif = kstest(abs(jc-0.5), lambda x : 2*x)
       jackindep = kstest([abs(jc[p]-jc[q]) for (p,q) in self.raw_pairs ], lambda x : x*(2-x))
       #P = jackunif[1] *jackindep[1]  ;jack = P * (1-log(P))
       jack = 2*jackunif[0] + jackindep[0]
       if jack<=0: return 1000.
       try: 
          co=self.cost_option 
          if co=='jack': return jack
          if co=='full': return full
       except:
           pass
       return max(jack, full)
   #@+node:gte.20180115110615.1: *4* 1-d version
       
   # This is a special 1-d case where the pairs for indep testing are obvious.
   # This is why a guide is part of Node.diagnostics.
   def test_1d_stuff(self):
       fullcol=2 ; jackcol=3
       qsa = self.guide_sorted_diagnostic_array
       fullunif = kstest(abs(qsa[:,1]-0.5), lambda x : 2*x)
       fullindep = kstest([abs(qsa[i,1]-qsa[i+1,1]) for i in range(0,self.length,2) ], lambda x : x*(2-x))
       P = fullunif[1] * fullindep[1]  ; full = P * (1-log(P))
       if full<=0: return 1000.
       jackunif = kstest(abs(qsa[:,2]-0.5), lambda x : 2*x)
       jackindep = kstest([abs(qsa[i,2]-qsa[i+1,2]) for i in range(0,self.length,2) ],  lambda x : x*(2-x))
       P = jackunif[1] * jackindep[1]  ; jack = P * (1-log(P))
       if jack<=0: return 1000.
       self.cost_components = [fullunif[0], fullindep[0], jackunif[0], jackindep[0]]
       return -log(min(jack, full))
       
       

   #@+at
   #    unif = kstest(abs(self.IRS_input[0]-0.5), lambda x : 2*x)
   #    # omit the jack extremes 0 and 1?
   #    indep = kstest(self.IRS_input[1], lambda x : x*(2-x) )
   #    product = unif[1]*indep[1]
   #    combine = product * (1-log(product))
   #    self.IRS_results = (combine, unif[1],unif[0], indep[1],indep[0])
   #    if combine<=0: return 1000
   #    cost = -log(combine)
   #    if isnan(cost): return 1000
   #    return cost
   #@+node:gte.20180102065719.1: *3* cost function

   # run -i ../nonpareil
   # agelen=loadtxt('AP_females',usecols=(10,11),skiprows=1,delimiter=",")
   # agelen[:,1]+=0.01*random_sample(shape(agelen)[0])
   # RPF = ProbField(agelen[:,0], agelen[:,1], subcauchy, linear_distance, [3, 1], True)
   # RPF.numerical_init([3,1000])
   # RPF.raw_pairs = find_1d_nearby_pairs(RPF.survey_array[:,0])
   #@+node:gte.20180102065826.1: *4* uniform
   def jack_prob(self,i):
       su = self.survey[i]
       return su.ogive[i]/(1-su.steps[i])

   def plausible_uniform_probs(self):
       probs = array([self.jack_prob(i) for i in self.forlen])
       return kstest(abs(probs[1:-1] - 0.5), lambda x : 2*x)

       # returns KS deviation and its probability

   #@+node:gte.20180102065859.1: *4* independent
   def jack_pair_prob_diff(self, i, k):
       # assert k>i
       si = self.survey[i]
       prob_i = si.ogive[i] / (1 - si.steps[i] - si.steps[k])
       sk = self.survey[k]
       prob_k = (sk.ogive[k] - sk.steps[i]) / (1 - sk.steps[i] - sk.steps[k])
       return abs(prob_i - prob_k)

   def plausible_indep_pairs(self, display=False):
       assert self.raw_pairs
       # no fix: making pairs is not the job of ProbField
       # gets called once and used repeatedly with different bandwidths
       # haven't decided yet where responsibility lies
       prob_pairs = [self.jack_pair_prob_diff(p[0], p[1]) \
                for p in self.raw_pairs] # if self.length>self.twindict[i]>i ]
       if display:
           G.plot_cdf(prob_pairs, srted=False)
           x=linspace(0,1,100)
           plot(x, x*(2-x))
       return kstest(prob_pairs, lambda x : x*(2-x) )

   #@+node:gte.20180102070007.1: *4* combined cost
   def plausible_unif_indep(self):
       product = self.plausible_uniform_probs()[1] * self.plausible_indep_pairs()[1]
      # print self.plausible_uniform_probs()[1]
     #  print self.plausible_indep_pairs()[1]
    #   print product
       return product * (1 - log(product))
   

   def silent_cost_function(self, bandwidths):
       self.numerical_init(bandwidths)
       pui = self.plausible_unif_indep()
       if pui<=0: return 1000
       cost = -log(pui)
       if isnan(cost): return 1000.
       return cost

   def cost_function(self, bandwidths):
       cost = self.silent_cost_function(bandwidths)
       print bandwidths, 'cost=', cost
       return cost

   import cma
   # 
   # cma.fmin(SF.cost_function, bandwidth, sigma0=0.2,
   #          options={
   #             'scaling_of_variables': [3.503,  37.341, 1.861],
   #             'tolfun': 1,
   #             'tolx': 0.3})
   # 
   #@+node:gte.20170328083829.1514: *3* resample
   #  for bootstrap confidence bounds 
   def resampled_survey_indices(self):
       rasa = random_sample(self.length)
       return [ bisect(s.ogive[1:-1], rasa[i]) for i,s in enumerate(self.survey) ]
       # list of integer indices 
       # chooses the index of each step in proportion to its height
       # in the cdf at the survey node

   def resampled_survey_observations(self, obs=None):
       if obs is None: obs=self.observed
       # the option of specifying obs is for double bootstrapping
       return obs[self.resampled_survey_indices()]
    
   #@+at
   # [bisect(s.ogive[1:-1], r) for s, r in zip(self.survey, random_sample(self.length) )]
   #@+node:gte.20170804110555.1: *3* IRS test
   #@+at
   # This is restricted to nearby pairs defined by adjacency of the last guide component:
   #     this is built into the use of qsa and abs(qsa[i,1]-qsa[i+1,1])
   #     and even here it will fail if self.length is odd.
   #@@c
   def IRS_cost_init(self, fj='worst', report=False):
       self.fj = fj ; self.report=report

   def IRS_cost(self, bandwidth):
       global fullunif, fullindep, jackunif, jackindep
       self.numerical_init(bandwidth)
       self.make_diagnostics()
       if self.report>2:
         print self.diagnostic_array[:,1]
         print self.diagnostic_array[:,2]
       qsa = self.guide_sorted_diagnostic_array
       if self.report>1: print '\nbandwidths: ', bandwidth
       if self.fj[:2]=='ja' and not self.report>1: full=1
       else:
           unif = kstest(abs(qsa[:,1]-0.5), lambda x : 2*x)
           indep = kstest([abs(qsa[i,1]-qsa[i+1,1]) for i in range(0,self.length,2) ], lambda x : x*(2-x))
           # restricted to the special case of 1-d nearby pairs and to self.length being even
           # indep statistic needs generalizing to deal with other ways of making pairs
           # what about omitting the automatic extremes 0 and 1?
           P = unif[1] * indep[1]  ; full = P * (1-log(P))
           # kstest returns [deviation, pvalue]
           if self.report>1: print 'full', array([P, full, -log(full), unif[0], unif[1], indep[0], indep[1] ])
           if full<=0: return 1000.
       if self.fj[:2]=='fu' and not self.report>1: jack=1
       else:
           unif = kstest(abs(qsa[:,2]-0.5), lambda x : 2*x)
           indep = kstest([abs(qsa[i,2]-qsa[i+1,2]) for i in range(0,self.length,2) ], lambda x : x*(2-x))
           P = unif[1] * indep[1]  ; jack = P * (1-log(P))
           if self.report>1: print 'jack', array([P, jack, -log(jack), unif[0], unif[1], indep[0], indep[1] ])
           if jack<=0: return 1000.
       if self.report>1: print '[%.3f %.3f] :  (%.3f, %.3f)' % (bandwidth[0], bandwidth[-1], jack, full)
       answer = -log(min(jack, full)) 
       if self.report: print bandwidth, answer  # OR delete this line and have cost and silent_cost
       return answer
   #@-others
#@+node:gte.20161013070538.1772: ** Independent Random Samples
#@+at
# Routines for computing the probability of survey observations in the cdfs at their respective nodes. 
# 
# ** probabilities routines find arrays of cumulative probabilities of observations in their respective cdfs, and differences of cumulative probabilities of nearby pairs.  Separate routines for full and jackknifed cdfs.
# 
# Routine to choose whether to jackknife
# 
# Utility routine to choose pair indices in the easy, 1-d guide case.  
# n-dimensional guides need bespoke routines.
# 
# Compute the KS probabilities of the two arrays and combine them into a single probability.
# Wrapper to express all this as a function of bandwidth only, for feeding to optimization routines.
#@+node:gte.20170305083705.1565: *3* full probabilities
def make_full_IRS_input(self):
    own_probs = [ su.ogive[i] + su.steps[i]/2.  \
      for i, su in enumerate(self.survey) ]
    pair_prob_diffs = [ abs(own_probs[p[0]]-own_probs[p[1]]) \
      for p in self.pair_indices ]
    self.IRS_input = (array(own_probs), array(pair_prob_diffs))   
#@+at
# def find_own_probs(self):
#     self.own_probs = array([su.ogive[i] + su.steps[i]/2. \
#        for i, su in enumerate(self.survey) ])
# 
# def own_prob(self,i):
#     su = self.survey[i]
#     return su.ogive[i] + su.steps[i]/2.
# 
# def pair_prob_diff(self, i, k):
#     self.prob_pairs = abs( self.own_probs[i] - self.own_probs[k] )
# 
#@+node:gte.20170305083705.1566: *3* jackknifed probabilities
def make_jack_IRS_input(self):
    own_probs = array([ su.ogive[i] / (1 - su.steps[i])  \
      for i, su in enumerate(self.survey) ])
    # only the function of (ogive,steps) is changed
    pair_prob_diffs = []
    for i,k in self.pair_indices:
       si, sk = self.survey[i], self.survey[k]
       prob_i = si.ogive[i] / (1 - si.steps[i] - si.steps[k])
       prob_k = (sk.ogive[k] - sk.steps[i]) / \
                (1 - sk.steps[i] - sk.steps[k])
       pair_prob_diffs.append(abs(prob_i-prob_k))
    pair_prob_diffs = array(pair_prob_diffs)  
    self.IRS_input = (own_probs, pair_prob_diffs)    
#@+node:gte.20170724143539.1: *3* choose whether to jackknife
def jack_IRS_type(self, jack):
    self.make_IRS_input = self.make_jack_IRS_input if jack else self.make_full_IRS_input
#@+node:gte.20170305083705.1564: *3* make linear pairs
def make_linear_pair_indices(self, inx=1):
   a = argsort(self.survey_array[:,inx])
   # recall column 0 is the observed random variable
   self.pair_indices = \
      [ [a[i], a[i+1]] for i in range(0,self.length-1, 2) ]
   # doesn't work when self.length is odd
   for p in self.pair_indices:
      if p[0]>p[1]: p[0],p[1]=p[1],p[0]    
#@+node:gte.20161013070538.1774: *3* single cost function
def test_IRS_input(self):
   unif = kstest(abs(self.IRS_input[0]-0.5), lambda x : 2*x)
   # omit the jack extremes 0 and 1?
   indep = kstest(self.IRS_input[1], lambda x : x*(2-x) )
   product = unif[1]*indep[1]
   combine = product * (1-log(product))
   self.IRS_results = (combine, unif[1],unif[0], indep[1],indep[0])
   if combine<=0: return 1000
   cost = -log(combine)
   if isnan(cost): return 1000
   return cost


def silent_cost_function(self, bandwidths):
    self.numerical_init(bandwidths)
    self.make_IRS_input()
    return self.test_IRS_input()

def cost_function(self, bandwidths):
    cost = self.silent_cost_function(bandwidths)
    print bandwidths, 'cost=', cost
    return cost

def silent_approx_cost(self, bandwidths):
    # (weighted) sum of KS deviations
    self.numerical_init(bandwidths)
    cost = 2 * self.plausible_uniform_probs()[0] + self.plausible_indep_pairs()[0]
    # crudely multiply unif by 2 because twice as many points
    #print bandwidths, cost
    return cost

def approx_cost_function(self, bandwidths):
    cost = self.silent_approx_cost(bandwidths)
    print bandwidths, 'cost=', cost
    return cost

def unif_cost_function(self, bandwidths):
    self.numerical_init(bandwidths)
    pui = self.plausible_uniform_probs()[0] # KS deviation
    if pui<=0: return 1000
    print bandwidths, 'cost=', pui
    return pui

# folded version of unif
# !!!  make probs an argument computed elsewhere
def display_IRS(self):
    unipro = self.plausible_uniform_probs()[1]
    indpai = self.plausible_indep_pairs()[1]
    absci = sort(abs(self.prob_unif - 0.5))
    plot(absci, ordi(self.length), 'k')
    plot([0, .5], [0,1], 'k:')
    plot(sort(self.prob_pairs), ordi(self.length//2), 'r')
    ar=linspace(0,1,100) ; plot(ar, (lambda x : x*(2-x))(ar), 'r:')
    text(.1,.8,'%.3f' % unipro)
    text(.6,.2,'%.3f' % indpai)

#@+node:gte.20161102080849.1419: *3* commands to estimate bandwidth
#@+at
# 
# from scipy.optimize import fmin_powell
#  
# fmin_powell(self.cost_function, bandwidths, xtol=0.03, ftol = 0.01)
# 
# 
# import cma
# 
# cma.fmin(SF.cost_function, bandwidth, sigma0=0.2,
#         options={
#         'scaling_of_variables': bandwidth,
#         'tolfun': 1,
#         'tolx': 0.3 })
#@+node:gte.20161114090241.1483: ** nearby pairs
# returns a list of pairs of indices

#@+others
#@+node:gte.20170305083705.1563: *3* general version
def find_nearby_pairs(guide, distance, extra_pairs=0, \
    temperature=1., chatty=False, ranseed=44):
    # use:  PF.raw_pairs = find_nearby_pairs(...)
    # temperature scaled to bandwidth
    global twin, indices, ncity
    seed(ranseed)
    ncity = len(guide)
    nover = 100*ncity ; nlimit = 10*ncity
    indices = ncity + ncity%2 + 2*extra_pairs

    d=zeros([indices, indices])
    for i in range(ncity):
        for k in range(i+1,ncity):
            d[i,k] = d[k,i] = distance(guide[i],guide[k])
    # (the distance to a noncity index from anywhere is zero)

    # initial try at pairs
    sumlength = 0.
    twin = [0]*indices ; length = [0.]*indices
    for i in range(0,indices,2):
        twin[i]=i+1 ; twin[i+1]=i
        length[i] = length[i+1] = d[i,i+1]
        # or do all the pairs of lengths at once?
        sumlength += length[i]

    if chatty:
        clf()
        print 'initial mean length %.3f'%(sumlength/(ncity//2))
        subplot(211)
        plotpairs(twin, guide, indices, ncity)

    trials=0
    for j in range(100): # temperature steps
      success=0 ; sumchange=0
      changes=[] ; decreases=0 ; accepted_changes=0
      for k in range(nover):
        trials += 1  
        g1 = randint(indices) ; g2 = randint(indices) # guides
        t1 = twin[g1] ; t2 = twin[g2] # and their twins
        if g2==g1 or t1==g1: continue  # require two different pairs
        dg = d[g1,g2] ; dt = d[t1,t2] # distances of new partners
        change =  dg + dt - length[g1] - length[g2]
        changes.append(change)
        if change<0: decreases += 1

        if metrop(change, temperature):     
            #  change partners
            accepted_changes += 1
            twin[g1]=g2 ; twin[g2]=g1
            length[g1] = length[g2] = dg
            twin[t1]=t2 ; twin[t2]=t1
            length[t2]=length[t1]= dt
            success +=1
            sumchange += change
        if success > nlimit:  
            break
      # end of trials at a constant temperature
      if chatty:
          print '%d trials' % k,
          print 'mean %.3f  sdev %.3f changes' % (mean(changes), sqrt(var(changes))),
          print '%d decreases;  %d accepted' % (decreases, accepted_changes)

      sumlength += sumchange
      # print '%.3f %d %.3f %.3f'%(temperature, success, sumchange, sumlength)
      sumchange=0
      temperature *= 0.9
      if success==0: break
    print 'final mean length %.3f  trials %d'%(sumlength/(ncity//2), trials)
    if chatty:
        subplot(212)
        plotpairs(twin, guide, indices, ncity)
    pairs = [ [g, twin[g]] for g in range(ncity) if g < twin[g] < ncity ]
    return pairs

def metrop(change, temperature):
    return change<0. or random_sample()<exp(-change/temperature)

def plotpairs(pair, nodes, indices, ncity, inx=[0,1]):
    # figure()
    i1=inx[1] ; i0=inx[0]
    for g in range(ncity):
        t = pair[g]
        if t>g:
            if t>=ncity:
               #pass 
               scatter(nodes[g][i1],nodes[g][i0])
            else:
              plot([nodes[g][i1], nodes[t][i1] ], \
                   [nodes[g][i0], nodes[t][i0] ])

#@-others
#@+node:gte.20170816171959.1: ** Domain class
#@+at 
# Nodes are the common currency of ProbFields; a Domain should have a list of nodes for the field to act on.
# The nodes come from a DataFrame, possibly after a set of manipulations.
# 
# UNDER CONSRUCTION
# Ambiguous whether the intermediate methods are to create new DataFrames or modify the original one.
# domain_object(DF) can be modified without changing the input DF because it's a class instance not an alias.
# But we want a way of retaining intermediate DFs somehow, like after winnowing and areas.
# What steps are irrevocalble?
# 
# Should guide_names be an argument to __init__, for the routine case where there is nothing to do except make_node_list?  Or as an optional argument?
#@@c
class Domain(object):
    def __init__(self, DF): # DataFrame
       self.F =  DF
       # can initialize from a .csv file: DOM=Domain(pd.read_csv('file.name'))

    #@+<< winnow >>
    #@+node:gte.20171218085522.1: *3* << winnow >>
    #@+<< associated distance functions >>
    #@+node:gte.20171218085812.1: *4* << associated distance functions >>

    def approx_earthdeg_sqdist(DF, i, k, cutoff, latscale=111.11, lonscale=77.28):  # 3Ps latitude
        di=DF.iloc[i] ; dk=DF.iloc[k]
        ati, oni = di[['lat','long']]  # comitted to a naming convention
        atk, onk = dk[['lat','long']]
        if atk < ati: return 'unsorted'
        if (atk-ati) * latscale > cutoff: return 'break here'
        atdiff=(ati-atk)*latscale ; ondiff=(oni-onk)*lonscale
        return atdiff*atdiff+ondiff*ondiff

    def euclid_sq_dist(DFa, i, k, cutoff):
        dlat = DFa.iloc[i,0] - DFa.iloc[k,0]
        dlon = DFa.iloc[i,1] - DFa.iloc[k,1]
        if dlat < 0: return 'unsorted'
        if dlat > cutoff: return 'break here'
        return dlat*dlat + dlon*dlon


    #
    #@-<< associated distance functions >>

    def winnow(self, cutoff, sqdist=approx_earthdeg_sqdist):
        lFD=len(self.F)
        sqcut=cutoff*cutoff 
        retain=[True]*lFD
        for i in range(lFD):
            if not retain[i]: continue
            for k in range(i+1,lFD):
                if not retain[k]: continue
                sqd = sqdist(self.F,i,k, cutoff)
                if sqd == 'unsorted':
                    print 'OOPS!  sort the input by increasing lat'
                    return
                if sqd == 'break here': break
                if sqd < sqcut: retain[k]=False
                #print i, k, F.index[i], F.index[k], sqd ?? what is index?
        self.F = self.F[retain]
    #@-<< winnow >>
    #@+<< areas around points >>
    #@+node:gte.20171218085940.1: *3* << areas around points >>

    import matplotlib.tri as deltri
    deg_rad = pi/180.

    #@+<<associated area functions>>
    #@+node:gte.20180107155607.1: *4* <<associated area functions>>
    deg_rad = pi/180.

    def make_deg_tri_area(lonlat, triangle):
       TV=array([lonlat[t] for t in triangle])
       # awaiting array formulation?

    def app_SphDist(a,b):
        ''' The distance between two points on the surface of the earth given (lat, long) coordinates in degrees.  
        An approximate formula good for separations less than a degree [and not too close to the poles?] '''
        earthrad = 6367  # km radius at 55 N (6366 for a 40,000km circumference sphere)
        return earthrad * deg_rad * \
               hypot( a[0]-b[0], (a[1]-b[1])*cos(deg_rad*(a[0]+b[0])/2) )   

    def earth_deg_sqdist(DFa,i,k):
        x = app_SphDist(DFa[i], DFa[k])
        return x*x

    def earth_deg_tri_area(latlon, triangle):
        TV=array([latlon[t] for t in triangle])
        da = app_SphDist(TV[0],TV[1])
        db = app_SphDist(TV[0],TV[2])
        dc = app_SphDist(TV[1],TV[2])
        ds = (da + db + dc)/2
        if not (ds>da and ds>db and ds>dc): 
           # print ds,da,db,dc
           return 0.
        return sqrt(ds * (ds-da) * (ds-db) * (ds-dc))

    #@-<<associated area functions>>

    def areas(self, xname, yname, outgroup=None, earthdeg=True, ret_triang=False):
        sketch_boundary=False
        self.F['x'] = self.F[xname] * (sin(deg_rad*self.F[yname]) if earthdeg else 1.)
        inxy = self.F[['x', yname]].values
        nVert = len(inxy) 
        if outgroup:
            outgroup['x'] = outgroup[xname] * (sin(deg_rad*outgroup[yname]) if earthdeg else 1.)
            outxy = outgroup[['x', yname]].values
        augmented = vstack([ inxy, outxy]) if outgroup else inxy
        triang = deltri.Triangulation(augmented[:,0], augmented[:,1])
        TR = triang.triangles # array of triples of vertex indices
        boundary = [False]*nVert
        indomain=[]
        for i,T in enumerate(TR):
          (a,b,c) = sorted(T) 
          if c<nVert: indomain.append(i)
          elif b<nVert:
             boundary[a]=boundary[b]=True
             if sketch_boundary: 
                plt.plot([self.F.iloc[a][1],self.F.iloc[b][1]],  [self.F.iloc[a][0],self.F.iloc[b][0]], 'k:') 
          elif a<nVert: boundary[a]=True
        domTR=TR[indomain]
        # is the difference obsolete once we have squashed longitudes?
        make_triangle_area =      earth_deg_tri_area if earthdeg else make_km_tri_area
        vareas=zeros(len(inxy))
        for T in domTR:
          A3 = make_triangle_area(inxy, T) /3.
          for t in T: vareas[t] += A3 
        self.F['area'] = vareas
        if ret_triang: return triang

    #@+at
    #@@c
    # plot triangulation
    def plotria(vertices, triangles): 
      for t in triangles:
        a,b,c=t[0],t[1],t[2]
        plot( [vertices[a][1],vertices[b][1]] , [vertices[a][0],vertices[b][0]] )
        plot( [vertices[a][1],vertices[c][1]] , [vertices[a][0],vertices[c][0]] )
        plot( [vertices[c][1],vertices[b][1]] , [vertices[c][0],vertices[b][0]] )
    #@-<< areas around points >>
    #@+<< select a subset >>
    #@+node:gte.20171218084412.1: *3* << select a subset >>
    def select(self, criterion):
        self.F = self.F[criterion(self.F)]
        # criterion is a function that returns True or False depending if the columns satisfy it
        # e.g. criterion = lambda x: x.lat<55 & x.lat >50
    #@-<< select a subset >>

    def make_node_list(self, guide_names):
        self.N = [ Node(x) for x in self.F[guide_names].values ]
        for (n, a) in zip(self.N, self.F.area): n.weight=a

#@+others
#@+node:gte.20180103154742.1: *3* triangulation
#@+at
# triangulation with concavities (bites out of the convex set)
# from matplotlib.delaunay import delaunay as mat_delaunay
# ## scipy.spatial.Delaunay(points, furthest_site=False, incremental=False, qhull_options=None)
# might be better?
# ## Or: The matplotlib.delaunay module was deprecated in version 1.4. Use matplotlib.tri.Triangulation instead.
# 
# degrad = 2*pi/360.
# def degcos(theta): return cos(theta*degrad)
# def my_delaunay(ar):  # changing the interface  OGMAP SPECIFIC
#     return mat_delaunay(ar[:,0],ar[:,1]*degcos(ar[:,0]))[2]
# 
# def concave_delaunay(vertices, nonvertices=None):
#     # vertices nparrays with lat and long the first columns
#     # a distance of 1 in long is a distance of xscale in lat
#     # so 1 degree longitude at 60N is 0.5 degrees latitude
#     if nonvertices==None: return my_delaunay(vertices)
#     nVer = len(vertices)
#     allvert = vstack([vertices, nonvertices])
#     dela = my_delaunay(allvert)
#     concave_tri_indices=[] 
#     for i,t in enumerate(dela):
#         for v in t:
#             if v>=nVer: break
#         else: concave_tri_indices.append(i)
#     return dela[concave_tri_indices]
#@+node:gte.20170918084502.1588: *3* prob field relevant to domain
def relevant_prob_field(DOM, PF, limit_area=3, chatty=False):
    global effective_survey_areas
    # 3 sq km is a bit under 1% of the area 'represented' by a 2J3KL Campelen survey set
    le = PF.length
    effective_survey_areas = zeros(le)
    for v in DOM.N:
        PF(v)
        effective_survey_areas += v.steps * v.weight
    rel_inx = [i for i,w in enumerate(effective_survey_areas) if w>limit_area]
    if chatty: print 'length of subsurvey %d   limit area %.3g' \
           % (len(rel_inx), limit_area)
    #prob_field.effective_survey_areas = effective_survey_areas
    #prob_field.rel_inx = rel_inx
    return ProbField(PF.survey_array[rel_inx], PF.REL)
    # assume redundancy adjustment

#@+at
# Find a better criterion than limit area: something related to how the effective area of the survey node is much smaller than the (median of?) the survey nodes that `matter'.  Needs sorting rather than judging each node in isolation.
# 
# Make use of a FieldinDomain object and its methods for selecting relevant indices.
# Based on the current (in a decreasing sorted) effective node area of the survey node being less than cutoff (5%?) of the median of all preceding ones.
#@+node:gte.20180107155321.1: *3* triangle area functions
deg_rad = pi/180.

def make_deg_tri_area(lonlat, triangle):
   TV=array([lonlat[t] for t in triangle])
   # awaiting array formulation?

def app_SphDist(a,b):
    ''' The distance between two points on the surface of the earth given (lat, long) coordinates in degrees.  
    An approximate formula good for separations less than a degree [and not too close to the poles?] '''
    earthrad = 6367  # km radius at 55 N (6366 for a 40,000km circumference sphere)
    return earthrad * deg_rad * \
           hypot( a[0]-b[0], (a[1]-b[1])*cos(deg_rad*(a[0]+b[0])/2) )   

def earth_deg_sqdist(DFa,i,k):
    x = app_SphDist(DFa[i], DFa[k])
    return x*x

def earth_deg_tri_area(latlon, triangle):
    TV=array([latlon[t] for t in triangle])
    da = app_SphDist(TV[0],TV[1])
    db = app_SphDist(TV[0],TV[2])
    dc = app_SphDist(TV[1],TV[2])
    ds = (da + db + dc)/2
    if not (ds>da and ds>db and ds>dc): 
       # print ds,da,db,dc
       return 0.
    return sqrt(ds * (ds-da) * (ds-db) * (ds-dc))

#@-others
#@+node:gte.20180110111357.1: ** ProbField in Domain class
class ProbFieldinDomain(object):
    def __init__(self, PF, DOM):
        self.DOM=DOM
        self.PF = PF
        self.survey_concs=PF.observed
        self.domain_areas=DOM.F.area.values
        self.weights=column_stack([N.steps for N in DOM.N])
        # steps will have to have been done beforehand or else in a try block
        # IMV, for example, will have made steps
                
        #  do all these whether we need them or not: quick matrix mults
        self.mean_domain_concs = dot(self.survey_concs, self.weights)
        self.domain_masses = self.mean_domain_concs * self.domain_areas
        self.survey_areas = dot(self.weights, self.domain_areas)

    def relevant_survey_indices(self, cutoff):
        '''The survey nodes with the largest effective areas in the domain, 
        sufficient that their sum exceeds the cutoff fraction (like 99%) 
        of the total domain area.'''
        suar = self.survey_areas
        rasuar = range(len(suar))
        if cutoff ==1: return rasuar
        ars = argsort(suar) ; ras = ars[::-1]
        csas = cumsum(suar[ras]) # ranked cumulative sum
        target = cutoff * csas[-1]
        for icut in rasuar:
            if csas[icut]>target: break
        return ras[range(icut+1)]

    #@+others
    #@+node:gte.20180202070401.1: *3* measures of unevenness

    def area_frac_for_mass_frac(self, mass_frac=0.5):
        '''The summed areas of the domain nodes with the highest mean concentrations, 
        sufficient to encompass the massfrac fraction (like 50%) of the integrated mean value 
        over the field over the domain.'''
        C=self.mean_domain_concs ;  rasC = argsort(C)[::-1] 
        A=self.domain_areas ;  M = self.domain_masses
        cm = cumsum(M[rasC])
        target = mass_frac * cm[-1]
        inx = bisect(cm, target)
        return sum(A[rasC][:(inx+1)]) / sum(A)
        # interp1d might be slightly more accurate and would allow mass_frac to be a list

    def mass_frac_for_area_frac(self, area_frac=0.25):
        C=self.mean_domain_concs ;  rasC = argsort(C)[::-1] 
        A=self.domain_areas ;  M = self.domain_masses
        ca = cumsum(A[rasC])
        target = area_frac * ca[-1]
        inx=bisect(ca, target)
        return sum(M[rasC][:(inx+1)]) / sum(M)

    def make_republic_array(self, selector):
        ''' for 'repubic plots' of mass in a region by stratum (the selector in that case) '''
        DOM = self.DOM.F
        strata = sorted(unique(DOM[selector]))
        DOM['mass'] = self.domain_masses
        rd=[]
        for s in strata: rd.append([s, sum(DOM[DOM[selector]==s].mass) ])
        self.republic_array = array(rd)  # doesn't save stratum number

    def make_republic_dictionary(self, selector):
        ''' for 'repubic plots' of mass in a region by stratum '''
        DOM = self.DOM.F
        strata = sorted(unique(DOM[selector]))
        DOM['mass'] = self.domain_masses
        rd={}
        for s in strata: rd[s] = sum(DOM[DOM[selector]==s].mass)
        self.republic_dictionary = rd
        
    def mass_frac_of_selected_subset(self, subset):
        thing = pd.DataFrame({'stratum':self.DOM.F.stratum, 'mass':self.domain_masses})
        subthing = thing[thing.stratum.isin(subset)] 
        return sum(subthing.mass) / sum(thing.mass)

    #@+node:gte.20180125075130.1562: *3* augmented cost function
    def augmented_jack_cost(self, bandwidth):
        PF=self.PF ; DOM=self.DOM
        PF.numerical_init(bandwidth)
        PF.make_diagnostics()
        for N in DOM.N: PF(N)
        c = PF.jack_cost()
        imv = WeightedSumofMeans(PF, DOM.N)
        dive = median(PF.diagnostic_array[:,2])
        print bandwidth, array([c,imv.base,dive])
        return c
    #@-others
#@+node:gte.20180129082827.1567: ** Uneven class
#@+at
# An Uneven object is initialized with a DataFrame with the areas, mean concentrations, and masses of a Domain, plus a column of labels of potential interest.  interpret_survey_series() will produce a list of such DFs, one for each year.
# 
# Basically this is a 'list' of attributes of the Nodes of a domain.
#@@c
from bisect import bisect

class Uneven(object):
    def __init__(self, PFiDO):
        AS = argsort(PFiDO.domain_areas)
        areas = PFiDO.domain_areas[AS]
        self.carea = cumsum(areas) ; self.totarea = self.carea[-1]
        self.mass = areas * PFiDO.mean_domain_concs[AS]
        self.cmass = cumsum(self.mass) ; self.totmass = self.cmass[-1]
        self.labels = PFiDO.DOM.F.stratum[AS]
    
    def area_with_mass_frac(self, fractions):
        answer=[]
        for f in fractions:
           inx = bisect(self.cmass, f*self.totmass)
           answer.append(self.carea[inx] / self.totarea)
        # interp1d might be slightly more accurate
    
    def mass_frac_of_selected_subset(self, subsets):
        # subset is a list of labels considered relevant
        # ?? need to have something in Uneven object that knows about labels.
        answer=[]
        for s in subsets:
           subDF = self.DF[self.DF.label.isin(subset)]
           answer.append(subDF)
        return sum(subDF.mass)
    
    def republic_list(self):
        return([ [label, sum(DF[DF.label==label].mass)] for label in self.labels])
        
    def gini_coeficient(self):
        pass
        # copy code from use_ogmap()
#@+node:gte.20170616083209.1534: ** Statistic class
#@+at
# The ability of ProbabField to compute can be realized in practice only at a finite set of nodes.  So a generic statistic consists of getting a list of nodes, computing their cdfs, and doing something with the result.
# In theory it is possible that additional nodes would be required depending on the results of partial calculations; but we won't go there (in any case that would only be a more complicated essence() ).
# 
# The calculations can be split into 'design' (__init__) and 'observation' (__call__) parts.  The design part deals with the step sizes of the cdf at the survey nodes; the observation part then takes account of the values observed by the survey.  The reason for doing this:  We are interested in Monte Carlo estimates of confidence bounds -- this is natural when we are claiming to have estimated the cdf everywhere.  If we have done the 'design' part once and saved the results, then it can sometimes be very quick to resample the cdf at the survey nodes and repeat the 'observation' part of the statistic with the resampled synthetic observations.
# 
# Remember, though, that the resampled observations will not be in increasing order (unlike Node.stepsat).  Statistics that depend on order, like quantiles, will need some thought.
# 
# A Statistic is a function; what we normally think of as the value of the statistic is stored in the attribute <.base> 
#@@c
from numpy import *

class Statistic(object):
    def __init__(self, field, nodes, extras=None, sdev=False):
        # allow for nodes and node weight arguments to come bundled:
        #if not isinstance(nodes[0], Node):
        #   (self.nodes, extras) = nodes
        self.nodes=nodes 
        self.field=field ; self.length=field.length ; self.observed = field.observed
        self.sdev=sdev
        for N in self.nodes: field(N)
        #self.cdfs = [field(N, retval=True) for N in nodes]
        self.has_variance=False  # OR class StatisticWithVariance(Statistic) ?
        self.essence_of_design(extras)
        self.base = self.T1 = self()
        
    def essence_of_design(self, extras):
        ''' Work to be done once with the design part of the statstic '''
        pass
    
    def essence_of_result(self, observed):
        ''' Work to be done with the observation part of the statistic:
            many times as the field is resampled. '''
        pass
    
    def __call__(self, observed=None):
        if observed is None:  observed=self.observed
        return self.essence_of_result(observed)
    
    #@+<<bootstrap confidence bounds>>
    #@+node:gte.20170616083209.1538: *3* <<bootstrap confidence bounds>>
    #@+at
    # There is a statistic T0 of a probability field F0.
    # There is a survey S1 of F0 from which we infer F1 and T1.
    # 
    # The ultimate aim if to produce two functions:
    #     confidence_bound(self, probability)
    #     P_exceed_T0(self, mass)
    # 'mass' here is a handy name for the randomly varying quantity we're interested in.
    # 
    # We simulate many resamples {S2} of F1, and make the corresponding {F2, T2}.
    # From the basic bootstrap conjecture:
    #      T2 differs from T1 in the same way that T1 differs from T0 
    # we form a set of candidate values for T0:  {CT0}.
    # 
    # If desired, we simulate many resamples {S3} from each F2, and compute how often the corresponding {CT1} exceeds the known T1.  This enables a modification of the basic bootstrap conjecture.
    # 
    # Need a generic call protocol for invoking it, or template.
    #@@c
    #@+<<resample>>
    #@+node:gte.20170616083209.1543: *4* <<resample>>
    def resample(self, n_resamples, double_samples=0):
       self.nRe=n_resamples; self.nDo=double_samples
       self.rRe=arange(self.nRe) ; self.rDo=arange(self.nDo)
       T2vector=[] ; T3matrix=[]
       for re in self.rRe:
          boot_obs = self.field.resampled_survey_observations()
          T2vector.append( self(boot_obs) )
          T3vector=[]
          if self.nDo:
            for du in self.rDo:
              reboot_obs = self.field.resampled_survey_observations(boot_obs)
              T3vector.append( self(reboot_obs) )
          T3matrix.append(T3vector)
       self.T2s=array(T2vector) ; self.T3s=array(T3matrix).T
       # double samples in columns, resamples in rows
       # the last axis of T3s and T2s match: CT0_method(T3s, T2s) will work
    #@+at
    #    self.qRe=quantiles(self.nRe)
    #    if self.nDo: self.qDo=quantiles(self.nDo)
    #@-<<resample>>
    #@+<<confidence functions>>
    #@+node:gte.20170616083209.1545: *4* <<confidence functions>>

    def make_confidence_functions(self, CT0_method, bdry='truncate'):
        # other boundary methods still need writing
        CT0s  = sort(CT0_method(self.T2s, self.base))
        # basic bootstrap conjecture
        self.basic_prob_exceed, self.basic_conf_bound = match_quantiles(CT0s, quantiles(self.nRe))
        if self.nDo:
            # double bootstrap conjecture
            CT1s =  CT0_method(self.T3s, self.T2s)
            greater = sum(CT1s>self.base, axis=0)
            # quantiles grow up from the bottom, 'greater' down from the top
            # count how often they overlap:
            frac_CT1_gt_T1 = array([sum(greater+i>=self.nDo) for i in range(self.nDo)])/self.nRe
            self.double_prob_exceed, self.double_conf_bound = match_quantiles(CT0s, frac_CT1_gt_T1)
            self.frac_CT1_gt_T1 = frac_CT1_gt_T1
        
    #@+at
    #     # low prob bound = max(conjectured_probs[0], self.prob_exceed(CT0s[0])  etc
    #     if self.nDo>self.nRe or self.nDo==0:
    #         msg = (CT0s[0], CT0s[-1], self.prob_exceed(CT0s[0]),  self.prob_exceed(CT0s[-1]) )
    #     else:
    #         msg = (self.conf_bound(probs[0]), self.conf_bound(probs[-1]), probs[0], probs[-1])
    #     print 'predictions for values outside (%.3f, %.3f) or probs outside (%.3f, %.3f) are' % msg
    #     print treatment[out[:2]]
    #@-<<confidence functions>>

    #@+others
    #@-others
    #@-<<bootstrap confidence bounds>>

#@+<<collection of discrepancy measures>>
#@+node:gte.20170709082420.1: *3* <<collection of discrepancy measures>>

def difference_discrepancy(T2, T1):
    return T2-T1

def difference_CT0(T2, T1):
    return 2*T1 - T2

def root_discrepancy(T2, T1):
    return sqrt(T2) - sqrt(T1)

def root_CT0(T2, T1):
    return T1 * (2 - sqrt(T2/T1))**2

def ratio_discrepancy(T2, T1):
    return T2/T1

def ratio_CT0(T2, T1):
    return T1*T1/T2

class interval_CT0(object):
    def __init__(self, lower=0, upper=1):
        self.L=lower ; self.U=upper
    def logist(self, x):
        return(x-self.L) / (self.U-x)
    def __call__(self, T2, T1):
        LT1 = self.logist(T1)
        A = LT1*LT1 / self.logist(T2)
        return (self.L + self.U*A) / (1+A)

def studentized_discrepancy(T2, T1):
    # the [1] component is the standard deviation
    return (T2[0] - T1[0]) / T2[1]

def studentized_CT0(T2, T1):
    return T1[0] + (T1[0] - T2[0]) * T1[1] / T2[1]

def other_pctile_CT0(T2, T1):
    return T2

def studentized_ratio_discrepancy(T2, T1):
    x=T2[0] ; s=0.5*sqrt(T2[1])
    a = x + s - sqrt(x*x+s*s)
    # From x-a to x+s-a is an interval of width s with the property that the multiples
    # of the subintervals from x-a to x and from x to x+s-a are the same.
    return (x-T1[0])/T1[0] * (x-a)/x
    return log(x) - log(T1[0]) / (log(x) - log(x-a))
    # need a measure that is proportional to x-a or its log
#@-<<collection of discrepancy measures>>
#@+<<match quantiles>>
#@+node:gte.20170706083737.1: *3* <<match quantiles>>
def match_quantiles(xvec, yvec, bdry='truncate', report_ends=False):
    treatment = {'ex': 'extrapolated', 'tr': '\ntruncated at boundary values', 'in': 'invalid'}
    #TP = TTP if bdry[:2]=='tr' else (ETP if bdry[:2]=='ex')
    if bdry[:2]=='tr': TP=TTP
    if bdry[:2]=='ex': TP=ETP
    if bdry[:2]== 'in': bdry=(nan, nan)
    lx=len(xvec) ; xquan=quantiles(lx)
    ly=len(yvec) ; yquan=quantiles(ly)
    y_of_xvec = TP(yquan, yvec, bounds_error=False)(xquan)
    x_to_y = TP(xvec, y_of_xvec,  bounds_error=False)
    x_of_yvec = TP(xquan, xvec,  bounds_error=False)(yquan)
    y_to_x = TP(yvec, x_of_yvec, bounds_error=False)
    
    # report boundaries and treatment
    if report_ends: 
        if lx>ly:
            bounds= tuple([decode(xvec, encode(xquan, yquan[i])) for i in [0,-1]] + [yvec[0],yvec[-1]])
        else:
            bounds= tuple([xvec[0],xvec[-1]] + [decode(yvec, encode(yquan, xquan[i])) for i in [0,-1]] )
        print 'x values outside (%.3f, %.3f) and y outside (%.3f, %.3f) are' % bounds ,
        print treatment[bdry[:2]]
    return x_to_y, y_to_x

def ETP(x, y, bounds_error):
    return interp1d(x, y, bounds_error=bounds_error, fill_value='extrapolate')

def TTP(x, y, bounds_error):
    return interp1d(x, y, bounds_error=bounds_error, fill_value=(y[0],y[-1]) )
#@+at
# def match_quantiles(xvec, yvec, bdry='extrapolate', report_ends=False):
#     treatment = {'ex': 'extrapolated', 'tr': '\ntruncated at boundary values', 'na': 'invalid'}
#     if bdry[:2]== 'na': bdry=(nan, nan)
#     lx=len(xvec) ; xquan=quantiles(lx)
#     ly=len(yvec) ; yquan=quantiles(ly)
#     y_of_xvec = TP(yquan, yvec, fill_value=bdry, bounds_error=False)(xquan)
#     x_to_y = TP(xvec, y_of_xvec, fill_value=bdry, bounds_error=False)
#     x_of_yvec = TP(xquan, xvec, fill_value=bdry, bounds_error=False)(yquan)
#     y_to_x = TP(yvec, x_of_yvec, fill_value=bdry, bounds_error=False)
#     
#     # report boundaries and treatment
#     if report_ends: 
#         if lx>ly:
#             bounds= tuple([decode(xvec, encode(xquan, yquan[i])) for i in [0,-1]] + [yvec[0],yvec[-1]])
#         else:
#             bounds= tuple([xvec[0],xvec[-1]] + [decode(yvec, encode(yquan, xquan[i])) for i in [0,-1]] )
#         print 'x values outside (%.3f, %.3f) and y outside (%.3f, %.3f) are' % bounds ,
#         print treatment[bdry[:2]]
#     return x_to_y, y_to_x
# 
#@@c
    
def encode(V, x, ends='extrapolate'):
    i = bisect(V, x)
    if i==0:
        e2=ends[:2]
        if e2=='ex': i=1 # extrapolate
        if e2=='tr': return (0, 0.) # truncate
        if e2=='na': return nan # invalid
    if i==len(V):
        e2=ends[:2]
        if e2=='ex': i=i-1 # extrapolate
        if e2=='tr': return (i-2, 1.) # truncate
        if e2=='na': return nan 
    vi=V[i-1]
    return (i-1, (x-vi) / (V[i]-vi) )
    
def decode(V, code):
    if code is nan: return nan
    (i, frac)=code ; vi=V[i]
    return vi + frac*(V[i+1]-vi)

def merge(x1, y1, x2, y2):
    ''' both xys are increasing pairs; so are merged series '''
    x12=[] ; y12=[]  ; i1=0 ; i2=0 ; l1=len(x1) ; l2=len(x2)
    while i1<l1 and i2<l2:
        if x1[i1]<=x2[i2]:
            assert(y1[i1]<=y2[i2])
            x12.append(x1[i1]) ; y12.append(y1[i1])
            i1+=1
        else:
            assert(y1[i1]>y2[i2])
            x12.append(x2[i2]) ; y12.append(y2[i2])
            i2+=1
    return x12, y12

            
                
#@-<<match quantiles>>
#@+<<weighted sum of means>>
#@+node:gte.20170622080520.1: *3* <<weighted sum of means>>
#@+at
# The mean of a cdf at a node is an average of survey observations, weighted by the step heights; hence, a linear combination.  So a weighted sum of means is also a linear combination of survey observations, whose coefficients depend only on the survey design, not its results.  This is convenient for resampling: the call can bypass the details of recomputing the cdfs at every statistic node, having done most of the work just once in essence_of_design.
# 
# If node_weights is None  the nodes are taken to have a weight attribute that can be extracted: for example DOM.N where DOM is a domain.
#@@c
class WeightedSumofMeans(Statistic):

    def essence_of_design(self, node_weights):
       if node_weights is None: node_weights = [n.weight for n in self.nodes]
       self.has_variance=True 
       steps = array([c.steps for c in self.nodes])
       self.survey_weights = inner(transpose(steps), node_weights)
       if self.sdev:  
         self.varmat = zeros([self.length, self.length])
         for n,N in enumerate(self.field.survey):
           st = N.steps 
           var = (diag(st) - outer(st,st)) 
           self.varmat += self.survey_weights[n]**2 * var

    def essence_of_result(self, observed):
        lcstat = dot(self.survey_weights, observed)
        if self.sdev:
           lcvar = self.vmv(observed, self.varmat)
           return (lcstat, sqrt(lcvar))
        # print observed, lcstat
        return lcstat #, lcvar

    @staticmethod
    def vmv(vector,matrix):
       # vector and matrix are ndarrays
       return inner(inner(vector,matrix),vector)
#@-<<weighted sum of means>>
#@+node:gte.20170908103523.1: ** Confidence class
#@+at
# Not sure where this was going or how much it's still needed
# 
#@@c
class Confidence(object):
    def __init__(self, statistic, CT0_method, bdry='truncate'):
        # other boundary methods still need writing
        # ?? where is bdry used?  would be in call to match_quantioles
        CT0s  = sort(CT0_method(statistic.T2s, statistic.base))
        # basic bootstrap conjecture
        self.basic_prob_exceed, self.basic_conf_bound = match_quantiles(CT0s, quantiles(statistic.nRe), bdry=bdry)
        # could do this with calls to TP
        self.nDo=statistic.nDo
        if statistic.nDo:
            # double bootstrap conjecture
            CT1s =  CT0_method(statistic.T3s, statistic.T2s)
            greater = sum(CT1s>statistic.base, axis=0)
            # quantiles grow up from the bottom, 'greater' down from the top
            # count how often they overlap:
            frac_CT1_gt_T1 = array([sum(greater+i>=statistic.nDo) for i in range(statistic.nDo)])/statistic.nRe
            self.double_prob_exceed, self.double_conf_bound = match_quantiles(CT0s, frac_CT1_gt_T1, bdry=bdry)
            self.frac_CT1_gt_T1 = frac_CT1_gt_T1
      
    def list_conf_bounds(self, probs):
        self.basic_conf_bounds_list = self.basic_conf_bound(probs)
        if self.nDo: self.double_conf_bounds_list = self.double_conf_bound(probs)
        
#@-others



#@-leo
