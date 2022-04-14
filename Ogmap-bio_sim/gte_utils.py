#@+leo-ver=5-thin
#@+node:gte.20140305074942.1353: * @file gte_utils.py
from math import *
from pylab import *
from numpy import *

#@+others
#@+node:gte.20140305074942.1354: ** amalgam
def locate(xx, # a list of numbers in increasing order
            y  #  a number
            ):
    '''Returns an integer index of the list, referring to the interval to the left of the number with the same index.
    A list of n elements (indexed from 0 to n-1) defines, and locate() considers, n+1 intervals (indexed from 0 to n) including two semi-infinite ones outside.
    Seems to work for cdf functions but not for interpolation search???
    '''
    left=-1 ; right=len(xx)
    while (right-left)>1:
        mid=(right+left)>>1  # integer part of division by 2
        if y>xx[mid]: left=mid
        else: right=mid
    return right

def order(L):
    ''' Given a list L, return a list giving the rank of each element of L (the position it would have in sorted(L)) , and the inverse of that list: namely the position in the original list of the element that has the rank corresponding to that number in the list (SAY THAT BETTER) '''
    rank=[None]*len(L) ; invrank=[None]*len(L)
    LL = sorted(enumerate(L), key=lambda l: l[1])
    # this can be replaced by a custom function passes an argument which is l[1] by default
    for i,x in enumerate(LL):
        rank[x[0]]=i ; invrank[i]=x[0]
    return rank, invrank

def logspace(a,b,c): return exp(linspace(log(a),log(b),c))

def prob_combine(probvec):
    prod = product(probvec) ;  logprod = log(prod)
    adjust = 1.
    for i in reversed(range(1, len(probvec))):
        adjust = 1. - logprod*adjust/float(i)
    return prod * adjust

def subarray(nd_arr, condition):
    return \
     nd_arr[ [irow for irow,row in enumerate(nd_arr) if condition(row)] ]
    #  only numpy arrays accept a list of indices for subsetting
    # condition() is a true or false statement about a row of the array

def lininterp(x0, y0, x1, y1, x):
    if x0==x1: return 0.5*(y0+y1)
    return (y0*(x1-x) + y1*(x-x0))/(x1-x0)

def sort_on_column(A, n):
    return A[lexsort((A[:,n],))]

def cdf_ordi(n):
    half = 0.5/n
    return linspace(half, 1-half, n)

def plot_cdf(A, srted=True):
    if not srted: A=sorted(A)
    plot(A, cdf_ordi(len(A)) )

def twopen(variable, dtype='float'):
    ''' Get data either from an array or from a file that genfromtxt can read to return an array. '''
    return array(variable) if isinstance(variable, (ndarray, list)) \
           else genfromtxt(variable, dtype=dtype)
#@+node:gte.20140305074942.1355: ** administer test
import os

class TestAdmin(object):
    def __init__(self, reference_file_name):
       #@+<<initialize>>
       #@+node:gte.20140305074942.1358: *3* <<initialize>>
       LABEL_PREFIX = '#['
       reference_file = open(reference_file_name)
       lines = self.ref_lines = reference_file.readlines()
       self.tokens = [ (i, lines[i]) 
               for i in range(len(lines)) if lines[i][:2]==LABEL_PREFIX ]
       #@-<<initialize>>
    def __call__(self, test):
        #@+<<initialize candidate object>>
        #@+node:gte.20140305074942.1359: *3* <<initialize candidate object>>
        self.candidate = []
        #@-<<initialize candidate object>>
        test(self)
        #@+<<complete candidate object>>
        #@+node:gte.20140305074942.1360: *3* <<complete candidate object>>
        #@+at
        # cand_file = open('candidate','w')
        # cand_file.writelines(self.candidate)
        # cand_file.close()
        #@@c
        self.lines_to_file(self.candidate, 'candidate')
        #@-<<complete candidate object>>
        #@+<<create reference object>>
        #@+node:gte.20140228160530.1294: *3* <<create reference object>>
        for (i,t) in enumerate(self.tokens):
            if t[1][2:-2] == self.label:
                # the test function being exercised is responsible for having set self.label
                startline=t[0]
                endline = self.tokens[i+1][0]
                break
        #print startline, endline
        reflines = self.ref_lines[startline+1:endline]
        self.lines_to_file(reflines, 'reference')
        #@+at
        # ref_file = open('reference','w')
        # ref_file.writelines(self.reflines)
        # ref_file.close()
        # 
        #@-<<create reference object>>
        #@+<<compare candidate and reference objects>>
        #@+node:gte.20140305074942.1362: *3* <<compare candidate and reference objects>>
        os.system('diff candidate reference')
        os.system('rm reference')
        #@-<<compare candidate and reference objects>>
    def plus(self, string):
        self.candidate.append(string+'\n')
    def set_label(self, string):
        self.label=string
    def lines_to_file(self, lines, filename):
       Fi = open(filename,'w')
       Fi.writelines(lines)
       Fi.close()

#@-others
#@-leo
