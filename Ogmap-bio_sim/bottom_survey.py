#@+leo-ver=5-thin
#@+node:gte.20170317064727.1505: * @file bottom_survey.py
#@@language python
#@+at
# WORK FLOW
# make_bs_file(setdetfile, columns, outfilename)
# Sort the det details file by year, then put header line at top
# preprocess_data(filename, n_var, tow_area = 0.02355) # makes tiebreak file
# find_median_bandwidth_ogmap(years, variable_name, default_bandwidth=[4., 42., 2.], \
#     SFA=None, NAFO=None, conf_bounds_at=[], \
#     guide_names=['lat','long','depth'], kernel=subcauchy, distance=trawl_distance)
# 
# ***********  In the top directory:
#     
# Pre-process data:
#     square root depth
#     tonnes bzw. thousands per sq km
#     break ties
# 
# Make subdirectories with smaller arenas as desired
# 
# 
# ***********  In a specific subarena directory
#     
# Master of ceremonies:  use_ogmap()
#     bring on the desired acts, in order
#     save the outputs in a self-documenting ascii file
# 
# Get guide and ranvar from data file:
#     specify the guide names in order
#     specify the variable names (could be several)
#     select data for a year (or is this in master of ceremonies?)
# 
# Construct the probability field for the arena:
#     default kernel and distance 
#     default but adjustable bandwidth
#     guide and a column of observed values
#     restrict survey to nodes relevant to the arena
#     redo the probability field
# 
# Adjust bandwidths
# 
# Compute integrated mean value:
#     confidence bounds
#     probability that the true value falls below reference levels
# 
# Compute Gini coefficient:
#     may want a different arena ("non-desert" strata)
#     confidence bounds would be time-consuming
# 
#@@c
from ogmap import *
import pandas as pd
import os
from random import *

#@+others
#@+node:gte.20150128073711.1609: ** Custom distance functions
(KLAT, KLON, KDEP, KSTRAT) = range(4)

class trawl_distance(object):
    def __call__(self, ga, gb):
        #ga = ga.guide ; gb = gb.guide
        return app_SphDist(ga,gb)/self.horscale \
               + abs(ga[KDEP]-gb[KDEP])/self.verscale
    def get_scale(self, bandwidth):
        self.horscale, self.verscale = bandwidth

def app_SphDist(a,b):
    ''' The distance between two points on the surface of the earth given latitude and longitude coordinates in degrees.  An approximate formula good for separations less than a degree [and not too close to the poles?]. '''
    deg_rad = pi/180.
    earthrad = 6367  # km radius at 55 N
    return earthrad * deg_rad * \
           hypot( a[KLAT]-b[KLAT], (a[KLON]-b[KLON])*cos(deg_rad*(a[KLAT]+b[KLAT])/2) )   

# the following two functions are for approximating the stratified random analysis of STRAP

class stratum_distance(Distance):
    def __call__(self, a, b):
        return 0 if a[KSTRAT]==b[KSTRAT] else 1
    def get_scale(self, bandwidth):
        pass

class stratum_kernel(Kernel):
    def get_shape(self, bandwidth):
        return bandwidth
    def __call__(self, distance):
        return 1. if distance<0.5 else 0 
#@+node:gte.20151214084701.1604: ** Preprocess survey data
#@+at
# The 'standard' measure of concentration in input data is numbers or kilograms;
# the `standard' measure of swept area is number of square kilometres.
# 
#@@c
# =============================================================================
# This is where you get the data ready that will be used in Ogmap. This is done in
# the parent directory, before you go into the subdirectory for the specific area
# and species that you will be working on. It is IMPORTANT that the starting headers 
# in your data align with the headers in the domain file (i.e. year, lat, long, 
# depth, NAFO, sfa, stratum, fisharea in the example data sets that I sent).                                                        
# =============================================================================                                  
def preprocess_data(filename, n_var, tow_area=0.02355, lengths=60, prefix='survey'):
    ''' Convert the variables to tonnes or thousands per square kilometre
       Add a column for sqare root depth
       Randomize the order for correctness tests '''
    df = pd.read_csv(filename, delim_whitespace=True)
    try:
        tow_area = 1
     # tow_area = input('tow area [CR for default Campelen 0.02355 sq km]? If standardized to swept area, enter 1 then CR ')
    except:
       tow_area = 0.02355
    tow_scale = 0.001/tow_area # --> thousands or tonnes per sq km
    for i in range(-n_var, 0):
       df.iloc[:,i] *= tow_scale
    df['rootdepth'] = sqrt(df['depth'])
    if n_var<lengths:  # if there are more than this they must be sets of numbers by length
       shu = range(len(df.index)) ; shuffle(shu)
       df_shu = df.iloc[shu]
       #df_shu.to_csv('randomized_survey.tpsk', sep=' ', float_format='%.5g', na_rep='nan', index=False)
       df.to_csv(prefix+'.tpsk', sep=' ', float_format='%.5g', na_rep='nan', index=False)
    else: 
       df.to_csv(prefix+'_length.tpsk', sep=' ', float_format='%.5g', na_rep='nan', index=False)

#This is the current program used for optimization in finding bandwidths
import cma
     
def find_median_bandwidth(years, variable_name, default_bandwidth=[4., 42., 2.], \
    arena='local.arena', \
    guide_names=['lat','long','rootdepth'], kernel=subcauchy, distance=trawl_distance):
    global data, year_data, PF, AR, IMV, RPFs, logfile
    #logfile = open('find_median_bw.bslog','w')
    #logfile.write(variable_name+'\n')
    print 'Starting bandwidths', default_bandwidth
    logfile.write( 'Starting bandwidths [ %.3g %.3g %.3g ]\n' % tuple(default_bandwidth) )
    try:
        survey_data
    except:
        data = assemble_data('survey.tpsk', guide_names, variable_name)
    ####  what happened to assemble_data() ???
    #@+<< interpret years if 'last' or 'all' >>
    #@+node:gte.20150820074125.1580: *3* << interpret years if 'last' or 'all' >>
    if years=='last': years = [int(data[-1,0])]
    elif years=='all': years = range(int(data[0,0]), int(data[-1,0])+1)
    else: assert isinstance(years, list)     
    #@-<< interpret years if 'last' or 'all' >>
    #if years=='last': years = [int(data[-1,0])]
    #elif years=='all': years = range(int(data[0,0]), int(data[-1,0])+1)
    #else: assert isinstance(years, list)
    #prob_field_for_arena(year, variable_name, bandwidth, guide_names, arena=AR, \
    #     kernel=subcauchy, distance=trawl_distance):

    RPFs={} ; bestbands=[]  
    for year in years:
      print year 
      logfile.write('year %d  ' % year)
      bandwidth = default_bandwidth # median crude for all NSRF survey w_fishable
      year_data = G.subarray(data, lambda x: x[0]==year and x[-1]>=0)
      # should be able to do this with bitwise operators &|~ 
      guide = year_data[:,1:] ; observed = year_data[:,-1]
      PF = ProbField(guide, observed, subcauchy, trawl_distance, bandwidth)
      RPF = prob_field_for_arena(PF, AR)
      RPF.raw_pairs = find_nearby_pairs(RPF.survey_array, RPF.distance)
      print 'pairs'
      RPF.cost_function(bandwidth)
      mincost = cma.fmin(RPF.silent_cost_function, bandwidth, sigma0=0.2,
         options = {'scaling_of_variables': [3.503,  37.341, 1.861],
                    'tolfun': .1, 'tolx': 0.01} )
      logfile.write(('best band [ '+3*'%.3g '+'] ') % tuple(mincost[0]))
      logfile.write( '  ;  -log(prob) %.3g\n' % mincost[1])
      bestbands.append(mincost[0])
      print mincost[0], mincost[1]
      logfile.flush() 
      RPFs[year]=RPF  
    beswtbands = array(bestbands)
    medians = median(bestbands,axis=0) 
    print medians
    logfile.write( 'median bandwidths  [ %.3g, %.3g, %.3g ] \n' % tuple(medians) )
    medfile = open(variable_name + '.median_bw','w')
    medfile.write( 'median bandwidths  [ %.3g, %.3g, %.3g ] \n' % tuple(medians) )
    medfile.close()
    return medians

#@+node:gte.20180119184202.1: ** BEST MAKE bwestim a separate function, with randomization applied only here
#@+at
# So bwopt=3 is not an option in interpret_..()
#@+node:gte.20180129082827.1566: ** EVEN NEWER logmap
# =============================================================================
# This section of code finds abundance or weight at length. It requires you to
# run preprocess data again on the data file containing length frequency information.
# =============================================================================
def logmap(n_var, surfile='survey_length.tpsk', domfile='local.domain'):
    DOM = Domain(pd.read_csv(domfile, delim_whitespace=True))
    guidenames=['lat','long','rootdepth']#,'stratum']
    bwlog = open('median.bandwidth') #uses bandwidth file
    bwlines = bwlog.readlines()
    bandwidth = eval(bwlines[-1][18:])    
    results={'PFiDOs':[]}
    DOM.make_node_list(guidenames)
    sur_df = pd.read_csv('../'+surfile, delim_whitespace=True)
    allyears = sorted(unique(sur_df.year))
    guides = sur_df[guidenames].values
    relevance = Relevance(subcauchy(), trawl_distance(), bandwidth)
    # initialize results dictionary
    lengthdist=[]
    for year in allyears:
        dfy = sur_df[sur_df.year==year]
        guides = dfy[guidenames].values
        nos_persk_permm = dfy.values[:, -n_var-1:-1]
        # numbers per sq km per mm length  ;  the last variable is rootdepth
        PF = ProbField(column_stack([arange(len(dfy)), guides]), relevance)
        # this next bit wants to be a function: for finding a relevant field and FID
        for N in DOM.N: PF(N)
        FID = ProbFieldinDomain(PF, DOM)
        rind = FID.relevant_survey_indices(1) #(0.99)
        # option for 0.99?  should be 1.0 if we know the survey is all in the domain
        RPF = ProbField(PF.survey_array[rind], PF.REL) 
        for N in DOM.N: RPF(N)
        PFiDO = ProbFieldinDomain(RPF, DOM)
        results['PFiDOs'].append(PFiDO)
        billions_permm = dot(PFiDO.survey_areas, nos_persk_permm[rind,:]) #/ 1.e6 #uncomment this division to have results as billions
        # numbers per mm length group
        lengthdist.append([year]+list(billions_permm))
        print year
    avg = mean(array(lengthdist)[:,1:], axis=0)
    lengthdist.append([nan]+list(avg))
    savetxt('length.distrib', array(lengthdist), fmt='%.1f')
    return lengthdist, results

# =============================================================================
# The following code is the code that finds the biomass or abundance of certain
# variables for the years and area that you wish. To define area, you will need
# to utilize a different folder containing the local domain of the area you
# are integrating over.
# 
# You can either manually enter the various requirements, or keep them in a 
# file called default.params. Using the file is much easier.
# =============================================================================
def use_ogmap(survey_file_name='survey.tpsk', domain_file_name='local.domain', \
              param_file_name='default.params', afterthoughts=False):
    global RPF, DOM, IMV, PF, logfile, FIDs, bestbands
    print 'ogmap analysis of trawl survey series'
    
    #@+<< assemble running instructions >>
    #@+node:gte.20171214091255.1603: *3* << assemble running instructions >>
    # Option of getting running instructions from a parameter file.
    try:
        param_lines = open(param_file_name).readlines()
        for line in param_lines: exec(line)
    except:
        pass

    # Define file of survey history
    varname = 'biomass'
    #varname = raw_input('Name of variable to analyse: ')
    dirtemp = os.getcwd()
    dirtemp1 = os.path.basename(dirtemp)
    logfile = open(varname + '_' + '_' + dirtemp1 +  '.log', 'w')
    #logfile = open(varname+'dat', 'w')
    logfile.write(varname + '     ' + os.getcwd() + '\n')

    try:
        guide_names
    except:
       try: guide_names = input('Guide names [CR for default]? ')
       except: guide_names = ['lat', 'long', 'rootdepth', 'stratum']

    try:
        years
    except:
        years = input("Years ('last', 'all' or array): ")

    try:
        bwopt
    except:
        print 'Bandwidth options:'
        print '1:  Specify the bandwidths'
        print '2:  Take the bandwidths from the median.bandwidth file'
        print '3:  Create the median.bandwidth file and then use the result'
        print 'If first run on species and area, need to do this, then rerun with option 2'
        bwopt = input('Bandwidth option? ')
    if bwopt==1:
        bw = input('Bandwidth array? ')
    elif bwopt==2:
        bwlog = open('median.bandwidth')
        bwlines = bwlog.readlines()
        bw = eval(bwlines[-1][18:])

    # Take it for granted we always want integrated mean value over the arena.

    try:
        cb
    except:
        cb = raw_input('Confidence bounds ([y]/n)? ')
    if len(cb) == 0 or cb[0] != 'n':  # if not ('n' in cb):
        try:
            cbmethod
        except:
           print 'CB method options:'
           print '1: ratio CT0\n2: difference CT0\n3: other percentile CT0'
           print '4: square root CT0 (default)'
           cbmeth = raw_input('Method option?')
           if cbmeth[:4] == 'diff': cbmethod=difference_CT0
           elif cbmeth[:4] == 'othe': cbmethod=other_pctile_CT0
           elif cbmeth[:4] == 'rati': cbmethod=ratio_CT0
           else: cbmethod=root_CT0
        #mirror = raw_input('Invert resurvey?')
        #if mirror[0]=='y':
        #    mirror = True
        #    cbmethod=other_pctile_CT0
        #else: mirror=False
        try:
            resamples
        except:
            resamples = input('Resamples: ') ; doublesamples=0
        if isinstance(resamples, tuple):
               (resamples, doublesamples) = resamples
        # resample method.  if invert:  cbmethod=other_pctile_CT0
        try:
            ref_probs
        except:
            try:
               ref_probs = input('Reference probs: ')
               if isinstance(ref_probs, (float,int)): ref_probs=[ref_probs]
            except: ref_probs=[]
        try:
            ref_vals
        except:
            try:
               ref_vals = input('Reference points: ')
               if isinstance(ref_vals, (float,int)): ref_vals=[ref_vals]
            except: ref_vals=[]

    try:
        gini
    except:
        gini = raw_input('Gini coefficient (y/[n])? ')
        gini = not (len(gini) == 0 or gini[0] != 'y')
    try:
        sanity
    except:
        sanity = raw_input('Sanity check: (y/[n])? ')
        sanity = not (len(sanity) == 0 or sanity[0] != 'y')
        
# =============================================================================
        #NOT CURRENTLY WORKING!!
#     try:
#         republic
#         do_republic = True
#     except:
#         do_republic = raw_input('Republic plot (y/[n])? ')
#         do_republic = not (len(do_republic) ==0 or do_republic[0] != 'y')
# =============================================================================
        #do_republic = False
    #@-<< assemble running instructions >>
    #@+<< assemble data for guide and observed >>
    #@+node:gte.20160201094014.1881: *3* << assemble data for guide and observed >>
    # try to allow for a true csv file as well
    try:
        df = pd.read_csv(survey_file_name, delim_whitespace=True)
    except:
        df = pd.read_csv('../'+survey_file_name, delim_whitespace=True)
    guide_names=['lat','long','rootdepth']    
    cols = ['year'] + guide_names + [varname]
    survey_data = df[cols].values
    # don't seem to need the .loc or iloc here

    allyears = sorted(unique(survey_data[:,0]))
    if isinstance(years, int): years = [years]
    elif years=='last': years = [allyears[-1]]
    elif years=='all': years = allyears
    else: assert isinstance(years, list)
    print years
       
    #@-<< assemble data for guide and observed >>
    #@+<< assemble domain >>
    #@+node:gte.20171214091255.1604: *3* << assemble domain >>
    DOM = Domain(pd.read_csv(domain_file_name, sep=' '))
    DOM.make_node_list(guide_names)
    #@-<< assemble domain >>
# =============================================================================
#     if do_republic:
#         republic_dict = {'year':[]}
#         for st in unique(DOM.F.stratum): republic_dict[st]=[]
#         republic_DF = pd.DataFrame(republic_dict)
# =============================================================================
    
    if bwopt==3: 
        bandwidth=bw =[3,20,1]
        bestbands=[]
    print 'Using bandwidths ', bw
    logfile.write( 'Using bandwidths  [ %.3g, %.3g, %.3g ] \n' % tuple(bw) )
    logfile.write('Year Estimate Mass50 LowCIprob UpCIprob LowCIval UpCIval LRP USR ProbLRP ProbUSR Gini\n')
    
    relevance = Relevance(subcauchy(), trawl_distance(), bw)
    FIDs={} ; RPFs={}

    for year in years:
        year_data = survey_data[ (survey_data[:,0]==year) & (survey_data[:,-1]>=0) ]
        guide = year_data[:,1:-1] ; observed = year_data[:,-1]/1000.
        # change units to millions bzw. kilotonnes per sq km
        PF = ProbField(column_stack([observed, guide]), relevance)        

        ###########
        # temporary for testing
        ################
        #RPF = relevant_prob_field(DOM, PF)
        #IMV = WeightedSumofMeans(RPF, DOM.N)
        #print IMV.base,
        for N in DOM.N: PF(N)
        FID = ProbFieldinDomain(PF, DOM)
        #PFiDO = ProbFieldinDomain(RPF, DOM)
        rind = FID.relevant_survey_indices(0.99)
        # option?  should be 1 if we know the survey is all in the domain
        RPF = ProbField(PF.survey_array[rind], PF.REL)
        #RIMV = WeightedSumofMeans(RRPF, DOM.N)
        #print RIMV.base
        #print PF.length, RPF.length, RRPF.length
        #continue
        ###############
        
        if bwopt==3: # only estimate bandwidths (and report median?)
          RPF.raw_pairs = find_nearby_pairs(RPF.survey_array, RPF.REL.distance)
          RPF.cost_function(bandwidth)
          mincost = cma.fmin(RPF.silent_cost_function, bandwidth, sigma0=0.2,
             options = {#'scaling_of_variables': bandwidth,
                        'tolfun': .1, 'tolx': 0.01} )
          logfile.write(('best band [ '+3*'%.3g '+'] ') % tuple(mincost[0]))
          logfile.write( '  ;  -log(prob) %.3g\n' % mincost[1])
          bestbands.append(mincost[0])
          print mincost[0], mincost[1]
          logfile.flush() 
          RPFs[year]=RPF 
          continue # ignore all of the other possible options 

        IMV = WeightedSumofMeans(RPF, DOM.N)
        print '\nYEAR %d  thousands of tonnes %.3g' % (year, IMV.base)
        # the next 5 lines are a Gini substitute
        FID = ProbFieldinDomain(RPF, DOM)
        #FID.make_mean_concs()
        area_frac = 100*FID.area_frac_for_mass_frac()
        FIDs[year] = FID
        print '%.1f percent of the area had 50 percent of the mass' % area_frac
        if sanity: # True or False
            mean_conc = mean(RPF.observed)
            total_area = sum(DOM.F.area)
            print 'sanity check area: %.3g  mean conc: %.3g  mass: %.3g' % \
                (total_area, mean_conc, total_area*mean_conc )
        logfile.write('%d %.3g ' % (year, IMV.base) )
        logfile.write('%.1f ' % area_frac)
        if len(cb) == 0 or cb[0] != 'n':  # if not ('n' in cb):
           IMV.resample(resamples, doublesamples)       
           #IMV.T2s.to_csv(r"C:\Shrimp\Ogmap\SFA4 PB NoClosed\testT2sA.csv", index=None, header=True)
           try: os.mkdir('CIResamples')
           except: pass
           yrtmp=int(year)
           outname=(r'CIResamples/'+varname+'_'+'resamples'+'_'+str(yrtmp)+'.csv')
           #outname2=(r'CIResamples/'+varname+'_'+'resamples2'+'_'+str(yrtmp)+'.csv')
           savetxt(outname, IMV.T2s, delimiter=',')
           #savetxt(outname2, IMV.T3s, delimiter=',')
           IMV.make_confidence_functions(cbmethod) 
 
        if cb == False:
           logfile.write('na na ')
        try:
            ref_probs
            p2 = list(IMV.basic_conf_bound(ref_probs))
            print 'Ref probs:', ref_probs, '  Values: ',
            print len(p2)*'%.3g ' % tuple(p2)
            #logfile.write('Ref probs [ ')
            for rp in ref_probs:
                logfile.write('%.3g ' % rp)
            #logfile.write(']   Values [ ')
            for rp in p2:
                logfile.write('%.3g ' % rp)
            #logfile.write(']\n')
        except: 
            logfile.write('na na ')
            pass
        
# =============================================================================
#         if do_republic:
#             PFiDO.make_republic_array('stratum')
#             republic_row = list(PFiDO.republic_array[:,1])+[year]
#             republic_DF.loc[len(republic_DF)] = republic_row
# =============================================================================
                            
        try:
            ref_vals
            pl = list(IMV.basic_prob_exceed(ref_vals))
            print 'Ref values: ', ref_vals,   'Probs of exceeding true', 
            print len(pl)*'%.3g ' % tuple(pl)
            if ref_vals<>[]:
                logfile.write((len(ref_vals)*'%.3g ') % tuple(ref_vals))
            #logfile.write('Ref values [ ')
            #logfile.write((len(ref_vals)*'%.3g ') % tuple(ref_vals) )
            #for rp in ref_vals:
            #    logfile.write('%.3g ' % rp)
            #logfile.write(']   Probs [ ')
            for rp in pl:
                logfile.write('%.3g ' % rp)
            #logfile.write(']\n')
            if ref_vals==[]:
                logfile.write('na na na na ')
        except: pass    
        if gini:
            gc = GiniCoefficient(RPF, DOM.N, DOM.F.area)()
            print 'Gini coeff %.3f' % gc
            logfile.write('%.3f\n' % gc)
        if gini == False:
                logfile.write('na \n')
        logfile.flush()
    #outside the years loop
    if bwopt==3: #write median bandwidth file
        bestbands = array(bestbands)
        medians = median(bestbands,axis=0) 
        print medians
        medfile = open('median.bandwidth','w')
        medfile.write( 'median bandwidths  [ %.3g, %.3g, %.3g ] \n' % tuple(medians) )
        medfile.close()
# =============================================================================
#     if do_republic:
#        republic_DF.to_csv(varname+'_republic.csv', float_format='%.2f', index=False)
# =============================================================================
    logfile.close()   

#@verbatim
#@+node:gte.20150723064139.1592: ** Gini coefficient
class GiniCoefficient(Statistic):
    ''' The degree of aggregation of the field of mean values of the probabilitty field '''
    # call GiniCoefficient(field, arena.nodes, arena.varea)
    def essence_of_design(self, varea):
        self.areas = varea
    
    def essence_of_result(self, observed):
        concs = [inner(N.steps, observed) for N in self.nodes]
        #concentrations = [field.meva(N, retval=True)[0] for N in self.nodes]
        conar = zip(concs, self.areas)
        conar.sort() ; conar=array(conar)
        return abstract_gini(conar[:,0], conar[:,1])

def abstract_gini(sorted_concentration_array, area_array):
    c=sorted_concentration_array ; a=area_array ; m=c*a # mass
    lorenz_integral = inner(a, cumsum(m) - m/2.) / sum(m) / sum(a)
    return 1-2*lorenz_integral

# write a Brown function based on inner(a[1:]-a[:-1],m[1:]+m[:-1])

def test_gini():
    pass
    # set the concentrations somehow, maybe by finagling what the cdfs have to be; maybe by breaking things up so that essence_of_result() can be fed with fake cdfs for concentration instead of computed ones.     
        
def lorenz_curve(sorted_concentrations, areas):
    c=array(sorted_concentrations) ; a=array(areas)
    m = a*c
    plot([0]+list(cumsum(a)),[0]+list(cumsum(m)))
    lorenz_integral = inner(a, cumsum(m) - m/2.) 
    gini = 1 - 2 * lorenz_integral / sum(m) / sum(a)
    print lorenz_integral, gini

def make_lorenz(a):
    def y(x): return  a*(a+1)/(x+a) -a
    absci = linspace(0,1,101)
    plot(absci, y(1-absci))

#@verbatim
#@+node:gte.20160808084718.1670: *3* exercise Gini
def quick_gini(datafile, years):
    AR=local_Arena()
    guide_names=['lat', 'long', 'rootdepth']
    varname='w_totals'
    bw=[ 3.46, 42.4, 2.37 ] 
    df=pd.read_csv(datafile, delim_whitespace=True)
    cols = ['year'] + guide_names + [varname]
    data = df[cols].values
    
    for year in years:
        year_data = data[ (data[:,0]==year) & (data[:,-1]>=0) ]
        guide = year_data[:,1:-1] ; observed = year_data[:,-1]/1000.
        # change units to millions bzw. kilotonnes per sq km
        PF = ProbField(guide, observed, subcauchy, trawl_distance, bw)
        # would need to offer option here for stratum_distance and _kernel
        RPF = prob_field_for_arena(PF, AR)
        gc = GiniCoefficient(RPF, AR.nodes, AR.varea)
        gc.resample(500,0)  
        gc.find_fraction_exceeding_T1(root_CT0)     
        print '%d %.3f [ %.3f %.3f %.3f ]' % tuple([year, gc.base] + \
           gc.statistic_value([.1,.5,.9] ) )
        gc.find_fraction_exceeding_T1(interval_CT0())     
        print '%d %.3f [ %.3f %.3f %.3f ]' % tuple([year, gc.base] + \
           gc.statistic_value([.1,.5,.9] ) )


#@+node:gte.20150817080109.1580: ** ProbField for Arena
#@+at
#@+at
#@+at
#@+at
# A limit area of 3 sq km is appropriate to Campelen survey density in 2J3KL (what about 3Ps?)
# but maybe not for small shellfish beds within 3Ps.
#@+node:gte.20160727081456.1669: ** Gini_size
def gini_size(year):
    global RPF, AR, Awt, Anos, Asize
    try:
       df = pd.read_csv('tiebreak', delim_whitespace=True)
    except:
       df = pd.read_csv('../tiebreak', delim_whitespace=True)
    cols = ['year', 'lat', 'long', 'rootdepth', 'w_totals', 'n_totals']
    data = df[cols][df.year==year].values
    
    bw = [4, 43, 2.33]
    AR = local_Arena()
    ref_probs = [.05, .1, .5, .9, .95]
    
    guide = data[:,1:] ; weight = data[:,-2]/1000.
    # change units to millions bzw. kilotonnes per sq km
    PF = ProbField(guide, weight, subcauchy, trawl_distance, bw)
    RPF = prob_field_for_arena(PF, AR)
    IMV = IntegratedMeanValue(RPF, AR)
    print '\nYEAR %d  thousands of tonnes %.3g' % (year, IMV.base)
    IMV.resample(599,0)
    IMV.find_fraction_exceeding_T1(root_CT0)
    print array(IMV.statistic_value(ref_probs))
    gc = GiniCoefficient(RPF, AR.nodes, AR.varea)()
    print 'Gini coeff %.3f' % gc
    
    num_obs=RPF.survey_array[:,-2]/1.e6 # to give avg size in grams 
    Awt=[]
    Anos=[]
    Asize=[]
    for NA in AR.nodes:
        RPF(NA)
        mass = inner(NA.steps, RPF.observed)
        number = inner(NA.steps, num_obs)
        Awt.append(mass)
        Anos.append(number)
        Asize.append(mass/number)

#@+node:gte.20160408081107.1638: ** STRAP compare
def strap_compare(param_file_name='default.params'):
    global PF, AR, IMV, STMV, logfile, data, df, STPF
    print 'ogmap analysis of trawl survey series'
    
    AR = local_Arena()
    
    # Option of getting running instructions from a parameter file.
    try:
        param_lines = open(param_file_name).readlines()
        for line in param_lines: exec(line)
    except:
        pass
    
    # Define file of survey history
    varname = raw_input('Name of variable to analyse: ')
    logfile = open(varname + '.log', 'w')
    logfile.write(varname+'\n')

    try:
        guide_names
    except:
       try: guide_names = input('Guide names [CR for default]? ')
       except: guide_names = ['lat', 'long', 'rootdepth']
 
    years = input("Years ('last', 'all' or array): ")
   
    try:
        bwopt
    except:
        print 'Bandwidth options:'
        print '1:  Specify the bandwidths'
        print '2:  Take the bandwidths from the .median_bw file'
        print '3:  Create the .median_bw file and then use the result'
        bwopt = input('Bandwidth option? ')
    if bwopt==1:
        bw = input('Bandwidth array? ')
    elif bwopt==2:
        bwlog = open(varname + '.median_bw')
        bwlines = bwlog.readlines()
        bw = eval(bwlines[-1][18:])

    # Take it for granted we always want integrated mean value over the arena.


    ############
    #  end of interactive input; processing starts next
    ############
    
    #@+<< assemble data for guide and observed >>
    #@+node:gte.20160408081107.1639: *3* << assemble data for guide and observed >>
    try:
        df = pd.read_csv('tiebreak', delim_whitespace=True)
    except:
        df = pd.read_csv('../tiebreak', delim_whitespace=True)
        
    cols = ['year'] + guide_names + [varname]
    data = df[cols].values

    if years=='last': years = [int(data[-1,0])]
    elif years=='all': years = range(int(data[0,0]), int(data[-1,0])+1)
    # data must be sorted by year for this to work
    else: assert isinstance(years, list)     
    #@-<< assemble data for guide and observed >>
    if bwopt==3: bw = find_median_bandwidth(years, varname)
    print 'Using bandwidths ', bw
    logfile.write( 'Using bandwidths  [ %.3g, %.3g, %.3g ] \n' % tuple(bw) )
    
    for year in years:
        year_data = data[ (data[:,0]==year) & (data[:,-1]>=0) ]
        guide = year_data[:,1:-1] ; observed = year_data[:,-1]/1000.
        # change units to millions bzw. kilotonnes per sq km
        PF = ProbField(guide, observed, subcauchy, trawl_distance, bw)
        # would need to offer option here for stratum_distance and _kernel
        #RPF = prob_field_for_arena(PF, AR)
        IMV = IntegratedMeanValue(PF, AR)
        print '\nYEAR %d  thousands of tonnes %.3g' % (year, IMV.base)
        logfile.write('%d  %.0f\n' % (year, IMV.base) )
        STF = ProbField(guide, observed, stratum_kernel, stratum_distance, [1,1,1])
        #STRPF = prob_field_for_arena(STF, AR)
        STMV = IntegratedMeanValue(STF, AR)
        scatter(IMV.survey_weights - STMV.survey_weights, PF.observed)
        logfile.write('\n')
        logfile.flush()

# =============================================================================
# EDIT this part of code to align with the headers in your data and arena file
# For example, you might have a header CMA instead of SFA. You can have as many
# columns as you wish in your domain file; just edit the original to assign
# different areas based on location (for example, I have the header fisharea).
# =============================================================================
def make_arena_directory(full_arename, sub_arename, \
        NAFO=None, SFA=None, strata=None, fisharea=None, omitstrata=None):
    ''' Run in the 'top' directory for all areas.
        The survey history file and full lookup arena file will probably live here.
        Create subarenas corresponding to, e.g. management areas.
    '''
    global LNshelf, arena
    try:
        os.mkdir(sub_arename)
    except:
        print 'arena exists'
    LNshelf = pd.read_csv(full_arename, sep=' ')
    if NAFO:
        if not isinstance(NAFO, list): NAFO=[NAFO] 
        arena = LNshelf[LNshelf.NAFO.isin(NAFO)]
    if SFA: arena = LNshelf[LNshelf.SFA == SFA]
    if strata: arena = LNshelf[LNshelf.stratum.isin(strata)]
    if fisharea: arena = LNshelf[LNshelf.fisharea == fisharea]    
    if omitstrata:
        arena = arena[~ arena.stratum.isin(omitstrata)]
    arena.to_csv(sub_arename+'/local.domain', sep=' ', float_format='%.4f',index=False)

#@+node:gte.20150129084044.1612: ** Arena
class local_Arena(object):  # should be Arena (Domain) !!!
    def __init__(self, arefile='local.arena'):
        arena = pd.read_csv(arefile, sep=' ')
        table = arena[['lat','long','rootdepth','stratum','area']].values #        table = arena[['lat','long','rootdepth','stratum','area']].values
        self.vertices = table[:,:-1]
        self.nodes = [Node(v) for v in self.vertices]
        self.varea = table[:,-1]
        self.nVer = len(self.varea)
        self.totarea = sum(self.varea)

def logmap_plot(filename, years, variable):
    global size, value
    results=open(filename)  # 'logmap.log')
    reslines = [l.split() for l in results]
    lv=len(variable)
    value={}
    for y in years:
        rightyear=False
        value[y]=zeros(81)
        for rl in reslines:
            if rl[0]=='year' and int(rl[1])==y:
                rightyear=True
                continue
            if rightyear and (variable in rl[0]):
                value[y][int(rl[0][lv:])] = float(rl[1])
            if rightyear and (rl[0]=='year'):
                break
        plot(1+arange(80), value[y][1:])
        title(str(y) + '  ' + variable)

#@+node:gte.20150420085109.1581: ** Read set details card
global setdet_dict
#@+node:gte.20150420085109.1609: *3* identifiers

def year(sd):
    # nominal year of survey
    year = int(sd[9:11]) ; month = int(sd[11:13])
    if month<3: year -=1
    year += 1900 if year>60 else 2000
    return year
#@+node:gte.20150420085109.1607: *3* typical guides
def latitude(sd):
    return float(sd[65:67]) + float(sd[67:70])/600.

def longitude(sd):
    return -float(sd[70:72]) - float(sd[72:75])/600.

def depth(sd):
    dep = sd[43:47]
    try: x=float(dep)
    except: x=0.
    return x

def rootdepth(sd):
    return sqrt(depth(sd))

def endlat(sd):
    return float(sd[97:99]) + float(sd[99:102])/600.

def endlong(sd):
    return -float(sd[102:104]) - float(sd[104:107])/600.

#@+node:gte.20150420085109.1610: *3* auxiliary guides
def temperature(sd):
    T = float(sd[63:65])*0.1
    return -T if sd[62]=='9' else T+float(sd[62])*10.

#@+node:gte.20150420120709.1586: *3* observations
#@+node:gte.20150420085109.1611: *3* dictionary of variables
#@+at
# The value of each dictionary entry is a function which in general will have been defined in an earlier node.# identifiers
#@@c
setdet_dict = {
   'unique': (lambda sd: sd[1:9]), # vessel-trip-set as a string
   'year': year,
   'NAFO': (lambda sd: sd[20:22]),
   
   # typical guides
   'lat': latitude,
   'long': longitude,
   'rootdepth': rootdepth,
   'depth': depth,
   'stratum': (lambda sd: int(sd[17:20]) ),
   
   # auxiliary guides
   'endlat': endlat,
   'endlong': endlong,
   'temp': temperature,
   'light': (lambda sd: int(sd[25:28]) ),
   'bottype': (lambda sd: sd[31]),
   
   # observations
   'number': (lambda sd: float(sd[84:90]) ),
   'weight': (lambda sd: float(sd[90:97])*0.01 )
   
}
#@+node:gte.20150421082927.1587: *3* make file for bottom_survey
def make_bottom_survey_file(setdetfile, columns, outfilename):
    outfile = open(outfilename,'w')
    outfile.write('#')
    for c in columns: outfile.write(c+' ')
    outfile.write('\n')
    setdet = open(setdetfile)
    for sd in setdet:
        for c in columns:
            outfile.write(str(setdet_dict[c](sd))+' ')
        outfile.write('\n') ; outfile.flush()
    outfile.close()
#@+node:gte.20150722095535.1578: ** map average size
def map_avg_size(years, variable_name, bandwidth, SFA=None, NAFO=None, \
    guide_names=['lat','long','depth'], kernel=subcauchy, distance=trawl_distance):
    alldata = genfromtxt('processed')
    labeldict = headers('processed')
    data = alldata[:, [labeldict[g] for g in ['year'] + guide_names  + \
       ['n_'+variable_name, 'w_'+variable_name] ] ]
    AR = local_Arena()
    RPFs={} 
    if years=='last': years = [int(data[-1,0])]
    elif years=='all': years = range(int(data[0,0]), int(data[-1,0])+1)
    else: assert isinstance(years, list)     
    for year in years:
      print year 
      #logfile.write('year %d\n' % year)
      year_data = G.subarray(data, lambda x: x[0]==year and not isnan(x[-1]) and x[-2]>0) 
      guide = year_data[:,1:-1] ; observed = year_data[:,-1] / year_data[:,-2]
      PF = ProbField(guide, observed, kernel, distance, bandwidth)
      relevant_indices = relevant_survey(AR, PF)
      RPF = ProbField(PF.survey_array[relevant_indices,:-1], \
            PF.survey_array[relevant_indices,-1], kernel, distance, bandwidth)
      RPFs[year]=RPF
    return RPFs 
 
def show_mean_map(PF, AR, empty_dots, extremes):
    z=extremes[:]
    for v in AR.nodes:
        PF(v)
        z.append(PF.meva(v, retval=True)[0])
    x = hstack( [empty_dots[0], AR.vertices[:,1] ]) 
    y = hstack( [empty_dots[1], AR.vertices[:,0] ])
    print max(z[2:]), min(z[2:]), extremes
    scatter(x,y,c=z,lw=0)

#@+node:gte.20150205090301.1622: ** utility functions
def explore():  # run in NSRF_B
    alldata=genfromtxt('tiebreak')
    all_13=G.subarray(alldata,lambda x: x[0]==2013)
    bdry=loadtxt('../sfa4_boundary.bln',skiprows=1)
    subplot(2,1,1)
    scatter(all_13[:,2],all_13[:,1], s=4*all_13[:,-3],hold=0,lw=0)
    plot(bdry[:,0],bdry[:,1],':k')
    text(-66,66,'2013')
    subplot(2,1,2)
    all_14=G.subarray(alldata,lambda x: x[0]==2014)
    scatter(all_14[:,2],all_14[:,1], s=4*all_14[:,-3],hold=0,lw=0)
    plot(bdry[:,0],bdry[:,1],':k')
    text(-66,66,'2014')

#@+at
# Use the methods defined in bottom_survey.py
#@@c
def explore_imv():
    global data, year_data, PF, AR, IMV, DMV
    AR = get_SFA_arena(4)
    data = assemble_data('tiebreak', ['lat','long','depth'], 'w_fishable')
    for year in [2013, 2014]:
      print year
      year_data = G.subarray(data, lambda x: x[0]==year and not isnan(x[-1])) 
      guide = year_data[:,1:-1] ; observed = year_data[:,-1]
      PF = ProbField(guide, observed, subcauchy, trawl_distance, [3.3, 40, 2])
      IMV = IntegratedMeanValue(PF, AR)
      print 'subcau', round(IMV.base/1000)
      DF = ProbField(guide, observed, logist, trawl_distance, [0.844, 36.1, .757])
      DMV = IntegratedMeanValue(DF, AR)
      print 'logist', round(DMV.base/1000)
      
#@+node:gte.20150130072611.1617: ** Find nearby pairs

def find_nearby_pairs(guide, distance, extra_pairs=0, temperature=1., chatty=False):
    # temperature scaled to bandwidth
    ncity = shape(guide)[0]
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
        g1 = randint(0,indices-1) ; g2 = randint(0,indices-1) # guides
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

def plotpairs(pair, nodes, indices, ncity):
    # figure()
    for g in range(ncity):
        t = pair[g]
        if t>g:
            if t>=ncity:
               #pass 
               scatter(nodes[g][1],nodes[g][0])
            else:
              plot([nodes[g][1], nodes[t][1]],[nodes[g][0], nodes[t][0] ])

#@+node:gte.20150723064139.1592: ** Gini coefficient
#@+at
#@+at
# lorenz_curve([2/3,4/3],[1/2,1/2])  # scaled to sum to 1
# 0.416666666667 0.166666666667
# 
# conar=zip(Awt,AR.varea) # mean weight, area it applies to
# conar.sort()
# figure()
# conar=array(conar)
# lorenz_curve(conar[:,0],conar[:,1])
# 
# 
#@+node:gte.20160808084718.1670: *3* exercise Gini
def quick_gini(datafile, years):
    AR=local_Arena()
    guide_names=['lat', 'long', 'rootdepth']
    varname='w_totals'
    bw=[ 3.46, 42.4, 2.37 ] 
    df=pd.read_csv(datafile, delim_whitespace=True)
    cols = ['year'] + guide_names + [varname]
    data = df[cols].values
    
    for year in years:
        year_data = data[ (data[:,0]==year) & (data[:,-1]>=0) ]
        guide = year_data[:,1:-1] ; observed = year_data[:,-1]/1000.
        # change units to millions bzw. kilotonnes per sq km
        PF = ProbField(guide, observed, subcauchy, trawl_distance, bw)
        # would need to offer option here for stratum_distance and _kernel
        RPF = prob_field_for_arena(PF, AR)
        gc = GiniCoefficient(RPF, AR.nodes, AR.varea)
        gc.resample(500,0)  
        gc.find_fraction_exceeding_T1(root_CT0)     
        print '%d %.3f [ %.3f %.3f %.3f ]' % tuple([year, gc.base] + \
           gc.statistic_value([.1,.5,.9] ) )
        gc.find_fraction_exceeding_T1(interval_CT0())     
        print '%d %.3f [ %.3f %.3f %.3f ]' % tuple([year, gc.base] + \
           gc.statistic_value([.1,.5,.9] ) )


  
#@+node:gte.20150227103526.1585: ** Year effects to normalize by

def make_promise():
    global promise
    # year = loadtxt('year')
    shrimp = genfromtxt('processed', usecols=[0,11,12]) # magic numbers for total nos & wt
    masspromise = zeros(len(shrimp))
    norms=[]
    for y in range(1996, 2015):
        yearshrimp = G.subarray(shrimp, lambda x: x[0]==y and not isnan(x[2]) )
        #nonzero = [(i,s) for (i,s) in enumerate(totmass) if year[i]==y and s>0]
        nonzero = G.subarray(yearshrimp, lambda x : x[2]>0)
        fracnonz = len(nonzero)/float(len(yearshrimp))
        mednonzero = median([x[2] for x in nonzero])
        meanall = mean(yearshrimp[:,2])
        mediall = median(yearshrimp[:,2])
        mnbyfracpos = mednonzero / fracnonz
        mnofraczero = mednonzero * (1-fracnonz)
        print '%d %d %.3f %.3f %.3f %.3f %.3f %.3f' % ( y, len(nonzero), fracnonz, mednonzero, meanall, mediall, mnbyfracpos, mnofraczero)
        norms.append([y, len(nonzero), fracnonz, mednonzero,meanall, mediall, mnbyfracpos, mnofraczero])
        # for x in nonzero: promise[x[0]] = x[1] / mednonzero
    #savetxt('promise', promise, fmt='%.3e')
    norms=array(norms)
    figure()
    plot(log(norms[:,3]))
    plot(log(norms[:,4]))
    plot(log(norms[:,5]))
    plot(log(norms[:,6]))
    plot(log(norms[:,7]))
    
    
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

#@+node:gte.20151231073302.1602: ** extras
#@+node:gte.20151231073302.1603: *3* effective areas of survey nodes
#@+at
# Plot the effective areas of survey nodes, adjusted to the average area,
# to demonstrate where the over- and under-influential sets are.
# excess = IMV.survey_weights - AR.totarea/RPF.length
# scatter(RPF.survey_array[:,1],RPF.survey_array[:,0],lw=0,color='r',s=-excess/5,hold=0)
# scatter(RPF.survey_array[:,1],RPF.survey_array[:,0],lw=0,color='g',s=excess/5)
#@-others
#@-leo
