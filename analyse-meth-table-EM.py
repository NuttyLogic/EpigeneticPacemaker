#!/usr/bin/python

#!/Library/Frameworks/Python.framework/Versions/Current/bin/python

#!/opt/python-2.5/bin/python



"""
This script tests fmin_slsqp using Example 14.4 from Numerical Methods for
Engineers by Steven Chapra and Raymond Canale.  This example maximizes the
function f(x) = 2*x0*x1 + 2*x0 - x0**2 - 2*x1**2, which has a maximum
at x0=2, x1=1.
"""

import string, sys, re, copy, os, time
import scipy.stats as stats
from scipy.optimize import fmin_slsqp
import numpy as np
## from numpy import random *
from  time import  *
from  tree_n import  *
from random import *
from methyl import *
from numpy import linalg as LA
from pace_maker_n import *

def times_n(line, n):
    line_n = []
    for x in line:
        line_n += [x] *n 
    print line_n
    return line_n
    exit(8)
       
#-------------------------------------------------------------------------------------#
def main (argv):

    fname = "meth-tab-hmn-n656-top-5000-cov-sites.txt"
    fname = "meth-tab-hmn-n656-top-1000-cov-sites.txt"
    fname = "meth-tab-hmn-n656-top-1000-var-sites.txt"
#    fname = "meth-tab-n50-top-50-var-sites-scl-100.txt"
#    fname = "meth-Horv-tab-n16-top-50-var-sites-scl-100.txt"
 #   fname = "meth-Horv-tab-n16-top-500-var-sites-scl-1.txt"
#    fname = "meth-tab-hmn-n500-top-1000-var-sites-scl-1.txt"
#    fname = "meth-tab-hmn-n200-top-500-var-sites-scl-1.txt"
#    fname = "meth-tab-hmn-n300-top-300-var-sites-scl-100.txt"
#    fname = "meth-tab-hmn-n300-top-300-var-sites-scl-1.txt"
#    fname = "meth-tab-hmn-n656-top-1000-abs-pcc-sites.txt"
#    fname = "meth-tab-hmn-n656-top-100-pcc-sites.txt"
#    fname = "meth-tab-kids-n78-top-1000-pcc-sites.txt"
#    fname = "meth-tab-kids-n78-top-100-pcc-sites.txt"
#    fname = "meth-tab-kids-n78-top-1000-abs-pcc-sites.txt"
    fname = "meth-tab-hmn-n656-top-1000-abs-pcc-sites.txt"
#    fname = "meth-tab-kids-n78-top-1000-abs-pcc-years.csv"
#    fname = "meth-tab-kids-n78-top-1000-var-sites.txt"
#    fname = "meth-tab-kids-n78-top-1000-cov-sites.txt"
#    fname = "meth-tab-kids-n78-top-1000-pcc-sites.txt"
#    fname = "meth-tab-kids-n78-top-1000-rnd-sites.txt"
    fname = "meth-tab-dog-n108-top-1000-abs-pcc-sites.txt"
    fname = "meth-tab-dog-n108-top-1000-cov-sites.txt"
    fname = "meth-tab-dog-n108-top-1000-var-sites.txt"
    fname = "meth-tab-GSE42861-n689-top-1000-abs-pcc-sites.txt"
    fname = "meth-tab-GSE87571-n732-top-1000-abs-pcc-sites.txt"
    fname = "meth-tab-GSE87571-n366-top-1000-abs-pcc-sites.txt"
    fname = "meth-tab-GSE87571_2of2-n366-top-1000-abs-pcc-sites.txt"
#    fname = "meth-tab-GSE78874-n259-top-5000-abs-pcc-sites.txt"
#    fname = "meth-tab-GSE78874-n259-top-3000-abs-pcc-sites.txt"
    fname = "meth-tab-GSE74193-n675-top-1000-abs-pcc-sites.txt"
    table, MC_times, sp_list, site_list = read_table(fname)
#    exit(8)

    n_s = len(table)
    n_t  = len(table[0])

#    sp_list_n = times_n(sp_list, 4)
#    MC_times_n = times_n(MC_times, 4)
#    table_n = []
#    for line in table:
#        table_n.append(times_n(line, 4))
    
#    r_rates_n, r_d_n, r_e_n = find_soln_n(A, y, MC_times, n_t, n_s)
#    print "\n\nreturned Rates:", "\n".join(["%f\t%f" % (r_rates_n[i], rates[i]) for i in range(len(rates))])
#    exit(8)
#    print "\n\nIntercept:", intercept
    start_pts = 1
#    n_s, n_t, A, y = preprocess (table)

#    MC_err, MC_rates, MC_d, PM_err, PM_times, PM_rates, PM_d, mhtn_itr  = mhtn_HC_rnd_pts(A, y, MC_times, n_t, n_s, table, start_pts)
    import datetime
    now = datetime.datetime.now()
    start_tm = "%s-%s-%s_%s:%s:%s" %(now.year, now.month, now.day, now.hour, now.minute, now.second)
    print "\n\n\nRun started at %s" % start_tm
#    MC_err, MC_rates, PM_err, PM_times, PM_rates, mhtn_itr  = mhtn_HC1( MC_times, n_t, n_s, table, start_pts)
#    MC_err, MC_rates, MC_d, PM_err, PM_times, PM_rates, PM_d = mhtn_correct_once(MC_times, table)
    MC_err, MC_rates, MC_d, PM_err, PM_times, PM_rates, PM_d = EM_MC_strt_pt(MC_times, table)
    now = datetime.datetime.now()
    end_tm = "%s-%s-%s_%s:%s:%s" %(now.year, now.month, now.day, now.hour, now.minute, now.second)
    print "\n\n\nRun ended at %s" % end_tm

   
    MC_rss = (MC_err**2)*n_s * n_t
#    x1 = np.random.randn(n_times) 
#    times1 = [times[i]*x1[i] for i in range(n_times)]
#    p_0 = times
#    exit(8)
#    result = fmin_slsqp(testfunc_new1, [p_0] ,  args=(A, y, n_t, n_s), \
#      #                  args=[table], \
#                        iprint=2, iter=1000000, full_output=1)


#    print "\n\n\nFinal result:", result
#    final_t = result[0]
#    PM_err = result[1][0]

#    PM_rss = (PM_err**2)*n_s * n_t
#    print "MHTN ended after %d iterations" %(mhtn_itr)
    mhtn_rss = (PM_err**2)*n_s * n_t
#    print    "\n\nMC err ", MC_err, "PM err ", PM_err, "MHTN Err:", mhtn_err, "MC RSS ", MC_rss, "PM RSS ", PM_rss
    print    "\n\nMC err ", MC_err,  "MHTN Err:", PM_err, "MC RSS ", MC_rss,   "mhtn_rss", mhtn_rss
   
    chi2 =  log(MC_rss/mhtn_rss)* n_t * n_s
    p_val = stats.chi2.cdf(x=chi2,  df=n_t)
    print "chi^2: %f, df: %d, p_val: %f"%( chi2,  n_t, p_val)
    fout = fname + ".CEM"
    output_results(fout, sp_list, site_list, MC_times, MC_rates, MC_d, PM_times, PM_rates, PM_d)

    return MC_rss, mhtn_rss

#--------------------------------------------------------------------------------------------#

if __name__ == "__main__":
   main(sys.argv[1:])
