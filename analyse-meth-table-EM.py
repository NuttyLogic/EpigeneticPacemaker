#!/usr/bin/python

#!/Library/Frameworks/Python.framework/Versions/Current/bin/python

#!/opt/python-2.5/bin/python



"""
This script tests fmin_slsqp using Example 14.4 from Numerical Methods for
Engineers by Steven Chapra and Raymond Canale.  This example maximizes the
function f(x) = 2*x0*x1 + 2*x0 - x0**2 - 2*x1**2, which has a maximum
at x0=2, x1=1.
"""

import scipy.stats as stats
from methyl import *
from pace_maker_n import *

def times_n(line, n):
    line_n = []
    for x in line:
        line_n += [x] *n 
    return line_n

#-------------------------------------------------------------------------------------#
def main (argv):

    fname = "/Users/colinfarrell/Documents/UPM_Project/EpigeneticClock_UniversalPaceMaker/tests/test_data/meth-tab-n100-m200.txt"
    table, MC_times, sp_list, site_list = read_table(fname)

    n_s = len(table)
    n_t = len(table[0])

    MC_err, MC_rates, MC_d, PM_err, PM_times, PM_rates, PM_d = EM_MC_strt_pt(MC_times, table)

    MC_rss = (MC_err**2)*n_s * n_t

    mhtn_rss = (PM_err**2)*n_s * n_t
    print("\n\nMC err ", MC_err,  "MHTN Err:", PM_err, "MC RSS ", MC_rss,   "mhtn_rss", mhtn_rss)
   
    chi2 = log(MC_rss/mhtn_rss)* n_t * n_s
    p_val = stats.chi2.cdf(x=chi2,  df=n_t)
    print("chi^2: %f, df: %d, p_val: %f"%( chi2,  n_t, p_val))
    fout = fname + ".CEM"
    output_results(fout, sp_list, site_list, MC_times, MC_rates, MC_d, PM_times, PM_rates, PM_d)

    return MC_rss, mhtn_rss

#--------------------------------------------------------------------------------------------#

if __name__ == "__main__":
   main(sys.argv[1:])
