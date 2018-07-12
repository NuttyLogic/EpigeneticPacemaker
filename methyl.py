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
from scipy.optimize import fmin_slsqp
import numpy as np
## from numpy import random *
from  time import  *
from  tree_n import  *
from random import *
from numpy import linalg as LA
from pace_maker_n import *
import matplotlib.pyplot
import pylab
import datetime


#--------------------------------------------------------------------------------------------#

def generate_log_grow_table(rates, times, d, var_t, var_s):

    n_t  = len(times)
    n_s = len(rates)

    table = [ [math.log(x,abs(b)) for x in times] for b in rates]

    return table, times
#--------------------------------------------------------------------------------------------#

def generate_table(rates, times, d, var_t, var_s):

    n_t  = len(times)
    n_s = len(rates)

    table = [ [0. for x in range(n_t)] for j in range(n_s)]
        #Here we generate the real times. We vary the chronological
        #time by var_t to obtain the PM effect. 
    t_n = [np.random.normal(times[j], var_t) for j in range (n_t)] 
    t_err = sqrt(sum([(times[j]-t_n[j])**2 for j in range (n_t)])/n_t)
    
    print "IN gen Table,  var_t %f, t_err %f, " %(var_t, t_err)
#    exit(8)
    s_err = 0.
    MC_err  = 0.
    for i in range(n_s):
        r_i = rates[i]
        d_i = d[i]
#        print "rate for site  ",i , " is ", r_i, "intercept", d_i
        for j in range(n_t):
            t_j = t_n[j]
            #  Here we generate the the methylation read for site i at
            #  individual j that is d_i + t_j*r_i. To generate noise,
            #  we add noise that is normal with var_s
            s_ij = (np.random.normal((d_i + t_j*r_i), var_s))
            s_err += (s_ij - (d_i + t_j*r_i))**2
            MC_err += (s_ij - (d_i + times[j]*r_i))**2
#            print "dist for time %f and rate %f and intercpt %f, is: %f" %(t_j, r_i, d_i, s_ij)
            table[i][j] = s_ij

    MC_err = sqrt(MC_err/(n_s*n_t))
    s_err = sqrt(s_err/(n_s*n_t))
    print "IN gen Table,  var_s %f, s_err %f, MC_err %f, %f" %(var_s, (s_err), MC_err, sqrt(s_err**2 + t_err**2) )
#    exit(8)
    return table, t_n

#--------------------------------------------------------------------------------------------#

def write_table(table, sp_list, sites, ages, fname):

    n_s = len(table)
    n_t  = len(table[0])

    if n_t != len(sp_list):
        ERR("Incompatible length n_t != len(times): %d, %d" %(n_t , len(sp_list)))

    if n_t != len(ages):
        ERR("Incompatible length n_t != len(ages): %d, %d" %(n_t , len(ages)))

    distTab = open(fname,'w')
    distTab.write("ID \t")
    distTab.write(" \t".join(["%s"%(sp_list[i]) for i in range(len(sp_list))]))
    distTab.write("\n")
    distTab.write("Age \t")
    distTab.write(" \t".join(["%f"%(ages[i]) for i in range(len(ages))]))
    
    distTab.write("\n")
    for i in range(n_s):
        distTab.write( "%s \t"%(sites[i] ))
        distTab.write(" \t".join(["%f"%(table[i][j]) for j in range(n_t)]))
        distTab.write("\n")

    distTab.close()        

 #--------------------------------------------------------------------------------------------#

def read_table(infile):
    f = open (infile)
    line = f.readline()
    line  = line.strip()
    print line
    cols = re.split('[\t\n\s,]+',line)
    print "\n\ncols", cols[:4], cols[-4:], len(cols)
#    exit(8)
    if  not re.match(r'^ID\s*',cols[0]):
        ERR ("Invalid first line %s"%cols[0])
    cols.pop(0)
    n_t = len (cols)
#    gr  = re.match(r'^age \(y\):\s*(\d+)',cols[13]).groups()
#    print "gr", gr, "len",len(gr)
#    exit(8)
    sp_list = [re.search(r'\s*(\w+)',a).groups()[0] for a in cols]
    print "sp list",  sp_list[:5], sp_list[-5:], len(sp_list)
    if len(sp_list) != n_t:
        ERR ("invalid len of names: %d vs %d"%(len(sp_list),n_t))

#    exit(8)
# Now the individual names
    line = f.readline()
    line  = line.strip()
    age_fac = 1.
    print line
    cols = re.split('[\t\n\s,]+',line)
    print "\n\ncols", cols[:4], len(cols)
    if  not re.match(r'^Age',cols[0]):
        ERR ("Invalid first line %s"%cols[0])
    if len(cols) != n_t+1:
        ERR ("invalid len of ages: %d vs %d"%(len(cols),n_t+1))
    ages = [float(a) * age_fac for a in cols[1:]]
    if len(ages) != n_t:
        ERR ("invalid len of ages %d vs %d"%(len(ages) , n_t))
    avg_age = sum(ages)/n_t
    print "ages ", ages[:5], ages[-5:], "avg %f"%(avg_age), "len %d"%( len(ages))
#    exit(8)
        
    table = []
    sites = []
    avg = []
    cov = []
    var = []
    cnt = 0
    #print "line:", line ,"\n\n\n\n\n\nkk\n"
    for line in  f:
        line  = line.strip()
#        print line
#        exit(8)
        cols = re.split('[\t\n\s,]+',line)[0:]
        l = len(cols)
        if (l != n_t +1):
            ERR("Invalid length %d instead of %d"% (l, n_t+1))
#        print l,  cols[:5]
        sites.append(re.search(r'\s*(\w+)',cols[0]).groups()[0])
        cols.pop(0)
        print "site is %s" %(sites[-1])
#        exit(8)
        row = [float(i) for i in cols]
        site_avg = sum(row)/n_t
        avg.append(site_avg)
        var.append ([cnt, sum([(x - site_avg)**2 for x in row])/n_t])
        cov.append ([cnt, sum([(row[i] - site_avg)*(ages[i] - avg_age) for i in range(n_t)])/n_t])
#        print "row", row[:5], avg[-1]
#        print "Site %d is %s, avg %f, var %f, cov %f "% (cnt, sites[-1], avg[-1], var[-1], cov[-1][1])

        if len(row) != n_t:
            ERR("Invalid length of row: %d instead of %d"% (len(row), n_t))
        table.append(row)
        cnt+=1
        if not cnt%1000:
            print "Processed %d lines"%cnt
        if cnt == 100000:
            break
#        exit(8)
    n_s = cnt
    print "cnt", cnt
    print avg[:10]
    print cov[:10]
#    if n_t != len(sex):
#        ERR("Invalid length for sex:  %d instead %d"% (n_t,  len(sex)))
#    if n_t != len(table):
#        ERR("Invalid length for sex:  %d instead %d"% (n_t,  len(table)))

#    exit(8)
    return table, ages, sp_list, sites

#--------------------------------------------------------------------------------------------#


#-----------------------------------------------------------------------------------_#

def find_soln_man(table, times):

    n_s = len(table)
    n_t  = len(table[0])
    Arr = []
    Arr2 = []
    ya = []
    n_cols = 2 * n_s
    sum_t = sum_t_sq = 0.
    for i in range(n_t):
        sum_t += times[i]
        sum_t_sq += (times[i])**2

        
    for i in range(n_cols):
        rowB = [0 for i1 in range(n_cols)]
        if i < n_s:
            rowB[i] = sum_t_sq
            rowB[i+n_s] = sum_t
        else:
            rowB[i - n_s] = sum_t
            rowB[i] = n_t
        Arr2.append(rowB)


                    
    for i in range(n_s):
        for j in range(n_t):
            t_j = times[j]
            row = [0 for i1 in range(n_cols)]   # all row entries except two will be zeros 
            row[i] = t_j
            row[i+n_s] = 1
            Arr.append(row)
            ya.append(table[i][j])

#    print len(A)
#    print A
    y = np.matrix(ya)
    B = np.matrix(Arr2)
    
    A = np.matrix(Arr)
    At = A.transpose()
    print A
    print "\n\ntranspose \n", At
#    exit(8)
    A_At = np.dot(At,A)
    Aa = np.asarray(A_At)
    Ba = np.asarray(B)
    for i in range (n_cols):
        print Aa[i], Ba[i]
        for j in range (n_cols):
            if abs(Aa[i][j] - Ba[i][j]) > .001:
                ERR ("Different vals at i %d, j %d: %f vs %f\n" % (i, j, Aa[i][j], Ba[i][j]))
    invB = LA.inv(B)
    invBa = np.asarray(invB)
    inv_A_At = LA.inv(A_At);
    invAa = np.asarray(inv_A_At)
    if len(invBa) != len(invAa):
        Err("unequal len")
    if len(invBa) != n_cols:
        Err("invalid  len")
    for i in range (n_cols):
        print invBa[i]
        for j in range (n_cols):
            if abs(invAa[i][j] - invBa[i][j]) > .001:
                ERR ("Different vals at i %d, j %d: %f vs %f\n" % (i, j, invAa[i][j], invBa[i][j]))
 

    print "lengths are OK", n_cols
    temp = np.dot(invB, At)
    tempa = np.asarray(temp)
    print "tempa", tempa
    temp1 = []
    print "invAa", invAa
    print "invBa", invBa
    print "At", At

    
    print "ns %d, n_t %d, n_cols %d" %(n_s, n_t, n_cols)
    for i in range (n_cols):
        print "tempa[i]", tempa[i]
        row = [0. for i1 in range(n_t*n_s)]   # all row entries except two will be zeros 
        for j in range(n_t):
            t_j = times[j]
            for k in range(n_s):
#                print "i %d, j %d, k %d, t_j %f, invBa[i][k] %f, invBa[i][k+n_s] %f ,  invBa[i][k] + invBa[i][k+n_s] %f,  tempa[i][n_s * j + k] %f "%\
#                    (i, j, k, t_j,  invBa[i][k], invBa[i][k+n_s],  invBa[i][k] + invBa[i][k+n_s], tempa[i][n_s * j + k])
                row[n_s * j + k] = t_j * invBa[i][k] + invBa[i][k+n_s]
        print "Row", row
        temp1.append(row)
            
    print len(tempa), len(tempa[0])
    print len(temp1), len(temp1[0])
#    exit(8)
    
    if n_s * n_t != (len(tempa[0])):
        ERR ("invalid ;ength %d vs %d\n" % ( n_s * n_t , (len(tempa))))
    print "invAa", invAa
    print "invBa", invBa
    print "At", At
    for i in range(n_cols):
        print "\n\n", i, tempa[i], '\n', temp1[i]
        for j in range(n_t*n_s):
            print i, j,  times[j], tempa[i][j], temp1[i][j]
            if abs(tempa[i][j] - temp1[i][j]) > .001:
                print 
                ERR ("Different vals at i %d, j %d, %f vs %f, diff %f\n" % (i, j, tempa[i][j], temp1[i][j], tempa[i][j] - temp1[i][j]))


    exit(8)
    yt = (y).transpose()
#    print "yt:", yt


    beta = np.dot(temp, yt)

    exit(8)

#    soln = LS_analytic (A, y)
    soln = LS_analytic_fast (A, y)
    if soln == None:
        print "this sample yielded a singular matrix. Try another sample"
        return None
#    print soln
#    exit(8)

    r_rates = soln[0][:n_s]
    r_d = soln[0][n_s:]
    r_err = sqrt(soln[1] / (n_s*n_t))
    return r_rates, r_d, r_err
    exit(8)
#--------------------------------------------------------------------------------------------#

def find_soln_direct (y, times, n_t, n_s):

    n_cols = 2 * n_s


    sumT = 0.
    sumTsqr = 0.
    for j in range(n_t):
            t_j = times[j]
            sumT += t_j
            sumTsqr += t_j**2

    B = []
    for i in range(n_s):
        row0 = (n_s * n_t) *[0.]
        for j in range(n_t):
            row0[i*n_t + j] = times[j]*(1/sumTsqr - sumT**2/((sumTsqr**2)*((sumT**2)/sumTsqr - n_t))) + sumT/(sumTsqr * (sumT**2/sumTsqr - n_t))
        B.append(row0)
    
 
    for i in range(n_s):
        row0 = (n_s * n_t) *[0.]
        for j in range(n_t):
            row0[i*n_t + j] = sumT * times[j]/(sumTsqr *( sumT**2/sumTsqr - n_t)) - 1/((sumT**2)/sumTsqr - n_t)
        B.append(row0)





    r_rates = []
    for i in  range(n_s):
        temp = 0.
        for j in range(n_t):
            temp += B[i][i*n_t + j] * y[i*n_t + j]
        r_rates.append(temp)

    r_d = []
    for i in  range(n_s):
        temp = 0.
        for j in range(n_t):
            temp += B[n_s + i][i*n_t + j] * y[i*n_t + j]
        r_d.append(temp)

    return r_rates, r_d
    exit(8)
#--------------------------------------------------------------------------------------------#

def find_soln(table, times):
    print "in find_soln"
    n_s = len(table)
    n_t  = len(table[0])
    print "n_s, n_t,", n_s, n_t
#    exit(8)
    A = []
    y = []
    n_cols = 2 * n_s  
    
    for i in range(n_s):
        for j in range(n_t):
            t_j = times[j]
            row = [0 for i1 in range(n_cols)]   # all row entries except two will be zeros 
            row[i] = t_j
            row[i+n_s] = 1
            A.append(row)
            y.append(table[i][j])

#    print len(A)
#    print A

#    soln = LS_analytic (A, y)
    soln = LS_analytic_fast (A, y)
    if soln == None:
        print "this sample yielded a singular matrix. Try another sample"
        return None
#    print soln
#    exit(8)

    r_rates = soln[0][:n_s]
    r_d = soln[0][n_s:]
    r_err = sqrt(soln[1] / (n_s*n_t))
    return r_rates, r_d, r_err
    exit(8)

    
#--------------------------------------------------------------------------------------------#

    
def preprocess1(table):
    n_s = len(table)
    n_t  = len(table[0])
    y = []

    y1 = [table[i][j]  for i in range(n_s) for j in range(n_t)]
    
    for i in range(n_s):
        for j in range(n_t):
             y.append(table[i][j])
             if y[-1] != y1[i*n_t + j]:
                 ERR("y[-1] != y1[i*n_t + ]: %6f, %6f"%( y[-1], y1[i*n_t + j]))

    
    return n_s, n_t, y
#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#

    
def preprocess(table):
    n_s = len(table)
    n_t  = len(table[0])
    n_cols = 2 * n_s  
    A = []
    y = []

    for i in range(n_s):
        for j in range(n_t):
            row = [0 for i1 in range(n_cols)]   # all row entries except two will be zeros 
            row[i+n_s] = 1  # this never changes
            A.append(row)
            y.append(table[i][j])

    return n_s, n_t, A, y
#--------------------------------------------------------------------------------------------#
def mhtn_1strt_pt_o (A, y, times, n_t, n_s, table):

    print "In mhtn_1strt_pt"
    imp = 1
    mhtn_itr = 0
    while imp > .00001:
        mhtn_itr += 1
        import datetime
        now = datetime.datetime.now()

        print "\n\n\nbefore find_soln_n", "%s-%s-%s_%0d:%0d:%02d" %(now.year, now.month, now.day, now.hour, now.minute, now.second)
        r_rates_n, r_d_n, r_e_n = find_soln_n(A, y, times, n_t, n_s)
        now = datetime.datetime.now()
        print "\n\n\nAfter find_soln_n, before find_soln_man1", "%s-%s-%s_%0d:%0d:%02d" %(now.year, now.month, now.day, now.hour, now.minute, now.second)
        r_rates1, r_d1,  = find_soln_direct(y, times, n_t, n_s)
        now = datetime.datetime.now()
        print "\n\n\nAfter  find_soln_man1", "%s-%s-%s_%0d:%0d:%02d" %(now.year, now.month, now.day, now.hour, now.minute, now.second)
        for i in range(n_s):
            if abs(r_rates_n[i] - r_rates1[i]) > .0000001:
                ERR ("Error in entry %d for rates: %6f vs %6f" %(i, r_rates_n[i] , r_rates1[i]))
        for i in range(n_s):
            if abs(r_d_n[i]  - r_d1[i] ) > .0000001:
                ERR ("Error in entry %d for r_ds: %6f vs %6f" %(i, r_d_n[i]  , r_d1[i] ))
#        print "\n\n\nsame vectors:\n"  , r_rates_n, r_d_n
#        exit(8)
                        
#            print "\n\nreturned Rates:", "\n".join(["%f\t%f" % (r_rates_n[i], r_d_n[i] ) for i in range(len(r_d_n))])
        prev_err = calc_error(table, times, r_rates_n, r_d_n)
        if abs(prev_err - r_e_n) > 0.00001:
            ERR("unequal errors MC_err !=  r_e_n")

        #    exit(8)
        if len(r_rates_n) != n_s:
            ERR("Invalid len len(r_rates_n) != n_s: %d vs %d" %(len(r_rates_n), n_s))

        if len(r_d_n) != n_s:
            ERR("Invalid len len(r_d_n) != n_s: %d vs %d" %(len(r_d_n), n_s))

#            sum_r = sum_r_sq = 0.
        sum_r_sq1 = sum([x**2 for x in r_rates_n])
        sum1 = sum([r_rates_n[j] * r_d_n[j] for j in range(n_s)])
            ## for i in range(n_s):
            ##     sum_r += r_rates_n[i]
            ##     sum_r_sq += (r_rates_n[i])**2

        sum2 = n_t * [0.]
        times_n = n_t * [0.]

#            print "sum2", sum2
        #    exit(8)
        for i in range(n_t):
            sum2[i] = sum([ r_rates_n[j] * table[j][i] for j in range(n_s)])
            times_n[i] = (sum2[i] - sum1)/ sum_r_sq1
#            times_n[i] = max((sum2[i] - sum1)/ sum_r_sq1, 0)

        new_err = calc_error(table, times_n, r_rates_n, r_d_n)
            
        imp = prev_err - new_err
        print "imp", imp
        if (imp < -0.00001 and prev_err > .00001):
            ERR("new_err > prev_err: %f vs %f" % (new_err, prev_err))
        now = datetime.datetime.now()
                
#        print "%s_%s:%s: prev err %f, new err %f, imp %f" % \
#            ( now.day, now.hour,now.minute,  prev_err,  new_err, imp)


        times = times_n

    return new_err, times_n, r_rates_n, r_d_n, mhtn_itr

#--------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
def PMEM (times, table, itr_limit):

    print "In mhtn_1strt_pt"
    n_s = len(table)
    n_t  = len(table[0])
    y = [table[i][j]  for i in range(n_s) for j in range(n_t)]
    imp = 1
    EM_itr = 0
    while imp > .00001:
        EM_itr += 1
        import datetime
        now = datetime.datetime.now()

        print "\n\n\nbefore find_soln_direct", "%s-%s-%s_%0d:%0d:%02d" %(now.year, now.month, now.day, now.hour, now.minute, now.second)
        r_rates1, r_d1,  = find_soln_direct(y, times, n_t, n_s)
        now = datetime.datetime.now()
        print "\n\n\nAfter  find_soln_man1", "%s-%s-%s_%0d:%0d:%02d" %(now.year, now.month, now.day, now.hour, now.minute, now.second)

                        
#            print "\n\nreturned Rates:", "\n".join(["%f\t%f" % (r_rates_n[i], r_d_n[i] ) for i in range(len(r_d_n))])
        prev_err = calc_error(table, times, r_rates1, r_d1)
        if EM_itr == 1:
            init_err = prev_err
            init_rates = r_rates1
            init_d = r_d1

        #    exit(8)
        if len(r_rates1) != n_s:
            ERR("Invalid len len(r_rates_n) != n_s: %d vs %d" %(len(r_rates_n), n_s))

        if len(r_d1) != n_s:
            ERR("Invalid len len(r_d_n) != n_s: %d vs %d" %(len(r_d_n), n_s))

#            sum_r = sum_r_sq = 0.
        sum_r_sq1 = sum([x**2 for x in r_rates1])
        sum1 = sum([r_rates1[j] * r_d1[j] for j in range(n_s)])
            ## for i in range(n_s):
            ##     sum_r += r_rates_n[i]
            ##     sum_r_sq += (r_rates_n[i])**2

        sum2 = n_t * [0.]
        times_n = n_t * [0.]

#            print "sum2", sum2
        #    exit(8)
        for i in range(n_t):
            sum2[i] = sum([ r_rates1[j] * table[j][i] for j in range(n_s)])
            times_n[i] = (sum2[i] - sum1)/ sum_r_sq1
#            times_n[i] = max((sum2[i] - sum1)/ sum_r_sq1, 0)

        new_err = calc_error(table, times_n, r_rates1, r_d1)
            
        imp = prev_err - new_err
        print "imp", imp
        if (imp < -0.00001 and prev_err > .00001):
            ERR("new_err > prev_err: %f vs %f" % (new_err, prev_err))
        now = datetime.datetime.now()
                
#        print "%s_%s:%s: prev err %f, new err %f, imp %f" % \
#            ( now.day, now.hour,now.minute,  prev_err,  new_err, imp)


        times = times_n
        if EM_itr == itr_limit: break
  
    return init_err, init_rates, init_d, new_err, times_n, r_rates1, r_d1, EM_itr

#--------------------------------------------------------------------------------------------#
#def ckSolnOK( times, r_rates)
#--------------------------------------------------------------------------------------------#
def mhtn_HC1(  MC_times, n_t, n_s, table, runs):

    now = datetime.datetime.now()
    start_tm = "%s-%s-%s_%s%s" %(now.year, now.month, now.day, now.hour, now.minute)
    print "\n\n%s_%s:%s: At mhtn_HC" %( now.day, now.hour,now.minute), 

#    runs = 5
    y = [table[i][j]  for i in range(n_s) for j in range(n_t)]

    MC_rates, MC_d = find_soln_man1(y,  MC_times, n_t, n_s)
    MC_err = calc_error(table, MC_times, MC_rates, MC_d)
    
    avg_age = sum(MC_times)/n_t
    print "avg age:", avg_age
    print "MC_ERR", MC_err
#    exit(4)
    sols = []
    avg_times = n_t *[0.]
    avg_rate = n_s *[0.]
    max_err = min_err = -1.
    tot_err = 0.
    for i in range(runs):
    
        rnd_times = [np.random.normal(avg_age, 3*sqrt(avg_age)) for j in range (n_t)] #
        print "Time tossed: %s"% " ".join(["%3.3f"%(rnd_times[j]) for j in range((10))])
        
        new_err, times_n, r_rates_n, r_d_n, mhtn_itr = mhtn_1strt_pt1(y,  rnd_times, n_t, n_s, table)
        tot_err += new_err
        if max_err == -1:
            max_err =  min_err = new_err
        if new_err > max_err:
            print " max_err updated from %f to %f" %( max_err, new_err)
            max_err = new_err
        if new_err < min_err:
            print " min_err updated from %f to %f" %( min_err, new_err)
            min_err = new_err

        for j in range(n_t):
            avg_times[j] += times_n[j]
        for j in range(n_s):
            avg_rate[j] +=  r_rates_n[j]
        sols.append([new_err, times_n, r_rates_n, r_d_n, mhtn_itr]) 

        print "Solns  %d returned from mhtn_1strt_pt, err %3.3f"% (i,  new_err), times_n[:5]

#    exit(8)            
    print "\n\n\nSolns returned from mhtn_HCt. min err %3.3f, max %3.3f" %( max_err, min_err)
    print " ".join(["(%3.3f,  %3.2f, %3.2f), "%(avg_times[k]/runs, MC_times[k], avg_times[k]/avg_times[0]) for k in range((runs))])

#    exit(8)
    return MC_err, MC_rates, new_err, [avg_times[j]/runs for j in  range(n_t)], [avg_rate[j]/runs for j in  range(n_s)], mhtn_itr
    exit(8)

#-------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------#
def mhtn_correct_once(MC_times, table):

    itr_limit = 1
    MC_err, MC_rates, MC_d, cor_once_err, cor_once_times, rcor_once_rates, cor_once_d, EM_itr = PMEM (MC_times, table, itr_limit)

    if EM_itr != 1:
        ERR( "Invalid EM_itr %d" % (EM_itr))


    return MC_err, MC_rates, MC_d, cor_once_err, cor_once_times, rcor_once_rates, cor_once_d

#--------------------------------------------------------------------------------------------#
def EM_MC_strt_pt(MC_times, table):

    itr_limit = 0
    MC_err, MC_rates, MC_d, PM_err, PM_times, PM_rates, PM_d, EM_itr = PMEM (MC_times, table, itr_limit)

    if EM_itr != 1:
        print "Iterated ", EM_itr, "times"


    return MC_err, MC_rates, MC_d, PM_err, PM_times, PM_rates, PM_d

#--------------------------------------------------------------------------------------------#
def mhtn_HC_rnd_pts(A, y,  MC_times, n_t, n_s, table, runs):

    now = datetime.datetime.now()
    start_tm = "%s-%s-%s_%s%s" %(now.year, now.month, now.day, now.hour, now.minute)
    print "\n\n%s_%s:%s: At mhtn_HC" %( now.day, now.hour,now.minute), 

#    runs = 5
    
    MC_rates, MC_d, MC_err = find_soln_n(A, y,  MC_times, n_t, n_s)
    avg_age = sum(MC_times)/n_t
    print "avg age:", avg_age
    print "MC_ERR", MC_err
#    exit(4)
    sols = []
    avg_times = n_t *[0.]
    avg_rate = n_s *[0.]
    avg_d = n_s *[0.]
    max_err = min_err = -1.
    tot_err = 0.
    for i in range(runs):
    
        rnd_times = [np.random.normal(avg_age, 3*sqrt(avg_age)) for j in range (n_t)] #
        print "Time tossed: %s"% " ".join(["%3.3f"%(rnd_times[j]) for j in range((10))])
        
        new_err, times_n, r_rates_n, r_d_n, mhtn_itr = mhtn_1strt_pt_o(A, y,  rnd_times, n_t, n_s, table)
        tot_err += new_err
        if max_err == -1:
            max_err =  min_err = new_err
        if new_err > max_err:
            print " max_err updated from %f to %f" %( max_err, new_err)
            max_err = new_err
        if new_err < min_err:
            print " min_err updated from %f to %f" %( min_err, new_err)
            min_err = new_err

        for j in range(n_t):
            avg_times[j] += times_n[j]
        for j in range(n_s):
            avg_rate[j] +=  r_rates_n[j]
            avg_d[j] +=  r_d_n[j]
        sols.append([new_err, times_n, r_rates_n, r_d_n, mhtn_itr]) 

        print "Solns  %d returned from mhtn_1strt_pt, err %3.3f"% (i,  new_err), times_n[:5]

#    exit(8)            
    print "\n\n\nSolns returned from mhtn_HCt. min err %3.3f, max %3.3f" %( max_err, min_err)
    print " ".join(["(%3.3f,  %3.2f, %3.2f), "%(avg_times[k]/runs, MC_times[k], avg_times[k]/avg_times[0]) for k in range((runs))])

#    exit(8)
    return MC_err, MC_rates, MC_d, new_err, [avg_times[j]/runs for j in  range(n_t)], [avg_rate[j]/runs for j in  range(n_s)], \
        [avg_d[j]/runs for j in  range(n_s)], mhtn_itr
    exit(8)

#-------------------------------------------------------------------------------------#

def find_soln_n(A, y, times, n_t, n_s):

#    n_cols = 2 * n_s  
#    print type(A), len(A), n_s, n_t
#    exit(8)
    for i in range(n_s):
        for j in range(n_t):
            t_j = times[j]
            A[i*n_t+j][i] = t_j

#    print len(A)
#    print A

#    soln = LS_analytic (A, y)
    soln = LS_analytic_fast (A, y)
    if soln == None:
        print "this sample yielded a singular matrix. Try another sample"
        return None
#    print soln
#    exit(8)

    r_rates = soln[0][:n_s]
    r_d = soln[0][n_s:]
    r_err = sqrt(soln[1] / (n_s*n_t))
    return r_rates, r_d, r_err
    exit(8)
#--------------------------------------------------------------------------------------------#

def calc_error_detailed(table, times, rates, d):
    n_s = len(table)
    n_t  = len(table[0])
    tot_err = 0.
    site_err = [0.] * n_s
    time_err = [0.] * n_t
    for i in range(n_s):
        
        for j in range(n_t):
            t_j = times[j]
#            print "predict is %f * %f + %f = %f" %( t_j, rates[i], d[i], t_j* rates[i] - d[i]), "Data = %f, err = (%f - %f)^2 = %f"%  (table[i][j],  table[i][j], t_j* rates[i] + d[i], (table[i][j] - t_j* rates[i] - d[i])**2  )
            err = (table[i][j] - t_j* rates[i] - d[i])**2
            site_err[i] += err
            time_err[j] += err
            tot_err += err

    return sqrt(tot_err / (n_s*n_t)), site_err, time_err


#--------------------------------------------------------------------------------------------#

def calc_error(table, times, rates, d):
    n_s = len(table)
    n_t  = len(table[0])
    tot_err = 0.
    
    for i in range(n_s):
        for j in range(n_t):
            t_j = times[j]
#            print "predict is %f * %f + %f = %f" %( t_j, rates[i], d[i], t_j* rates[i] - d[i]), "Data = %f, err = (%f - %f)^2 = %f"%  (table[i][j],  table[i][j], t_j* rates[i] + d[i], (table[i][j] - t_j* rates[i] - d[i])**2  )
            err = (table[i][j] - t_j* rates[i] - d[i])**2
            tot_err += err

    return sqrt(tot_err / (n_s*n_t))



#--------------------------------------------------------------------------------------------#

def output_results(run_id, sp_list, site_list, MC_times, MC_rates, MC_d, PM_times, PM_rates, PM_d):
    n_t = len(MC_times)
    n_s = len(MC_rates)
    if n_t != len(sp_list):
        ERR ("Invalid len %d"%(n_t))
    if n_s != len(site_list):
            ERR ("Invalid len %d"%(n_s))

    ftimesrep = open("meth-"+run_id+"-times-MCvsPM.csv",'w')
    ftimesrep.write("SampleID, MC-age, PM-age, ratio\n")
    time_srt = sorted([ [MC_times[j], PM_times[j]] for j in range(n_t)], key=lambda (k,v): (k))
    print "\n".join(["%3.3f, %3.3f"%(time_srt[j][0],  time_srt[j][1]) for j in range((n_t))])
#    exit(8)
    
    for j in range (n_t):
         ftimesrep.write("%s, %f, %f, %f\n"%(sp_list[j], MC_times[j], PM_times[j], MC_times[j]/PM_times[j]))

    ftimesrep.close()

    
    fratesrep = open("meth-"+run_id+"-rates-MCvsPM.csv",'w')
    fratesrep.write("Site, MC-rate, PM-rate, rate ratio, MC-d, PM-d\n")
    for j in range (n_s):
        print "%s, %f, %f, %f, %f, %f, %f\n"% \
           (site_list[j], MC_rates[j], PM_rates[j], MC_rates[j]/PM_rates[j],  MC_d[j], PM_d[j],  MC_d[j]/PM_d[j])

        fratesrep.write("%s, %f, %f, %f, %f, %f, %f\n"% \
                        (site_list[j], MC_rates[j], PM_rates[j], MC_rates[j]/PM_rates[j],  MC_d[j], PM_d[j],  MC_d[j]/PM_d[j]))

    fratesrep.close()

         
#--------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------#



def simulate_run (n_sites, n_times, var_t, var_s, trend):

    new_run = 1
#    new_run = 0
    seed(os.urandom(10))

    fname = "meth-tab-n%d-m%d.txt"%(n_sites, n_times)
#    fname = "meth-tab-n656-top-50-cov-sites.txt"
    if new_run:

        rates =  np.random.randn(n_sites)
#        rates = [uniform(2,10) for i in range (n_sites)]
        print "Rates", (rates)
#        exit(8)
    
#        times = [exp(2+i) for i in np.random.randn(n_times)]
        MC_times = [ uniform(0,100) for i in range (n_times)]
#        print "\n\ntimes", MC_times
#        exit(8)
    
        PM_trend_times = [x**trend * 100/100**trend   for x in MC_times]

        
        intercept =  [uniform(0,1) for i in range (n_sites)]
        table, PM_actual_times = generate_table(rates, PM_trend_times, intercept, var_t, var_s)
#        t_err = sqrt(sum([(PM_times[j] - times[j])**2 for j in range(n_times)])/n_times)
#        print "t err", t_err

        

        sp_list = ["sp%d" %(x) for x in range(n_times)]

        sites = ["st%d" %(x) for x in range(n_sites)]

        write_table(table, sp_list, sites, MC_times, fname)
#        exit(8)
    else:
        table, MC_times, sp_list, sites = read_table(fname)
#    exit(8)

#    n_s, n_t, A, y = preprocess (table)
    
#    n_s, n_t, y = preprocess1 (table)

#    MC_err, MC_rates, PM_err, PM_times, PM_rates, mhtn_itr  = mhtn_HC(A, y, times, n_t, n_s, table, mhtn_strts)
    MC_err, MC_rates, MC_d, PM_err, PM_ret_times, PM_ret_rates, PM_ret_d = EM_MC_strt_pt(MC_times, table)

#    MC_err, MC_rates, PM_err, PM_times, PM_rates, mhtn_itr  = mhtn_HC1( y, times, n_t, n_s, table, mhtn_strts)
#    if abs(MC_err - MC_err1) > .00000001:
#        ERR("difference in MC_err errors: %10f va %10f"% (MC_err ,MC_err1))
#    if abs(PM_err - PM_err1) > .00000001:
#        ERR("difference in PM_err errors: %10f va %10f"% (PM_err ,PM_err1))

    #    mhtn_err, mhtn_itr = mhtn_HC(A, y, times, n_t, n_s, table)
#    r_rates, r_d, r_e = find_soln_man(table, times)
#    r_rates_n, r_d_n, r_e_n = find_soln_n(A, y, times, n_t, n_s)
#    print "\n\nreturned Rates:", "\n".join(["%f\t%f" % (r_rates_n[i], rates[i]) for i in range(len(rates))])
#    exit(8)
#    print "\n\nIntercept:", intercept

#    MC_err = calc_error(table, times, r_rates_n, r_d_n)
#    print "MC err ", MC_err, r_e_n
#    if abs(MC_err - r_e_n) > 0.00001:
#        ERR("unequal errors MC_err !=  r_e_n")

    MC_rss = (MC_err**2)*n_sites * n_times
##     x1 = np.random.randn(n_times) 
##     times1 = [times[i]*x1[i] for i in range(n_times)]
##     p_0 = times1
## #    exit(8)
##     result = fmin_slsqp(testfunc_new1, [p_0] ,  args=(A, y, n_t, n_s), \
##       #                  args=[table], \
##                         iprint=2, iter=1000000, full_output=1)


##     print "\n\n\nFinal result:", result
##     final_t = result[0]
##     PM_err = result[1][0]

    PM_rss = (PM_err**2)*n_sites * n_times
#    print "MHTN ended after %d iterations" %(mhtn_itr)
    
    print    "\n\nMC err ", MC_err, "PM err ", PM_err, "MC RSS ", MC_rss, "PM RSS ", PM_rss
    return MC_rss, PM_rss
#-----------------------------------------------------------------------------------_#

def testfunc_new1(p, *args):
#    exit(8)
#    table = args[0]
#    n_s = len(table)
#    n_t  = len(table[0])
    A, y, n_t, n_s = args[:] 
#    print "A, y, n_t, n_s", A, y, n_t, n_s
#    exit(8)
    
    times = p # The parameter is the ST edges the procedure  -
    if n_t != len(times):
        ERR("n_t != len(times): %d, %d" % (n_t, len(times)))
    # the supertree edges
    #    print "supertree edges (betas):", betas
    #    print "gene specific rates:", rs
    for t1 in times:
        if t1 <= 0:
            return np.inf 
    

    r_rates, r_d, r_e =  find_soln_n(A, y, times, n_t, n_s)
    print "at testfunc_new1, err", r_e

#    exit(8)

        
#    err = calc_error (table, times, r_rates, r_d)
#    if err == None:
#        ERR("unidentified error: None")

#    if not err:
#        ERR("unidentified error: None")

#    import datetime
#    now = datetime.datetime.now()
#    start_tm = "%s-%s-%s_%s%s" %(now.year, now.month, now.day, now.hour, now.minute)
#    print "%s_%s:%s" %( now.day, now.hour,now.minute), \
#          "final total err: %f  "% (r_e)
#    exit(8)
    return r_e



#-----------------------------------------------------------------------------------_#
def ieqcons1(x, *args):
#    print "in test_ieqcons x. x", x
    A, y, n_t, n_s = args[:] 
#    n_t = args[1]
#    b = x[:-n_t]
#    r = x[-n_t:]
#    const = []
#    const.append(b[0]+b[2] - b[3])
#    const.append(b[3]+b[2] - b[1])
    times = x # The parameter is the ST edges the procedure  -
    if n_t != len(times):
        ERR("n_t != len(times): %d, %d" % (n_t, len(times)))
    const = [i for i in x]
#    print const + const1
#    exit(8)
#    return array(const + const1)
    return array( const)

#-----------------------------------------------------------------------------------_#

def testfunc(p, *args):
#    print "at testfunc_new"
#    exit(8)
    table = args[0]
    n_s = len(table)
    n_t  = len(table[0])
#    A, y, n_t, n_s = args[:]
#    print "A, y, n_t, n_s", A, y, n_t, n_s
#    exit(8)
    
    times = p # The parameter is the ST edges the procedure  -
    if n_t != len(times):
        ERR("n_t != len(times): %d, %d" % (n_t, len(times)))

#    for t1 in times:
#        if t1 <= 0:
#            return np.inf 
    
    r_rates, r_d, r_e = find_soln(table, times)
#    r_rates, r_d, r_e = find_soln_n(A, y, times, n_t, n_s)

#    exit(8)

        
#    err = calc_error (table, times, r_rates, r_d)
#    if err == None:
#        ERR("unidentified error: None")

#    if not err:
#        ERR("unidentified error: None")

#    import datetime
#    now = datetime.datetime.now()
#    start_tm = "%s-%s-%s_%s%s" %(now.year, now.month, now.day, now.hour, now.minute)
#    print "%s_%s:%s" %( now.day, now.hour,now.minute), \
#          "final total err: %f  "% (r_e)
#    exit(8)
    return r_e

