#!/usr/bin/python


###!/opt/python-2.5/bin/python
###!/Library/Frameworks/Python.framework/Versions/Current/bin/python


from phylo import *
from methyl import *

#exit(8)

def prepareData(infile, top_times, top_sites, crit):
    f = open (infile)
    line = f.readline()
    line  = line.strip()
    print line
    cols = re.split('\s+\"',line)
    print "\n\ncols", cols[:4], cols[-4:], len(cols)
#    exit(8)
    if  not re.match(r'^!Sample_characteristics_ch1\s*',cols[0]):
        ERR ("Invalid first line %s"%cols[0])
    print "cols0", cols[0]
    cols.pop(0)
    print "cols0", cols[0]
    n_t = len (cols)
#    gr  = re.match(r'^age \(y\):\s*(\d+)',cols[13]).groups()
    gr  = re.match(r'^age:\s*(-?\d+(\.\d+)?)',cols[13]).groups()
    print "gr", gr, "len",len(gr)
#    exit(8)
#    ages = [float(re.match(r'^age \(y\):\s*(\d+)',a).groups()[0]) for a in cols]
    ages = []
    missing_age = []
    none_cnt = 0
    sum_ages = 0
    for cnt in range(n_t):
        a = cols[cnt]
        print a
        if re.match(r'^age:\s*(-?\d+(\.\d+)?)',a):
            age = float(re.match(r'^age:\s*(-?\d+(\.\d+)?)',a).groups()[0])
            print "OK age", age
            sum_ages += age
            ages.append(age)
        if re.match(r'^age:\s*NA',a):
            print "NOT OK age"
            none_cnt +=1
            missing_age.append(cnt)
        if re.match(r'^age:\s*(-?\d+(\.\d+)?)',a):
            print "OK age", a
        else:
            ERR ( "NOT OK age %s"%( a))

#        print re.match(r'^age:\s*(\d+)',a).groups()
#    exit(4)
    print "none_cnt", none_cnt
    print "missing age", missing_age
#    ages = [float(re.match(r'^age:\s*(\d+)',a).groups()[0]) for a in cols]
    real_nt =  n_t - none_cnt
    if len(ages) != real_nt:
        ERR ("invalid len of ages %d vs %d"%(len(ages) , real_nt))
    avg_age = sum(ages)/len(ages)
    std_age = sqrt(sum([(ages[i] - avg_age)**2 for i in range(len(ages))]))
    print "ages ", ages[:5], ages[-5:], len(ages), ", cnt:", cnt
#    exit(8)
# Now the individual names
    line = f.readline()
    line  = line.strip()
    print line
    cols = re.split('\"\s+\"',line)
    print "\n\ncols", cols[:4], len(cols)
    if  not re.match(r'^\"ID_REF',cols[0]):
        ERR ("Invalid first line %s"%cols[0])
    if len(cols) != n_t+1:
        ERR ("invalid len of names: %d vs %d"%(len(cols),n_t+1))
    sp_list = [re.search(r'\s*(\w+)',a).groups()[0] for a in cols[1:]]
    print "sp list",  sp_list[:5], sp_list[-5:], len(sp_list)
    if len(sp_list) != n_t:
        ERR ("invalid len of names: %d vs %d"%(len(sp_list),n_t))
    for x in reversed(missing_age):
        print "no age for ", x, sp_list[x], ". removing "
        sp_list.pop(x)
#    exit(8)
        
    table = []
    sites = []
    avg = []
    cov = []
    var = []
    pcc = []
    cnt = 0
    #print "line:", line ,"\n\n\n\n\n\nkk\n"
    cnt_empty = 0
    discarded = 0
    for line in  f:
        line  = line.strip()
#        print line
#        exit(8)
        cols = re.split('[\t\n\s,]+',line)[0:]
        l = len(cols)
        if (l != n_t +1):
            print "line:", line
            print ("Invalid length fr row %d:  %d instead of %d"% (cnt, l, n_t+1))
            print "col[0] %s, col[1] %s,  col[2] %s,  col[3] %s,last %s" %(cols[0], cols[1], cols[2], cols[3], cols[-1])
            continue
            ERR("Invalid length fr row %d:  %d instead of %d"% (cnt, l, n_t+1))
#        print l,  cols[:5]
        sites.append(re.search(r'\s*(\w+)',cols[0]).groups()[0])
        cols.pop(0)
#        exit(8)
        try :
            row = [float(i) for i in cols]
        except ValueError:
            print "Invalid val was found at line ", cnt, "site ", sites[-1]
            discarded += 1
            continue
        for x in reversed(missing_age):
            row.pop(x)
        if len(row) != real_nt:
                ERR("Invalid length for row %d:  %d instead of %d"% (len(row) != real_nt))
##         row = []
##         for x in cols:
##                 if  re.match(r'^\s*$',x):
##                     print "Empty val:", x, "at line ", cnt,"Len:", l
##                     row.append(0.0)
##                     cnt_empty +=1
##                 else:
## #                    print "Not empty ", x
##                     row.append(float (x))
#        row = [ float(i) for i in cols]
        avg_site = sum(row)/real_nt
        std_site = sqrt(sum([(row[i] - avg_site)**2 for i in range(real_nt)]))
        avg.append(avg_site)
        var.append ([cnt, sum([(x - avg_site)**2 for x in row])/real_nt])
        cov.append ([cnt, sum([(row[i] - avg_site)*(ages[i] - avg_age) for i in range(real_nt)])/real_nt])
        pcc.append([cnt,  abs(sum([(row[i] - avg_site)*(ages[i] - avg_age) for i in range(real_nt)]))/(std_site*std_age)])
#        print "row", row[:5], avg[-1]
#        print "Site %d is %s, avg %f, var %f, cov %f "% (cnt, sites[-1], avg[-1], var[-1], cov[-1][1])

        if len(row) != real_nt:
            ERR("Invalid length of row %d: %d instead of %d"% (cnt, len(row), real_nt))
        table.append(row)
        cnt+=1
        if not cnt%1000:
            print "Processed %d lines, discarded %d"% (cnt, discarded)
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

    if n_s < top_sites:
        top_sites = n_s
        
    if n_t < top_times:
        top_times = n_t

    var.sort(key=lambda tup: tup[1], reverse=True) # sort in place by the 1st element
    cov.sort(key=lambda tup: tup[1], reverse=True) # sort in place by the 1st element
    pcc.sort(key=lambda tup: tup[1], reverse=True) # sort in place by the 1st element
    print "first 5 and last 5 var:",  var[:5], var[-5:]
#    print "first 5 and last 5 cov:",  cov[:5], cov[-5:]
    for j in range(top_sites):
        print "Var", var[j], "Cov", cov[j]
#        print cov[j]
#    exit(8)
    
    if crit == "var":
        top_site_names = [sites[var[j][0]] for j in range(top_sites)]
        distTabfn = "meth-tab-hmn-n%d-top-%d-var-sites.txt"%(top_times, top_sites)
        tableTop =  [ table[var[j][0]][:top_times] for j in range(top_sites)]

    elif crit == "cov":
        top_site_names = [sites[cov[j][0]] for j in range(top_sites)]
        distTabfn = "meth-tab-hmn-n%d-top-%d-cov-sites.txt"%(top_times, top_sites)
        tableTop =  [ table[cov[j][0]][:top_times] for j in range(top_sites)]
    elif crit == "pcc":
        top_site_names = [sites[pcc[j][0]] for j in range(top_sites)]
        distTabfn = "meth-tab-hmn-n%d-top-%d-abs-pcc-sites.txt"%(top_times, top_sites)
        distTabfn = "meth-tab-GSE42861-n%d-top-%d-abs-pcc-sites.txt"%(top_times, top_sites)
        distTabfn = "meth-tab-GSE87571-n%d-top-%d-abs-pcc-sites.txt"%(top_times, top_sites)
        tableTop =  [ table[pcc[j][0]][:top_times] for j in range(top_sites)]
    else:
        ERR("UNrecognized criterion %s"%(crit))
       
    write_table(tableTop, sp_list[:top_times], top_site_names, ages[:top_times], distTabfn)

    exit(8)
    print "Now transposing the table to have sites as lines and times (persons) as colunms"
    tableT =  [ [0. for x in range(top_times)] for j in range(top_sites)]    


    for i in range(top_times):
        for j in range(top_sites):
            tableT[j][i] = table[i][var[j][0]]


    return tableT, ages, sp_list, top_site_names

#exit(8)
#--------------------------------------------------------------------------------------------#

         
#--------------------------------------------------------------------------------------------#
def main(argv):
    print "in main"
    if not argv[0]:
        ERR("No input file was supplied")
    infile = argv[0]
    top_sites = 1000
    top_times = 1000
#    crit = "var"
    crit = "pcc"
#    crit = "cov"

    table, MC_times, sp_list, site_list = prepareData(infile, top_times, top_sites, crit)
    if len(MC_times) != len(table[0]):
        ERR ("Incompatible time: %d vs %d" %(len(MC_times) , len(table[0])))
    n_t = len(MC_times)
    n_s = len(site_list)
    print "after prepareData, before MC soln, table size",  len(table), len(table[0])

#-----------------------------------------------------------------------------------_#


if __name__ == "__main__":
   main(sys.argv[1:])
#   main (["yyy1-f10000"])
#   main (["/Users/ssagi/Downloads/fff1.txt"])
   main (["Users/sagisnir/Downloads/GSE87571_additional_sample_chararcteristics.csv"])
