#!/Library/Frameworks/Python.framework/Versions/Current/bin/python

#!/opt/python-2.5/bin/python


##!/usr/bin/python



import string, sys, re, copy, os, time, random, math
from scipy.optimize import fmin_slsqp
import numpy  as np
from  time import  *
from  tree_n import  *
from  phylo import  *
import time as tm
from  random import *
from numpy import linalg as LA
#-----------------------------------------------------------------------------------_#

def guess_start_vals(gtrees, n_e):
    n_t = len(gtrees)
    print "n_t:", n_t
    orig_n_t = n_t 
    if n_t > 100:
        n_t = 100
    else:
        #        n_t = (n_t+1)/2
        pass
        
    start_vals = None
    count = 0
    while start_vals == None:
        gene_sample = (sample(range(orig_n_t), n_t))

        gene_list = [gtrees[i] for i in gene_sample]
        print gene_list

        start_vals = solve_sample_LS(gene_list, n_e)
        if start_vals == None:
            # ERR("failed failed")
            pass
        count +=1
        if count > 4:
            ERR("could not find start values for more than 20 times")
          #        exit(8)
        
    
    print "\n\nAfter loop of %d samples: start values:" %(count),  start_vals
#    exit(8)
#    print "\n\nexp:", [math.exp(x1) for x1 in start_vals[:n_t - 1]]
    
    rs_sample = [1.] + [math.exp(x1) for x1 in start_vals[:n_t - 1]] # we get the rates only for the genes in the sample
    betas =  [math.exp(x1) for x1 in start_vals[n_t-1 :]]
#    print "betas:\n", betas
    
#    exit(8)
    if len(rs_sample) != n_t:
        ERR("length of rs and n_t are different %d vs %d" %(len(rs), n_t))
      
    if len(betas) != n_e:
        ERR("length of betas and n_e are different %d vs %d" %(len(betas), n_e))
        
      
      
        #print "\n\nSampled rs (%d)" % len(rs_sample), rs_sample
        #print "\n\nbetas (%d)" % len(betas), betas
        #test_err = calc_error (gene_list, betas, rs_sample)
    
        #    print "test Err", test_err
#    exit(8)

    rs = obtain_rates_from_edges(betas, gtrees)  # Now we get the rates of all the genes.

      
    print "\n\nrs (%d)" % len(rs), rs
    print "\n\nbetas (%d)" % len(betas), betas

#    exit(8)
    return rs, betas

#----------------------------------------------------------------------

#-----------------------------------------------------------------------------------_#

def solve_sample_LS(sample_trees, n_e):

    n_t = len(sample_trees)
        #    n_e = len(gtrees[0])

    #    We here build the matrix to the LS analytic solution
    #    First n_t-1 columns are for the 2nd-n_t'th gene rate
    #    last n_e columns, are for ST edges
    #    We have a row for every GT value. Here we are interested only in 1:1 mappings
    #    where a GT edge maps to a ST edge

    A = [] # our matrix for the LS solution
    n_cols = n_t -1 + n_e # We assign arbitrarily 1 to gene 1's \
    # rate, so it is not represented in the
    # matrix 
    e_lens = []
    edge_appeared = [0]* n_e

    print "Number of gene trees: ", n_t
#    exit(8)
    for t1 in range(n_t):
        print "Now analysing tree ", t1, sample_trees[t1]
        gt_edges = sample_trees[t1]
        print "#edge in the tree %d" % (len(gt_edges))
        for e1 in range(len(gt_edges)):
            [GT_e_len, mapped_path] = gt_edges[e1]
            if len(mapped_path) == 1:
                print "1:1 mapping to edge ", mapped_path
                st_edge = list(mapped_path)[0]
                if GT_e_len <= 0.:
                    print "Non positive val for GT_e_len: %f" % (GT_e_len)
                    #                    break
                    #                    ERR("Non positive val for GT_e_len: %f" % (GT_e_len))
                    GT_e_len = .0001
                e_lens.append(math.log(GT_e_len))
                row = [0 for i1 in range(n_cols)]   # all row entries except two will be zeros 
                row[n_t - 1 + st_edge] = 1
                edge_appeared[st_edge] = 1
                if t1 > 0: # for gene no 0, indication is added as it is assigned value 1 arbitrarily
                    row[t1-1 ] = 1
                print "row ", row
                A.append(row)

            else:
                print "edge %d with len %3.3f mapped to few edges %s" \
                %(e1, GT_e_len, str(mapped_path))
                # this is JUST 
                WARN( "This is NOT A PERFECT TREE!!!!!\n\n\n" \
                     "edge %d with len %3.3f mapped to few edges %s" \
                %(e1, GT_e_len, str(mapped_path))) 

    if sum(edge_appeared) < len(edge_appeared):
        WARN  ("Not all edges touched, only %d out of %d\n" %\
            (sum(edge_appeared), len(edge_appeared)))
        return None
    else:
        WARN( "All edges touched,  %d out of %d\n" % \
          (sum(edge_appeared), len(edge_appeared)))
          #        exit(9)
    if len(A) != len(e_lens):
        ERR ("Unequal lengths of edge lengths and matrix")

    print "edge lens:", e_lens

    soln = LS_analytic (A, e_lens)
    if soln == None:
        print "this sample yielded a singular matrix. Try another sample"
        return None
    
    return soln

#----------------------------------------------------------------------
#----------------------------------------------------------------------

def LS_analytic_fast (Arr, b):

    A = np.matrix(Arr)
    y = np.matrix( b)
#    print type(A), type(y)
    
 #   print ("In  LS_analytic_fast, before   lstsq")
#    exit(8)
    x = np.linalg.lstsq(Arr,b)
#    print ("In  LS_analytic_fast, after   lstsq")
    
    return x
    exit(8)

    
#----------------------------------------------------------------------

def LS_analytic (Arr, e_lens):


    A = np.matrix(Arr)
    y = np.matrix( e_lens)
#    fmat = open("dro-mat-LS.dat", 'w')
#    print >> fmat , (A)
#    print >> fmat , "\n\n", (y)
#    print >> fmat , "\n\n", [math.exp(x) for x in e_lens]

#    fmat.close()
    
    n_rows, n_cols = A.shape

    u, s, v = np.linalg.svd(A)
    rank = sum(s > 1e-10)
    if rank < n_cols:
        WARN("Matrix rank %d is smaller than #cols %d and hence singular." % \
            (rank, n_cols))
        return None

#    print "Rank ", rank, "\nv", v, "\ns", s


#    print "The A matrix, rank (%d):\n" % (rank ), A
    #    exit(8)
    At = A.transpose()
#    print "\n\ntranspose \n", At
#    exit(8)
    A_At = np.dot(At,A)
#    print "multiply \n", A_At

    inv_A_At = LA.inv(A_At);
    #B = matrix ('2,5,7').transpose()
    #print "B", B
    #C = dot (A, B) 

    #print "C", C
    #C1 = matrix ('33,32,104, 60, 55, 59, 170').transpose()
    #exit(8)

#    print "After inversion\n", inv_A_At
    temp = np.dot(inv_A_At, At)

#    print "e_lens", e_lens, 
#    y = np.matrix( e_lens)
    # y = matrix( [1, 2, 4, 5])   # test 
    #y = [log( x) for x in [2, 4, 1, 2]]
    rows, cols = y.shape
#    print "y: ", y, "\n\ny[0]" , "\n\nexp:\n", [math.exp(y[0,i]) for i  in range(cols)]
#    exit(8)
    yt = (y).transpose()
#    print "yt:", yt


    beta = np.dot(temp, yt)

#    print "sol, beta:\n", beta, [math.exp(x) for x in beta]
#    exit(8)
    return beta




#-----------------------------------------------------------------------------------_#
def process_gtrees(GTs, ST_leaves, ST_newick):

    ST_no_len = re.sub(r':[-+]?\d*(.\d+)?',r'',ST_newick) #remove length
    gt_edges = []
    GTs_sum_edges = 0
    ST_sp_list = ST_leaves.keys()
    for gt in GTs:
        print "Gene tree read %s\n" % (gt)
        gt = re.sub(r'_\w+',r'',gt) # remove name suffix for Drosophila
        print "Gene tree read %s\n" % (gt)
        if not re.match(r'\s*(([(),]+[\w:\.]+)+\)+;?)',gt):
           ERR("Invalid tree:"+line) 
        gr =  re.match(r'\s*(([(),]+[\w:\.]+)+\)+;?)',gt).groups()
        #    print "\n\nTree:", line,"\n" , gr
        #    for i in range (len(gr)):
        #        print "\n\n",i, gr[i],"\n\n\n"
        #        pass
        #    exit(8)

        print gr
        GT = re.sub(r'\)\d+',r')', gr[0]) #remove bootstrap value
        print GT
#        exit(8)
        temp =  re.sub(r':[-+]?\d*(.\d+)?',r'',GT) #remove length
        temp =  re.sub(r';',r'',temp) #remove ;
        GT_sp_list = re.findall('[^,()]+',temp)
        edges = 2*len(GT_sp_list) - 3
        GTs_sum_edges += edges
        if set(GT_sp_list) > set(ST_sp_list):
            ERR("Not all species in GT are in ST %s" %(str(sorted(set(GT_sp_list)-set(ST_sp_list)))))

        print "GT:", GT

        GT_leaves = {};
        GT =  re.sub(r';',r'',GT) #remove ;
        TreeNode.count = 0
        GT_root = insertSubTree_n(GT, GT_leaves)
#        rescale_edge_len(GT_root) 
        if len(GT_root.children) < 3:
            ERR("GT %s is binary" % (GT))
            unroot_tree(GT_root)
            print "back from unroot"
            if not (GT_root.son.brother.brother):
                ERR("Tree remaind binary after unrooting");
            GT = printTree(GT_root,"")

    #    exit(8)
        GT_ck =  printDetailedTree(GT_root, "")

    #    dist = obtain_RF(GT_ck, GT)
    #    if dist != 0:
    #        ERR("NOT a zero dist between trees")

        GT_sp_list =  re.findall('[^,()]+',re.sub(r':[-+]?\d*(.\d+)?',r'',GT))
        super_sub_tree = obtain_subtree(ST_no_len, GT_sp_list)

        #        dist = obtain_RF(super_sub_tree, re.sub(r':[-+]?\d*(.\d+)?',r'',GT))
        dist  =0 # debug
        if dist != 0:
            ERR("NOT a zero dist between trees")
#        exit(8)

    #    r_i =   math.exp(np.random.normal())
    #    r_i = i+1.   # 
    #    print "original gene specific rate", r_i
    #    assign_gene_rate(con_GT_root, r_i, 0.5)
    #    assign_gene_rate(con_GT_root, 1, 0) # just for test
    #    GT_array.append((GT_root,GT_leaves))

        map_edges(GT_root, ST_leaves)
        ck_edge_mapping(ST_leaves, GT_leaves)
        # For the Drosophila we sheck that all ST edges were mapped to
        st_me = []
        mapped_edges = []
        getMappedEdges(GT_root, mapped_edges)
        print "GT edges", mapped_edges
#        exit(8)
        for m_e in mapped_edges:
            #print "m_e:", m_e[1]
            st_me += m_e[1]
            #print st_me
        if set (st_me) != set(range(len(st_me))):
            pass
#            ERR("Not all edges in ST were maped")

        else:
            print "All edges in ST were mapped"  
        gt_edges.append(mapped_edges)
        

    print "\n\n\n\nFinished generating Trees:\n\n\n%s\n\n\n" \
        %(". \n".join([str(x) for x in gt_edges]))
#    exit(8)
    return gt_edges
    
#-----------------------------------------------------------------------------------_#

def exceptional_edges(gtrees, betas, rs):
    
    n_t = len(gtrees)
    print "\n\nin Clac Err. n_t\n\n\n ", n_t,"betas:", betas, \
          "\nr's:", rs
#    exit(8)
    tot_err = 0.
    tot_edges = 0
    for i in range(n_t):
        #        print "\n\nAnalysing tree %d with rate %f"%(i, rs[i]),"\ntree:", gtrees[i]
        if tot_err == np.inf:
            ERR("tot_err arrived at inf")
        tree_err = 0.
        gt_edges = gtrees[i]
        for j in range(len(gt_edges)):
            GT_e_len = gt_edges[j][0]
            mapped_path = gt_edges[j][1]
            path_len = sum([betas[i1] for i1 in mapped_path])
            mapped_path_s = [[str(i1), str(betas[i1]), str(math.log(betas[i1]))] for i1 in mapped_path]
#            print "edge %d in tree %d has length %3.3f and mapped edges %s with tot len %f" % \
#                (j, i, GT_e_len, str(mapped_path), path_len)
            
#            exit(8)
#            print "Analsing edge ", j, ": ",gtrees[i][j] 
            if GT_e_len <= 0:
                #                print "Invalid value for GT_e_len", GT_e_len
                GT_e_len = 0.00001
                tot_err += (np.inf)
                break
#                exit(8)
            if path_len <= 0:
                print "Invalid value for path_len", path_len
                tot_err += (np.inf)
                break
                exit(8)
            if rs[i] <= 0:
                print "Invalid value for rs[i] ", rs[i]
                tot_err += (np.inf)
                break
                exit(8)
#            print "finding err for edge %d: %f, with beta %f and rate %f: %f, (w/o log: %f)" %\
#                  (j, GT_e_len, path_len, rs[i], math.log(GT_e_len/(path_len * rs[i])) , \
#                  GT_e_len/(path_len * rs[i]))

            expect_len = path_len * rs[i]
            p_val = .05
            cdf =  math.exp(-GT_e_len/expect_len)
            print "p-val for : path_len %f, rate %f, expct %f, obsrv %f, cdf %f \n\n" % \
                    (path_len , rs[i],  path_len * rs[i], GT_e_len,  cdf)
            if cdf < p_val:
                print "EXCEPTIONAL edge len: path_len %f, rate %f, expct %f, obsrv %f, p-val %f\n\n" % \
                    (path_len , rs[i],  path_len * rs[i], GT_e_len,  cdf)
#                exit(8)


#----------------------------------------------------------------------


def calc_error(gtrees, betas, rs):
    
    n_t = len(gtrees)
#    print "\n\nin Clac Err. n_t\n\n\n ", n_t,"betas:", betas, \
#          "\nr's:", rs
    #    exit(8)
    tot_err = 0.
    tot_edges = 0
    for i in range(n_t):
        #        print "\n\nAnalysing tree %d with rate %f"%(i, rs[i]),"\ntree:", gtrees[i]
        if tot_err == np.inf:
            ERR("tot_err arrived at inf")
        tree_err = 0.
        gt_edges = gtrees[i]
        for j in range(len(gt_edges)):
            GT_e_len = gt_edges[j][0]
            mapped_path = gt_edges[j][1]
            path_len = sum([betas[i1] for i1 in mapped_path])
            mapped_path_s = [[str(i1), str(betas[i1]), str(math.log(betas[i1]))] for i1 in mapped_path]
#            print "edge %d in tree %d has length %3.3f and mapped edges %s with tot len %f" % \
#                (j, i, GT_e_len, str(mapped_path), path_len)
            
#            exit(8)
#            print "Analsing edge ", j, ": ",gtrees[i][j] 
            if GT_e_len <= 0:
                #                print "Invalid value for GT_e_len", GT_e_len
                GT_e_len = 0.00001
                tot_err += (np.inf)
                break
#                exit(8)
            if path_len <= 0:
                print "Invalid value for path_len", path_len
                tot_err += (np.inf)
                break
                exit(8)
            if rs[i] <= 0:
                print "Invalid value for rs[i] ", rs[i]
                tot_err += (np.inf)
                break
                exit(8)
#            print "finding err for edge %d: %f, with beta %f and rate %f: %f, (w/o log: %f)" %\
#                  (j, GT_e_len, path_len, rs[i], math.log(GT_e_len/(path_len * rs[i])) , \
#                  GT_e_len/(path_len * rs[i]))

            residue = (math.log(GT_e_len/(path_len * rs[i])))**2
#            print "mapped_path_s: ", mapped_path_s
#            print "Error on path (%d) with len %f (gene len %f): %f. Exp: %f. New len %f. Adding %f\n" % \
#              (len(mapped_path), math.log(path_len), math.log(GT_e_len), math.log(GT_e_len/(path_len * rs[i])), \
#               (GT_e_len/(path_len * rs[i])), GT_e_len, residue)

            tree_err+= residue
            tot_edges +=1
            #        print "Error for tree %f\n" %(tree_err)
        tot_err += tree_err

#    print "upon returning from calc error. tot_err: %f, avg err %f" % ( tot_err, tot_err/ tot_edges  )
    if tot_err != tot_err: # a check for NaN
        ERR ("Value NaN for init_err returned from clac_error %f" %(tot_err))
    if tot_edges != tot_edges: # a check for NaN
        ERR ("Value NaN for tot_g_edges %f returned from clac_error" %(tot_edges))
      
    return tot_err, tot_edges     

#----------------------------------------------------------------------


def  obtain_rates_from_edges(betas, gtrees):

    n_t = len(gtrees)
    rs =  [0] * (n_t)
#    print "rs (%d):" %(len(rs)), rs
    #    exit(8)
    for i in range(0,n_t):
        #        print "\n\nAnalysing tree ", gtrees[i]
        gt_edges = gtrees[i]
        sum_err = 0.
        n_gt_edges = len(gt_edges)
        for j in range(n_gt_edges):
            GT_e_len = gt_edges[j][0]
            mapped_path = gt_edges[j][1]
            path_len = sum([betas[i1] for i1 in mapped_path])
            if (path_len) <= 0.:
                err += (np.inf)
                return err
                ERR("Path len: %f"% (path_len))
            #            print "edge %d in tree %d has length %3.3f and mapped edges %s with tot len %f" % \
            #                (j, i, GT_e_len, str(mapped_path), path_len)
            sum_err += (math.log(GT_e_len) - math.log(path_len))
        rs[i] = math.exp(sum_err/n_gt_edges)

		#    print "returning rs:", rs
    return rs


#----------------------------------------------------------------------



def findconso(n, consts):
    print "in findconso, node id", n.id

    if n.son != None:
        #        print "Node ", n.id, "has son ",n.son.id,". continue recursing"
        if n.name != None:
            print "ERROR ERROR ERROR. Node ", n.id, "has son AND a name."
            exit(8)
        findconso(n.son, consts)
            #        print "Back at internal node ", n.id, " after visiting all sons"
        setattr(n,'height',n.son.height + n.son.edge_len)    # Set len up to the node
    else:
        #        print "Node ", n.id, "has NO son. Must be a leaf"
        if n.name == None:
            print "ERROR ERROR ERROR. Node ", n.id, "has NO son AND NO name."
            exit(8)
            #        print "\nAt a leaf", n.id, " with name ", n.name
        setattr(n,'height',0)    # Set len up to the node
            
    if n.brother != None:
        #        print "Node ", n.id, "has brother. ",n.brother.id,".continue recursing. "
        findconso(n.brother, consts)
        #        print "Back at node ", n.id, " at height ",  n.height, \
        #            " after returning from brother ", n.brother.id, " at height ",  n.brother.height
        #        print "\n\nadding constraint n.height+n.edge_len  - n.brother ", \
        #            n.height+n.edge_len - (n.brother.height + n.brother.edge_len)
        consts.append( n.height+n.edge_len - (n.brother.height + n.brother.edge_len))
    print "Returning from findconso:", consts
    return consts

#--------------------------------------------------------------------------------------------#

def findConstraints(n, consts):

    if n.son != None:
        #        print "Node ", n.id, "has son ",n.son.id,". continue recursing"
        print "in findConstraint for internal, node id", n.id, consts
        if n.name != None:
            print "ERROR ERROR ERROR. Node ", n.id, "has son AND a name."
            exit(8)
        findConstraints(n.son, consts)
#        print "Back at internal node ", n.id, " after visiting all sons"
        setattr(n,'height',n.son.height + n.son.edge_len)    # Set len up to the node
        if n.son.brother == None:
            ERR("A single son to node %d" %(n.id))
        n = n.son
        while n.brother:
            print "Now finding constraints to brother ", n.brother.id
            findConstraints(n.brother, consts)
            print "Adding constraint for %d and its brother %d: %f" % \
                  ( n.id, n.brother.id, n.height + n.edge_len - \
                    (n.brother.height + n.brother.edge_len))
            consts.append( n.height + n.edge_len - \
                    (n.brother.height + n.brother.edge_len))
            n = n.brother
    else:
        #        print "Node ", n.id, "has NO son. Must be a leaf"
        print "in findConstraint for a leaf, node id %d, name %s" %(n.id, n.name),  consts
        if n.name == None:
            print "ERROR ERROR ERROR. Node ", n.id, "has NO son AND NO name."
            exit(8)
#        print "\nAt a leaf", n.id, " with name ", n.name
        setattr(n,'height',0)    # Set len up to the node

    print "Returning from findConstraints at node %d with height %f and e_len %f:" % \
          (n.id, n.height, n.edge_len), consts
            

#--------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------#

def findConstraints_n(n, consts):

    if n.children:
        #        print "Node ", n.id, "has son ",n.son.id,". continue recursing"
        print "in findConstraint_n for internal, node id", n.id, consts
        if n.name:
            print "ERROR ERROR ERROR. Node ", n.id, "has son AND a name."
            exit(8)
        if len( n.children < 2):
            ERR("Num of children <2")
        cld0 = n.children[0]
        findConstraints_n(cld0, consts)
        print "Back at internal node ", n.id, " after visiting all sons"
        setattr(n,'height',cld0.height + cld0.edge_len)    # Set len up to the node
        for cld1 in (n.children[1:]):
            print "Now finding constraints to brother ", cld1.id
            findConstraints_n(cld1, consts)
            print "Adding constraint for %d and its brother %d: %f" % \
                  ( cld0.id, cld1.id, cld0.height + cld0.edge_len - \
                    (cld1.height + cld1.edge_len))
            consts.append( cld0.height + cld0.edge_len - \
                    (cld1.height + cld1.edge_len ))
    else:
        #        print "Node ", n.id, "has NO son. Must be a leaf"
        print "in findConstraints_n for a leaf, node id %d, name %s" %(n.id, n.name),  consts
        if n.name == None:
            print "ERROR ERROR ERROR. Node ", n.id, "has NO son AND NO name."
            exit(8)
#        print "\nAt a leaf", n.id, " with name ", n.name
        setattr(n,'height',0)    # Set len up to the node

    print "Returning from findConstraints_n at node %d with height %f and e_len %f:" % \
          (n.id, n.height, n.edge_len), consts
            

#--------------------------------------------------------------------------------------------#


def check_binary_tree(root1):


    if (root1.father):
        ERR("Root %d has a father" % root1.id)
    if (root1.edge_len):
        ERR("Root %d has a edge len %f" % root1.id, root1.edge_len) 


    pstack = [root1]
        
    while pstack:
        n = pstack.pop()
#            print "Now at node id %d " % (n.id)
        if n.children:
            #                print "Node %d is internal" % n.id
            if n.name:
                ERR("Node %d has son %d and a name %s simultaneously" % \
                            (n.id, n.children[0].id, n.name))
            if len( n.children < 2):
                ERR("Num of children <2")

            if len( n.children > 2):
                ERR("Num of children > 2")

            for n1 in (n.children):
                if n1.father != n:
                    ERR("Father %d has son %d that does not point on it" % \
                            (n.id, n1.id))
                pstack.append(n1)

        else:
            print "arrived at leaf %s" %n.name
            if not n.name:
                ERR("Leaf %d has no name " %n.id)





#--------------------------------------------------------------------------------------------#
def ckTree(root1):

    if (root1.father):
        ERR("Root %d has a father" % root1.id)
    if (root1.edge_len):
        ERR("Root %d has a edge len %f" % root1.id, root1.edge_len) 


    pstack = [root1]
        
    while pstack:
        n = pstack.pop()
#            print "Now at node id %d " % (n.id)
        if n.son:
            #                print "Node %d is internal" % n.id
            if n.name:
                ERR("Node %d has son %d and a name %s simultaneously" % \
                            (n.id, n.son.id, n.name))
            n1 = n.son
            while n1:
                if n1.father != n:
                    ERR("Father %d has son %d that does not point on it" % \
                            (n.id, n.son.id))
                pstack.append(n1)
                n1 = n1.brother

            else:
#                print "arrived at leaf %s" %n.name
                if not n.name:
                    ERR("Leaf %d has no name " %n.id)


#--------------------------------------------------------------------------------------------#

def find_MC_constraints(root1, scale):

        pstack = [root1]
        cnt = 0
        while pstack:
            n = pstack.pop()
            print "in numberTreeEdges, node id", n.id
            if n.father:
                print "rescaling edge for node ", n.id, " : ", cnt
                if not n.edge_len:
                    n.edge_len = 0.000013
                    
                if ( n.edge_len < 0.000013):
                    n.edge_len = 0.000013

                n.edge_len *= scale
            
            if n.son:
                print "Node %d is internal" % n.id
                print "rout at leaf ", n.name, "set route to". n.root_from_leaf
                n1 = n.son
                setattr(n,'root_from_leaf', n1,root_from_leaf + set(n1.edge_num))    # Set route to emptyset

                while n1:
                    pstack.append(n1)
                    n1 = n1.brother

            else:
                setattr(n,'root_from_leaf',set([]))    # Set route to emptyset
                print "arrived at leaf %s" %n.name
                print "rout at leaf ", n.name, "set route to". n.root_from_leaf


#----------------------------------------------------------------------
 



def numberTreeEdges(root1):

        pstack = [root1]
        cnt = 0
        while pstack:
            n = pstack.pop()
            print "in numberTreeEdges, node id", n.id
            if n.father:
                print "Numbering edge for node ", n.id, " : ", cnt
                setattr(n,'edge_num',cnt)
                cnt += 1
            if n.children:
                print "Node %d is internal" % n.id
                for n1 in  n.children:
                    pstack.append(n1)

            else:
                print "arrived at leaf %s" %n.name

#        exit(8)

#--------------------------------------------------------------------------------------------#

def getMappedEdges(n, edges):
    print "in getMappedEdges, node id", n.id
    if n.father:
        if not n.mapped_edges:
            ERR("No edges are mapped to edge at node %d" % n.id)

        print "getting mapped edges for node ", n.id, " : ", n.mapped_edges
        if n.edge_len <= 0.0:
            #            ERR("Invalid edge len %f for gene" % (n.edge_len))
            edges.append([0.00001, n.mapped_edges])
        else:
            edges.append([n.edge_len, n.mapped_edges])
        
    if n.son:
        print "going to son ", n.son.id
        getMappedEdges(n.son, edges)

    else:
        print "At a leaf ", n.name
        
    if n.brother:
        print "going to brother ", n.brother.id
        getMappedEdges(n.brother, edges)
#--------------------------------------------------------------------------------------------#



def getTreeEdges(root1):

    pstack = [root1]
    edges = []
    ind = 0
    while pstack:
        n = pstack.pop()
        print "in getTreeEdges, node id", n.id
        if n.father:
            edges.append(n.edge_len)
            if n.edge_num != ind:  # see if we track the same path when edges length was assigned
                ERR("Inconsistency in edge numbers. " \
                    "Expected %d but found %d" % \
                    (ind,  n.edge_num))
            ind += 1
        if n.son:
            print "Node %d is internal" % n.id
            n1 = n.son
            while n1:
                pstack.append(n1)
                n1 = n1.brother

        else:
            print "arrived at leaf %s" %n.name
            pass
    
    return edges


#--------------------------------------------------------------------------------------------#


def assignTreeEdges(root1, edges):
#    print "\n\nin assignTreeEdges, root id", root1.id

    pstack = [root1]
    cnt = 0
    while pstack:
        n = pstack.pop()
#        print "in assignTreeEdges, node id", n.id
        if n.father:
            if n.edge_num != cnt:  # see if we track the same path when edges length was assigned
                ERR("Inconsistency in edge numbers. " \
                    "Expected %d but found %d" % \
                    (cnt,  n.edge_num))
            n.edge_len = edges[cnt]
            cnt += 1
        if n.son:
#            print "Node %d is internal" % n.id
            n1 = n.son
            while n1:
                pstack.append(n1)
                n1 = n1.brother

        else:
#            print "arrived at leaf %s" %n.name
            pass
    
    if cnt != len(edges):
        ERR("Error in assignTreeEdges. Not all edges assigned. Only %d out of %d" % \
                (cnt,len(edges)) )
 
#--------------------------------------------------------------------------------------------#


def get_mapped_route(leaves, x, y):
    if not x in leaves.keys():
        ERR ("Leaf "+x+" is not in the leaves of the tree")

    if not y in leaves.keys():
        ERR ("Leaf "+y+" is not in the leaves of the tree")


    print "\n\n\nfinding mapped edges between "+x+" and "+y


    lca  = LCA(leaves[x],leaves[y])
    if not (lca):
        ERR("Could not find LCA from %s to %s"%(x,y))

    n = leaves[x]
    print "augmenting path, starting from %s, id %d toward lca %d"% (n.name, n.id, lca.id)
    mapped_lst1 = []
    while (n != lca): 
        print "\nat node id %d,  " % (n.id, ) ,
        mapped_lst1 += n.mapped_edges
        if not (n.father):
            ERR("Arrived at root %d but not to LCA %d" %(n.id, lca.id))
        n = n.father
        

    print "Stopped at the lca %d, with mapped edges " % n.id, mapped_lst1

    n = leaves[y]
    print "Now augmenting path, starting from %s, id %d toward lca %d"% (n.name, n.id, lca.id)
    mapped_lst2 = []
    while (n != lca): 
        print "\nat node id %d,  " % (n.id) ,
        mapped_lst2 += n.mapped_edges
        if not (n.father):
            ERR("Arrived at root %d but not to LCA %d" %(n.id, lca.id))
        n = n.father
        

    print "Stopped at the lca %d, with  mapped edges " % n.id, mapped_lst2






    if set(mapped_lst1) & set(mapped_lst2) != set([]):
        print "INtersection is not empty ",mapped_lst1, mapped_lst2
        ERR("INtersection is not empty at node %d" %n.id)

    print "reached LCA %d. " % lca.id, ". Path is ", mapped_lst1+ mapped_lst2
    return  mapped_lst1 +  mapped_lst2




#-----------------------------------------------------------------------------------_#
def ck_edge_mapping(super_leaves, sub_leaves):

    sub_sp_list = sub_leaves.keys()
    for i in range(len(sub_sp_list)):
        x = sub_sp_list[i]
        for j in range(i+1,len(sub_sp_list)):
            y = sub_sp_list[j]
            print "Checking rout between leafs %s and %s" % (x,y)
            lca  = LCA(sub_leaves[x],sub_leaves[y])
            n = sub_leaves[x]
            sub_edge_lst1 = []
            while (n != lca): 
#                print "\nat node id %d,  edge num %d" % (n.id, n.edge_num) ,
                sub_edge_lst1 += (n.mapped_edges)
                if not (n.father):
                    ERR("Arrived at root %d but not to LCA %d" %(n.id, lca.id))
                n = n.father
            
            n = sub_leaves[y]
            sub_edge_lst2 = []
            while (n != lca): 
#                print "\nat node id %d,  edge num %d" % (n.id, n.edge_num) ,
                sub_edge_lst2 += (n.mapped_edges)
                if not (n.father):
                    ERR("Arrived at root %d but not to LCA %d" %(n.id, lca.id))
                n = n.father
 
            if set(sub_edge_lst1) & set(sub_edge_lst2) != set([]):
                print "INtersection is not empty ",sub_edge_lst1, sub_edge_lst2
                ERR("INtersection is not empty at node %d" %n.id)
            sub_edge_lst =  sub_edge_lst1 + sub_edge_lst2
            sup_edge_lst = get_edge_list(super_leaves, x, y)
            print "edge list between %s and %s in ST is:\n %s and in GT:\n %s\n"% \
                (x, y, str(sup_edge_lst), str (sub_edge_lst))
            if  set(sup_edge_lst) != set(sub_edge_lst):
                ERR("Error in edge maping between %s and %s" %(x, y))



#-----------------------------------------------------------------------------------_#

def get_edge_list(leaves, x, y):
    if not x in leaves.keys():
        ERR ("Leaf "+x+" is not in the leaves of the tree")

    if not y in leaves.keys():
        ERR ("Leaf "+y+" is not in the leaves of the tree")


    print "\n\n\nfinding edges between "+x+" and "+y, leaves[x].edge_num


    lca  = LCA(leaves[x],leaves[y])
    if not (lca):
        ERR("Could not find LCA from %s to %s"%(x,y))

    n = leaves[x]
    print "augmenting path, starting from %s, id %d toward lca %d"% (n.name, n.id, lca.id)
    edge_lst1 = []
    while (n != lca): 
        print "\nat node id %d,  edge num %d" % (n.id, n.edge_num) ,
        edge_lst1.append(n.edge_num)
        if not (n.father):
            ERR("Arrived at root %d but not to LCA %d" %(n.id, lca.id))
        n = n.father
        

    print "Stopped at the lca %d, with path " % n.id, edge_lst1

    n = leaves[y]
    print "Now augmenting path, starting from %s, id %d toward lca %d"% (n.name, n.id, lca.id)
    edge_lst2 = []
    while (n != lca): 
        print "\nat node id %d,  edge num %d" % (n.id, n.edge_num) ,
        edge_lst2.append(n.edge_num)
        if not (n.father):
            ERR("Arrived at root %d but not to LCA %d" %(n.id, lca.id))
        n = n.father
        

    print "Stopped at the lca %d, with path " % n.id, edge_lst2






    if set(edge_lst1) & set(edge_lst2) != set([]):
        print "INtersection is not empty ",edge_lst1, edge_lst2
        ERR("INtersection is not empty at node %d" %n.id)

    print "reached LCA %d. " % lca.id, ". Path is ", edge_lst1 + edge_lst2
    return edge_lst1 + edge_lst2


#-----------------------------------------------------------------------------------_#
def map_edges(n, ST_leaves):
    print "at MAP EDGES in node ",n.id
    if (n.children):
        print "Node ", n.id, "has son ",n.son.id," and is internal."
        if (n.name):
            ERR("ERROR ERROR ERROR. Node %d has children AND a name."%(n.id))
            exit(8)
        if  (n.father): 
            print "\n\nArrived at an internal node %d. Looking for a defining quartet\n" % \
                    n.id
            n1, n2, n3, n4 = get_defining_quartet(n)
            print "\n\nAfter get quartet. The defining leaves are %s, %s, %s, %s \n\n\n"%\
                (n1.name, n2.name, n3.name, n4.name)
            print "before paths"
#            exit(8)
            n1_n2_edges = get_edge_list(ST_leaves, n1.name, n2.name)

            print "\n\n%s - %s: " % (n1.name, n2.name), n1_n2_edges
#            exit(8)
            n1_n3_edges = get_edge_list(ST_leaves, n1.name, n3.name)
            n1_n4_edges = get_edge_list(ST_leaves, n1.name, n4.name)
            n2_n3_edges = get_edge_list(ST_leaves, n2.name, n3.name)
            n2_n4_edges = get_edge_list(ST_leaves, n2.name, n4.name)
            n3_n4_edges = get_edge_list(ST_leaves, n3.name, n4.name)
            
            print "\nThe qartet is %s,%s|%s,%s\n"% ( n1.name, n2.name, n3.name, n4.name)
            print "Paths are:\n"
            print "%s - %s" % (n1.name, n2.name), sorted(n1_n2_edges)
            print "%s - %s" % (n1.name, n3.name), sorted(n1_n3_edges)
            print "%s - %s" % (n1.name, n4.name), sorted(n1_n4_edges)
            print "%s - %s" % (n2.name, n3.name), sorted(n2_n3_edges)
            print "%s - %s" % (n2.name, n4.name), sorted(n2_n4_edges)
            print "%s - %s" % (n3.name, n4.name), sorted(n3_n4_edges)

                
            if (set(n1_n3_edges) | set(n2_n4_edges)) !=  (set(n1_n4_edges) | set(n2_n3_edges)):
                print (set(n1_n3_edges) | set(n2_n4_edges)), " versus\n", \
                    (set(n1_n4_edges) | set(n2_n3_edges))
                ERR("Problem in edge sets")
            if (set(n1_n3_edges) | set(n2_n4_edges)) - \
                    (set(n1_n4_edges) & set(n2_n3_edges)) != \
                    (set(n1_n2_edges) | set(n3_n4_edges)):
                print "ALL edges:" , sorted((set(n1_n3_edges) | set(n2_n4_edges)))
                print " central edges:", sorted(set(n1_n4_edges) & set(n2_n3_edges))
                print "Leaf edges:", sorted((set(n1_n2_edges) | set(n3_n4_edges)))
                ERR("Problem in edge sets")
                
                #            exit(8)
            mapped_edges = set(n1_n4_edges) & set(n2_n3_edges)
            print "mapped edges ", mapped_edges
            if len(mapped_edges) != 1:
                if n.father.father:
                    print "father is NOT the root"
#                    print "at node ", n.id , " with father ", n.father.id ," with father ", n.father.father.id, \
#                          " with father ", n.father.father.father.id 
                else:
                    print "father IS  the root"
                lca  = LCA(ST_leaves[n2.name],ST_leaves[n3.name])
                if lca.father:
                    print "LCA of ", n2.name, n3.name, "is NOT the root\n"
 #                   ERR("There should be a 1:1 mapping between edges, but found 1:%d" % (len(mapped_edges)))
                else:
                    print "LCA of ", n2.name, n3.name, " IS the root\n"
                    
                #            exit(8)

            setattr(n,'mapped_edges',mapped_edges) 
        else:
            print "Node %d is the root" % n.id

        map_edges(n.son, ST_leaves)
    else:
        print "Node ", n.id, "has NO son. Must be a leaf"
        if not (n.name):
            ERR ("ERROR ERROR ERROR. Node %dhas NO son AND NO name." %  n.id)

        print "\nAt a leaf %s (node id %d)" %( n.name, n.id)
        print "Now finding defining triplet\n"

        n1, n2 = get_defining_triplet(n)
        if (not n1) or (not n2):
            ERR("Get triplet for node "+str(n.id)+" returned NULL node")
        print "Nodes returned from get Triplets are", n1.id, n2.id
        if not (n1.name):
            ERR("node "+str(n1.id)+" has no name")
        if not (n2.name):
            ERR("node "+str(n2.id)+" has no name")
        print "Names of leaves retrurned from get triplet are ", n1.name, n2.name
        n_n1_edges = get_edge_list(ST_leaves, n.name, n1.name)
        n_n2_edges = get_edge_list(ST_leaves, n.name, n2.name)
        n1_n2_edges = get_edge_list(ST_leaves, n1.name, n2.name)
        print "\nThe triplet is %s,%s|%s\n"% ( n1.name, n2.name, n.name)
        print "Paths are:\n"
        print "%s - %s" % (n.name, n1.name), n_n1_edges
        print "%s - %s" % (n.name, n2.name), n_n2_edges
        print "%s - %s" % (n1.name, n2.name), n1_n2_edges
#        exit(8)

        if (set(n_n1_edges) ^ set(n_n2_edges)) !=  (set(n1_n2_edges) ):
            print set(n_n1_edges) ^ set(n_n2_edges), "versus:\n", (set(n1_n2_edges))
            ERR("Problem in edge sets")

        if (set(n_n1_edges) ^ set(n1_n2_edges)) !=  (set(n_n2_edges) ):
            print set(n_n1_edges) ^ set(n1_n2_edges) , "versus:\n", set(n_n2_edges)
            ERR("Problem in edge sets")
        
        setattr(n,'mapped_edges',set(n_n1_edges) & set(n_n2_edges))
 
#        exit(8)
           
    if (n.brother):
#        print "Node ", n.id, "has brother. ",n.brother.id,".continue recursing. "
        map_edges(n.brother, ST_leaves)

#-----------------------------------------------------------------------------------_#


def test_ieqcons(p,*args):
    """ Lefthandside of the inequality constraint """
    print "at test ieqcons"
    exit(8)
    gtrees = args[0]
    superT = args[1]
    n_e =  args[2] # no. of edges in the supertree
    remained_son = args[3]
    print len(gtrees), gtrees[:]
    n_t = len(gtrees) # no. of gene trees
    print   "n_t %d, n_e %d\n" % (n_t, n_e)
    print "remained_son", remained_son.id
#    exit(8)


    betas = p # The parameter is the ST edges the procedure  -
    if n_e != len(betas):
        ERR("n_e != len(betas): %d, %d" % (n_e, len(betas)))
    # the supertree edges

#    print "supertree edges (betas):", betas

    assignTreeEdges(superT, [x for x in betas])
    #    exit(8)
    
    consts = []
    n = superT.son
    if n == None:
        ERR("No son to root %d"% (superT.id))
        
    root_height = 0
    while n:
        print "Adding constraints at son of root:", n.id
        findConstraints(n,consts)

        print "\n\nNode %d is MC at height: %f, len to root %f "%(n.id, n.height, n.height+n.edge_len)
        if n != remained_son:
            if root_height:
                consts.append(root_height - (n.height+n.edge_len))
            else:
                print "\n\nsetting height %f\n" %(n.height+n.edge_len)
                root_height = n.height+n.edge_len
        n = n.brother

    print "remained_son %d" %( remained_son.id)
    print "\nRoot height %f remained_son.height %f, remained_son.edge_len %f\n" % \
        (root_height, remained_son.height, remained_son.edge_len)

    # We ignore all equality constrains and only output the two for
    # the root
    ieconsts = []
    ieconsts.append (root_height - remained_son.height + remained_son.edge_len)
    ieconsts.append (remained_son.height  - root_height + remained_son.edge_len)
    print "\n\nInequality Constraints:", ieconsts
    exit(8)
    return  np.array(ieconsts)



#-----------------------------------------------------------------------------------_#


def test_eqcons_new(p,*args):
    """ Lefthandside of the equality constraint """
    print "at test eqcons"
    gtrees = args[0]
    superT = args[1]
    n_e =  args[2] # no. of edges in the supertree
    remained_son = args[3]
    print len(gtrees), gtrees[:]
    n_t = len(gtrees) # no. of gene trees
    print   "n_t %d, n_e %d\n" % (n_t, n_e)
    print "remained_son", remained_son.id
#    exit(8)


    betas = p # The parameter is the ST edges the procedure  -
    if n_e != len(betas):
        ERR("n_e != len(betas): %d, %d" % (n_e, len(betas)))
    # the supertree edges

#    print "supertree edges (betas):", betas

    assignTreeEdges(superT, [x for x in betas])
    #    exit(8)
    
    consts = []
    n = superT.son
    if n == None:
        ERR("No son to root %d"% (superT.id))
        
    root_height = 0
    while n:
        print "Adding constraints at son of root:", n.id
        findConstraints(n,consts)

        print "\n\nNode %d is MC at height: %f, len to root %f "%(n.id, n.height, n.height+n.edge_len)
        if n != remained_son:
            if root_height:
                consts.append(root_height - (n.height+n.edge_len))
            else:
                print "\n\nsetting height %f\n" %(n.height+n.edge_len)
                root_height = n.height+n.edge_len
        n = n.brother

    print "remained_son %d" %( remained_son.id)
    print "\nRoot height %f remained_son.height %f, remained_son.edge_len %f\n" % \
        (root_height, remained_son.height, remained_son.edge_len)

    if root_height >= remained_son.height - remained_son.edge_len and \
        root_height <= remained_son.height + remained_son.edge_len :
        print "\n\n\nRemaind son is in OK height\n";
#        consts.append(0.0000000001)

    else:
        print "\n\n\nRemianed son is NOT in good height\n";
        if root_height < remained_son.height - remained_son.edge_len:
            print "remained_son is too high: %f for root height %f" % \
                      (remained_son.height - remained_son.edge_len, root_height)
            consts.append(remained_son.height - remained_son.edge_len - root_height)
        else:
            if  root_height > remained_son.height + remained_son.edge_len :
                print "remained_son is too low: %f for root height %f" % \
                      (remained_son.height + remained_son.edge_len, root_height)
                consts.append(root_height - remained_son.height - remained_son.edge_len)
            else:
                ERR("SMTH wrong with heights")

    t1 = printTree(superT, "")     
    print "\n\n\n\nBefore returning from constraints. betas:", betas

    print "\nTree:\n",t1
    print "\n\nConstraints:", consts
#    exit(8)
    return  np.array(consts)

#    return array([0])



#-----------------------------------------------------------------------------------_#


def test_eqcons(p,*args):
    """ Lefthandside of the equality constraint """
    print "at test eqcons"
    gtrees = args[0]
    superT = args[1]
    n_e =  args[2] # no. of edges in the supertree
    remained_son = args[3]
    print len(gtrees), gtrees[:]
    n_t = len(gtrees) # no. of gene trees
    print   "n_t %d, n_e %d\n" % (n_t, n_e)
    print "remained_son", remained_son.id
#    exit(8)


    betas = p # The parameter is the ST edges the procedure  -
    if n_e != len(betas):
        ERR("n_e != len(betas): %d, %d" % (n_e, len(betas)))
    # the supertree edges

#    print "supertree edges (betas):", betas

    assignTreeEdges(superT, [x for x in betas])
    #    exit(8)
    
    consts = findconso(superT,[])


    t1 = printTree(superT, "")     
    print "\n\n\n\nBefore returning from constraints. betas:", betas

    print "\nTree:\n",t1
    print "\n\nConstraints:", consts
#    exit(8)
    return  np.array(consts)

#    return array([0])



#----------------------------------------------------------------------

def ck_validity(superT, super_leaves, GT_array, betas):
    print "In ck_validity", len(betas)
    assignTreeEdges(superT, betas)

#    STN = printDetailedTree(superT,"")
#    stree_f = open ("stree.tre",'w')
#    stree_f.write(STN)
#    stree_f.close()

    sup_sp_list = super_leaves.keys()
    n = len (sup_sp_list)
    sp_dist = {}
    for i in range(n):
        for j in range(i+1, n):
            x = sup_sp_list[i]
            y = sup_sp_list[j]
            d_st = find_dist_between_leaves(x,y , super_leaves)
            sp_dist[x, y] = d_st
            print "In ck_validity, ST dist between %s and %s is %f3\n" % \
                (x,y,d_st)


    for [GT_root, GT_leaves] in GT_array:
        
        sub_sp_list = GT_leaves.keys()
        for i in range(len(sub_sp_list)):
            x = sub_sp_list[i]
            for j in range(i+1,len(sub_sp_list)):
                y = sub_sp_list[j]
                d_st = sp_dist[x, y]
                print "Checking rout between leafs %s and %s" % (x,y)
                mapped_route = get_mapped_route(GT_leaves, x, y)
                d_mapped = sum([betas[i1] for i1 in mapped_route])
                if abs(d_mapped - d_st) > .00000001:
                    ERR("inconsistent length between leaves %s and %s: %f vs %f, diff %f" %\
                            (x, y, d_mapped, d_st,d_mapped - d_st  ))
                else:
                    print "dist  between leaves %s and %s: %f vs %f, diff %f is OK" %\
                            (x, y, d_mapped, d_st,d_mapped - d_st  )


    print "\n\n\n\nLeaving ck_validity. All dists are OK!!!\n\n\n\n"
#    exit(8)
#-----------------------------------------------------------------------------------_#
    
#----------------------------------------------------------------------


def testfunc_new(p, *args):
  
    gtrees = args[0]
    superT = args[1]
    n_e =  args[2] # no. of edges in the supertree
    remained_son = args[3]
    
#    print len(gtrees), gtrees[:]
    n_t = len(gtrees) # no. of gene trees
#    print   "n_t %d, n_e %d\n" % (n_t, n_e)
    #    exit(8)


    betas = p # The parameter is the ST edges the procedure  -
    if n_e != len(betas):
        ERR("n_e != len(betas): %d, %d" % (n_e, len(betas)))
    # the supertree edges
    #    print "supertree edges (betas):", betas
    #    print "gene specific rates:", rs
    for b in betas:
        if b <= 0:
            return np.inf 
    

    #    exit(8)

    rs = obtain_rates_from_edges(betas, gtrees)
    if not rs:
        err = (np.inf)
        return err
#    print "Infered rs:", rs
    #    exit(8)
        
    err, n_gn_edges = calc_error(gtrees, betas, rs) 
    if err == None:
        ERR("unidentified error: None")

    if not err:
        ERR("unidentified error: None")

    import datetime
    now = datetime.datetime.now()
    start_tm = "%s-%s-%s_%s%s" %(now.year, now.month, now.day, now.hour, now.minute)
    print "%s_%s:%s" %( now.day, now.hour,now.minute), \
          "final total err: %f relative: %f "% (err, err/n_gn_edges)
#    exit(8)
    return err



#-----------------------------------------------------------------------------------_#


def test_ieqcons(p,*args):
    return p


#-----------------------------------------------------------------------------------_#



def testfunc(p, *args):
  
    gtrees = args[0]
    superT = args[1]
    n_e =  args[2] # no. of edges in the supertree
#    print len(gtrees), gtrees[:]
    n_t = len(gtrees) # no. of gene trees

#    print   "n_t %d, n_e %d\n" % (n_t, n_e)
#    exit(8)
    if len(p)  != n_t+ n_e:
        print "ERROR ERROR: unequal length of trees or r's"
        exit(8)
        
    betas = p[:n_e] # First n_e parameters in the input, are betas -
    # the supertree edges
    rs = p[n_e:]  # The last parameters in the input, are the specific
#    print "supertree edges (betas):", betas
#    print "gene specific rates:", rs

#    exit(8)



    err = 0.
    for i in range(n_t):
#        print "\n\nAnalysing tree ", gtrees[i]
        gt_edges = gtrees[i]
        for j in range(len(gt_edges)):
            GT_e_len = gt_edges[j][0]
            mapped_path = gt_edges[j][1]
            path_len = sum([betas[i1] for i1 in mapped_path])
#            print "edge %d in tree %d has length %3.3f and mapped edges %s with tot len %f" % \
#                (j, i, GT_e_len, str(mapped_path), path_len)
            
#            exit(8)
#            print "Analsing edge ", j, ": ",gtrees[i][j] 
            if GT_e_len <= 0:
#                print "Invalid value for GT_e_len", GT_e_len
                GT_e_len = 0.00001
#                err += (np.inf)
#                break
#                exit(8)
            if path_len <= 0:
#                print "Invalid value for path_len", path_len
                err += (np.inf)
                break
                exit(8)
            if rs[i] <= 0:
                print "Invalid value for rs[i] ", rs[i]
                err += (np.inf)
                break
                exit(8)
#            print "finding err for edge %d: %f, with beta %f and rate %f: %f, (w/o log: %f)" %\
#                (j, gtrees[i][j], betas[j], rs[i], log(gtrees[i][j]/(betas[j] * rs[i])) , 
#                 (gtrees[i][j]/(betas[j] * rs[i])))
            residue = (math.log(GT_e_len/(path_len * rs[i])))**2
            
            err+= residue
            
#    print "finall err:", err
    return err



#---------------------------------------------------------------------------------------#
def calc_error_old(orig_rs, result, n_trees):
    final_rho  = result[-n_trees:]
    betas = result[:-n_trees]
	#    print "final betas", betas
    print "final rho ",final_rho,  "\n\n\n", 
    print "orig rho", orig_rs

    sum_err = 0
    cnt_pairs = 0
    for i in range(n_trees):
        for j in range(i+1,n_trees):
            cnt_pairs +=1
            if final_rho[i] > final_rho[j]:
                print "final_rho(%d,%d)=%3.3f/%3.3f =%f, orig: rho(%d,%d)=%3.3f/%3.3f =%3.3f, err: %3.3f\n" % \
                    (i,j, final_rho[i], final_rho[j], final_rho[i]/final_rho[j], i,j, \
                         orig_rs[i],orig_rs[j], orig_rs[i]/orig_rs[j],\
                         abs(((final_rho[i]/final_rho[j]) - \
                                (orig_rs[i]/orig_rs[j]))/ (final_rho[i]/final_rho[j])))
            else: 
                print "final_rho(%d,%d) =%3.3f/%3.3f = %f, orig: rho(%d,%d)=%3.3f/%3.3f = %3.3f, err: %3.3f\n" % \
                    (j,i,  final_rho[j], final_rho[i],  final_rho[j]/final_rho[i], j,i, \
                         orig_rs[j],orig_rs[i], float(orig_rs[j])/orig_rs[i],\
                         abs(((final_rho[j]/final_rho[i]) - \
                                (orig_rs[j]/orig_rs[i]))/ (final_rho[j]/final_rho[i]))
    )

            sum_err += abs(((final_rho[j]/final_rho[i]) - \
                                (orig_rs[j]/orig_rs[i]))/ (final_rho[j]/final_rho[i]))


    print "\n\n\navg err %f\n\n" % (sum_err/(float(cnt_pairs)))
    
    return (sum_err/(float(cnt_pairs)))

#---------------------------------------------------------------------------------#
