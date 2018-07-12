#!/usr/bin/python

#!/opt/python-2.5/bin/python

##!/Library/Frameworks/Python.framework/Versions/Current/bin/python





import string, sys, re, copy, os,   getopt, time
from numpy import *
from scipy import *
from phylo import genAllTrees, addTaxonAtAllPositions
from tree import obtainMatchedParen,ERR, WARN
from time import *
import numpy
from subprocess import call

#
#-----------------------------------------------------------------------------------_#

class TreeNode_n:
    count = 0
    def __init__(self):
        self.id = TreeNode_n.count
        self.father = None
        self.children = None
        self.name = None
        self.edge_len = None
        self.btstrp = None
        TreeNode_n.count +=1
#        print '(Initialized node: %d)' % self.id
	
    def tell(self):
        '''Tell my details.'''
        print 'Name:"%s" Age:"%s"' % (self.name, self.age),

#
#-----------------------------------------------------------------------------------_#
def buildTrees(k):
    all_trees = [] 
    genAllTrees("(1,2,3)", 4, k, all_trees)
    roots_a = []
    leaves_a = []
    for tree_s in all_trees:
        print tree_s
        exit(8)
        leaves = {};
        root = insertSubTree_n(tree_s, leaves)
        roots_a.append(root)
        leaves_a.append(leaves)

    return roots_a, leaves_a


#-----------------------------------------------------------------------------------_#
"""
def prepareStack(n, roots_stk, leaves_stk, nodes_stk):
#    print "in prepareStack, node id", n.id
    if n.children:
#        print "Node ", n.id, "has children. continue recursing"
        if n.name:
            ERR( "ERROR ERROR ERROR. Node "+ n.id+ "has son AND a name.")
        for cld in n.children:
            if cld.children:
#                print "cld %d is internal. visit it"%(cld.id)
                k =  len(cld.children)
                if k > 2:
                    print "child %d of current node %d has %d sons. Opening an entry for it" %\
                        (cld.id, n.id, k)
                    roots_a, leaves_a = buildTrees(k)
                    
                    roots_stk.append(roots_a)
                    leaves_stk.append(leaves_a)
                    nodes_stk.append(cld)
                prepareStack(cld, roots_stk, leaves_stk, nodes_stk)
            else:
#               print "cld %d is a leaf. No need to visit it"%(cld.id)
                pass
 

#-----------------------------------------------------------------------------------_#
def genAllRefinements(level, n_levels, root, r_stk, l_stk, n_stk):
    n = n_stk[level]
#    print "In genAllRefinements, level %d, final level %d node %d" %      (level, n_levels,n.id)
    if level == n_levels:
        treen = printDetailedTree_n(root, "")
#        print "\n\nArrived at final level. Tree is: treen"
        return
    root_arr = r_stk[level]
    leaves_arr = l_stk[level]
    deg = len(root_arr)
    if n.father:
        f = n.father
        print "Node handled %d is not the root. Replacing it in its father %d" %\
          ( n.id, f.id)
        if not n in f.children:
            ERR("child %d is not in his father %d children list" %\
                (n.id, f.id))
    
        f.children.remove(n)
        if  n in f.children:
            ERR("child %d is still in his father %d children list" %\
                (n.id, f.id))
        for i in range(deg):
            t_leaves = leaves_arr[i]
            t_root = leaves_arr[i]
            t_s = printDetailedTree_n(t_root, "")
            print "Tree to replace:", t_s
            exit(8)

"""            

#-----------------------------------------------------------------------------------_#

def printTree_n(n, outS, detailed = False):
#    print "in printTree_n, node id", n.id
    if n.children:
#        print "Node ", n.id, "has son ",n.son.id,". continue recursing"
        if n.name:
            ERR( "ERROR ERROR ERROR. Node "+ n.id+ "has son AND a name.")
        print '(',
        outS += '('
        outS = printTree_n(n.children[0], outS, detailed)
        for cld in (n.children[1:]):
            outS += ','
            outS = printTree_n(cld, outS, detailed)
        else:
            outS += ')'
        if detailed:
             outS += "[%d]"%n.id
#        outS += ')'+"[%d]"%n.id
        if n.btstrp :
            print "Found bootstrap %f" % (n.btstrp)
            outS += "%5.8f"% (n.btstrp)
        if n.edge_len != None:
            outS += ":%5.8f"% (n.edge_len)
        print ')',
    else:
#        print "Node ", n.id, "has NO children. Must be a leaf"
        if n.name == None:
            print "ERROR ERROR ERROR. Node ", n.id, "has NO son AND NO name."
            exit(8)
#        print "\nAt a leaf", n.name
        outS += n.name
        if detailed:
             outS += "[%d]"%n.id
#        outS += n.name+"[%d]"%n.id
        if n.btstrp :
            ERR ("Found a btstrp %f. A leaf cannot have btstrp" %  (n.btstrp))
        if n.edge_len != None:
            outS += ":%5.8f"% (n.edge_len)
        print  n.name,
            
    
    return outS

#-----------------------------------------------------------------------------------_#
def contract_tree(n):
    print "at contract tree in node ",n.id
    if (n.children):
        print "Node ", n.id, "has son ",n.children[0].id," and is internal."
        if (n.name):
            ERR( "ERROR ERROR ERROR1. Node %d has son AND a name." % (n.id))
        children_n =[]
        for cld in  n.children:
           children_n.append(contract_tree(cld))
           
        n.children = children_n
        if len (n.children)  == 1:
            s = n.children[0]
            print "\nSingle Son, Single Son. node %d  has only one son %d  "\
                %(n.id,  s.id)
            if (s.edge_len):
                if (n.edge_len):
                    s.edge_len += n.edge_len
            else: 
                if (n.edge_len):
                    s.edge_len =  n.edge_len
            return s
    else:
        if not (n.name):
            ERR("Node %d has no son and also no name" %n.id)
        print "At a leaf %d with name %s" %(n.id, n.name)

    print "returning from %d"% n.id
    return n

#-----------------------------------------------------------------------------------_#



#-----------------------------------------------------------------------------------_#
def contract_low_bootstrap(n, thresh):
#    print "at contract low boootstrap in node ",n.id
    if (n.children):
#        print "Node ", n.id, "has son  and is internal."
        if (n.name):
            ERR( "ERROR ERROR ERROR1. Node %d has son AND a name." % (n.id))

        stk = []
        for cld in n.children:
            stk.append(cld)
        for cld in stk:
            contract_low_bootstrap(cld, thresh)

#        print "Back at contract low boootstrap, node %d"%n.id

        if n.btstrp and n.btstrp < thresh:
            print "Contracting a node %d with low bootstrap  %f"%(n.id, n.btstrp)
            if not n.father:
                ERR("root node should have bootstrap")

            f = n.father
            f.children.remove(n)
            f.children +=  n.children
            for cld in n.children:
                cld.father = f
        
        else:
            if n.btstrp:
#                print "Node %d has ok bootstrap %d" %(n.id, int(n.btstrp))
                pass
            else:
#                print "Node %d has no bootstrap" %(n.id)
                pass
    else:
        if n.btstrp:
            ERR("A leaf %d cannot have bootstrap value %f" %(n.id, n.btstrp))
#        print "Node %d is a leaf, we don't check bootstrap"% (n.id)

#-----------------------------------------------------------------------------------_#

#-----------------------------------------------------------------------------------_#

def printDetailedTree_n(n, outS):
#    print "in printTree_n, node id", n.id
    if n.children:
#        print "Node ", n.id, "has son ",n.son.id,". continue recursing"
        if n.name:
            ERR( "ERROR ERROR ERROR. Node "+ n.id+ "has son AND a name.")
        print '(',
        outS += '('
        outS = printDetailedTree_n(n.children[0], outS)
        for cld in (n.children[1:]):
            outS += ','
            outS = printDetailedTree_n(cld, outS)
        else:
            outS += ')'+"[%d]"%n.id
        if n.btstrp :
            print "Found bootstrap %f" % (n.btstrp)
            outS += "%5.8f"% (n.btstrp)
        if n.edge_len != None:
            outS += ":%5.8f"% (n.edge_len)
        print ')',
    else:
#        print "Node ", n.id, "has NO children. Must be a leaf"
        if n.name == None:
            print "ERROR ERROR ERROR. Node ", n.id, "has NO son AND NO name."
            exit(8)
#        print "\nAt a leaf", n.name
#        outS += n.name
        outS += n.name+"[%d]"%n.id
        if n.btstrp :
            ERR ("Found a btstrp %f. A leaf cannot have btstrp" %  (n.btstrp))
        if n.edge_len != None:
            outS += ":%5.8f"% (n.edge_len)
        print  n.name,
            
    
    return outS

#-----------------------------------------------------------------------------------_#
#-----------------------------------------------------------------------------------_#

def addSon_n(father, son):
#    print "in addSon_n for father %d and a new son %d\n" % \    (father.id,  son.id)
    if father.name != None:
        ERR("ERROR ERROR. Trying to insert a son to a leaf\n")
    if father.children:
        father.children.append(son)
    else:
        father.children = [son]
    son.father = father



#--------------------------------------------------------------------------------------------------------------


def insertSubTree_n(tree, leaves):
#    print "\nIn insertSubTree_n. Tree:"# , tree
    father =  TreeNode_n()
    if tree[0] != '(' or tree[-1] != ')' :
        ERR("tree does not start and end with ( and ):"+tree)
    ind = 1
#    print "entered insert subtree. Allocated a father with id ", father.id
#    print "subtree is ", tree
    while tree[ind] != ')':
#        print "looking at char %d: %c\n" % (ind,  tree[ind] )
        if tree[ind] == '(':
#            print "Found a ( at %d, inserting a new son \n"  % ind
            matched = obtainMatchedParen(tree, ind)
#            print "the matched \')\' is at ", matched
#            print "resulted subtree:"+tree[ind:matched+1]
            son =  insertSubTree_n(tree[ind : matched+1], leaves)
            ind = matched +1
#            print "returnd to subtree loop after finding son %d." \
 #               " char now is %c\n" % (son.id, tree[ind] )
            if  re.match(r'[-+]?(\d+(\.\d*)?|\.\d+)',tree[ind:]):
                btstrp = (re.match(r'[-+]?(\d+(\.\d*)?|\.\d+)',tree[ind:])).group(0)
#                print "Found bootstrap ", btstrp
                ind += len(btstrp) 
                son.btstrp = float(btstrp)
            if tree[ind] == ':':
#                print "Found edge length at pos %d" %(ind)
                if not re.match(r'[-+]?(\d+(\.\d*)?|\.\d+)',tree[ind+1:]):
                    ERR("Found \':\' after "+name+" but could not parse len at "+ \
                            tree[ind:])   
                edge_len = (re.match(r'[-+]?(\d+(\.\d*)?|\.\d+)',tree[ind+1:])).group(0)
#                print "Found edge len ", edge_len
                ind += len(edge_len) +1
                son.edge_len = float(edge_len)

        elif tree[ind]  == ',':
#            print "found , inserting son with root id %d " \
#	     " under father %d\n" % (son.id,  father.id)
            ind +=1
            addSon_n(father, son)
        elif re.match(r'(\w+)',tree[ind:]):
            name = re.match(r'(\w+)',tree[ind:]).group(0)
#            print "found a name ", name
            son = TreeNode_n()
            son.name =  name
            if name in leaves.keys():
                ERR ("name "+name + " appears twice in the leaves")
            leaves[name] = son
            ind += len(name)
            if tree[ind] == ':':
                if not re.match(r'[-+]?(\d+(\.\d*)?|\.\d+)',tree[ind+1:]):
                    ERR("Found \':\' after "+name+" but could not parse len at "+ \
                            tree[ind:])   
                edge_len = (re.match(r'[-+]?(\d+(\.\d*)?|\.\d+)',tree[ind+1:])).group(0)
#                print "Found edge len ", edge_len
                ind += len(edge_len) +1
                son.edge_len = float(edge_len)
        else:
            ERR ("Invalid char \" "+tree[ind]+ "\" at "+str(ind))
    
#    print "exited from insertSubtreeWithBtstrp loop. char is now %c\n" % tree[ind]
 
            
    addSon_n(father, son)
    return father


#-----------------------------------------------------------------------------------_#
def find_num_resolutions(root1, upBound, max_cld):
#     print "At find_num_resolutions node %d" % root1.id

     pstack = [root1]
     numRes = 1
     while pstack:
         n = pstack.pop()
#         print "in  find_num_resolutions , node id", n.id
         if not (n.children):
             ERR("node has no children")
         cnt = len(n.children)
         if cnt < 2:
             ERR("node has only one son")

         for cld in n.children:
             if cld.children:
                 pstack.append(cld)
             else:
#                 print "no sons" 
                 pass
#         print ("Node has %d sons"% cnt)
         if cnt > 2:
             if cnt > max_cld:
                 print "too many children to do the calculations: %d"  % (cnt)
                 return 100* upBound
             if n.father or cnt > 3:
                 localNumRes = factorial(2*cnt-3)/(2**(cnt-2)*factorial(cnt-2))
                 numRes *= localNumRes
                 print "number of sons %d, local resolutions %f, accumul #resolution %f"% \
                     (cnt, localNumRes, numRes)
                 if numRes > upBound:
                     print "arrived the upper bound: %d" % numRes
                     return numRes
             else:
                 print "This is the root node %d that has %d children" %( n.id, cnt)
                 pass

     return numRes


#-----------------------------------------------------------------------------------_#

def insert_element2queue (que, e_id, t):

    for i in range(len(que)):
#        print "Now at element ", i, " val ",  que[i]
        if t < que[i][1]:
#            print "inserting  the element before ",  que[i][1]
            que.insert(i, [e_id,t])
            break
    else:        
#        print "Inserting element at the end "
        que.append([e_id,t])

#    for i in range(len(que)):
#       print "Now at element ", i, " val ",  que[i]



#-----------------------------------------------------------------------------------_#

def yule_tree(n, rate):
    TreeNode_n.count = 0
    leaves = {}
    events = []
    root =  TreeNode_n()
    son1 = TreeNode_n()
    son1.father = root
    root.children = []
    son2 = TreeNode_n()
    son2.father = root
    root.children.append(son1)
    root.children.append(son2)

    next  = random.exponential(rate)  
    events.append([0,next])
    next  = random.exponential(rate)
    insert_element2queue(events,1,next)

    print events
#    exit(8)
    son1.edge_len = son2.edge_len = 0
    cnt = 2
    son1.name = str(0)
    leaves[son1.name] = son1
    son2.name =  str(1)
    leaves[son2.name] = son2
    print "name1 %s, name2 %s \n" % (son1.name , son2.name )
    now = 0

    for i in range(2,n):
        print "Now creating node leaf ", i, "now", now, "\nevents", events[i-2:]      
        [split_id, split_t] = events[i-2]
        delta = split_t - now
        print "Now is ", now, "next split in leaf ",split_id, "  at ", split_t, "delta ", delta
        lv_keys = leaves.keys()
        for k in lv_keys:
            leaves[k].edge_len += delta
            
        now = split_t
        split_n = leaves[str(split_id)]
        print "in Yule loop . now inserting ", i, \
            " spliting node", split_n.id, "at time ", split_t

#        This will be the new leaf node
        new_n = TreeNode_n()
        new_n.edge_len = 0
        new_n.name =  str(i)
        leaves[new_n.name] = (new_n)
        print "leaves size ", len(leaves)

#        This will be the new internal  node father to the new leaf 
        new_int = TreeNode_n()
        print "split node ", split_n.id, ", father ", split_n.father.id
        new_int.edge_len = split_n.edge_len
        grandpa = split_n.father 
        new_int.father = grandpa 
        if not split_n in grandpa.children:
            ERR("Could not find split node %d in granpa %id children" %(split_n.id, grandpa.id))
        grandpa.children.remove(split_n)
        grandpa.children.append(new_int)

#        if grandpa.son == split_n:
#            grandpa.son = new_int
#        else:
#            if grandpa.son.brother != split_n:
#                ERR("no link from father %d to son %d" %(grandpa.id, split_n.id))
#            else:
#                grandpa.son.brother = new_int

        new_int.children = [split_n]
        split_n.edge_len = 0
        split_n.father = new_int
        new_int.children.append(new_n)
        new_n.father = new_int 
        new_split = numpy.random.exponential(rate) 
        print "The bifurcating leaf ", split_id, " will split at time ", now + new_split
        insert_element2queue(events, split_id, now + new_split)

        new_split = numpy.random.exponential(rate) 
        print "The new leaf ", i, " will split at time ", now + new_split
        insert_element2queue(events, i, now + new_split)
        print "new father ", new_int.id, "new sons:", new_int.children[0].name, ", ", \
            new_int.children[1].name

        treen = printTree_n(root, "", detailed = True)
        print "\n\n\n", treen
    
    else:
        print "finished adding leaves. Now adding final length. now", now, "events", events[i-1:]      
        [split_id, split_t] = events[i-1]
        delta = split_t - now
        print "final delta ", delta
        lv_keys = leaves.keys()
        for k in lv_keys:
            leaves[k].edge_len += delta


    treen = printTree_n(root, "", detailed = True)

    print "\n\n\n", treen
    return root, leaves
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------_#
def unroot_tree(root):
    print "At unroot tree for node %d" % root.id
    if not (root.children):
        ERR("Root has no son")
    if len (root.children) < 2:
        ERR("Root as only one son")
    if len (root.children) > 2:
        print ("Root as more two sons. no need to unroot")
        return
    [son1,son2] = root.children
    print "son is %d, brother is %d" %( son1.id, son2.id)

    if  (son2.children):
        print "Son2 %d has sons %d and %d" %  (son2.id, son2.children[0].id, son2.children[1].id)
        if  len (son2.children) <2:
            ERR("Only one son to brother")
# We lift all son2's children as the root children

        if son2.edge_len:
            son1.edge_len += son2.edge_len
        for cld in son2.children:
            print "moving grandson %d to be root's son" %(cld.id)
            root.children.append(cld)
            cld.father = root
        root.children.remove(son2)
    else:
        print "Brother has NO sons"
        if  not (son1.children):
            print "Son1 has NO sons"
            return
        print "Son %d has sons"%son1.id
        if  len (son1.children)< 2:
            ERR("Only one son to son %d" % son1.id)

        if son1.edge_len:
            son2.edge_len += son1.edge_len
        for cld in son1.children:
            root.children.append(cld)
            cld.father = root
        root.children.remove(son1)

    print "Root sons are now",
    for cld in root.children:
        print cld.id,

#-----------------------------------------------------------------------------------_#


def check_MC_unrooted_tree(root, remained_son):
    root_height = 0
    for n in root.children:
        if not checkMCTree(n):
            print "\nNot MC at son %d"% (n.id)
            return 0
            break
        print "\n\nNode %d is MC at height: %f, len to root %f "%(n.id, n.height, n.height+n.edge_len)
        if n != remained_son:
            if root_height:
                if abs(root_height - (n.height+n.edge_len)) > 0.000001:
                    print "NON MC at node %d. at height: %f." % (n.id, n.height)
                    print "root_height %f, node height+ node.edge_len: %f, diff: %012.10f\n" % \
                          (root_height, n.height+n.edge_len, root_height - (n.height+n.edge_len))
                    return 0
                    break
            else:
                print "\n\nsetting height %f\n" %(n.height+n.edge_len)
                root_height = n.height+n.edge_len

    print "remained_son %d" %( remained_son.id)
    print "\nRoot height %f remained_son.height %f, remained_son.edge_len %f\n" % \
        (root_height, remained_son.height, remained_son.edge_len)

    if root_height >= remained_son.height - remained_son.edge_len and \
        root_height <= remained_son.height + remained_son.edge_len :
        print "\n\n\nMC MC MC\n";
        return 1
    else:
        print "\n\n\nNO MC, NO MC, NO MC\n";
        
        return 0

#-----------------------------------------------------------------------------------_#



def checkMCTree(n):
    print "in checkMCTree, node id", n.id
    if n.children:
        print "Node ", n.id, "has son ", n.children[0].id,". continue recursing"
        if n.name != None:
            print "ERROR ERROR ERROR. Node ", n.id, "has son AND a name."
            exit(8)
        print "Check MC at children"
        for son in n.children:
            if not checkMCTree(son):
                print "Son ", son.id, " is not MC"
                return 0
        print "All children are MC. Now check MC at node"
        son = n.children[0]
        setattr(n,'height',son.height + son.edge_len)  
        for  n0 in n.children:
            if abs(n0.height + n0.edge_len - n.height) > 0.0001:
                print "\n\n\nNon MC tree at node ", n.id, " with height %3.4f" % \
                    n.height, "and son with height ",   n0.height + n0.edge_len
                return 0
        print "Node ", n.id, " at height  %3.4f" % n.height, " is MC"
    else:
        print "Node ", n.id, "has NO son. Must be a leaf"
        if n.name == None:
            print "ERROR ERROR ERROR. Node ", n.id, "has NO son AND NO name."
            exit(8)
#        print "\nAt a leaf", n.name
        setattr(n,'height',0)    # Set height 0 for a leaf
    
    return 1



#-----------------------------------------------------------------------------------_#
def root_at_edge(n):

    print "In  root_at_edge for edge above node id %d" %(n.id)

    f = n.father  
    f.children.remove(n)
    newRoot = TreeNode_n()
    newRoot.children.append(n)
    newRoot.father = f
    n.father = newRoot
    push_to_root(newRoot)
    newRoot.father = None
    return newRoot


#-----------------------------------------------------------------------------------_#

def push_to_root(n):
#   print "In push to root for node id %d" %(n.id)


   f = n.father
   f.children.remove(n)
   while f:
#       print "In swap between son %d and father %d " %           (n.id,  f.id)
       g = f.father
       if g:
#           print "removing father %d from his father %d"%(f.id, g.id )
           g.children.remove(f)
#           print "After Remove son"
#           print "tree after remove son", printDetailedTree_n(root, "")
#            exit(8)
       else:
#           print "No father to %d, probably the root and we should stop" %f.id
           pass
           
#       print "Inserting the new son. Sons now are:", [x for x in n.children ] if n.children
       if n.children:
           n.children.append(f)
       else:
           n.children = [f]
       f.father = n
#       print "After inserting the new son. Sons now are:", [x.id for x in n.children]
       n = f
       f = g
       if f:
#           print "at end of iteration. F now is %d, n is %d" %(f.id, n.id)
           pass
       else:
#           print "last iteration - f is NULL"
           if len(n.children) == 1:
               print ("only one son to last node %d"%n.id), n.children[0].id, n.children[0].name
               n.father.children.remove(n)
               n.children[0].father = n.father
               n.father.children += n.children
               #exit(8)
        



#-----------------------------------------------------------------------------------_#

def genAllRootedTrees(nleaves):
    ncurrent = 3
    trees =  genAllTrees("(1,2,3)",ncurrent, nleaves)

    print "trees (%d)"%(len(trees)), trees
    rooted_trees = []
#exit(8)
    for tree in trees:
        rooted_trees += genAllRootings(tree)
    return rooted_trees


#-----------------------------------------------------------------------------------_#

def genAllRootings(tree):

    print "In genAllRootings for tree ", tree
    tree =  re.sub(r';',r'',tree ) #remove ;
    temp_t = tree
    temp_t = re.sub(r'\)[-+]?\d*(\.\d+)?\s*',r')', temp_t) #remove bootstrap value
    temp_t = re.sub(r':[-+]?\d*(\.\d+)?\s*$',r'',temp_t) #remove length at the end
    temp_t =  re.sub(r':[-+]?\d*(\.\d+)?(e-\d+)?',r'',temp_t) #remove length

    sp_list  = re.findall('[^,()]+',temp_t)
    sp_n = len(sp_list)

    rootings = []
    A = addTaxonAtAllPositions(tree,"sagiOut")
    if len(A) != 2*sp_n -3:
        print "ERROR in addTaxonAtAllPositions:"
        print "Input tree:", temp_t
        print "All insertions:", "\n\n".join(str(i) +")\t"+ A[i] for i in range(len(A)))
#    print "A:", A
    for t1 in A:
#        print "\nNow working with tree ", t1
        TreeNode_n.count = 0
        leaves = {};
        root = insertSubTree_n(t1, leaves)
#        check_binary_tree(root)
#        print "\nTree %s is binary and OK"%(t1)
        if sp_n != len(leaves)-1:
            print "\nleaves.keys():", leaves.keys()
            print "\n\nsp_list:", sp_list
            ERR("problem with  sp_n != len(leaves) -1 : %d != %d" %( sp_n, len(leaves)-1))
        n =  leaves["sagiOut"]
        if not (n):
            ERR("sagiOut was not found")
        push_to_root(n)
#        print "\nafter push:",  printTree_n(n,"")
        if len(n.children) != 1:
            ERR("strange number of children after push to root %d"%( len(n.children) ))
        newRoot = n.children[0]
        newRoot.father = None
        treen = printTree_n(newRoot, "")
#        print "\n\nTree before rooting ", (re.sub(r':[-+]?\d*(\.\d+)?(e-\d+)?',r'',(re.sub(r'\[\d+\]','', t1))))
#        print "\n\n\nTree after rooting ", treen, "\n",(re.sub(r':[-+]?\d*(\.\d+)?(e-\d+)?','',(re.sub(r'\[\d+\]','', treen))))
        check_binary_tree(newRoot)
#        exit(8)
        rootings.append(treen)

#    if len(rootings) != (2*sp_n - 3):
#        print "\nrootings (%d):" %(len(rootings)), # rootings
#        print "\n\nsp_list (%d):"%len(sp_list),  sp_list
#        WARN("Problem:  len(rootings) != (2*sp_n - 3): %d != %d"% ( len(rootings),  (2*sp_n - 3)))

    return rootings

#-----------------------------------------------------------------------------------_#

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
            if len(n.children) < 2:
                ERR("Only one sone to node %d" %(n.id))

            if len(n.children) > 2:
                if n.father:
                    ERR("More than two sons to a node %d that is NOT the root" %(n.id))
                else:
                    if len(n.children) > 3:
                        ERR("Root node %d has more than 3 children: %d"% (n.id, len(n.children) ))

            for cld in n.children:
                if cld.father != n:
                    ERR("Father %d has son %d that does not point on it" % \
                            (n.id, cld.id))
                pstack.append(cld)
        else:
#            print "arrived at leaf %s" %n.name
            if not n.name:
                ERR("Leaf %d has no name " %n.id)


#-----------------------------------------------------------------------------------_#
def LCA(n1, n2):
    # First stage: delete attr
    print "now starting from ", n1.id
    n = n1
    while n.father !=  None: 
#        print "Now at node ",n.id
        if hasattr(n,'lca'):
            delattr(n,'lca')
        n = n.father

#    print "stopped at node ", n.id
    if hasattr(n,'lca'):
        delattr(n,'lca')
        
#    print "now starting from ", n2.id
    n = n2
    while n.father != None: 
        setattr(n,'lca',1)
        n = n.father

#    print "stopped at node ", n.id
    setattr(n,'lca',1)

#    print " starting again from ", n1.id
    n = n1
    while n.father != None: 
#        print "Now at node ",n.id
        if hasattr(n,'lca'):
#            print "stopped at ", n.id
            return n
        n = n.father

#    print "stopped at node ", n.id
    if not hasattr(n,'lca'):
        ERR("At LCA. Arrived at a node w/o father (root?) but with no lca-mark")
    return n


#-----------------------------------------------------------------------------------_#
def get_defining_triplet(n):
    print "\nAt get_defining_triplet for leaf %s\n" % n.name
    if not (n.father):
        ERR("Node id "+str(n.id)+" has no father. Can't obtain triplet")
    f = n.father
# We first check if father has three sons        
    if len (f.children) < 2:
        ERR("Father %d has a single son %d" % (f.id, f.son.id))
    if len (f.children) > 2:
#        print "Three sons at node %d (possibly the root)\n" % f.id
#        exit(8)
        brothers = []
        for b in f.children:
            if b != n:
                print "found a son %d that is not me. Adding to brothers" % \
                    b.id
                brothers.append(b)
                if len(brothers) == 2:
#                    print "Two brother %s found. leaving" % str(brothers)
                    break
            else:
                print "Found myself: %d" %b.id
        if len(brothers) < 2:
            ERR("Did not find two bothers. brothers: %s" %str(brothers))
                
        print "getting descendants of both brothers %d, %d" % \
            (brothers[0].id, brothers[1].id )
        d1 = get_descendant(brothers[0])
        d2 = get_descendant(brothers[1])
        return d1, d2
            

    print "Father %d has only two sons: %d, %d" %     \
        (f.id, f.children[0].id, f.children[1].id)
    if f.children[0] != n and  f.children[1] != n:
        ERR("Non of father sons %d and %d is me %d" %  \
                (f.children[0].id, f.children[1].id, n.id))
    if f.children[0] == n:
        b  = f.children[1]
    else :
        b = f.children[0]

    if not (f.father):
        print "father ",f.id," has no father. Is a root"
        if not (b.children):
            ERR("No father to father (root?) %d and no son to brother %d", f.id, )
        g1 = b.children[0]
        if len (b.children) < 2: 
            ERR("Grandson "+str(g1.id)+" must have a brother ")
        g2 = b.children[1]
        print "Found two grandsons of root", g1.id, g2.id
        d1 = get_descendant(g1)
        d2 = get_descendant(g2)
        return d1, d2
    else:
        u =  get_brother(f)
        if u == None:
            ERR("No bother to father "+str(f.id))
        
        d1 = get_descendant(b)
        d2 = get_descendant(u)
        return d1, d2


#-----------------------------------------------------------------------------------_#

#-----------------------------------------------------------------------------------_#
def get_defining_quartet(n):
    print "\nAt get defining quartet for node %d" % n.id
    # first get the two descendants under the node
    if not (n.father):
        ERR("Node id "+str(n.id)+" has no father. Can't obtain triplet")
    f = n.father
    if not (n.children):
            ERR("Internal starting node "+str(n.id)+" has no son")
    if  len(f.children) < 2:
        ERR("Father %d has a single son %d" % (f.id, f.children[0].id))
    s1 = n.children[0]
    s2 = n.children[1]
    if s2 == None:
        ERR("No brother to son "+str(s1.id)+" of base node "+str(n.id))
    d1 = get_descendant(s1)
    d2 = get_descendant(s2)

# We first check if father has three sons        
    if len(f.children) > 2:
        print "Three sons at node %d (possibly the root)\n" % f.id
#        exit(8)
        brothers = []
        print "searching for brothers"
        for b in f.children:
            if b != n:
                print "found a son %d that is not me. Adding to brothers" % \
                    b.id
                brothers.append(b)
                if len(brothers) == 2:
                    print "Two brother %s found. leaving" % str(brothers)
                    break
            else:
                print "Found myself: %d" %b.id
        if len(brothers) < 2:
            ERR("Did not find two bothers. brothers: %s" %str(brothers))
                
        print "getting descendants of both brothers %d, %d" % \
            (brothers[0].id, brothers[1].id )
        d3 = get_descendant(brothers[0])
        d4 = get_descendant(brothers[1])
        return d1, d2, d3, d4
            

    print "Father %d has only two sons: %d, %d" %     \
        (f.id, f.children[0].id, f.children[1].id)
    if f.children[0] != n and  f.children[1] != n:
        ERR("Non of father sons %d and %d is me %d" %  \
                (f.children[0].id, f.children[1].id, n.id))
    if f.children[0] == n:
        b  = f.children[1]
    else :
        b = f.children[0]

    if not (f.father):
        print "father ",f.id," has no father. Is a root"
        if not (b.children):
            ERR("No father to father (root?) %d and no son to brother %d", f.id, )
        g1 = b.children[0]
        if  len(b.children) < 2: 
            ERR("Grandson "+str(g1.id)+" must have a brother ")
        g2 = b.children[1]
        print "Found two grandsons of root", g1.id, g2.id
        d3 = get_descendant(g1)
        d4 = get_descendant(g2)
        return d1, d2, d3, d4
    else:
        u =  get_brother(f)
        if u == None:
            ERR("No bother to father "+str(f.id))
        
        d3 = get_descendant(b)
        d4 = get_descendant(u)
        return d1, d2, d3, d4


#-----------------------------------------------------------------------------------_#
def get_descendant(n):
    print "At get descendant for node ",n.id
    if (n.name):
#        print "starting node", n.id, " is a leaf", n.name
#    while hasattr(n,'son'):
        pass
    while n.children:
#        print "descending at ",n.id
        if n.name:
            ERR("Node "+str(n.id)+" has son "+str(n.children[0].id)+" and a name "+n.name)
        n = n.children[0]
    
#    if not hasattr(n,'name'):
    if not n.name:
        ERR("Node "+str(n.id)+" has NO son and also NO name ")
#    print "ended descending at ", n.id, " with name ", n.name
#    exit(8)
    return n

#-----------------------------------------------------------------------------------_#
#-----------------------------------------------------------------------------------_#
def get_brothers(n):
    if not (n.father):
        ERR("No father to node "+str(n.id)+" and hence also no brother")
    f=n.father
    if not (f.children):
        ERR("No son to father node "+str(f.id))
    if len(f.children) < 2:
        ERR("Father %d has a single son %d" % (f.id,
                                               f.children[0].id))
    bros = []
    for s in f.children:
        print "at son ", s.id
        if s == n:
            print "this son is me ", n.id
        else:
            print "arrived at a son ",s.id ," that is not me"
            bros.append(s)
    bros2 = [s for s in f.children if s != n]
    if set(bros) != set(bros2):
        print "Brothers are different:\n", str([x.id for x in bros]),"\nversus:\n", str([x.id for x in bros2])
        exit(7)
        
    return bros


#-----------------------------------------------------------------------------------_#



def main(argv):
    print "in main"
    tree = "(((((((((((((((a,b)))),c))),d))),((((e:.2):.2),((((g:.1):.1):.1):.1):.1):.1):.1):.1):.1)):.1)"

#    genAllRootings("(a,(b,(c,d)))")
#    exit(8)
    t_leaves = {};
    TreeNode_n.count = 0
    print tree
    t_root = insertSubTree_n(tree, t_leaves)
    t_s = printTree_n(t_root,"")
    print "\n\n\nt_s:", t_s

    btstrp_thresh = 80
    rn = contract_tree(t_root)
#    contract_low_bootstrap(t_root, btstrp_thresh)
    t_s = printTree_n(rn,"")
    print "\n\n\nTree after contract :", t_s
  


#-----------------------------------------------------------------------------------_#


if __name__ == "__main__":
   main(sys.argv[1:])
