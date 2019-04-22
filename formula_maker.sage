
# coding: utf-8

# In[1]:


from itertools import product, combinations_with_replacement

def preimages(part, el):
    """Given an (ordered) partition `part`, say A_0,...,A_k of A,
    and an element `el` returns a generator that yields tuples of the form
    ((x,el)) or ((x,el),(y,el)) where x,y belong to A and at least one of them belongs to A_k"""
    last = part[-1]
    for x, y in combinations_with_replacement(last, 2):
        yield ((x, el), (y, el))
    prev = sum(part[:-1], [])
    for x, y in product(last, prev):
        yield ((x, el), (y, el))

def dags_for_partition(part):
    """Given an (ordered) partition `part`, say A_0,...,A_k of a set S, returns a generator that yields
    all DAGs with set of nodes S such that:
    (1): Every node in A_0 has no incoming arc.
    (2): Every node in A_i (i>0) has 1 or 2 incoming arcs,
    (2a): all of them with source in A_j (j<i) and
    (2b): at least one of them with source in A_{i-1}"""
    list_of_maps = []
    for i in range(1, len(part)):
        for l in part[i]:
            list_of_maps.append(preimages(part[:i], l))
    for f in product(*list_of_maps):
        dag = DiGraph(multiedges=True)
        # dag = defaultdict(list)
        for l in f:
            #for pair in l:
            #    dag[pair[0]].append(pair[1])
            dag.add_edges(l)
        yield dag

def dags_for_partition_mod_iso(part):
    found = []
    for dag in dags_for_partition(part):
        for other in found:
            if dag.is_isomorphic(other):
                break
        else:
            found.append(dag)
            yield(dag)

def restricted_partitions(num_nodes):
    comps = Compositions(num_nodes-1)
    for comp in comps:
        l = [[0]]
        last = 1
        for i in comp:
            l.append(list(range(last,last+i)))
            last += i
        #print comp,"->",l
        yield l
        
cached_dags = {}        
def structural_dags_mod_iso(num_nodes):
    if num_nodes in cached_dags:
        print "cached"
        for dag in cached_dags[num_nodes]:
            yield dag
        return
    print "non cached"
    found = []
    numparts = 0
    for part in restricted_partitions(num_nodes):
        numparts += 1
        print "partition ", numparts , part
        for dag in dags_for_partition_mod_iso(part):
            found.append(dag)
            print(len(found))
            yield dag
    cached_dags[num_nodes] = found


# In[2]:


from sage.symbolic.function_factory import function
L = function("L")
labelings = function("labelings")
num_trees = function("num_trees")
num_nodes = function("num_nodes")
from collections import Counter
import pickle

cached_symbolic = {}
def counting_formula_symbolic(num_hybrid):
    try:
        return cached_symbolic[num_hybrid]
    except:
        pass
    dags = structural_dags_mod_iso(num_hybrid+1)
    symbols = list(var('n%d' % i) for i in range(num_hybrid+1))
    total = 0
    for dag in dags:
        subproduct = 1 / dag.automorphism_group().order()
        for v in dag.vertices():
            c = Counter(Counter(dag.outgoing_edges(v)).values())
            if c[1] + c[2] > 0:
                subproduct *= L(symbols[v],c[1],c[2])
        total += subproduct
    cached_symbolic[num_hybrid]=total
    return total

def t(ni):
    return (2*ni-3).multifactorial(2)
def s(ni):
    return 2*ni-1
def lab(ni,ai,bi):
    return factorial(s(ni)-1+ai+2*bi) / (factorial(s(ni)-1) * 2^bi)

def count_BTC(num_leaves,num_hybrid):
    if num_hybrid == 0:
        return t(num_leaves)
    symbolic = counting_formula_symbolic(num_hybrid)
    symbols = list(var('n%d' % i) for i in range(num_hybrid+1))
    # print symbols
    num_nodes = num_hybrid+1
    iterators = num_nodes * [range(1,num_leaves-num_hybrid+1)]
    elements = cartesian_product(iterators)
    total = 0
    for el in elements:
        if sum(el) != num_leaves:
            continue
        substitutions = {symbols[i]: el[i] for i in range(num_nodes)}
        symbolic_part = symbolic.subs(substitutions)
        factor = symbolic_part.substitute_function(L,lab)
        # print map(t,el)
        factor *= prod(map(t,el))
        factor *= multinomial(*el)
        total += factor
    return total


# In[ ]:

for n in range(7,8):
    print count_BTC(n,n-1)


# In[ ]:


# print cached_symbolic

