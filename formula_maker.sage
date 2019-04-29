from itertools import product, combinations_with_replacement
from sage.symbolic.function_factory import function
from collections import Counter

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
    all multiDAGs with set of nodes S such that:
    (1): Every node in A_0 has no incoming arc.
    (2): Every node in A_i (i>0) has 2 incoming arcs,
    (2a): all of them with source in A_j (j<i) and
    (2b): at least one of them with source in A_{i-1}"""
    list_of_maps = []
    for i in range(1, len(part)):
        for l in part[i]:
            list_of_maps.append(preimages(part[:i], l))
    for f in product(*list_of_maps):
        dag = DiGraph(multiedges=True)
        for l in f:
            dag.add_edges(l)
        yield dag

def dags_for_partition_mod_iso(part):
    """Same as dags_for_partition but yields only non-isomorphic ones.
    """
    found = []
    for dag in dags_for_partition(part):
        for other in found:
            if dag.is_isomorphic(other):
                break
        else:
            found.append(dag)
            yield(dag)

def restricted_partitions(num_nodes):
    """Returns a generator that yields all ordered partitions of the set {0,...,num_nodes-1}
    of the form [[0],[1,...,k1],[k1+1,...,k2],...,[kl+1,...,num_nodes-1]]
    """
    comps = Compositions(num_nodes-1)
    for comp in comps:
        l = [[0]]
        last = 1
        for i in comp:
            l.append(list(range(last,last+i)))
            last += i
            yield l
        
cached_dags = {}        
def structural_dags_mod_iso(num_nodes):
    """Returns a generator that yields all structural multiDAGs with `num_nodes` nodes
    modulo isomorphism
    """
    if num_nodes in cached_dags:
        for dag in cached_dags[num_nodes]:
            yield dag
        return
    found = []
    numparts = 0
    for part in restricted_partitions(num_nodes):
        numparts += 1
        for dag in dags_for_partition_mod_iso(part):
            found.append(dag)
            yield dag
    cached_dags[num_nodes] = found


L = function("L")

cached_symbolic = {}
def counting_formula_symbolic(num_hybrid):
    """Computes the symbolic formula for the number of BTC networks with given
    number of hybrid nodes.
    """
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
    """Returns the number of BTC networks with given number of leaves and hybrid nodes.
    """
    if num_hybrid == 0:
        return t(num_leaves)
    symbolic = counting_formula_symbolic(num_hybrid)
    symbols = list(var('n%d' % i) for i in range(num_hybrid+1))
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
        factor *= prod(map(t,el))
        factor *= multinomial(*el)
        total += factor
    return total


for n in range(1,5):
    print "Formula: F(%d)=%s" % (n,counting_formula_symbolic(n))
    for h in range(n):
        print "Count: |BTC_{%d,%d}|=%d" % (n,h,count_BTC(n,h))

