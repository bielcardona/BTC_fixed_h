from itertools import permutations, combinations_with_replacement #, product
from collections import defaultdict, Counter
import networkx as nx


def simple_product(it1, it2):
    if not it1 or not it2:
        return
    it1x = iter(it1)
    it2x = iter(it2)
    gone = []
    x = it1x.__next__()
    for y in it2x:
        gone.append(y)
        yield x, y
    for x in it1x:
        for y in gone:
            yield x, y


def product(*iterables):
    """
    https://stackoverflow.com/questions/55882454/is-there-a-way-to-efficiently-compute-the-product-of-two-or-more-iterators
    """
    if len(iterables) == 0:
        yield []
        return
    if len(iterables) == 1:
        for e in iterables[0]:
            yield [e]
        return
    if len(iterables) == 2:
        yield from simple_product(*iterables)
        return
    it1, *rest = iterables
    gone = []
    x = it1.__next__()
    for t in product(*rest):
        gone.append(t)
        yield (x,) + t
    for x in it1:
        for t in gone:
            yield (x,) + t


def partitions(l):
    """Returns a generator that yields all partitions of the list `l`."""
    if len(l) == 0:
        yield []
        return

    if len(l) == 1:
        yield [l]
        return

    first = l[0]
    for smaller in partitions(l[1:]):
        # insert `first` in each of the subpartition's subsets
        for n, subset in enumerate(smaller):
            yield smaller[:n] + [[first] + subset] + smaller[n + 1:]
        # put `first` in its own subset
        yield [[first]] + smaller


def restricted_ordered_partitions(l):
    """Returns a generator that yields all ordered partitions of `l`
    such the first subset has only one element.
    """
    for i in range(len(l)):
        first = l[i]
        others = l[:i] + l[i + 1:]
        for partition in partitions(others):
            for reordering in permutations(partition):
                yield [[first]] + list(reordering)


def k_partitions(l, k):
    """Returns a generator that yields all partitions of `l` with `k` blocks.
    """
    if not l:
        return
    if k == 1:
        yield [l]
        return
    last = l[-1]
    others = l[:-1]
    for previous in k_partitions(others, k-1):
        yield previous + [[last]]
    for previous in k_partitions(others, k):
        for i in range(k):
            yield previous[:i] + [previous[i] + [last]] + previous[i+1:]


def k_partitions_tuple(l, k):
    """
    Same as k_partitions but now the blocks of each partition is a tuple
    """
    for part in k_partitions(l, k):
        yield [tuple(p) for p in part]


def permutations_with_repetitions(c):
    """
    Returns a generator that yields oll permutations with repetitions with
    the elements of the counter `c`.
    For instance permutations_with_repetitions(Counter('aab') yields
    the permutations aab, aba baa
    """
    unique = [x for x in c.keys() if c[x] > 0]
    if len(unique) == 1 and c[unique[0]] == 1:
        yield unique
        return
    for first in unique:
        cant = c.copy()
        cant[first] -= 1
        for p in permutations_with_repetitions(cant):
            yield [first] + p


def join(l):
    """Returns the concatenation of all the sublists in `l`."""
    return sum(l, [])


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
        dag = nx.MultiDiGraph()
        for l in f:
            dag.add_edges_from(l)
        yield dag


def structural_dags(nodes):
    """Returns a generator that yields all structural multiDAGs having as set of nodes `nodes`
    """
    for partition in restricted_ordered_partitions(nodes):
        for dag in dags_for_partition(partition):
            yield dag


def dags_for_partition_as_dict(part):
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
        dag = defaultdict(list)
        for l in f:
            for pair in l:
                dag[pair[0]].append(pair[1])
        yield dag


def structural_dags_as_dict(nodes):
    """Returns a generator that yields all structural multiDAGs having as set of nodes `nodes`
    """
    for partition in restricted_ordered_partitions(nodes):
        for dag in dags_for_partition_as_dict(partition):
            yield dag

# Tree generator


def add_elementary(t, u, kind=None, label=None):
    """Adds to the tree `t` an elementary node above `u` and labels
    it with an attribute `kind` with value `label`.
    Returns the identifier of this new node.
    """
    parents = t.predecessors(u)
    uprime = t.number_of_nodes()
    for parent in parents:
        t.add_edge(parent, uprime)
        t.remove_edge(parent, u)
    t.add_edge(uprime, u)
    if label:
        labels = nx.get_node_attributes(t, kind)
        labels[uprime] = label
        nx.set_node_attributes(t, kind, labels)
    return uprime


def add_leaf(t, uprime, taxon):
    """Adds to the tree `t` a leaf hanging from `uprime` and labels it with `taxon`.
    """
    new_leaf = t.number_of_nodes()
    t.add_edge(uprime, new_leaf)
    labels = nx.get_node_attributes(t, 'label')
    labels[new_leaf] = taxon
    nx.set_node_attributes(t, 'label', labels)


def tree_generator(taxa):
    """
    Returns a generator that yields all binary phylogenetic trees over `taxa`.
    """
    if len(taxa) == 1:
        t = nx.DiGraph()
        t.add_nodes_from([0])
        labels = {0: taxa[0]}
        nx.set_node_attributes(t, 'label', labels)
        yield t
        return
    taxon = taxa[-1]
    for tant in tree_generator(taxa[:-1]):
        for u in tant.nodes():
            t = tant.copy()
            uprime = add_elementary(t, u)
            add_leaf(t, uprime, taxon)
            yield t


def split_permutation(perm, sep='|'):
    """
    Splits the list `perm` into sublists, indicated by the separator `sep`.
    Example: ['a','|','b','c','|','|','d','|'] -> [['a'],['b','c'],[],['d'],[]]
    :param perm:
    :param sep:
    :return:
    """
    output = []
    actual = []
    for x in perm:
        if x == sep:
            output.append(actual)
            actual = []
        else:
            actual.append(x)
    output.append(actual)
    return output


def tree_annotations_generator(tree, symbols):
    if len(symbols) == 0:
        t = tree.copy()
        yield t
        return
    num_nodes = tree.number_of_nodes()
    all_symbols = symbols.copy()
    all_symbols['|'] = num_nodes-1
    perms = permutations_with_repetitions(all_symbols)
    for p in perms:
        words = list(split_permutation(p))
        t = tree.copy()
        for (u, w) in zip(tree.nodes(), words):
            for letter in w:
                add_elementary(t, u, 'T', letter)
        yield t


def annotated_tree_generator(taxa, symbols):
    for tree in tree_generator(taxa):
        for at in tree_annotations_generator(tree, symbols):
            root = [u for u in at.nodes() if at.in_degree(u) == 0][0]
            new_root = add_elementary(at, root)
            nx.set_node_attributes(at, 'H', {new_root: taxa})
            yield at


def annotated_forest_generator(part, dag):
    tgs = [annotated_tree_generator(block, Counter(dag[block])) for block in part]
    fg = product(*tgs)
    for f in fg:
        forest = nx.disjoint_union_all(f)
        Hlabels = nx.get_node_attributes(forest,'H')
        hybrids = {v: k for k, v in Hlabels.items()}
        Tlabels = nx.get_node_attributes(forest,'T')
        for u, Tlab in Tlabels.items():
            forest.add_edge(u, hybrids[Tlab])
        root = [u for u in forest.nodes() if forest.in_degree(u) == 0][0]
        forest.remove_node(root)
        yield forest


def BTC_generator(taxa, h):
    for part in k_partitions_tuple(taxa, h+1):
        for dag in structural_dags_as_dict(part):
            for f in annotated_forest_generator(part,dag):
                yield f
