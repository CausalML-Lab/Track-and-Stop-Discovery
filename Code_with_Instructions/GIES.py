import networkx as nx
import itertools as itr
import random
from typing import Union, List
import numpy as np
from graph_utils import fill_vstructures, direct_chordal_graph
from scipy.special import binom
import pyAgrum as gum

def tree_of_cliques(
        degree: int,
        min_clique_size: int,
        max_clique_size: int,
        ngraphs: int = 1,
        nnodes: int = None
):
    if ngraphs == 1:
        counter = random.randint(min_clique_size, max_clique_size)
        source_clique = list(range(counter))
        previous_layer_cliques = [source_clique]
        current_layer_cliques = []
        arcs = set(itr.combinations(source_clique, 2))

        while counter < nnodes:
            for parent_clique in previous_layer_cliques:
                for d in range(degree):
                    if counter < nnodes:
                        clique_size = random.randint(min_clique_size, max_clique_size)
                        intersection_size = min(len(parent_clique)-1, random.randint(int(clique_size/2), clique_size-1))
                        num_new_nodes = clique_size - intersection_size
                        print(intersection_size)

                        indices = set(random.sample(parent_clique, intersection_size))
                        intersection = [
                            parent_clique[ix] for ix in range(len(parent_clique))
                            if ix in indices
                        ]
                        new_clique = intersection + list(range(counter, counter+num_new_nodes))
                        current_layer_cliques.append(new_clique)
                        arcs.update(set(itr.combinations(new_clique, 2)))
                        counter += num_new_nodes
            previous_layer_cliques = current_layer_cliques.copy()
        g = nx.DiGraph()
        g.add_edges_from(arcs)
        if not nx.is_connected(g.to_undirected()):
            raise RuntimeError
        return g
    else:
        return [tree_of_cliques(degree, min_clique_size, max_clique_size, nnodes=nnodes) for _ in range(ngraphs)]


def hairball(degree: int, num_layers: int = None, nnodes: int = None) -> nx.DiGraph:
    """
    Return a `degree`-ary tree with `nnodes` nodes or `num_layers` layers.
    """
    nnodes = sum((degree**i for i in range(num_layers))) if nnodes is None else nnodes
    return nx.full_rary_tree(degree, nnodes, create_using=nx.DiGraph)


def hairball_plus(
        degree: int,
        e_min: int = None,
        e_max: int = None,
        ngraphs: int = 1,
        num_layers: int = None,
        nnodes: int = None,
        edge_prob: float = None,
        nontree_factor: float = None
):
    """
    Generate a "hairball", then add a random number of edges between `e_min` and `e_max`, then
    triangulate the graph to make it chordal.
    """
    if ngraphs == 1:
        g = hairball(degree, num_layers=num_layers, nnodes=nnodes)
        order = list(nx.topological_sort(g))
        if e_min is not None:
            num_extra_edges = random.randint(e_min, e_max)
        elif edge_prob is not None:
            nnodes = g.number_of_nodes()
            num_missing_edges = binom(nnodes, 2) - (nnodes - 1)
            num_extra_edges = np.random.binomial(num_missing_edges, edge_prob)
        else:
            num_extra_edges = int(nontree_factor*(nnodes - 1))
        missing_edges = [(i, j) for i, j in itr.combinations(order, 2) if not g.has_edge(i, j)]
        extra_edges = random.sample(missing_edges, num_extra_edges)
        g.add_edges_from(extra_edges)
        fill_vstructures(g, order)
        return g
    else:
        return [hairball_plus(degree, e_min, e_max, num_layers=num_layers, nnodes=nnodes) for _ in range(ngraphs)]


def random_directed_tree(nnodes: int):
    """
    Generate a random undirected tree, then pick a random root to make it directed.
    """
    g = nx.random_tree(nnodes)
    root = random.randint(0, nnodes-1)
    d = nx.DiGraph()

    queue = [root]
    while queue:
        current_node = queue.pop()
        nbrs = list(g.neighbors(current_node))
        d.add_edges_from([(current_node, nbr) for nbr in nbrs])
        queue += nbrs
        g.remove_node(current_node)
    return d


def random_chordal_graph(nnodes, p=.1, ngraphs=1) -> Union[nx.DiGraph, List[nx.DiGraph]]:
    if ngraphs == 1:
        g = nx.erdos_renyi_graph(nnodes, p)
        perm = random.sample(set(range(nnodes)), nnodes)
        for node in perm:
            nbrs = set(g.neighbors(node)) | {node}
            g.add_edges_from((i, j) for i, j in itr.combinations(nbrs, 2))

        d = nx.DiGraph()
        for node in perm:
            d.add_edges_from([(node, nbr) for nbr in g.neighbors(node)])
        return d
    else:
        return [random_chordal_graph(nnodes, p=p) for _ in range(ngraphs)]


def shanmugam_random_chordal(nnodes, density):
    while True:
        d = nx.DiGraph()
        d.add_nodes_from(set(range(nnodes)))
        order = list(range(1, nnodes))
        for i in order:
            num_parents_i = max(1, np.random.binomial(i, density))
            parents_i = random.sample(list(range(i)), num_parents_i)
            d.add_edges_from({(p, i) for p in parents_i})
        for i in reversed(order):
            for j, k in itr.combinations(d.predecessors(i), 2):
                d.add_edge(min(j, k), max(j, k))

        perm = np.random.permutation(list(range(nnodes)))
        d = nx.relabel.relabel_nodes(d, dict(enumerate(perm)))

        return d


def tree_plus(nnodes: int, e_min: int, e_max: int, ngraphs: int = 1):
    if ngraphs == 1:
        g = random_directed_tree(nnodes)
        order = list(nx.topological_sort(g))
        num_extra_edges = random.randint(e_min, e_max)
        extra_edges = random.sample(list(itr.combinations(order, 2)), num_extra_edges)
        g.add_edges_from(extra_edges)
        fill_vstructures(g, order)
        return g
    else:
        return [tree_plus(nnodes, e_min, e_max) for _ in range(ngraphs)]


def random_chordal_graph2(nnodes: int, k: int, ngraphs: int = 1, ensure_connected=True) -> Union[List[nx.DiGraph], nx.DiGraph]:
    if ngraphs == 1:
        for ix in itr.count():
            if ix > 100:
                raise ValueError("100 iterations without a connected graph, please change parameters")
            t = nx.random_tree(nnodes)

            subtrees = []
            for i in range(nnodes):
                x = random.randint(0, nnodes-1)
                t_i = nx.Graph()
                t_i.add_node(x)
                t_i_nbrs = {x: set(t.neighbors(x))}
                k_i = random.randint(1, 2*k-1)
                for j in range(k_i-1):
                    y = random.sample(t_i_nbrs.keys(), 1)[0]
                    z = random.sample(t_i_nbrs[y], 1)[0]
                    t_i.add_edge(y, z)
                    t_i_nbrs[y] -= {z}
                    t_i_nbrs[z] = set(t.neighbors(z)) - {y}
                    if not t_i_nbrs[y]:
                        del t_i_nbrs[y]
                    if not t_i_nbrs[z]:
                        del t_i_nbrs[z]
                subtrees.append(t_i)

            g = nx.Graph()
            g.add_nodes_from(range(nnodes))
            for (i, t_i), (j, t_j) in itr.combinations(enumerate(subtrees), 2):
                if t_i.nodes & t_j.nodes: g.add_edge(i, j)

            if not ensure_connected or nx.is_connected(g):
                return direct_chordal_graph(g)
    else:
        return [random_chordal_graph2(nnodes, k) for _ in range(ngraphs)]


if __name__ == '__main__':
    import causaldag as cd

    d_nx = shanmugam_random_chordal(10, .1)
    print(d_nx.number_of_nodes())
    print(nx.chordal_graph_cliques(d_nx.to_undirected()))
    print(nx.is_connected(d_nx.to_undirected()))
    d = cd.DAG.from_nx(d_nx)
    print(d.vstructures())



def convert(adjacency_list):    
    graph_dict = {}
    for u, v in adjacency_list:
        if u in graph_dict:
            graph_dict[u].append(v)
        else:
            graph_dict[u] = [v]
        if v in graph_dict:
            graph_dict[v].append(u)
        else:
            graph_dict[v] = [u]
    return graph_dict 

# Python3 program to implement greedy
# algorithm for graph coloring
 
def addEdge(adj, v, w):
     
    adj[v].append(w)
     
    # Note: the graph is undirected
    adj[w].append(v) 
    return adj
 
# Assigns colors (starting from 0) to all
# vertices and prints the assignment of colors
def greedyColoring(adj, V):
     
    result = [-1] * V
 
    # Assign the first color to first vertex
    result[0] = 0;
 
 
    # A temporary array to store the available colors.
    # True value of available[cr] would mean that the
    # color cr is assigned to one of its adjacent vertices
    available = [False] * V
 
    # Assign colors to remaining V-1 vertices
    for u in range(1, V):
         
        # Process all adjacent vertices and
        # flag their colors as unavailable
        for i in adj[u]:
            if (result[i] != -1):
                available[result[i]] = True
 
        # Find the first available color
        cr = 0
        while cr < V:
            if (available[cr] == False):
                break
             
            cr += 1
             
        # Assign the found color
        result[u] = cr
 
        # Reset the values back to false
        # for the next iteration
        for i in adj[u]:
            if (result[i] != -1):
                available[result[i]] = False
 
    # Print the result
    #for u in range(V):
        #print("Vertex", u, " --->  Color", result[u])
    return result

def indices_of_elements(input_list):
    unique_elements = set(input_list)
    result = [[] for _ in range(max(unique_elements) + 1)]

    for index, element in enumerate(input_list):
        result[element].append(index)

    return result

def adj_list_to_string_with_vertices(adj_list):
    vertices = sorted(set(v for edge in adj_list for v in edge))
    vertices_str = ';'.join(map(str, vertices))

    adj_str = ';'.join([f"{edge[0]}->{edge[1]}" for edge in adj_list])

    return f"{vertices_str};{adj_str}"
    
    
    
def block_sample_intervention(bn,intv,s):
    import numpy as np

    bn1  = gum.BayesNet(bn)
    import numpy as np
    for j in intv:
        for i in bn1.parents(j):
            bn1.eraseArc(gum.Arc(i,j))
       
       
        bn1.cpt(j)[:] = [0.5 ,0.5]

           
                                        
            
    df = gum.generateSample(bn1, n=s, name_out=None, show_progress=False, with_labels=True, random_order=False)
    aa = df[0]
    return aa  



def Sample_and_update_dist(bb,index,Ps,Sample_No):
     
    def binary_vector_to_decimal(binary_vector):
        decimal = 0
        for digit in binary_vector:
            decimal = decimal * 2 + digit
        return decimal

    
    bb = bb[Sample_No+1]
    bb = binary_vector_to_decimal(bb)
    Ps[index][bb] =   Ps[index][bb] + 1
    return    






# Copyright 2022 Olga Kolotuhina, Juan L. Gamella

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:

# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.

# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.

# 3. Neither the name of the copyright holder nor the names of its
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

"""The main module, containing the implementation of GIES, including
the logic for the insert, delete and turn operators. The
implementation is directly based on the 2012 GIES paper by Hauser &
Bühlmann, "Characterization and Greedy Learning of Interventional
Markov Equivalence Classes of Directed Acyclic Graphs" -
https://www.jmlr.org/papers/volume13/hauser12a/hauser12a.pdf

Further credit is given where due.

Additional modules / packages:

  - gies.utils contains auxiliary functions and the logic to transform
    a PDAG into a CPDAG, used after each application of an operator.
  - gies.scores contains the modules with the score classes:
      - gies.scores.decomposable_score contains the base class for
        decomposable score classes (see that module for more details).
      - gies.scores.gauss_obs_l0_pen contains a cached implementation
        of the gaussian likelihood BIC score used in the original GES
        paper.
   - gies.test contains the modules with the unit tests and tests
     comparing against the algorithm's implementation in the R package
     'pcalg'.

"""

import numpy as np
import utils1 as utils



def fit_bic(
    data, I, A0=None, phases=["forward", "backward", "turning"], iterate=True, debug=0,bn=None
):
    """Run GIES on the given data, using the Gaussian BIC score
    (l0-penalized Gaussian Likelihood). The data is not assumed to be
    centered, i.e. an intercept is fitted.

    To use a custom score, see gies.fit.

    Parameters
    ----------
    data : list of numpy.ndarray
        Every matrix in the list corresponds to an environment,
        the n x p matrix containing the observations, where columns
        correspond to variables and rows to observations.
    I: a list of lists
        The family of intervention targets, with each list being the
        targets in the corresponding environment.
    A0 : numpy.ndarray, optional
        The initial I-essential graph on which GIES will run, where where `A0[i,j]
        != 0` implies the edge `i -> j` and `A[i,j] != 0 & A[j,i] !=
        0` implies the edge `i - j`. Defaults to the empty graph
        (i.e. matrix of zeros).
    phases : [{'forward', 'backward', 'turning'}*], optional
        Which phases of the GIES procedure are run, and in which
        order. Defaults to `['forward', 'backward', 'turning']`.
    iterate : bool, default=True
        Indicates whether the given phases should be iterated more
        than once.
    debug : int, optional
        If larger than 0, debug are traces printed. Higher values
        correspond to increased verbosity.

    Returns
    -------
    estimate : numpy.ndarray
        The adjacency matrix of the estimated I-essential graph
    total_score : float
        The score of the estimate.

    Raises
    ------
    TypeError:
        If the type of some of the parameters was not expected,
        e.g. if data is not a numpy array.
    ValueError:
        If the value of some of the parameters is not appropriate,
        e.g. a wrong phase is specified.

    Example
    -------

    Data from a linear-gaussian SCM (generated using
    `sempler <https://github.com/juangamella/sempler>`__)

    >>> import numpy as np
    >>> data = [np.array([[3.23125779, 3.24950062, 13.430682, 24.67939513],
    ...                  [1.90913354, -0.06843781, 6.93957057, 16.10164608],
    ...                  [2.68547149, 1.88351553, 8.78711076, 17.18557716],
    ...                  [0.16850822, 1.48067393, 5.35871419, 11.82895779],
    ...                  [0.07355872, 1.06857039, 2.05006096, 3.07611922]])]
    >>> interv = [[]]

    Run GIES using the gaussian BIC score:

    >>> import gies
    >>> gies.fit_bic(data, interv)
    (array([[0., 1., 1., 0.],
           [0., 0., 0., 0.],
           [1., 1., 0., 1.],
           [0., 1., 1., 0.]]), -15.641090039220082)

    """
    # Initialize Gaussian BIC score (precomputes scatter matrices, sets up cache)
    cache = GaussIntL0Pen(data, I,bn=bn)
    # Unless indicated otherwise, initialize to the empty graph
    A0 = np.zeros((cache.p, cache.p)) if A0 is None else A0
    return fit(cache, A0, phases, iterate, debug)


def fit(
    score_class,
    A0=None,
    phases=["forward", "backward", "turning"],
    iterate=True,
    debug=0,
):
    """
    Run GIES using a user defined score.

    Parameters
    ----------
    score_class : gies.DecomposableScore
        An instance of a class which inherits from
        gies.decomposable_score.DecomposableScore (or defines a
        local_score function and a p attribute, see
        gies.decomposable_score for more info).
    A0 : np.array, optional
        the initial I-essential graph on which GIES will run, where where A0[i,j]
        != 0 implies i -> j and A[i,j] != 0 & A[j,i] != 0 implies i -
        j. Defaults to the empty graph.
    phases : [{'forward', 'backward', 'turning'}*], optional
        which phases of the GIES procedure are run, and in which
        order. Defaults to ['forward', 'backward', 'turning'].
    iterate : bool, default=True
        Indicates whether the given phases should be iterated more
        than once.
    debug : int, optional
        if larger than 0, debug are traces printed. Higher values
        correspond to increased verbosity.

    Returns
    -------
    estimate : np.array
        the adjacency matrix of the estimated CPDAG
    total_score : float
        the score of the estimate

    """
    # Select the desired phases
    if len(phases) == 0:
        raise ValueError("Must specify at least one phase")
    # Unless indicated otherwise, initialize to the empty graph
    A0 = np.zeros((score_class.p, score_class.p)) if A0 is None else A0
    # GES procedure
    #total_score = score_class.full_score(A0)
    total_score = 0
    A, score_change = A0, np.Inf
    # Run each phase
    while True:
        last_total_score = total_score
        for phase in phases:
            if phase == "forward":
                fun = forward_step
            elif phase == "backward":
                fun = backward_step
            elif phase == "turning":
                fun = turning_step
            else:
                raise ValueError('Invalid phase "%s" specified' % phase)
            print("\nGES %s phase start" % phase) if debug else None
            print("-------------------------") if debug else None
            while True:
                score_change, new_A = fun(A, score_class, max(0, debug - 1))
                # TODO
                # Was > 0, but might be stuck in loop for the turn operator
                if score_change > 0.01:
                    # Transforming the partial I-essential graph into an I-essential graph
                    A = utils.replace_unprotected(new_A, score_class.interv)
                    total_score += score_change
                else:
                    break
            print("-----------------------") if debug else None
            print("GES %s phase end" % phase) if debug else None
            print("Total score: %0.10f" % total_score) if debug else None
            [print(row) for row in A] if debug else None
        if total_score <= last_total_score or not iterate:
            break
    return A, total_score


def forward_step(A, cache, debug=0):
    """
    Scores all valid insert operators that can be applied to the current
    I-essential graph A, and applies the highest scoring one.

    Parameters
    ----------
    A : np.array
        the adjacency matrix of an I-essential graph, where A[i,j] != 0 => i -> j
        and A[i,j] != 0 & A[j,i] != 0 => i - j.
    cache : DecomposableScore
        an instance of the score class, which computes the change in
        score and acts as cache.
    debug : int, optional
        if larger than 0, debug are traces printed. Higher values
        correspond to increased verbosity

    Returns
    -------
    score : float
        the change in score resulting from applying the highest
        scoring operator (note, can be smaller than 0).
    new_A : np.array
        the adjacency matrix of the partial I-essential graph resulting from applying the
        operator (not yet a I-essential graph).

    """
    # Construct edge candidates (i.e. edges between non-adjacent nodes)
    fro, to = np.where((A + A.T + np.eye(len(A))) == 0)
    edge_candidates = list(zip(fro, to))
    # For each edge, enumerate and score all valid operators
    valid_operators = []
    print("  %d candidate edges" % len(edge_candidates)) if debug > 1 else None
    for (x, y) in edge_candidates:
        valid_operators += score_valid_insert_operators(
            x, y, A, cache, debug=max(0, debug - 1)
        )
    # Pick the edge/operator with the highest score
    if len(valid_operators) == 0:
        print("  No valid insert operators remain") if debug else None
        return 0, A
    else:
        scores = [op[0] for op in valid_operators]
        score, x, y, T = valid_operators[np.argmax(scores)]
        # Apply operator
        new_A = insert(x, y, T, A, cache.interv)
        print(
            "  Best operator: insert(%d, %d, %s) -> (%0.10f)" % (x, y, T, score)
        ) if debug else None
        print(new_A) if debug else None
        return score, new_A


def backward_step(A, cache, debug=0):
    """
    Scores all valid delete operators that can be applied to the current
    I-essential graph A, and applies the highest scoring one.

    Parameters
    ----------
    A : np.array
        the adjacency matrix of a I-essential graph, where A[i,j] != 0 => i -> j
        and A[i,j] != 0 & A[j,i] != 0 => i - j.
    cache : DecomposableScore
        an instance of the score class, which computes the change in
        score and acts as cache.
    debug : int
        if larger than 0, debug are traces printed. Higher values
        correspond to increased verbosity

    Returns
    -------
    score : float
        the change in score resulting from applying the highest
        scoring operator (note, can be smaller than 0).
    new_A : np.array
        the adjacency matrix of the partial I-essential graph resulting from applying the
        operator (not yet an I-essential graph).

    """
    # Construct edge candidates:
    #   - directed edges
    #   - undirected edges, counted only once
    fro, to = np.where(utils.only_directed(A))
    directed_edges = zip(fro, to)
    fro, to = np.where(utils.only_undirected(A))
    undirected_edges = filter(lambda e: e[0] > e[1], zip(fro, to))  # zip(fro,to)
    edge_candidates = list(directed_edges) + list(undirected_edges)
    assert len(edge_candidates) == utils.skeleton(A).sum() / 2
    # For each edge, enumerate and score all valid operators
    valid_operators = []
    print("  %d candidate edges" % len(edge_candidates)) if debug > 1 else None
    for (x, y) in edge_candidates:
        valid_operators += score_valid_delete_operators(
            x, y, A, cache, debug=max(0, debug - 1)
        )
    # Pick the edge/operator with the highest score
    if len(valid_operators) == 0:
        print("  No valid delete operators remain") if debug else None
        return 0, A
    else:
        scores = [op[0] for op in valid_operators]
        score, x, y, H = valid_operators[np.argmax(scores)]
        # Apply operator
        new_A = delete(x, y, H, A, cache.interv)
        print(
            "  Best operator: delete(%d, %d, %s) -> (%0.4f)" % (x, y, H, score)
        ) if debug else None
        print(new_A) if debug else None
        return score, new_A


def turning_step(A, cache, debug=0):
    """
    Scores all valid turn operators that can be applied to the current
    I-essential graph A, and applies the highest scoring one.

    Parameters
    ----------
    A : np.array
        the adjacency matrix of a I-essential graph, where A[i,j] != 0 => i -> j
        and A[i,j] != 0 & A[j,i] != 0 => i - j.
    cache : DecomposableScore
        an instance of the score class, which computes the change in
        score and acts as cache.
    debug : int
        if larger than 0, debug are traces printed. Higher values
        correspond to increased verbosity

    Returns
    -------
    score : float
        the change in score resulting from applying the highest
        scoring operator (note, can be smaller than 0).
    new_A : np.array
        the adjacency matrix of the partial I-essential graph resulting from applying the
        operator (not yet an I-essential graph).

    """
    # Construct edge candidates:
    #   - directed edges, reversed
    #   - undirected edges
    fro, to = np.where(A != 0)
    edge_candidates = list(zip(to, fro))
    # For each edge, enumerate and score all valid operators
    valid_operators = []
    print("  %d candidate edges" % len(edge_candidates)) if debug > 1 else None
    for (x, y) in edge_candidates:
        valid_operators += score_valid_turn_operators(
            x, y, A, cache, debug=max(0, debug - 1)
        )
    # Pick the edge/operator with the highest score
    if len(valid_operators) == 0:
        print("  No valid turn operators remain") if debug else None
        return 0, A
    else:
        scores = [op[0] for op in valid_operators]
        score, x, y, C = valid_operators[np.argmax(scores)]
        # Apply operator
        new_A = turn(x, y, C, A, cache.interv)
        print(
            "  Best operator: turn(%d, %d, %s) -> (%0.15f)" % (x, y, C, score)
        ) if debug else None
        print(new_A) if debug else None
        return score, new_A


# --------------------------------------------------------------------
# Insert operator
#    1. definition in function insert
#    2. enumeration logic (to enumerate and score only valid
#    operators) function in score_valid_insert_operators


def insert(x, y, T, A, I):
    """
    Applies the insert operator:
      1) Orients all edges of the chain component of y according to a perfect elimination ordering,
        such that for all t in T the previously undirected edge becomes t -> y
      2) adds the edge x -> y

    Parameters
    ----------
    x : int
        the origin node (i.e. x -> y)
    y : int
        the target node
    T : iterable of ints
        a subset of the neighbors of y which are not adjacent to x
    A : np.array
        the current adjacency matrix
    I : list of lists
        list of interventions

    Returns
    -------
    new_A : np.array
        the adjacency matrix resulting from applying the operator

    """
    # Check inputs
    T = sorted(T)
    if A[x, y] != 0 or A[y, x] != 0:
        raise ValueError("x=%d and y=%d are already connected" % (x, y))
    if len(T) == 0:
        pass
    elif not (A[T, y].all() and A[y, T].all()):
        raise ValueError("Not all nodes in T=%s are neighbors of y=%d" % (T, y))
    elif A[T, x].any() or A[x, T].any():
        raise ValueError("Some nodes in T=%s are adjacent to x=%d" % (T, x))

    new_A = A.copy()
    A_undir = utils.only_undirected(A)

    # Finding the nodes of the chain component of y
    chain_comp_y = utils.chain_component(y, A)
    # Getting set C = T u NA_xy
    C = utils.na(y, x, A) | set(T)
    # Getting all the nodes of the chain component in the right order
    nodes = list(C) + [y] + list(chain_comp_y - {y} - C)
    # Computing the perfect elimination ordering according to the above order
    ordering = utils.maximum_cardinality_search(A_undir, nodes)
    # Orienting the edges of the chain component according to the perfect elimination ordering
    new_A = utils.orient_edges(A, ordering)
    # Adding edge x -> y
    new_A[x, y] = 1
    return new_A


def score_valid_insert_operators(x, y, A, cache, debug=0):
    """Generate and score all valid insert(x,y,T) operators involving the edge
    x-> y, and all possible subsets T of neighbors of y which
    are NOT adjacent to x.

    Parameters
    ----------
    x : int
        the origin node (i.e. x -> y)
    y : int
        the target node
    A : np.array
        the current adjacency matrix
    cache : instance of gies.scores.DecomposableScore
        the score cache to compute the score of the
        operators that are valid
    debug : int
        if larger than 0, debug are traces printed. Higher values
        correspond to increased verbosity

    Returns
    -------
    valid_operators : list of tuples
        a list of tubles, each containing a valid operator, its score
        and the resulting connectivity matrix

    """
    p = len(A)
    if A[x, y] != 0 or A[y, x] != 0:
        raise ValueError("x=%d and y=%d are already connected" % (x, y))
    # One-hot encode all subsets of T0, plus one extra column to mark
    # if they pass validity condition 2 (see below)
    T0 = sorted(utils.neighbors(y, A) - utils.adj(x, A))
    if len(T0) == 0:
        subsets = np.zeros((1, p + 1), dtype=np.bool)
    else:
        subsets = np.zeros((2 ** len(T0), p + 1), dtype=np.bool)
        subsets[:, T0] = utils.cartesian(
            [np.array([False, True])] * len(T0), dtype=np.bool
        )
    valid_operators = []
    print("    insert(%d,%d) T0=" % (x, y), set(T0)) if debug > 1 else None
    while len(subsets) > 0:
        print(
            "      len(subsets)=%d, len(valid_operators)=%d"
            % (len(subsets), len(valid_operators))
        ) if debug > 1 else None
        # Access the next subset
        T = np.where(subsets[0, :-1])[0]
        passed_cond_2 = subsets[0, -1]
        subsets = subsets[1:]
        # Check that the validity conditions hold for T
        na_yxT = utils.na(y, x, A) | set(T)
        # Condition 1: Test that NA_yx U T is a clique
        cond_1 = utils.is_clique(na_yxT, A)
        if not cond_1:
            # Remove from consideration all other sets T' which
            # contain T, as the clique condition will also not hold
            supersets = subsets[:, T].all(axis=1)
            subsets = utils.delete(subsets, supersets, axis=0)
        # Condition 2: Test that all semi-directed paths from y to x contain a
        # member from NA_yx U T
        if passed_cond_2:
            # If a subset of T satisfied condition 2, so does T
            cond_2 = True
        else:
            # Check condition 2
            cond_2 = True
            for path in utils.semi_directed_paths(y, x, A):
                if len(na_yxT & set(path)) == 0:
                    cond_2 = False
                    break
            if cond_2:
                # If condition 2 holds for NA_yx U T, then it holds for all supersets of T
                supersets = subsets[:, T].all(axis=1)
                subsets[supersets, -1] = True
        print(
            "      insert(%d,%d,%s)" % (x, y, T),
            "na_yx U T = ",
            na_yxT,
            "validity:",
            cond_1,
            cond_2,
        ) if debug > 1 else None
        # If both conditions hold, apply operator and compute its score
        if cond_1 and cond_2:
            # Compute the change in score
            aux = na_yxT | utils.pa(y, A)
            old_score = cache.local_score(y, aux)
            new_score = cache.local_score(y, aux | {x})
            print(
                "        new: s(%d, %s) = %0.6f old: s(%d, %s) = %0.6f"
                % (y, aux | {x}, new_score, y, aux, old_score)
            ) if debug > 1 else None
            # Add to the list of valid operators
            valid_operators.append((new_score - old_score, x, y, T))
            print(
                "    insert(%d,%d,%s) -> %0.16f" % (x, y, T, new_score - old_score)
            ) if debug else None
    # Return all the valid operators
    return valid_operators


# --------------------------------------------------------------------
# Delete operator
#    1. definition in function delete
#    2. enumeration logic (to enumerate and score only valid
#    operators) function in valid_delete_operators


def delete(x, y, H, A, I):
    """
    Applies the delete operator:
      1) Orients all edges of the chain component of y according to a perfect elimination ordering,
        such that for every node h in H:
           * orients the edge y -> h
           * if the edge with x is undirected, orients it as x -> h
      2) deletes the edge x -> y or x - y

    Note that H must be a subset of the neighbors of y which are
    adjacent to x. A ValueError exception is thrown otherwise.

    Parameters
    ----------
    x : int
        the "origin" node (i.e. x -> y or x - y)
    y : int
        the "target" node
    H : iterable of ints
        a subset of the neighbors of y which are adjacent to x
    A : np.array
        the current adjacency matrix
    I : list of lists
        list of interventions

    Returns
    -------
    new_A : np.array
        the adjacency matrix resulting from applying the operator

    """
    H = set(H)
    # Check inputs
    if A[x, y] == 0:
        raise ValueError("There is no (un)directed edge from x=%d to y=%d" % (x, y))
    # neighbors of y which are adjacent to x
    na_yx = utils.na(y, x, A)
    if not H <= na_yx:
        raise ValueError(
            "The given set H is not valid, H=%s is not a subset of NA_yx=%s"
            % (H, na_yx)
        )

    # Apply operator
    new_A = A.copy()
    A_undir = utils.only_undirected(A)
    # Finding the nodes of the chain component of y
    chain_comp_y = utils.chain_component(y, A)
    # Getting the set C = NA_yx \ H
    C = na_yx - H
    # Case 1: x is a neighbor of y
    if x in utils.neighbors(y, A):
        # Getting all the nodes of the chain component in the order: C, x, y, ...
        nodes = list(C) + [x] + [y] + list(chain_comp_y - {y} - C)
        # Computing the perfect elimination ordering according to the above order
        ordering = utils.maximum_cardinality_search(A_undir, nodes)
    # Case 2: x is not a neighbor of y
    else:
        # Getting all the nodes of the chain component in the order: C, y, ...
        nodes = list(C) + [y] + list(chain_comp_y - {y} - C)
        # Computing the perfect elimination ordering according to the above order
        ordering = utils.maximum_cardinality_search(A_undir, nodes)
    # Orienting the edges of the chain component according to the perfect elimination ordering
    new_A = utils.orient_edges(A, ordering)
    # delete the edge between x and y
    new_A[x, y], new_A[y, x] = 0, 0
    return new_A


def score_valid_delete_operators(x, y, A, cache, debug=0):
    """Generate and score all valid delete(x,y,H) operators involving the edge
    x -> y or x - y, and all possible subsets H of neighbors of y which
    are adjacent to x.

    Parameters
    ----------
    x : int
        the "origin" node (i.e. x -> y or x - y)
    y : int
        the "target" node
    A : np.array
        the current adjacency matrix
    cache : instance of gies.scores.DecomposableScore
        the score cache to compute the score of the
        operators that are valid
    debug : int
        if larger than 0, debug are traces printed. Higher values
        correspond to increased verbosity

    Returns
    -------
    valid_operators : list of tuples
        a list of tubles, each containing a valid operator, its score
        and the resulting connectivity matrix

    """
    # Check inputs
    if A[x, y] == 0:
        raise ValueError("There is no (un)directed edge from x=%d to y=%d" % (x, y))
    # One-hot encode all subsets of H0, plus one column to mark if
    # they have already passed the validity condition
    na_yx = utils.na(y, x, A)
    H0 = sorted(na_yx)
    p = len(A)
    if len(H0) == 0:
        subsets = np.zeros((1, (p + 1)), dtype=np.bool)
    else:
        subsets = np.zeros((2 ** len(H0), (p + 1)), dtype=np.bool)
        subsets[:, H0] = utils.cartesian(
            [np.array([False, True])] * len(H0), dtype=np.bool
        )
    valid_operators = []
    print("    delete(%d,%d) H0=" % (x, y), set(H0)) if debug > 1 else None
    while len(subsets) > 0:
        print(
            "      len(subsets)=%d, len(valid_operators)=%d"
            % (len(subsets), len(valid_operators))
        ) if debug > 1 else None
        # Access the next subset
        H = np.where(subsets[0, :-1])[0]
        cond_1 = subsets[0, -1]
        subsets = subsets[1:]
        # Check if the validity condition holds for H, i.e. that
        # NA_yx \ H is a clique.
        # If it has not been tested previously for a subset of H,
        # check it now
        if not cond_1 and utils.is_clique(na_yx - set(H), A):
            cond_1 = True
            # For all supersets H' of H, the validity condition will also hold
            supersets = subsets[:, H].all(axis=1)
            subsets[supersets, -1] = True
        # If the validity condition holds, apply operator and compute its score
        print(
            "      delete(%d,%d,%s)" % (x, y, H),
            "na_yx - H = ",
            na_yx - set(H),
            "validity:",
            cond_1,
        ) if debug > 1 else None
        if cond_1:
            # Compute the change in score
            aux = (na_yx - set(H)) | utils.pa(y, A) | {x}
            # print(x,y,H,"na_yx:",na_yx,"old:",aux,"new:", aux - {x})
            old_score = cache.local_score(y, aux)
            new_score = cache.local_score(y, aux - {x})
            print(
                "        new: s(%d, %s) = %0.6f old: s(%d, %s) = %0.6f"
                % (y, aux - {x}, new_score, y, aux, old_score)
            ) if debug > 1 else None
            # Add to the list of valid operators
            valid_operators.append((new_score - old_score, x, y, H))
            print(
                "    delete(%d,%d,%s) -> %0.16f" % (x, y, H, new_score - old_score)
            ) if debug else None
    # Return all the valid operators
    return valid_operators


# --------------------------------------------------------------------
# Turn operator
#    1. definition in function turn
#    2. enumeration logic (to enumerate and score only valid
#    operators) function in score_valid_insert_operators


def turn(x, y, C, A, I):
    """
    Applies the turning operator: For an edge x - y or x <- y,
      1) Orients all edges of the chain component of x according to a perfect elimination ordering if x <- y
      2) Orients all edges of the chain component of y according to a perfect elimination ordering,
        such that for all c in C, the previously undirected edge c -> y
      3) orients the edge as x -> y

    Parameters
    ----------
    x : int
        the origin node (i.e. x -> y)
    y : int
        the target node
    C : iterable of ints
        a subset of the neighbors of y
    A : np.array
        the current adjacency matrix
    I : list of lists
        list of interventions


    Returns
    -------
    new_A : np.array
        the adjacency matrix resulting from applying the operator

    """
    # Check inputs
    if A[x, y] != 0 and A[y, x] == 0:
        raise ValueError("The edge %d -> %d is already exists" % (x, y))
    if A[x, y] == 0 and A[y, x] == 0:
        raise ValueError("x=%d and y=%d are not connected" % (x, y))
    if not C <= utils.neighbors(y, A):
        raise ValueError("Not all nodes in C=%s are neighbors of y=%d" % (C, y))
    if len({x, y} & C) > 0:
        raise ValueError("C should not contain x or y")
    # Apply operator
    C = set(C)
    new_A = A.copy()
    A_undir = utils.only_undirected(A)
    # Case 1: directed edge y -> x
    if A[x, y] == 0 and A[y, x] != 0:
        # Finding the nodes of the chain component of x
        chain_comp_x = utils.chain_component(x, A)
        # Getting all the nodes of the chain component in the order: x, ...
        nodes_x = [x] + list(chain_comp_x - {x})
        # Computing the perfect elimination ordering according to the above order
        ordering_x = utils.maximum_cardinality_search(A_undir, nodes_x)
        # Orienting the edges of the chain component according to the perfect elimination ordering of x
        new_A = utils.orient_edges(A, ordering_x)

        # Finding the nodes of the chain component of y
        chain_comp_y = utils.chain_component(y, A)
        # Getting all the nodes of the chain component in the order: C, y ...
        nodes_y = list(C) + [y] + list(chain_comp_y - {y} - C)
        # Computing the perfect elimination ordering according to the above order
        ordering_y = utils.maximum_cardinality_search(A_undir, nodes_y)
    # Case 2: undirected edge y - x
    else:
        # Finding the nodes of the chain component of y
        chain_comp_y = utils.chain_component(y, A)
        # Getting all the nodes of the chain component in the order: C, y, x ...
        nodes_y = list(C) + [y] + [x] + list(chain_comp_y - {y} - C - {x})
        # Computing the perfect elimination ordering according to the above order
        ordering_y = utils.maximum_cardinality_search(A_undir, nodes_y)
    # Orienting the edges of the chain component according to the perfect elimination ordering of y
    new_A = utils.orient_edges(new_A, ordering_y)
    # Turn edge x -> y
    new_A[y, x] = 0
    new_A[x, y] = 1
    return new_A


def score_valid_turn_operators(x, y, A, cache, debug=0):
    """Generate and score all valid turn(x,y,C) operators that can be
    applied to the edge x <- y or x - y, iterating through the valid
    subsets C of neighbors of y.

    Parameters
    ----------
    x : int
        the origin node (i.e. orient x -> y)
    y : int
        the target node
    A : np.array
        the current adjacency matrix
    cache : instance of gies.scores.DecomposableScore
        the score cache to compute the score of the
        operators that are valid
    debug : int
        if larger than 0, debug are traces printed. Higher values
        correspond to increased verbosity

    Returns
    -------
    valid_operators : list of tuples
        a list of tubles, each containing a valid operator, its score
        and the resulting connectivity matrix

    """
    if A[x, y] != 0 and A[y, x] == 0:
        raise ValueError("The edge %d -> %d already exists" % (x, y))
    if A[x, y] == 0 and A[y, x] == 0:
        raise ValueError("x=%d and y=%d are not connected" % (x, y))
    # Different validation/scoring logic when the edge to be turned is
    # essential (x <- x) or not (x - y)
    if A[x, y] != 0 and A[y, x] != 0:
        return score_valid_turn_operators_undir(x, y, A, cache, debug=debug)
    else:
        return score_valid_turn_operators_dir(x, y, A, cache, debug=debug)


def score_valid_turn_operators_dir(x, y, A, cache, debug=0):
    """Logic for finding and scoring the valid turn operators that can be
    applied to the edge x <- y.

    Parameters
    ----------
    x : int
        the origin node (i.e. x -> y)
    y : int
        the target node
    A : np.array
        the current adjacency matrix
    cache : instance of gies.scores.DecomposableScore
        the score cache to compute the score of the
        operators that are valid
    debug : bool or string
        if debug traces should be printed (True/False). If a non-empty
        string is passed, traces are printed with the given string as
        prefix (useful for indenting the prints from a calling
        function)

    Returns
    -------
    valid_operators : list of tuples
        a list of tubles, each containing a valid operator, its score
        and the resulting connectivity matrix

    """
    # One-hot encode all subsets of T0, plus one extra column to mark
    # if they pass validity condition 2 (see below). The set C passed
    # to the turn operator will be C = NAyx U T.
    p = len(A)
    T0 = sorted(utils.neighbors(y, A) - utils.adj(x, A))
    if len(T0) == 0:
        subsets = np.zeros((1, p + 1), dtype=np.bool)
    else:
        subsets = np.zeros((2 ** len(T0), p + 1), dtype=np.bool)
        subsets[:, T0] = utils.cartesian(
            [np.array([False, True])] * len(T0), dtype=np.bool
        )
    valid_operators = []
    print("    turn(%d,%d) T0=" % (x, y), set(T0)) if debug > 1 else None
    while len(subsets) > 0:
        print(
            "      len(subsets)=%d, len(valid_operators)=%d"
            % (len(subsets), len(valid_operators))
        ) if debug > 1 else None
        # Access the next subset
        T = np.where(subsets[0, :-1])[0]
        passed_cond_2 = subsets[0, -1]
        subsets = subsets[1:]  # update the list of remaining subsets
        # Check that the validity conditions hold for T
        C = utils.na(y, x, A) | set(T)
        # Condition 1: Test that C = NA_yx U T is a clique
        cond_1 = utils.is_clique(C, A)
        if not cond_1:
            # Remove from consideration all other sets T' which
            # contain T, as the clique condition will also not hold
            supersets = subsets[:, T].all(axis=1)
            subsets = utils.delete(subsets, supersets, axis=0)
        # Condition 2: Test that all semi-directed paths from y to x contain a
        # member from C U neighbors(x)
        if passed_cond_2:
            # If a subset of T satisfied condition 2, so does T
            cond_2 = True
        else:
            # otherwise, check condition 2
            cond_2 = True
            for path in utils.semi_directed_paths(y, x, A):
                if path == [y, x]:
                    pass
                elif len((C | utils.neighbors(x, A)) & set(path)) == 0:
                    cond_2 = False
                    break
            if cond_2:
                # If condition 2 holds for C U neighbors(x), that is,
                # for C = NAyx U T U neighbors(x), then it holds for
                # all supersets of T
                supersets = subsets[:, T].all(axis=1)
                subsets[supersets, -1] = True
        # If both conditions hold, apply operator and compute its score
        print(
            "      turn(%d,%d,%s)" % (x, y, C),
            "na_yx =",
            utils.na(y, x, A),
            "T =",
            T,
            "validity:",
            cond_1,
            cond_2,
        ) if debug > 1 else None
        if cond_1 and cond_2:
            # Compute the change in score
            new_score = cache.local_score(
                y, utils.pa(y, A) | C | {x}
            ) + cache.local_score(x, utils.pa(x, A) - {y})
            old_score = cache.local_score(y, utils.pa(y, A) | C) + cache.local_score(
                x, utils.pa(x, A)
            )
            print(
                "        new score = %0.6f, old score = %0.6f, y=%d, C=%s"
                % (new_score, old_score, y, C)
            ) if debug > 1 else None
            # Add to the list of valid operators
            valid_operators.append((new_score - old_score, x, y, C))
            print(
                "    turn(%d,%d,%s) -> %0.16f" % (x, y, C, new_score - old_score)
            ) if debug else None
    # Return all the valid operators
    return valid_operators


def score_valid_turn_operators_undir(x, y, A, cache, debug=0):
    """Logic for finding and scoring the valid turn operators that can be
    applied to the edge x - y.

    Parameters
    ----------
    x : int
        the origin node (i.e. x -> y)
    y : int
        the target node
    A : np.array
        the current adjacency matrix
    cache : instance of gies.scores.DecomposableScore
        the score cache to compute the score of the
        operators that are valid
    debug : bool or string
        if debug traces should be printed (True/False). If a non-empty
        string is passed, traces are printed with the given string as
        prefix (useful for indenting the prints from a calling
        function)

    Returns
    -------
    valid_operators : list of tuples
        a list of tubles, each containing a valid operator, its score
        and the resulting connectivity matrix

    """
    # Proposition 31, condition (ii) in GIES paper (Hauser & Bühlmann
    # 2012) is violated if:
    #   1. all neighbors of y are adjacent to x, or
    #   2. y has no neighbors (besides u)
    # then there are no valid operators.
    non_adjacents = list(utils.neighbors(y, A) - utils.adj(x, A) - {x})
    if len(non_adjacents) == 0:
        print(
            "    turn(%d,%d) : ne(y) \\ adj(x) = Ø => stopping" % (x, y)
        ) if debug > 1 else None
        return []
    # Otherwise, construct all the possible subsets which will satisfy
    # condition (ii), i.e. all subsets of neighbors of y with at least
    # one which is not adjacent to x
    p = len(A)
    C0 = sorted(utils.neighbors(y, A) - {x})
    subsets = np.zeros((2 ** len(C0), p + 1), dtype=np.bool)
    subsets[:, C0] = utils.cartesian([np.array([False, True])] * len(C0), dtype=np.bool)
    # Remove all subsets which do not contain at least one non-adjacent node to x
    to_remove = (subsets[:, non_adjacents] == False).all(axis=1)
    subsets = utils.delete(subsets, to_remove, axis=0)
    # With condition (ii) guaranteed, we now check conditions (i,iii)
    # for each subset
    valid_operators = []
    print("    turn(%d,%d) C0=" % (x, y), set(C0)) if debug > 1 else None
    while len(subsets) > 0:
        print(
            "      len(subsets)=%d, len(valid_operators)=%d"
            % (len(subsets), len(valid_operators))
        ) if debug > 1 else None
        # Access the next subset
        C = set(np.where(subsets[0, :])[0])
        subsets = subsets[1:]
        # Condition (i): C is a clique in the subgraph induced by the
        # chain component of y. Because C is composed of neighbors of
        # y, this is equivalent to C being a clique in A. NOTE: This
        # is also how it is described in Alg. 5 of the paper
        cond_1 = utils.is_clique(C, A)
        if not cond_1:
            # Remove from consideration all other sets C' which
            # contain C, as the clique condition will also not hold
            supersets = subsets[:, list(C)].all(axis=1)
            subsets = utils.delete(subsets, supersets, axis=0)
            continue
        # Condition (iii): Note that condition (iii) from proposition
        # 31 appears to be wrong in the GIES paper; instead we use the
        # definition of condition (iii) from Alg. 5 of the paper:
        # Let na_yx (N in the GIES paper) be the neighbors of Y which
        # are adjacent to X. Then, {x,y} must separate C and na_yx \ C
        # in the subgraph induced by the chain component of y,
        # i.e. all the simple paths from one set to the other contain
        # a node in {x,y}.
        subgraph = utils.induced_subgraph(utils.chain_component(y, A), A)
        na_yx = utils.na(y, x, A)
        if not utils.separates({x, y}, C, na_yx - C, subgraph):
            continue
        # At this point C passes both conditions
        #   Compute the change in score
        new_score = cache.local_score(y, utils.pa(y, A) | C | {x}) + cache.local_score(
            x, utils.pa(x, A) | (C & na_yx)
        )
        old_score = cache.local_score(y, utils.pa(y, A) | C) + cache.local_score(
            x, utils.pa(x, A) | (C & na_yx) | {y}
        )
        print(
            "        new score = %0.6f, old score = %0.6f, y=%d, C=%s"
            % (new_score, old_score, y, C)
        ) if debug > 1 else None
        #   Add to the list of valid operators
        valid_operators.append((new_score - old_score, x, y, C))
        print(
            "    turn(%d,%d,%s) -> %0.16f" % (x, y, C, new_score - old_score)
        ) if debug else None
    # Return all valid operators
    return valid_operators






# Copyright 2022 Olga Kolotuhina, Juan L. Gamella

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:

# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.

# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.

# 3. Neither the name of the copyright holder nor the names of its
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

"""
"""

import numpy as np
from decomposable_score import DecomposableScore

# --------------------------------------------------------------------
# l0-penalized Gaussian log-likelihood score for a sample from a single
# (observational) environment


class GaussIntL0Pen(DecomposableScore):
    """
    Implements a cached l0-penalized gaussian likelihood score for the GIES setting.

    """

    def __init__(self, data, interv, lmbda=None, cache=True, debug=0,bn=None):
        """Creates a new instance of the class.

        Parameters
        ----------
        data : list of numpy.ndarray
            every matrix in the list corresponds to an environment,
            the nxp matrix containing the observations of each
            variable (each column corresponds to a variable).
        interv: a list of lists
            a list of the interventions sets which
            corresponds to the environments in data
        lmbda : float or NoneType, optional
            the regularization parameter. If None, defaults to the BIC
            score, i.e. lmbda = 1/2 * log(n), where n is the number of
            observations.
        cache : bool, optional
           if computations of the local score should be cached for
           future calls. Defaults to True.
        debug : int, optional
            if larger than 0, debug are traces printed. Higher values
            correspond to increased verbosity.

        """
        super().__init__(data, interv, cache=cache, debug=debug)
        self.p = self._data[0].shape[1]
        self.n_obs = np.array([len(env) for env in self._data])
        # Computing the sample covariances
        self._data = [sample - sample.mean(axis=0) for sample in self._data]
        self.sample_cov = np.array(
            [1 / self.n_obs[ind] * env.T @ env for (ind, env) in enumerate(self._data)]
        )
        # Discarded options for computing the sample covariances
        #  a) This is different to how it is computed in the PCALG pacakge
        # self.sample_cov = np.array([np.cov(env, rowvar=False, ddof=0) for env in self._data])
        #  c) This also sometimes yields outputs which are different to what is computed in PCALG
        # sample_cov = []
        # for i, env in enumerate(self._data):
        #     mean = np.mean(env, axis=0)
        #     aux = env - mean
        #     sample_cov.append(1 / self.n_obs[i] * aux.T @ aux)
        # self.sample_cov = np.array(sample_cov)
        self.bn = bn 
        self.N = sum(self.n_obs)
        self.lmbda = 0.5 * np.log(self.N) if lmbda is None else lmbda
        self.num_not_interv = np.zeros(self.p)
        self.part_sample_cov = np.zeros((self.p, self.p, self.p))

        # Check that the interventions form a conservative family of targets
        for j in range(self.p):
            if sum(i.count(j) for i in self.interv) == len(self._data):
                raise ValueError("The family of targets is not conservative")

        # Computing the numbers of non-interventions of a variable and the corresponding partial covariance matrix
        for k in range(self.p):
            for (i, n) in enumerate(self.n_obs):
                if k not in set(self.interv[i]):
                    self.num_not_interv[k] += n
                    self.part_sample_cov[k] += self.sample_cov[i] * n
            self.part_sample_cov[k] = self.part_sample_cov[k] / self.num_not_interv[k]

    def full_score(self, A):
        """
        Given a DAG adjacency A, return the l0-penalized log-likelihood of
        a sample from a single environment, by finding the maximum
        likelihood estimates of the corresponding connectivity matrix
        (weights) and noise term variances.

        Parameters
        ----------
        A : np.array
            The adjacency matrix of a DAG, where A[i,j] != 0 => i -> j.

        Returns
        -------
        score : float
            the penalized log-likelihood score.

        """
        print('hi')
        # Compute MLE
        B, omegas = self._mle_full(A)
        likelihood = 0
        for j, sigma in enumerate(self.part_sample_cov):
            gamma = 1 / omegas[j]
            likelihood += self.num_not_interv[j] * (
                np.log(gamma)
                - gamma
                * (np.eye(self.p) - B)[:, j]
                @ sigma
                @ (np.eye(self.p) - B)[:, j].T
            )
        l0_term = self.lmbda * (np.sum(A != 0) + self.p)
        score = 0.5 * likelihood - l0_term
        return score

    # Note: self.local_score(...), with cache logic, already defined
    # in parent class DecomposableScore.

    def _compute_local_score(self, k, pa):
        #print('*')
        """
        Given a node and its parents, return the local l0-penalized
        log-likelihood of a sample from a single environment, by finding
        the maximum likelihood estimates of the weights and noise term
        variances.

        Parameters
        ----------
        x : int
            a node.
        pa : set of ints
            the node's parents.

        Returns
        -------
        score : float
            the penalized log-likelihood score.

        """
        #print('hi***  ', k)
        #print(pa)
        #print(self.bn.parents(int(k)))
        #print(self.interv)
        #print('***********')
        pa = list(pa)
        # Compute MLE
        b, sigma = self._mle_local(k, pa)
        # Compute log-likelihood (without log(2π) term)
        n = self.num_not_interv[k]
        likelihood = -0.5 * n * (1 + np.log(sigma))
        #  Note: the number of parameters is the number of parents (one
        #  weight for each) + the marginal variance of x
        l0_term = self.lmbda * (len(pa) + 1)
        score = likelihood - l0_term
        #return score
        arcs = set([])
        Pa_k = set([])
        for v in self.interv:
           tmp = self.interv
           arcs = arcs | Learn_cut(v,self.bn.arcs(),self._data[tmp.index(v)])
        for a in arcs:
            if (a[1] == k ):
                             Pa_k.add(a[0]) 
        return (-len(set(pa)-Pa_k) - len(Pa_k-set(pa)))
    # --------------------------------------------------------------------
    #  Functions for the maximum likelihood estimation of the
    #  weights/variances

    def _mle_full(self, A):
        """
        Finds the maximum likelihood estimate for the whole graph,
        specified by the adjacency A.

        Parameters
        ----------
        A : np.array
            The adjacency matrix of a DAG, where A[i,j] != 0 => i -> j.

        Returns
        -------
        B : np.array
            the connectivity (weights) matrix, which respects the
            adjacency in A.
        omegas : np.array
            the estimated noise-term variances of the observed
            variables.

        """
        B = np.zeros(A.shape)
        omegas = np.zeros(self.p)
        for j in range(self.p):
            parents = np.where(A[:, j] != 0)[0]
            B[:, j], omegas[j] = self._mle_local(j, parents)
        return B, omegas

    def _mle_local(self, k, pa):
        """Finds the maximum likelihood estimate of the local model
        between a node and its parents.

        Parameters
        ----------
        x : int
            a node.
        pa : set of ints
            the node's parents.

        Returns
        -------
        b : np.array
            an array of size p, with the estimated weights from the
            parents to the node, and zeros for non-parent variables.
        sigma : float
            the estimate noise-term variance of variable x

        """
        pa = list(pa)
        b = np.zeros(self.p)
        S_k = self.part_sample_cov[k]
        S_kk = S_k[k, k]
        S_pa_k = S_k[pa, :][:, k]
        b[pa] = _regress(k, pa, S_k)
        sigma = S_kk - b[pa] @ S_pa_k
        return b, sigma


def _regress(j, pa, cov):
    # compute the regression coefficients from the
    # empirical covariance (scatter) matrix i.e. b =
    # Σ_{j,pa(j)} @ Σ_{pa(j), pa(j)}^-1
    return np.linalg.solve(cov[pa, :][:, pa], cov[j, pa])





def Learn_cut(X,Edges,data):
     import numpy as np
     #data = block_sample_intervention(bn,X,SAMPLES)
     #data = np.array(data,dtype=int)
     #print(data.shape)
     arcs = set([])
     Cut_X = []
     for e in Edges:
          if (e[0] in X and e[1] not in X):
                                        Cut_X.append((e[0],e[1]))
          if (e[0] not in X and e[1] in X):
                                       Cut_X.append((e[1],e[0]))
     index = 100               
     if (index >= 5):
          #print(j)
          for e in Cut_X: 
              #print(index)
              #print
              from causallearn.utils.cit import CIT
              chisq_obj = CIT(data, "chisq") # construct a CIT instance with data and method name
              pValue = chisq_obj(e[0], e[1])
              #tmp = mutual_info_score(data[0:index,e[0]],data[0:index,e[1]] )
              if (pValue >= 0.01):
                                narc = (e[1],e[0])  
              else:                     
                                narc = (e[0],e[1])  

              arcs.add(narc)
              
     #print(arcs) 
    
     return arcs  





def GIES_discovery(Num_Grphs,nodes,degree,Max_samples,Gap):
    
            import numpy as np
            import pandas as pd
            import gies 
            import numpy as np
            np.bool = np.bool_
            data_save = []
            # Assuming you have already defined and trained your Bayesian network with interventions

            
            import pyAgrum as gum
            import pyAgrum.causal as c
            
            import numpy as np
            for grphs_no in range(0,Num_Grphs, 1):
                        print('Graph No',grphs_no)
                        a = shanmugam_random_chordal(nodes,degree)
                        adjacency_list = list(a.edges)
                        graph_dict = convert(adjacency_list)
                        r = greedyColoring(graph_dict,len(graph_dict))
                        I = indices_of_elements(r)
                        #print(I)
                        tmp = adj_list_to_string_with_vertices(adjacency_list)
                        
                        
                        bn = gum.fastBN(tmp)
                        bn.generateCPTs()
                        Pv =1
                        for n in bn.nodes():
                           Pv = Pv * bn.cpt(n)
                        
                 
                           
                        
                        import graphical_models as gm
                        dag = gm.DAG()
                        
                        dag.add_arcs_from(bn.arcs(), check_acyclic=True)
                        bn.arcs()
            
            
            
            
            
                        Truedag = bn.arcs();
                        
                        #print('*')
                        ao = np.zeros((len(bn.nodes()),len(bn.nodes())))
                        for e in Truedag: 
                            ao[e[0]][e[1]] = 1
                            ao[e[1]][e[0]] = 1
                        
                        
                        I = indices_of_elements(r)    
                        for iii in range(1):
                          
                            #bn.generateCPTs()
                            List_O =[]
                            for i in I:
                              d1= block_sample_intervention(bn,i,100000)
                              d1 = np.array(d1,dtype=int)
                              List_O.append(d1)
                        
                              
                            
                        
                            
                            shd = np.zeros(Max_samples)
                            samples  =0
                            for jj in range(Gap,Max_samples,Gap):
                                    List =[]
                                    for i in I:
                                        d1 = List_O[I.index(i)]
                                        List.append(d1[0:jj,:])
                                        
                                    data = np.stack(List)
                                    estimate, score = fit_bic(data, I,bn=bn)
                                    est = set([])
                                    for i in range(len(bn.nodes())):
                                          for j in range(len(bn.nodes())):
                                              if estimate[i,j] == 1 and ( (i,j) in  Truedag or (j,i) in  Truedag ):
                                                                    est.add((i,j)) 
                        
                                    #print(jj)
                                    #print(len(Truedag-est)+len(est-Truedag))
                                    #print(shd)
                                    shd[samples:samples + Gap*len(I)] = (len(Truedag-est)+len(est-Truedag))/2 
                                    samples = samples + Gap*len(I) 
                                    if (samples >= (Max_samples+100)):
                                                                   break 
                                    if (len(Truedag-est)+len(est-Truedag)  == 0 and (jj>= 2000)):
                                                                                  break
                            #print(jj,'***',len(I))
                            shd[0:200] = len(bn.arcs())
                            data_save.append(shd)
                            #print(estimate)
                            #print(est)
                            #print((len(Truedag-est)+len(est-Truedag))/2)
                            #print(iii)
                            #print(shd)        
                            
                         
                        
                            #print('*************')
            return data_save    
                        
                        
                        
                                 
