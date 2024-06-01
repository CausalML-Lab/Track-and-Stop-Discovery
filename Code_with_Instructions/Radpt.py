import pyAgrum as gum


def is_chordal(graph):
    """
    Check if a graph is chordal (clique tree).
    """
    for node in graph:
        neighbors = graph[node]
        for i in range(len(neighbors)):
            for j in range(i + 1, len(neighbors)):
                if neighbors[j] not in neighbors[i] and neighbors[i] not in neighbors[j]:
                    return False
    return True

def perfect_elimination_ordering(graph):
   

    peo = []
    remaining_nodes = list(graph.keys())

    while remaining_nodes:
        min_vertex = min(remaining_nodes, key=lambda v: len(set(remaining_nodes) & set(graph[v])))
        peo.append(min_vertex)
        remaining_nodes.remove(min_vertex)
        # Update the graph by removing min_vertex and making it a clique with its neighbors.
        for neighbor in graph[min_vertex]:
            graph[neighbor] = list(set(graph[neighbor]) | set(graph[min_vertex]))

    return peo

# Example usage:
if __name__ == "__main__":
    # Example adjacency list for a chordal graph.
    graph = {
        0: [1],
        1: [0],
        
    }

    peo = perfect_elimination_ordering(graph)
    
from causaldag import DAG
import random
import networkx as nx
import numpy as np

from collections import defaultdict
from itertools import combinations
from networkx.algorithms import chordal_graph_cliques

import math
import sys
sys.path.insert(0, './PADS')


'''
Get (directed) adjacency list from networkx graph
'''
def get_adj_list(nx_G):
    adj_list = dict()
    for node, nbrdict in nx_G.adjacency():
        adj_list[node] = [nbr for nbr in nbrdict.keys()]
    return adj_list

'''
Compute clique tree of a chordal graph
'''
def compute_clique_tree(nx_G):
 

    # Compute clique graph H
    # Each node is a maximal clique of G
    # Each edge has weight equal to size of intersection between two maximal cliques
    maximal_cliques = list(chordal_graph_cliques(nx_G))
    H = nx.Graph()
    for i in range(len(maximal_cliques)):
        H.add_node(i, mc=maximal_cliques[i])
    for u,v in combinations(H.nodes(), 2):
        wt = len(H.nodes[u]['mc'].intersection(H.nodes[v]['mc']))
        H.add_edge(u, v, weight=wt)

    # Compute clique tree by taking maximum spanning tree of H
    T = nx.maximum_spanning_tree(H)

    return T

'''
Compute maximum independent set via PEO in O(n + m) time.

References:
Fănică Gavril. Algorithms for minimum coloring, maximum clique, minimum covering by cliques, and maximum independent set of a chordal graph. SIAM Journal on Computing, 1972
Jospeh Y-T Leung. Fast algorithms for generating all maximal independent sets of interval, circular-arc and chordal graphs. Journal of Algorithms, 1984
'''
def maximum_independent_set(nx_G):
    

    adj_list = get_adj_list(nx_G)
    G = dict()
    for v in nx_G.nodes():
        G[v] = adj_list[v]
    peo = perfect_elimination_ordering(G)
    #print(peo)
    # Sanity check: For every vertex v, v and its later neighbors form a clique

    for idx in range(len(peo)):
        v = peo[idx]
        later_nbrs = []
        for later_idx in range(idx+1, len(peo)):
            w = peo[later_idx]
            if w in adj_list[v]:
                later_nbrs.append(w)
        

    actual_to_peo = dict()
    for idx in range(len(peo)):
        actual_to_peo[peo[idx]] = idx

    S = set()
    for idx in range(len(peo)):
        v = peo[idx]
        no_earlier_nbr_in_S = True
        for w in adj_list[v]:
            if actual_to_peo[w] < idx and w in S:
                no_earlier_nbr_in_S = False
                break
        if no_earlier_nbr_in_S:
            S.add(v)

    # Sanity check: S is an independent set
    for u in S:
        for v in S:
            if u == v:
                continue
          

    # Sanity check: S is a maximal independent set
    for v in nx_G.nodes():
        if v not in S:
            has_neighbor_in_S = False
            for w in adj_list[v]:
                if w in S:
                    has_neighbor_in_S = True
                    break
     
    return S

'''
Compute minimum vertex cover by removing maximum independent set from vertex set.
'''
def minimum_vertex_cover(nx_G):
   
    maxIS = maximum_independent_set(nx_G)
    minVC = set(nx_G.nodes()).difference(maxIS)

    # Sanity check: minVC is vertex cover
 

    # Sanity check: minVC is minimal
    required = set()
    for u,v in nx_G.edges():
        if u in minVC and v not in minVC:
            required.add(u)
        if u not in minVC and v in minVC:
            required.add(v)
   

    return minVC

'''
Refactored from [CSB22]

The following few lines are copied from the proof of Lemma 1
Let n = p_d a^d + r_d and n = p_{d-1} a^{d-1} + r_{d-1}
1) Repeat 0 a^{d-1} times, repeat the next integer 1 a^{d-1} times and so on circularly from {0,1,...,a-1} till p_d * a^d.
2) After that, repeat 0 ceil(r_d/a) times followed by 1 ceil(r_d/a) times till we reach the nth position. Clearly, n-th integer in the sequence would not exceed a-1.
3) Every integer occurring after the position a^{d-1} p_{d-1} is increased by 1.
'''
def intervention_subroutine(A: set, k: int):
    #print(A,k)

    if k == 1:
        return [frozenset({v}) for v in A]
    else:
        # Setup parameters. Note that [SKDV15] use n and x+1 instead of h and L
        h = len(A)
        k_prime = min(k, h/2)
        a = math.ceil(h/k_prime)
        
        #print(h,a)  
        
        L = math.ceil(math.log(h,a))
    

        # Execute labelling scheme
        S = defaultdict(set)
        for d in range(1, L+1):
            a_d = pow(a,d)
            r_d = h % a_d
            p_d = h // a_d
            a_dminus1 = pow(a,d-1)
            r_dminus1 = h % a_dminus1 # Unused
            p_dminus1 = h // a_dminus1
            
            B = list(A)
            for i in range(1, h+1):
                node = B[i-1]
                if i <= p_d * a_d:
                    val = (i % a_d) // a_dminus1
                else:
                    val = (i - p_d * a_d) // math.ceil(r_d / a)
                if i > a_dminus1 * p_dminus1:
                    val += 1
                S[(d,val)].add(node)
        return [frozenset(x) for x in S.values()]

'''
Balanced partitioning on trees
'''
def tree_balanced_partitioning(nx_G, L: int):
  

    adj_list = get_adj_list(nx_G)
    V = set(nx_G.nodes())
    n = len(V)
    A = set()

    while len(V) > np.ceil(n/(L+1)):
        # Consider the subgraph on remaining vertices
        H = nx_G.subgraph(V)
       

        # Root arbitrarily
        root = list(V)[0]

        # Perform DFS
        dfs_preorder = list(nx.dfs_preorder_nodes(H, source=root))
        dfs_index = dict()
        for v in dfs_preorder:
            dfs_index[v] = len(dfs_index)
        children = defaultdict(set)
        for v in V:
            for w in adj_list[v]:
                if w in V and dfs_index[w] > dfs_index[v]:
                    children[v].add(w)

        # Compute sizes of subtrees T_u for each node u
        # Traverse in reverse DFS ordering
        subtree_size = dict()
        for v in dfs_preorder[::-1]:
            subtree_size[v] = 1
            for w in children[v]:
                subtree_size[v] += subtree_size[w]

        # if-else cases
        u = None
        for v in V:
            # Find a subtree T_u such that |T_u| = 1 + ceil(n/(L+1))
            if subtree_size[v] == 1 + np.ceil(n/(L+1)):
                u = v
                break
        if u is None:
            # Find a subtree T_u such that |T_u| > 1 + ceil(n/(L+1)) and all children w have |T_w| <= ceil(n/(L+1))
            # Start from root and recurse (see proof)
            done = False
            v = root
            while not done:
                # Check if all children w have |T_w| <= ceil(n/(L+1))
                ok = True
                for w in children[v]:
                    if subtree_size[w] > np.ceil(n/(L+1)):
                        v = w
                        ok = False
                        break
                if ok:
                    u = v
                    done = True

        # Add u to A
        
        A.add(u)

        # Remove T_u from V
        old_V_size = len(V)
        def remove_subtree_from_V(node):
            V.remove(node)
            for w in children[node]:
                remove_subtree_from_V(w)
        remove_subtree_from_V(u)
        new_V_size = len(V)

        # Sanity check: |V| drops by subtree_size[u]
   
    # Sanity check: |A| <= L and subtrees in G[V\A] have size <= ceil(n/(L+1))
    G_check = nx_G.subgraph(set(nx_G.nodes()).difference(A))
   

    return A

'''
Modified from [CSB22]'s separator_policy
'''
def r_adaptive_policy(dag: DAG, r: int, k: int, verbose: bool = False,bn= None,t_samples = None,data = None) -> set:
    I_count = 0

    intervened_nodes = set()
    current_cpdag = dag.cpdag()
    if verbose: print(f"Remaining edges: {current_cpdag.num_edges}")

    # If we are given r > log_2(n), we will use the additional adaptivity to check whether we should skip interventions
    # The current implementation of [CSB22] is actually n-adaptive since it ALWAYS performs such checks before intervening
    n = dag.nnodes
    L = np.ceil(np.power(n, 1/r))
    checking_budget = max(0, r - int(np.ceil(np.log2(n))))
    r -= checking_budget

    # Subroutine to extract essential graph and then remove oriented arcs
    def get_essential_graph_without_oriented_arcs():
        nonlocal current_cpdag

        G = nx.Graph()
        G.add_nodes_from(current_cpdag.nodes)
        G.add_edges_from(current_cpdag.edges)
        return G

    # Subroutine to perform interventions
    def perform_one_round_of_interventions(intervention_set): 
        nonlocal intervened_nodes
        nonlocal current_cpdag
        nonlocal checking_budget

        for intervention in intervention_set:
            
            relevant_vertices = [v for v in intervention]
            if checking_budget > 0:
                G = get_essential_graph_without_oriented_arcs()
                relevant_vertices = [v for v in intervention if G.degree[v] > 0]
            checking_budget -= 1
            if len(relevant_vertices) > 0:
                intervened_nodes.add(frozenset(relevant_vertices))
                node = intervention
                #print('I am intervened ', node)
                ARCS = Learn_cut(bn,node,bn.arcs(),t_samples,data[list(node)[0]])
        
                
                dag1 = dag.copy()
                dag1.add_arcs_from(ARCS, check_acyclic=False)
                
                        
                        
                current_cpdag = current_cpdag.interventional_cpdag(dag1, intervention)

    # Loop for (up to) r-1 rounds
    while r > 1:
        r -= 1
        if current_cpdag.num_arcs != dag.num_arcs:
            G = get_essential_graph_without_oriented_arcs()
            intervention_set = set()
            for cc_nodes in nx.connected_components(G):
                if len(cc_nodes) > 1:
                    H = G.subgraph(cc_nodes).copy()

                    if H.size() == len(cc_nodes) * (len(cc_nodes) - 1) / 2:
                        # If H is a clique, add all vertices to intervention set
                        intervention_set.update(cc_nodes)
                    else:
                        # Compute clique tree T_H
                        T_H = compute_clique_tree(H)

                        # Compute L-balanced partitioning on T_H
                        for mc_node in tree_balanced_partitioning(T_H, L):
                            # Add all vertices involved to the maximal clique
                            intervention_set.update(T_H.nodes()[mc_node]["mc"])

            # Non-adaptive alternative (treat as final round, even though we may have more adaptivity)
            # Compute G-separating system on remaining relevant vertices
            # For atomic interventions, this corresponds to the minimum vertex cover
            relevant_nodes = set()
            for cc_nodes in nx.connected_components(G):
                if len(cc_nodes) > 1:
                    relevant_nodes.update(cc_nodes)
            minVC = minimum_vertex_cover(G.subgraph(relevant_nodes))
            

            # Perform a round of intervention (choose the cheaper option)
            if len(minVC) <= len(intervention_set):
                # Treating as final round and assign all remaining budget for checking
                checking_budget += r
                r = 0
                if (len(minVC) ==0):
                    break
                perform_one_round_of_interventions(intervention_subroutine(minVC, k))
                I_count = I_count + len(intervention_subroutine(minVC, k))/2
            else:
                if (len(intervention_set) ==0):
                    break
                perform_one_round_of_interventions(intervention_subroutine(intervention_set, k))
                I_count = I_count + len(intervention_subroutine(minVC, k))/2
    # If still not fully oriented, perform a single round of non-adaptive interventions

    if current_cpdag.num_arcs != dag.num_arcs:

        G = get_essential_graph_without_oriented_arcs()

        # Compute G-separating system on remaining relevant vertices
        # For atomic interventions, this corresponds to the minimum vertex cover
        relevant_nodes = set()
        for cc_nodes in nx.connected_components(G):
            # Sanity check: Number of maximal cliques in each connected component at the final non-adaptive round is at most L = ceil(n^(1/r))
            T = compute_clique_tree(G.subgraph(cc_nodes))
            
            if len(cc_nodes) > 1:
                relevant_nodes.update(cc_nodes)
        minVC = minimum_vertex_cover(G.subgraph(relevant_nodes))
 
        if(len(minVC) != 0 ): 
                 
        # Perform a round of intervention
                 perform_one_round_of_interventions(intervention_subroutine(minVC, k))
                 I_count = I_count + len(intervention_subroutine(minVC, k))

    import math 
    return math.ceil(I_count),current_cpdag




def Learn_cut(bn,X,Edges,SAMPLES,data):
     import numpy as np
     #data = block_sample_intervention(bn,X,SAMPLES)
     data = np.array(data,dtype=int)
     #print(data.shape)
     arcs = set([])
     Cut_X = []
     for e in Edges:
          if (e[0] in X and e[1] not in X):
                                        Cut_X.append((e[0],e[1]))
          if (e[0] not in X and e[1] in X):
                                       Cut_X.append((e[1],e[0]))
     index = SAMPLES               
     if (index >= 5):
          #print(j)
          for e in Cut_X: 
              #print(index)
              #print
              from causallearn.utils.cit import CIT
              chisq_obj = CIT(data[0:index,:], "chisq") # construct a CIT instance with data and method name
              pValue = chisq_obj(e[0], e[1])
              #tmp = mutual_info_score(data[0:index,e[0]],data[0:index,e[1]] )
              if (pValue >= 0.005):
                                narc = (e[1],e[0])  
              else:                     
                                narc = (e[0],e[1])  

              arcs.add(narc)
              
     #print(arcs) 
    
     return arcs     
     
     
     
     
     
def block_sample_intervention(bn,intv,SAMPLES):
    import numpy as np

    bn1  = gum.BayesNet(bn)
    import numpy as np
    for j in intv:
        for i in bn1.parents(j):
            bn1.eraseArc(gum.Arc(i,j))
       
       
        bn1.cpt(j)[:] = [0.5 ,0.5]

           
                                        
            
    df = gum.generateSample(bn1, n=SAMPLES, name_out=None, show_progress=False, with_labels=True, random_order=False)
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
    
    
    
    
    
import networkx as nx
import itertools as itr
import random
from typing import Union, List
import numpy as np
from graph_utils import fill_vstructures, direct_chordal_graph
from scipy.special import binom


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
    
    
    
    
    
def r_apative(Num_Grphs,nodes,degree,Max_samples,Gap):

    
            import pyAgrum as gum
            import pyAgrum.causal as c
            
            import numpy as np
            Data_save=[]
            samples_list = [] 
    
            
            for grphs_no in range(0,Num_Grphs, 1):
                    print('Garph No.',grphs_no)
                    a = shanmugam_random_chordal(nodes,degree)
                    adjacency_list = list(a.edges)
                    graph_dict = convert(adjacency_list)
                    r = greedyColoring(graph_dict,len(graph_dict))
                    I = indices_of_elements(r)
                    tmp = adj_list_to_string_with_vertices(adjacency_list)
                    
                    bn = gum.fastBN(tmp)
                    
                    bn.generateCPTs()
                    Pv =1
                    for n in bn.nodes():
                       Pv = Pv * bn.cpt(n)
                    
                    import graphical_models as gm
                    dag = gm.DAG()
                    dag.add_arcs_from(bn.arcs(), check_acyclic=True)
                    len(bn.arcs())
                   
                    import graphical_models as gm
                    dag = gm.DAG()
                    
                    dag.add_arcs_from(bn.arcs(), check_acyclic=True)
                    
                    for j in range(1):
                              shd = np.zeros(int(Max_samples))
                              random.seed(897986274)
                              bn.generateCPTs()
                              data =[]
                              for x in bn.nodes():
                                 data.append(block_sample_intervention(bn,{x},100000))
                              t_index =0    
                              for t in range(Gap,Max_samples,Gap):
                                        t_samples = t
                                        
                                        I_count, cpdag = r_adaptive_policy(dag, 4, 1, verbose =  False,bn = bn,t_samples = t_samples,data=data)   
                                      
                                        #print(cpdag.arcs)
                                        #print(bn.arcs())
                                        
                                        #print(t)
                                        #print(shd)
                                        #print(len(cpdag.arcs - bn.arcs())) 
                                        shd[t_index:t_index+Gap*I_count] = len(cpdag.arcs - bn.arcs())
                                      
                                      
                                        t_index = t_index + Gap*I_count 
                                  
                                        #print(t,len(cpdag.arcs - bn.arcs()),'asd',   I_count)
                                    
                              #print(t)
                              #print(len(cpdag.arcs - bn.arcs()))  
                              
                              shd[0:200] = len(bn.arcs())
                              #print(shd[0:t_index - Gap*I_count])
                              #print('***********')
                              Data_save.append(shd)
                              samples_list.append(t_index - Gap*I_count+1) 
                  
            return Data_save     