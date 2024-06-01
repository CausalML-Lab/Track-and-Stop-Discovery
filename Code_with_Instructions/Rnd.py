import pyAgrum as gum
class MPDAG:
  def __init__(self, mpdag):
           
           self.nodes = mpdag.nodes()
           self.Arcs =  set([])
           self.Edges =  mpdag.arcs()
           self.apply_meek_rules()
        
  def removeedge(self,s,t):
      self.Edges.remove((s,t))
    
  def removearc(self,s,t):
      self.Arcs.remove((s,t))  
    
  def addedge(self,s,t):
      self.Edges.add((s,t))
    
  def addarc(self,s,t):
     
      if (s,t) in self.Edges:
             self.Edges.remove((s,t))
      elif(t,s) in self.Edges:
             self.Edges.remove((t,s))
      else:
              return 
        
      self.Arcs.add((s,t))
      self.apply_meek_rules()

  def edges(self):
      return self.Edges

  def arcs(self):
      return self.Arcs
 
  def size(self):
       return len(self.nodes)
  
  def ancestors(self,Y,Anc=[],X=[]):
       for v in Y:
            if (v not in Anc):
                Anc.append(v)
            for e in self.Arcs:    
                if ((e[1] == v) & (e[0] not in Anc)  & (e[0] not in X)):
                        Anc.append(e[0])
                        self.ancestors([e[0]],Anc,X)
                        
       return Anc 
  
  def undirected_neighbours(self,I):
      ln = [];
      for e in self.Edges: 
           if (e[0] in I and e[1] not in I and e[1] not in ln ):
                 ln.append(e[1])
           if (e[0] not in I and e[1]  in I and e[0] not in ln):
                 ln.append(e[0])  
      return ln      
            
  def Parents(self,Y,Pa=[]):
       for v in Y:
            
            for e in self.Arcs:    
                if ((e[1] == v) & (e[0] not in Pa) & (e not in Y) ):
                        Pa.append(e[0])
                        
                        
       return Pa
    
  def apply_meek_rules(self):
        import itertools

        def is_adjacent(a, b):
            return (a, b) in self.Edges or (b, a) in self.Edges or (a, b) in self.Arcs or (b, a) in self.Arcs

        def is_nonadjacent(a, b):
            return not is_adjacent(a, b)
        flag = False;
           
        # Rule 1: Orient b - c into b -> c whenever there is an arrow a -> b such that a and c are nonadjacent.
        for e in self.Arcs:
            a = e[0]
            b = e[1]
            for c in self.nodes:
                if c != b and is_nonadjacent(a, c)  and c!=a and  ((b, c) in self.Edges or (c, b) in self.Edges ) :
                    self.addarc(b, c)
                    flag = True;
                    break
            if (flag):
                   break;  
        # Rule 2: Orient a - b into a -> b whenever there is a chain a -> c -> b.
        for a, c, b in itertools.permutations(self.nodes, 3):
            if (a, c) in self.Arcs and (c, b) in self.Arcs and  ( (a, b) in self.Edges or (b, a) in self.Edges )  :
                self.addarc(a, b)
                break;
        # Rule 3: Orient a - b into a -> b whenever there are two chains a - k -> b and a - l -> b such that k and l are nonadjacent.
        for a, b, k, l in itertools.permutations(self.nodes, 4):
            if k != l and is_nonadjacent(k, l) and  ( (a, k) in self.Edges or (k, a) in self.Edges ) and ( (a, l) in self.Edges or (l, a) in self.Edges )  and (k, b) in self.Arcs and (l, b) in self.Arcs  and  ( (a, b) in self.Edges or (b, a) in self.Edges ) :
                self.addarc(a, b)
                break
        # Rule 4: Orient a - b into a -> b whenever there is an edge a - k and chain k -> l -> b such that k and b are nonadjacent.
        for a, k, l, b in itertools.permutations(self.nodes, 4):
            if k != b and is_nonadjacent(k, b) and (k, l) in self.Arcs and (l, b) in self.Arcs and ( (a, k) in self.Edges or (k, a) in self.Edges )  and    ( (a, b) in self.Edges or (b, a) in self.Edges ) :
                self.addarc(a, b)  
                break
    




def block_sample_intervention(bn,intv):
    import numpy as np

    bn1  = gum.BayesNet(bn)
    import numpy as np
    for j in intv:
        for i in bn1.parents(j):
            bn1.eraseArc(gum.Arc(i,j))
       
       
        bn1.cpt(j)[:] = [0.5 ,0.5]

           
                                        
            
    df = gum.generateSample(bn1, n=70000, name_out=None, show_progress=False, with_labels=True, random_order=False)
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
    
    
    
def Learn_cut(bn,X,arcs,data,index,Edges):
     import numpy as np
     data = np.array(data,dtype=int)
     #print(data.shape)
     Cut_X = []
     for e in Edges:
          if (e[0] in X and e[1] not in X):
                                        Cut_X.append((e[0],e[1]))
          if (e[0] not in X and e[1] in X):
                                       Cut_X.append((e[1],e[0]))
                    
     if (index >= 5):
          #print(j)
          for e in Cut_X: 
              #print(index)
              #print
              from causallearn.utils.cit import CIT
              chisq_obj = CIT(data[0:index,:], "chisq") # construct a CIT instance with data and method name
              pValue = chisq_obj(e[0], e[1])
              #tmp = mutual_info_score(data[0:index,e[0]],data[0:index,e[1]] )
              if (pValue >= 0.05):
                                narc = (e[1],e[0])  
              else:                     
                                narc = (e[0],e[1])  

              arcs.add(narc)
              if ((narc[1],narc[0]) in arcs):
                                               arcs.remove((narc[1],narc[0]))
     #print(arcs) 
    
     return arcs  


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





def Random_Interventions(Num_Grphs,nodes,degree,Max_samples,Gap):
        import numpy as np
        import math
        # Import the necessary modules
        import pyAgrum as gum
        import pyAgrum.causal as c
        
        Data_save=[]
        samples_list  = []
        for grphs_no in range(0,Num_Grphs, 1):
                print('Graph No.',grphs_no) 
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
                Truedag = bn.arcs();
        
                Edges = bn.arcs() 
    
                for j in range(1):
                    bn.generateCPTs()
                    arcs = set([])
                 
                    
                    I = indices_of_elements(r)
                    sz_i = len(I)
                    categories = range(sz_i)
                    probabilities = list(np.ones(sz_i)/sz_i)
                    avg_int_size  = 0
                    for i in I:
                        avg_int_size = avg_int_size + len(i)
                    avg_int_size = avg_int_size/sz_i
                    
                    shd = np.zeros(int(Max_samples))
                    sample_size = 1
                    index = list((np.zeros(sz_i)))
                    index_arr = np.zeros((200001,sz_i))
                    data =[]
                    for i in I:
                           data.append(block_sample_intervention(bn,i))
                            
                    for t in range(200001):
                       #print(t)
                       acttt= np.random.choice(categories, size=sample_size, p=probabilities)   
                       index[int(acttt[0])] =  index[int(acttt[0])] + 1
                       index_arr[t,:] = index     
                    samples = 0  
                    #print(I)
                    for t in range(Gap,Max_samples,Gap):
                        #print(t)
                         for i in I:  
                                
                            arcs =  Learn_cut(bn,i,arcs,data[I.index(i)],int(index_arr[t,I.index(i)]),Edges)
                         #print(len(Truedag -arcs) ,' safr', avg_int_size )
                      
                         shd[samples:samples + math.floor(avg_int_size*Gap)] = len(Truedag-arcs)
                         samples = samples +   math.floor(avg_int_size*Gap)
                         
                     
                  
                    shd[0:500] = len(bn.arcs())
                    Data_save.append(shd)
                    samples_list.append(samples-  math.floor(avg_int_size*Gap)) 
                    #print(samples_list)
            
        return Data_save