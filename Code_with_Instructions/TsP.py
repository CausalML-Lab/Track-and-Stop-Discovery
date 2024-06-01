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
    


def block_sample_intervention(bn,intv,T,Config):
    
    bn1  = gum.BayesNet(bn)
    import numpy as np
    for j in intv:
        for i in bn1.parents(j):
            bn1.eraseArc(gum.Arc(i,j))
       
        if (Config[intv.index(j)] == 1):
                                        bn1.cpt(j)[:] = [0 ,1]

        else:    
                                        bn1.cpt(j)[:] = [1 ,0]
            
    df = gum.generateSample(bn1, n=100000, name_out=None, show_progress=False, with_labels=True, random_order=False)
    aa = df[0]
    bb =np.array(aa[T],dtype=int)
    return bb  



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

def find_trees(mpdag):

    
    def dfs(v, visited,tree):
        visited[v] = True
        for neighbor in graph[v]:
            if not visited[neighbor]:
                tree.append(neighbor)
                dfs(neighbor, visited,tree)
        return tree,visited
    graph = []
    for i in range(mpdag.size()):
      graph.append([i])
    for e in mpdag.edges():
      graph[e[0]].append(e[1])
      graph[e[1]].append(e[0])
        
    n = len(graph)  # Number of vertices in the graph
    visited = [False] * n
    trees = []

    for v in range(n):
        if not visited[v]:
            tree = [v]
            tree, visited =  dfs(v, visited,tree)
            trees.append(tree)

    return trees



def find_trees_d(mpdag):

    
    def dfs(v, visited,tree):
        visited[v] = True
        for neighbor in graph[v]:
            if not visited[neighbor]:
                tree.append(neighbor)
                dfs(neighbor, visited,tree)
        return tree,visited
    graph = []
    for i in range(mpdag.size()):
      graph.append([i])
    for e in mpdag.arcs():
      graph[e[0]].append(e[1])
              
    n = len(graph)  # Number of vertices in the graph
    visited = [False] * n
    trees = []

    for v in range(n):
        if not visited[v]:
            tree = [v]
            tree, visited =  dfs(v, visited,tree)
            trees.append(tree)

    return trees


def PCO(mpdag,D):
    CC= find_trees(mpdag)
    
    PCO =[];
    from itertools import chain
    def flatten_chain(matrix):
       return list(chain.from_iterable(matrix))

    def intersection(lst1, lst2):
        return list(set(lst1) & set(lst2))
    i = -1;
    while len(CC)>0:
       i = (i+1)% len(CC)
       C = CC[i]
       #print('*')
       #print(CC)
       #print(C)
       Cbar = list(CC)
       Cbar.remove(C)
       Cbar = flatten_chain(Cbar)
       flag = True
       for e in mpdag.arcs():
          if ((e[0] in C) & (e[1] in Cbar)):
             flag = False
       for e in mpdag.edges():
          if ( ((e[0] in C) | (e[1] in C)) &  ((e[0] in Cbar) | (e[1] in Cbar))):
             flag = False         
       if(flag):
          #print('yes')
          CC.remove(C)
          tmp = intersection(C,D)  
          if (len(tmp) > 0):
            PCO.append(tmp)   

    PCO.reverse()
    return PCO
    
    
def IdentifyCausaLEffectinMPDAG(mpdag,y,x,Pv):
    D = mpdag.ancestors(y,Anc=[],X=x)
    B_d = PCO(mpdag,D)
    pygdo_x = 1
    for v in B_d:
         Pav = mpdag.Parents(v,Pa=[])
         tmp1 = mpdag.nodes - ( set(v) | set(Pav))
         tmp1 = list(map(str, tmp1))
         tmp2 = mpdag.nodes - (set(Pav))
         tmp2= list(map(str, tmp2))

         if (len(tmp1) > 0):    
               Ptmp = Pv.margSumOut(tmp1)/Pv.margSumOut(tmp2)
         else:
               Ptmp = Pv/Pv.margSumOut(tmp2)
         
         pygdo_x = pygdo_x * Ptmp  

    tmp =  set(D) - set(y)
    tmp = list(map(str, tmp))
    if (len(tmp)>0):
       pygdo_x =pygdo_x.margSumOut(tmp)
    a= list(pygdo_x.names)
    a.reverse()
    #print(a)
    return a,pygdo_x 


def enumerate_causaleffects(mpdag_arr,V,I,Pv,Config):
  def flatten_list_recursive(lst):
    flattened_list = []
    for item in lst:
        if isinstance(item, list):
            flattened_list.extend(flatten_list_recursive(item))
        else:
            flattened_list.append(item)
    return flattened_list
    
    
  le =[]
  for mm in mpdag_arr:
     #print(mpdag_arr.index(mm))   
     (a,pygdo_x) = IdentifyCausaLEffectinMPDAG(mm,V,I,Pv)
     Cause = []
     Int = []
     new_config=[]
     for var in a:
        if int(var) in I:
                     Int.append(var)
                     new_config.append(Config[I.index(int(var))]) 
        else:
                     Cause.append(int(var))   
     sz = len(Int)
     Cause= sorted(Cause)
     Cs =[]   
     for c in Cause:
        Cs.append(str(c))
     new_a = Int + Cs
     new_a.reverse()
     pygdo_x = pygdo_x.reorganize(new_a)
     bb = pygdo_x.toarray()
     bb = np.round(bb,decimals = 7)   
     bb = bb.tolist()
 
     j = 0
     for i in range(sz):
       bb = bb[new_config[j]]
       j = j+1
    
     le.append(flatten_list_recursive(bb))
  return le,Cs 


def OrientCut_and_Enumeratmpdags(mpdag,X):
    import numpy as np
    import copy
    def decimal_to_binary_vector(decimal_number, vector_length):
        binary_string = bin(decimal_number)[2:]  # Remove the '0b' prefix from the binary string
        binary_string = binary_string.zfill(vector_length)  # Pad with leading zeros to the desired length
        binary_vector = [int(bit) for bit in binary_string]  # Convert each character to an integer and store in a list
        return binary_vector
    MPDAG_LIST = []
    Edge_list =[]
    Cut_X = []
    M =[]
    for e in mpdag.edges():
      if (e[0] in X and e[1] not in X)  or (e[0] not in X and e[1] in X) :
             Cut_X.append(e)
    if (len(Cut_X) > 0):            
        for i in range(2**(len(Cut_X))):    
           M.append(decimal_to_binary_vector(i,len(Cut_X)))


        M = np.array(M)
        sz = M.shape
        r = sz[0]
        c = sz[1]

        for rows in range(r):
            #print(rows)
            mp = copy.deepcopy(mpdag)
            for col in range(c):
                e = Cut_X[col]     
                if (M[rows,col] == 0):
                                      new_arc_h = e[0]
                                      new_arc_t = e[1]
                else:

                                      new_arc_h = e[1]
                                      new_arc_t = e[0]
                tmp =(new_arc_h,new_arc_t)
                mp.addarc(new_arc_h,new_arc_t)
            #print(mp.edges())
           #print(mp.arcs())
            if (mp.arcs() not in Edge_list):
                                            MPDAG_LIST.append(mp)
                                            Edge_list.append(mp.arcs())
            #print('*************')
    else:
           mp = copy.deepcopy(mpdag)
           MPDAG_LIST.append(mp)     
    return MPDAG_LIST,Cut_X


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
    
    
    
def Track_and_stop(Num_Grphs,nodes,degree,Max_samples):
    Data_save = []
    def euclidean_distance(row1, row2):
        return np.sqrt(np.sum((row1 - row2) ** 2))
    
    def min_row_distance(matrix):
        min_distance = float('inf')
    
        for i in range(len(matrix)):
            for j in range(i + 1, len(matrix)):
                distance = euclidean_distance(matrix[i], matrix[j])
                min_distance = min(min_distance, distance)
    
        return min_distance
    
    def int_to_binary_list(number, bits):
            binary_format = "{:0" + str(bits) + "b}"
            binary_representation = binary_format.format(number)
            binary_list = [int(bit) for bit in binary_representation]
            return binary_list

    import pyAgrum as gum
    import pyAgrum.causal as c
    
    import numpy as np
    for grphs_no in range(0,Num_Grphs, 5):
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
        
            mpdg = MPDAG(bn)
        
        
            mpdg.Edges = bn.arcs()
            mpdg.Arcs =  set({})
        
            Cuts = []    
            MPDAG_LIST = []
             
            I = indices_of_elements(r)
            for i in I:
                tmp,Cts = OrientCut_and_Enumeratmpdags(mpdg,i)
                MPDAG_LIST.append(tmp)
                Cuts.append(Cts)
                
            for iii in range(5):
                if (grphs_no+iii >= Num_Grphs):
                                               break  
                print('Garph No.',grphs_no+iii)
                Truedag = bn.arcs();
                
              
                bn.generateCPTs()
                Pv =1
                for n in bn.nodes():
                   Pv = Pv * bn.cpt(n)
                import copy
                from scipy.optimize import linprog
                import numpy as np
                import math
                from scipy.special import rel_entr
                mpdag = copy.deepcopy(mpdg)
             
                VARS =[]
                PD =[]
                Ps =[]
                Cnts =[]
                data = []
                
            
                delta = 0.00001
                
                golden_intv = np.zeros(len(I))
                for i in I:
                    V = list(mpdg.nodes - set(i))
                    p_search_min = np.zeros(2**len(i))
                    for v in range(2**len(i)):
                       config = int_to_binary_list(v,len(i)) 
                       tmp,tmp1 = enumerate_causaleffects(MPDAG_LIST[I.index(i)],V,i,Pv,config)
                       tmp = np.array(tmp)
                       p_search_min[v] =  min_row_distance(tmp)   
                    golden_intv[I.index(i)] = np.argmax(p_search_min)          
                    
                for i in I:
                   I_index = I.index(i) 
                   V = list(mpdg.nodes - set(i))
                   Ps.append([0]*(2**len(V)))
                   #config = np.random.randint(0, 2**((len(i)))-1, size=1, dtype=int)[0]
                   config = int_to_binary_list(int(golden_intv[I.index(i)]),len(i)) 
                   tmp,tmp1 = enumerate_causaleffects(MPDAG_LIST[I_index],V,i,Pv,config)
                   PD.append(tmp)
                   Cnts.append([0]*(2**len(V)))
                   VARS.append(tmp1)
                   data.append(block_sample_intervention(bn,i,tmp1,config))
            
                sz_i = len(I)
                sz_d = len(MPDAG_LIST)
            
                Zt =0
                Nt = np.zeros(sz_i,dtype= int)
                sZt_vec = np.zeros(sz_d,dtype= float)
                aub = np.zeros((sz_d,sz_i+1),dtype= float)
                aub[:,0] = -1*np.ones(sz_d)
                aeq = np.ones((1,sz_i+1))
                aeq[0][0] =0
                bub = np.zeros(sz_d)
                optim_res = np.zeros(sz_d,dtype= float)
                bnds = [(0,None)]*(sz_i+1)
                tmp = np.zeros(sz_i+1)
                tmp[0] = 1
                d_star = np.zeros(sz_i)
                alpha_star = np.zeros(sz_i)
                t = 1
                
                shd = np.zeros(Max_samples)
                for i in  I:
            
                       Sample_and_update_dist(data[I.index(i)],I.index(i),Cnts,Nt[I.index(i)])
                       Nt[I.index(i)]+=1
                       Ps[I.index(i)] = (np.array(Cnts[I.index(i)])/Nt[I.index(i)]).tolist() 
                samples = 0
                while (True):
                  if ( samples > Max_samples):
                             break        
                  I_index = 0
                  Dstar_E = set({})
                  for i in I:
                               
                               KL_vector = np.zeros(len(MPDAG_LIST[I_index]))
                               dindex = 0 
                               for m in  MPDAG_LIST[I_index]:
                                  
                            
                            
                                   KL_vector[dindex] = sum( rel_entr(Ps[I_index],PD[I_index][dindex], out=None) )/math.log(2) 
                                   if math.isinf( KL_vector[dindex]):
                                                                        KL_vector[dindex]  = 10^5    
                               
                                   dindex = dindex+1 
                               alpha_star[I_index] = 1/KL_vector.min() 
                               d_star[I_index] = np.argmin(KL_vector)
                               ntmp = MPDAG_LIST[I_index][int(d_star[I_index])].arcs()  
                               ntmp1 = ntmp.copy()
                               #for e in ntmp1:
                                    #if( (e[0],e[1]) not in Cuts[I_index] and (e[1],e[0]) not in Cuts[I_index] ):
                                                                                                                  #ntmp.remove(e) 
                               Dstar_E =  Dstar_E | ntmp
                               I_index =   I_index + 1
                              
                  
                  
                              
                                   
                  alpha_star = alpha_star/(sum(alpha_star)) 
                  if (min(Nt) < 25*math.sqrt(t)):                                                         
                                                  act = np.argmin(Nt)                           
                  else:
                                                  act = np.argmax(t*alpha_star - Nt) 
                  
            
                  Sample_and_update_dist(data[act],act,Cnts,Nt[act])
                  Nt[act]+=1
                  Ps[act] = (np.array(Cnts[act])/Nt[act]).tolist()      
                  shd[samples:samples+len(I[act])] = len(Truedag - Dstar_E) +len(Dstar_E-Truedag)
                  #print(len(Truedag - Dstar_E) +len(Dstar_E-Truedag))  
                  t = t + 1
                  #print(Zt,' ',math.log(t/delta),'  ',t,'  ',shd[t-2])
                  if (len(Truedag - Dstar_E) +len(Dstar_E-Truedag) == 0 and t >= 1000 and len(Truedag) == len(Dstar_E) ):
                                                     break
                  samples = samples + len(I[act])
                  #print(len(Dstar_E-Truedag)+len(Truedag-Dstar_E),t,len(Dstar_E),len(Truedag))
                #print('Gold ',golden_intv)
                #print(samples)
                #print(len(Dstar_E-Truedag)+len(Truedag-Dstar_E))
                #print(len(Dstar_E))
                #print("***")
                #print(Dstar_E) 
             
                
                Data_save.append(shd.tolist())
                #print(iii)

    return Data_save
            