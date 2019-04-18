import numpy as np
from numpy import linalg as LA
import itertools
import profile
#import xlwt
#import bronkerbosch
#import re
from sage.all import *
#from sage.graphs.connectivity import connected_components


def construct_margulis_graph(n):
    vertices = itertools.product(range(n),range(n))

    dict = {}
    for v in vertices:
        n1 = ((v[0] + 2*v[1]) % n, v[1])
        n2 = (v[0], (2*v[0]+v[1]) % n)
        n3 = ((1+v[0] + 2*v[1]) % n, v[1])
        n4 = (v[0], (1+2*v[0]+v[1]) % n)
        neighbors = [n1,n2,n3,n4]
        print neighbors
        dict[v] = neighbors

    G = DiGraph(dict)
    G = Graph(G, multiedges=True)
    return G


def construct_chorded_cycles(n):
    vertices = range(n)
    dict ={}
    for v in vertices:
        n1 = (v - 1) %n
        n2 = (v + 1) %n
        if v != 0:
            n3 = inverse_mod(v,n)
        else:
            n3 = 0
        neighbors = [n1,n2,n3]
        dict[v] = neighbors
    G = Graph(dict)
    G.allow_multiple_edges(False)
    G.allow_loops(False)
    return G

def construct_other_cycle(n):
    #constructs the cycle defined by connecting k to k+n/2 mod n
    vertices = range(n)
    dict = {}
    for v in vertices:
        n1 = (v - 1) %n
        n2 = (v + 1) %n
        n3 = (v+ ceil(n/2.0) ) % n
        neighbors = [n1,n2,n3]
        dict[v] = neighbors
    G = Graph(dict)
    G.allow_multiple_edges(False)
    G.allow_loops(False)
    return G

def calculate_separation(G):
    #given a graph G, compute the ratio of a/t.
    a = len(G.independent_set())
    t = G.lovasz_theta()
    return a/t

def print_table_of_separation_values(max):
    best_sep = 1
    best_n = 0
    data = {}
    print "n  spectralGap"
    for n in range(5,max):
        if n in Primes():
            #G = construct_chorded_cycles(n)
            G = construct_other_cycle(n)
            Gc = G.complement()
            # LoG = G.line_graph()
            # LoGc = Gc.line_graph()
            sepG = calculate_separation(G)
            sepGc = calculate_separation(Gc)
            #sepLoG = calculate_separation(LoG)
            #sepLoGc = calculate_separation(LoGc)
            #print n, [sepG, sepGc, sepLoG, sepLoGc]
            print n, [sepG, sepGc]
            #eigenvalues = np.linalg.eigh(G.adjacency_matrix())[0]
            #print n, "  ",eigenvalues[-1] - eigenvalues[-2]
#print_table_of_separation_values(150)

def rotation_map_from_dict(dict):
    return(lambda x,y: [dict[x][y], dict[dict[x][y]].index(x)])

def dict_from_rotation(r, vertex_iterator, degree_iterator):
    dict = {}
    print r((0,1),(0,0))
    deg = list(degree_iterator)
    for v in vertex_iterator:
        print v
        neighbors=[]
        for d in deg:
            neighbors.append(r(v,d)[0])
        print neighbors
        dict[v] = neighbors
    return dict

class int_counter:
    #produces a integer each time it's called.
    i = 0
    def tick(self):
        self.i+=1
        return self.i -1
def canonize_dict(dict):
    #returns a new dictionary where each vertex is replaced with an integer.
    new_dict = {}
    replacement_mapping = {}
    counter = int_counter()
    def replace(vertex, replacement_mapping, counter):
        if replacement_mapping.has_key(vertex):
            return replacement_mapping[vertex]
        else:
            next_number = counter.tick()
            replacement_mapping[vertex] = next_number
            return next_number

    for key in dict.keys():
        #print key
        new_key = replace(key, replacement_mapping, counter)
        #print counter.tick()
        #print new_key
        new_neighbors = [replace(k, replacement_mapping,counter) for k in dict[key] ]
        new_dict[new_key] = new_neighbors
    return new_dict
#
def zig_zag(G,H):
    #Takes two rotation maps, returns the zig-zag rotation map.
    def to_return(vertex,edge):
        v,a = vertex
        i,j = edge

        ap, ip = H(a,i)
        w, bp = G(v,ap)
        b, jp = H(bp, j)
        return [(w,b), (jp,ip)]
    return to_return
    #assumes that G,H
#
# G = {0:[1,2,3], 1:[0,2,3], 2:[0,1,3], 3:[0,1,2]}
# g = rotation_map_from_dict(G)
# H = {0: [1,2], 1: [2,0], 2: [0,1]}
# h = rotation_map_from_dict(H)
# ziggy = zig_zag(g,h)
# print ziggy((0,1), (0,1))
# print dict_from_rotation(ziggy, itertools.product(range(4),range(3)), itertools.product(range(2), range(2)))
# print canonize_dict(dict_from_rotation(ziggy, itertools.product(range(4),range(3)), itertools.product(range(2), range(2))))


def build_ortho_graph(vects): #puts edges between vectors when they are orthogonal.
    epsilon = 10**-4
    vectors = [v[0] for v in vects] #just the vectors
    labels = [v[1] for v in vects] #just the labels
    adjacency_list = {}
    for v in range(len(vectors)):
        current_list = []
        for b in range(len(vectors)):
            if v!=b and abs( np.inner(vectors[v],vectors[b]) )<epsilon:
                current_list.append(labels[b])
        adjacency_list[labels[v]]=current_list
    return adjacency_list

def adjacency_list_from_matrix(M):
    dict = {}
    for n in range(len(M)):
        dict[n] = [m for m in range(len(M)) if M[m][n]==1]
    return dict
def zig_zag_product_example():
    #Tries to compute some constants for P z G,
    #But indep. set of G and L. theta of Gc time out.
    G = graphs.BidiakisCube()
    P = Graph([GF(25), lambda i,j: i!=j and (i-j).is_square()])
    g = adjacency_list_from_matrix(np.array(G.adjacency_matrix()).tolist())
    p = adjacency_list_from_matrix(np.array(P.adjacency_matrix()).tolist())
    g_rot = rotation_map_from_dict(g)
    p_rot = rotation_map_from_dict(p)
    print p
    print g
    p_zag_g = zig_zag(p_rot,g_rot)

    G = canonize_dict(dict_from_rotation(p_zag_g, itertools.product(range(25),range(12)), itertools.product(range(3),range(3))))
    G = Graph(G)
    # print "G indep. set"
    # print len(G.independent_set())
    # print "G L. theta"
    # print G.lovasz_theta()
    Gc = G.complement()
    print "Gc indep. set"
    print len(Gc.independent_set())
    print "Gc L. theta"
    print Gc.lovasz_theta()

    print "spectral gap"
    print np.linalg.eigh(G.adjacency_matrix())[0]

    g = G.plot()
    save(g,'/tmp/dom.png',axes=False,aspect_ratio=True)
    os.system('display /tmp/dom.png')

def sum_of_four_squares(n):
    #lists all ways that n can be written as a sum of 4 squares.
    if n ==0:
        return([[0,0,0,0]])
    if n ==1:
        return([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1],
               [-1,0,0,0],[0,-1,0,0],[0,0,-1,0],[0,0,0,-1]])
    if n%4==0:
        sub_solution = sum_of_four_squares(n/2)
        new_solution = [ [(x[0]- x[1])/2, (x[0]+x[1])/2, (x[2]-x[3])/2, (x[2]-x[3])/2] for x in sub_solution]
        return new_solution
    if n%4==2:
        sub_solution = sum_of_four_squares(n/2)
        new_solution = [ [x[0]- x[1], x[0]+x[1], x[2]-x[3], x[2]-x[3]] for x in sub_solution]


#zig_zag_product_example()
# cab_vects =   [[1,0,0,0],
#               [0,1,0,0],
#               [0,0,1,1],
#               [0,0,1,-1],
#               [1,-1,0,0],
#               [1,1,-1,-1],
#               [1,1,1,1],
#               [1,-1,1,-1],
#               [1,0,-1,0],
#               [0,1,0,-1],
#               [1,0,1,0],
#               [1,1,-1,1],
#               [-1,1,1,1],
#               [1,1,1,-1],
#               [1,0,0,1],
#               [0,1,-1,0],
#               [0,1,1,0],
#               [0,0,0,1]]
# cab_vects = map(lambda x: np.array(x)/np.linalg.norm(x), cab_vects)
# cannonized_cab_vects = canonize_dict(build_ortho_graph(zip(cab_vects, map(lambda x: str(x), cab_vects))))


#z = rotation_map_from_dict({0: [1,2], 1: [2,0], 2: [1,0]})
#print z(1,0)
#print dict_from_rotation(z, 3,2)

#print canonize_dict({'a': ['b','c'], 'b': ['c','a'], 'c': ['b','a']})
#print_table_of_separation_values(1000)
#
# G = construct_chorded_cycles(37)
#
# print G.is_regular()
# print G.is_directed()
# print "G indep. set"
# print len(G.independent_set())
# print "G L. theta"
# print G.lovasz_theta()
# Gc = G.complement()
# print "Gc indep. set"
# print len(Gc.independent_set())
# print "Gc L. theta"
# print Gc.lovasz_theta()
#
# g = G.plot()
# save(g,'/tmp/dom.png',axes=False,aspect_ratio=True)
# os.system('display /tmp/dom.png')

def add_signs(x):
    to_return = []
    if x == []:
        return [[]]
    if x[0] >0:
        to_return.extend( [[x[0]] + h for h in add_signs(x[1:])])
        to_return.extend( [[-x[0]] + h for h in add_signs(x[1:])])
    else:
        to_return.extend( [[x[0]] + h for h in add_signs(x[1:])])
    return to_return
def sum_of_four_squares(n):
    list = []
    top = int(math.ceil(n**0.5) + 1)
    for i in range(top):
        for j in range(top):
            for k in range(top):
                for h in range(top):
                    if i**2 + j**2 +k**2 +h**2 == n:
                        new_elements = add_signs([i,j,k,h])
                        list.extend(new_elements)
    return list

def sum_divisors(n):
    #sums the divisors of n
    current_sum = 0
    for i in range(1,n+1):
        if n%i==0:
            current_sum += i
    return current_sum
def test_sum_of_four_squares():
    for i in range(1,30):
        if i%2==1:
            sos = sum_of_four_squares(i)
            if not len(sos) == 8*sum_divisors(i):
                print "failed sum of squares test"
                print i
                print sos
                print len(sos)
                print sum_divisors(i)
                return
            for s in sos:
                if s[0]**2+s[1]**2+s[2]**2+s[3]**2 !=i:
                    print "failed sos test"
                    return
    print "passed test"
    return

#test_sum_of_four_squares()
#print sum_of_four_squares(26)
class glement():#implements the general linear group on Z
    def __init__(self,q,matrix):
        self.q = q
        self.matrix=[[matrix[0][0]%q, matrix[0][1]%q], [matrix[1][0]%q, matrix[1][1]%q]] #make sure the entries are modulo q.
    def determinant(self):
        return (self.matrix[0][0]*self.matrix[1][1] - self.matrix[1][0]*self.matrix[0][1])%self.q
    def inv(self):
        a = self.matrix[0][0]
        b = self.matrix[0][1]
        c = self.matrix[1][0]
        d = self.matrix[1][1]
        q = self.q
        det = (a*d-b*c)%q
        i = inverse_mod(det,q)
        return pslement(q,[[(i*d)%q, (i*-b)%q],[(i*-c)%q,(i*a)%q]])
    def string_name(self):
        return str(self.matrix)
    def __mul__(self, another_matrix):
        if self.q != another_matrix.q:
            print "can't multiply these matricies!"
        a = self.matrix
        b = another_matrix.matrix
        q = self.q
        e1= (a[0][0]*b[0][0]+a[0][1]*b[1][0]) %q
        e2= (a[0][0]*b[0][1]+a[0][1]*b[1][1]) %q
        e3= (a[1][0]*b[0][0]+a[1][1]*b[1][0]) %q
        e4= (a[1][0]*b[0][1]+a[1][1]*b[1][1]) %q
        new_matrix = [[ e1,e2 ],[ e3,e4 ]]
        return pslement(self.q, new_matrix)


class pslement(glement):
        #pselement is exactly the same as glement, except that we identify multiplication by [[-1,0][0,-1]]
        #We also need to check that the determinant is 1. This is currently a weakness in this code.
    def __init__(self, q, matrix):
        glement.__init__(self, q, matrix)
    def __eq__(self, another_matrix):
        if self.q != another_matrix.q:
            return False
        a = self.matrix
        b = another_matrix.matrix
        q = self.q
        if (a[0][0] - b[0][0]) %q ==0 and (a[0][1]-b[0][1])%q==0 and (a[1][0]-b[1][0])%q==0 and (a[1][1]-b[1][1])%q==0:
            return True
        elif (a[0][0] + b[0][0]) %q ==0 and (a[0][1]+b[0][1])%q==0 and (a[1][0]+b[1][0])%q==0 and (a[1][1]+b[1][1])%q==0:
            return True
        return False

class pglement(glement):
    def __init(self, q, matrix):
        glement.__init__(self, q, matrix)
    def __eq__(self, another_matrix):
        if self.q != another_matrix.q:
            return False
        a = self.matrix
        b = another_matrix.matrix
        q = self.q
        for r in range(q):
            if (a[0][0] - r*b[0][0]) %q ==0 and (a[0][1] - r*b[0][1])%q==0 and (a[1][0] - r*b[1][0])%q==0 and (a[1][1] - r*b[1][1])%q==0:
                return True
        return False

def list_psl(q):
    #returns a list of all the elements of psl.
    list = []
    for a,b,c,d in itertools.product(range(q), repeat = 4):
        mat = pslement(q,[[a,b],[c,d]])
        if mat.determinant()%q==1 and not mat in list:
            list.append(mat)
    return list

def list_gsl(q):
    list = []
    for a,b,c,d in itertools.product(range(q), repeat = 4):
        mat = pglement(q,[[a,b],[c,d]])
        if mat.determinant()!=0 and not mat in list:
            list.append(mat)
    return list

def list_special_generators(p,q):
    #These come from the number of ways to write p as the sum of 4 squares.
    answers = sum_of_four_squares(p)
    i=0
    S=[]
    for x in range(q):
        if x**2 % q == q-1:
            i=x
            break
    if kronecker(p,q)==1:
        for c in range(q):
            if (c**2 -p) %q ==0:
                normalizer = inverse_mod(c,q)
                break
        for a in answers:
            if a[0]>0 and a[0]%2==1 and a[1]%2==0 and a[2]%2==0 and a[3]%2==0:
                el = pslement(q,[ [ normalizer*(a[0]+i*a[1]), normalizer*(a[2]+i*a[3])], [ normalizer*(-a[2]+i*a[3]), normalizer*(a[0]-i*a[1])]] )
                S.append(el)
        return S
    elif kronecker(p,q)==-1:
        for a in answers:
            if a[0]>0 and a[0]%2==1 and a[1]%2==0 and a[2]%2==0 and a[3]%2==0:
                el = pslement(q,[ [ (a[0]+i*a[1]), (a[2]+i*a[3])], [ (-a[2]+i*a[3]), (a[0]-i*a[1])]] )
                S.append(el)
        return S


def createRamanujanGraph(p,q):
    #assume p and q are primes congruent to 1 mod 4
    if kronecker(p,q)==1:
        total_nodes = list_psl(q)
    elif kronecker(p,q)==-1:
        total_nodes = list_gsl(q)
    S = list_special_generators(p,q) #this set defines the edges in the graph
    adjacency_list = {}
    for node1 in total_nodes:
        neighbors = []
        standard_neighbors = [] #chooses a standard representative element.
        for s in S:
            neighbors.append(node1*s)
        for node2 in total_nodes:
            if node2 in neighbors:
                standard_neighbors.append(node2.string_name())
        adjacency_list[node1.string_name()]= standard_neighbors
    return adjacency_list

def test_list_psl(q):
    list = list_psl(q)
    for q1, q2 in itertools.combinations(list,2):
        if q1 == q2:
            print "messed up"
    print "test complete"

class depth_first_search:
    def __init__(self,adj_list, vertex):
        self.adj_list = adj_list
        self.vertex = vertex
        self.visited = []
        self.partially_visited=[] #visited, but not each of its children.
    def search(self):
        while True:
            cur_vertex = self.vertex
            self.partially_visited.append(self.vertex)
            self.visited.append(self.vertex)
            while self.vertex in self.partially_visited:
                try: self.vertex = (set(self.adj_list[self.vertex])-set(self.visited)).pop()
                except KeyError:
                    self.partially_visited.remove(self.vertex)
                    if len(self.partially_visited)==0:
                        raise StopIteration
                    else:
                        self.vertex = self.partially_visited[0]

            yield self.vertex
def induced_adjacency_list(adj_list, vertices):
    #Calculates the adjacency matrix for the induced graph on the vertices.
    vertices = list(vertices)
    new_list = {}
    for v in vertices:
        new_list[v]=list(set(adj_list[v]).intersection(set(vertices)))
    return new_list

# def profile():
#     %prun _=createRamanujanGraph(5,13)
# profile()

#print len(list_special_generators(13,17))
# l=list_psl(17)
# for l in l:
#     if l.matrix is None:
#         print "aweful"
r= createRamanujanGraph(5,13)
r = canonize_dict(r)
#print "r[0] is"
#print r[0]
#r = {0:[1,2], 1:[0,3,2,4], 2:[0,1,3], 3:[1,2,4], 4:[5,3], 5:[4,1], 6:[7], 7:[6]}

#d = depth_first_search(r,0)
#i = induced_adjacency_list(r,d.search())
G = Graph(r)
print i
if G.is_connected():
     print "connected!"
 if G.is_regular():
     print "regular "
# for vertex in d.search():
#     print vertex
# print r
#G = Graph(r)
# g = G.plot()
#print G.is_regular()
#print G.is_connected()
# for v in G.vertices():
#     print v
#connected_component = G.breadth_first_search(G.vertices()[0])
#print G.is_directed()
#rp = {}
# for v in connected_component:
#     print v
#     rp[v]=[x for x in r[v] if x in connected_component]
# print rp
#print r
#subgraph = G.subgraph(connected_component)
#adj = subgraph.adjacency_list()
#print adj
# A = G.adjacency_matrix()
# print G.is_regular()
# eigen =  np.linalg.eigh(A)[0]
# print eigen[0]
# print eigen[-2]
# print eigen[-1]
#
# Gc = subgraph.complement()
# print "independent set"
# print len(Gc.independent_set())
# print "l theta"
# print Gc.lovasz_theta()
# g = G.plot()
# save(g,'/tmp/dom.png',axes=False,aspect_ratio=True)
# os.system('display /tmp/dom.png')
# for i in range(6):
#      print list_special_generators(13,17)[i].matrix
# print "S"
# #createRamanujanGraph(5,13
# item = list_psl(13)[6]
# S = list_special_generators(5,13)
# s = S[4]
# print s.matrix
# print s.inv().matrix
# print (s*(s.inv())).matrix
# if s.inv() in S:
#     print "yipee"
# print s.matrix
# print (item*s).matrix
# res = pslement(13,[[11, 1], [0, 4]])
# print res.matrix
#print item*(res.inv()).matrix
# print list_psl(13)[5].matrix
# print list_psl(13)[5].inv().matrix
# print (list_psl(13)[5]*list_psl(13)[5].inv()).matrix
#test_list_psl(13)
#print map(lambda x: x.matrix,list_psl(5))
