#These calculations will test if the punctured Hadamard code might produce vectors which exhibit state independent contextuality.
import numpy as np
from sage.all import *
from sage.combinat.matrices.hadamard_matrix import symmetric_conference_matrix
import itertools
epsilon = 0.001
def hadamard_matrix(n):
    #returns a hadamard matrix of size 2^(n-1) using the construction H H
    #                                                             H-H
    if n == 1:
        return np.array([[1]])
    H = hadamard_matrix(n-1)
    return np.block([[H, H], [H, -H]])

def puncture_hadamard(n,removal):
    #assumes removal is a list of the indices to be removed.
    H = hadamard_matrix(n).tolist()
    H = [ [element for index, element in enumerate(row) if (not index in removal)] for row in H]
    return H

def build_ortho_graph(vects): #puts edges between vectors when they are orthogonal.
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

def build_ortho_graph_from_gram_matrix(gram_matrix):
    #Takes a gram matrix, adds an edge whenever there is a 0 between entries i and j.
    #Make sure the gram matrix is psd!
    s = np.shape(gram_matrix)[0]#in the case of cab_vects, this is 18.
    labels = [str(i) for i in range(s)]
    adjacency_list = {}
    for v in range(s):
        current_list = []
        for b in range(s):
            if abs(gram_matrix[v][b]) < epsilon:
                current_list.append(labels[b])
        adjacency_list[labels[v]]=current_list
    return adjacency_list

def check_conference_matrix():
    removal = [16,17,18,19,20,21,22,23,24,25]
    C=symmetric_conference_matrix(26)

    punctured = [ [element for index, element in enumerate(row) if (not index in removal)] for row in C]
    #print C
    #print punctured
    ortho_list = build_ortho_graph(zip(punctured,map(lambda x: str(x), punctured)))
    G = Graph(ortho_list)
    G.allow_multiple_edges(False)
    g = G.plot()
    #print "lovasz_theta"
    #print G.lovasz_theta()
    #print "independence number"
    #print len(G.independent_set())
    connected_vertices = max(G.connected_components(),key=lambda x:len(x))
    #print len(connected_vertices)
    S = G.subgraph(connected_vertices)
    s = S.plot()
    print connected_vertices
    print "S Lovasz Theta"
    print S.lovasz_theta()
    print len(S.vertices())
    print len(S.independent_set())

    for v in connected_vertices:
        print v
    save(g,'/tmp/dom.png',axes=False,aspect_ratio=True)
    os.system('display /tmp/dom.png')
#print hadamard_matrix(4)
#print np.array(puncture_hadamard(5,[8,9,10,11,12,13,14,15]))

def check_conference_matrix2():
    removal = [16,17,18,19,20,21,22,23,24,25]
    C=symmetric_conference_matrix(26)
    punctured = [ [element for index, element in enumerate(row) if (index in removal)] for row in C]
    ortho_list = build_ortho_graph(zip(punctured,map(lambda x: str(x), punctured)))
    G = Graph(ortho_list)
    G.allow_multiple_edges(False)
    g = G.plot()
    connected_vertices = max(G.connected_components(),key=lambda x:len(x))
    S = G.subgraph(connected_vertices)
    s = S.plot()
    vectors = map(lambda x: np.array(x), punctured)
    handle = np.array(np.random.rand(10))
    handle = handle/np.linalg.norm(handle)
    print reduce(lambda x,y: x+y, [np.dot(handle,v)**2 for v in vectors],0)
    print G.lovasz_theta()
    print len(G.independent_set())
    print len(G.vertices())
    print np.linalg.norm(vectors[0])
#check_conference_matrix2()

def check_pauli_matrices():
    #Pauli Matrices
    #Wait a second... it doesn't make sense to look at pauli matrices.
    X=np.array([[0,1],[1,0]])
    Y=np.array([[0,-1j],[1j,0]])
    Z=np.array([[1,0],[0,-1]])
    I=np.identity(2)

    M = np.kron(Y,Y)
    print M
class string_counter:
    #produces a different string each time it's called.
    i = 0
    def produce_string(self):
        self.i+=1
        return str(self.i)

def check_fourier_transform(n):
    #See the unitary matrix on the wikipedia page https://en.wikipedia.org/wiki/Quantum_Fourier_transform
    N = 2**n
    o = math.e**(2*math.pi*1j/N)
    F= [[o**(row*col) for col in range(N)] for row in range(N)]
    removal = [j for j in range(N) if j%64==0]
    punctured = [ [element for index, element in enumerate(row) if (not index in removal)] for row in F]
    names=string_counter()
    ortho_list = build_ortho_graph(zip(punctured,[names.produce_string() for p in punctured]))
    G = Graph(ortho_list)
    G.allow_multiple_edges(False)
    print len(G.edges())
    print len(G.vertices())
    print G.is_perfect()
    #if not G.is_perfect():
    #    print j
    #g = G.plot()

    #save(g,'/tmp/dom.png',axes=False,aspect_ratio=True)
    #os.system('display /tmp/dom.png')
#check_fourier_transform(7)


cab_vects =   [[1,0,0,0],
              [0,1,0,0],
              [0,0,1,1],
              [0,0,1,-1],
              [1,-1,0,0],
              [1,1,-1,-1],
              [1,1,1,1],
              [1,-1,1,-1],
              [1,0,-1,0],
              [0,1,0,-1],
              [1,0,1,0],
              [1,1,-1,1],
              [-1,1,1,1],
              [1,1,1,-1],
              [1,0,0,1],
              [0,1,-1,0],
              [0,1,1,0],
              [0,0,0,1]]

def sic_set_from_gram_matrix(simple_cab_vects):
    cab_vects = map(lambda x: np.array(x)/np.linalg.norm(x), simple_cab_vects )
    gram_matrix = [ [np.dot(v1,v2) for v1 in cab_vects] for v2 in cab_vects]
    gram_matrix = np.array(gram_matrix)
    #cab_gram = np.linalg.matrix_rank(gram_matrix)
    dual_psd = 4.5*np.identity(18) - gram_matrix
    print gram_matrix
    print dual_psd
    print np.linalg.matrix_rank(dual_psd)


def use_cab_vects_as_eig_vects(simple_cab_vects):
    #This shows that the columns of the matrix V are the eigenvectors of its gram matrix.
    cab_vects = map(lambda x: np.array(x)/np.linalg.norm(x), simple_cab_vects )
    gram_matrix = [ [np.dot(v1,v2) for v1 in cab_vects] for v2 in cab_vects]
    cab_vects = np.array(cab_vects)
    turned_vectors = cab_vects.T
    new_matrix = np.zeros((18,18))
    for i in range(4):
        print turned_vectors[i]
        print np.outer(turned_vectors[i],turned_vectors[i])
        new_matrix += np.outer(turned_vectors[i],turned_vectors[i])
    print new_matrix
    graph = build_ortho_graph_from_gram_matrix(new_matrix)
    G = Graph(graph)
    g = G.plot()
    save(g,'/tmp/dom.png',axes=False,aspect_ratio=True)
    os.system('display /tmp/dom.png')

    graph2 = build_ortho_graph(zip(cab_vects,map(lambda x: str(x), cab_vects)))
    H = Graph(graph2)
    h = H.plot()
    save(h,'/tmp/dom2.png',axes=False,aspect_ratio=True)
    os.system('display /tmp/dom2.png')

    vector_gram_matrix = np.ravel(np.array([[np.dot(c1,c2) for c1 in cab_vects] for c2 in cab_vects]))
    vector_new_matrix = np.ravel(new_matrix)
    M = np.array([vector_gram_matrix, vector_new_matrix ])
    print M
    print np.linalg.matrix_rank(M)

    print G.is_isomorphic(H)

def check_knesser_for_sic(n,k):
    vectors = []
    new_vector = 0
    for r in itertools.combinations(range(n),k):
        new_vector = [1 if b in r else 0 for b in range(n)]
        vectors.append(new_vector)
    print vectors
    vectors = map(lambda x: np.array(x)/np.linalg.norm(x), vectors)
    vectors = np.array(vectors)
    gram_matrix = [[np.dot(v1,v2) for v1 in vectors.T] for v2 in vectors.T]
    outer_matrix = reduce(lambda x,y: x+y, [np.outer(v1,v1) for v1 in vectors], np.zeros((len(vectors[0]),len(vectors[0]))) )
    print outer_matrix
    print np.linalg.eigh(outer_matrix)
#check_knesser_for_sic(7,3)

def eigenvalues_of_cabello_graph(simple_cab_vects):
    cab_vects = map(lambda x: np.array(x)/np.linalg.norm(x), simple_cab_vects )
    gram_matrix = [ [np.dot(v1,v2) for v1 in cab_vects] for v2 in cab_vects]
    ortho_graph = build_ortho_graph(zip(cab_vects,map(lambda x: str(x), cab_vects)))
    G = Graph(ortho_graph)
    print np.linalg.eigh(G.adjacency_matrix())
eigenvalues_of_cabello_graph(cab_vects)


#use_cab_vects_as_eig_vects(cab_vects)
#sic_set_from_gram_matrix(cab_vects)


# for s in powerset(range(2**(5-1))):
#     if len(s)%2==0:
#         p = np.array(puncture_hadamard(5,s)).tolist()
#         labels = map(lambda x: str(x), p)
#         ortho_list = build_ortho_graph(zip(p,labels))
#         G = Graph(ortho_list)
#         G.allow_multiple_edges(False)
#         if not G.is_perfect():
#             print "perfection found!"
#             print s
# print "found nothing"
#print "is it perfect?"
#print G.is_perfect()
#
# g = G.plot()
# save(g,'/tmp/dom.png',axes=False,aspect_ratio=True)
# os.system('display /tmp/dom.png')
