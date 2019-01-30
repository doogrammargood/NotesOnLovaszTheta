import numpy as np
from numpy import linalg as LA
import itertools
import bronkerbosch
import re
from sage.all import *

USE_LARGEST_CLIQUES = True;
epsilon = 0.001
#Pauli Matrices
X=np.array([[0,1],[1,0]])
Y=np.array([[0,-1j],[1j,0]])
Z=np.array([[1,0],[0,-1]])
I=np.identity(2)

def preprocessObservables(observables): #checks that the observables are valid. Returns their dimension.
    dimension = 0
    for O in observables:
        if not np.array_equal(O.conj().T,O): #check that the matrices are hermitian
            return -1
        if dimension == 0:
            dimension = O.shape[0]
        else:
            if O.shape[0] != dimension or O.shape[1] != dimension:
                return -1
    return dimension

def checkCommutation(A,B): #returns true if A and B commute. Assumes they are Hermitian!!!
    if np.array_equal(np.matmul(A,B), np.matmul(B,A)):
        return True
    else:
        return False

def build_graph(observables): #This puts edges between observables which commute. The contexts will be cliques in this graph.
    adjacency_list = {}
    for a in range(len(observables)):
        current_list = []
        for b in range(len(observables)):
            if a!=b and checkCommutation(observables[a], observables[b]):
                current_list.append(b)
        adjacency_list[a] = current_list
    return adjacency_list

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

def complete_3d_vectors(vects): #returns a new list of vectors. If these vectors are 3d, then we complete every pair of orthogonal vectors to a vertex.
    vectors = [v[0] for v in vects] #strip the labels
    hidden_vectors = [] #vectors added to complete cliques
    if len(vectors[0])!=3:
        print "Only use complete_3d_vectors when the vectors are 3d!" #otherwise, cross product is not defined.
    for v in range(len(vectors)):
        for b in range(v,len(vects)):
            if v!=b and abs( np.inner(vectors[v],vectors[b]) )<epsilon:
                hidden_vectors.append([np.cross(vectors[v],vectors[b]), 'u' + vects[v][1] +vects[b][1] + '*'])
    #labels = ['u' + str(i) for i in range(len(hidden_vectors))]
    #print map(list, zip(hidden_vectors, dummy_labels))
    #hidden_vectors = removeRedundantVectors(map(list,zip(hidden_vectors, dummy_labels)))
    total_vectors = removeRedundantVectors(hidden_vectors + vects)
    print total_vectors
    print len(vects)
    print len(total_vectors)
    return total_vectors


def simultaneouslyDiagonalize(observables): #Assumes observables have eigenvalues {-1, 1}.
    d=observables[0]
    for i in range(len(observables)):
        if i != 0:
            d = d + (2.**i)*observables[i]
    a, b = LA.eigh(d)
    return np.transpose(b)

def createLabels(vectors, dimension): #takes in a list of vectors, creates a list of (vector, label) tuples.
    labeled_vectors = []
    for counter in range(len(vectors)):
        label = "C-" + str(counter//dimension) + "-" + str(counter % dimension)
        labeled_vectors.append((vectors[counter], label))
    return labeled_vectors

def clean_yu_oh_labels(vectors): #cleans the labels for the completed yu-oh graph. Not for general use.
    pattern = re.compile(r"(\*,y)|(\*,h) |(\*,z)")
    for v in vectors:
        if pattern.match(v[1]):
            v[1] = v[1].split(',')[-1]

def removeRedundantVectors(vectors): #assumes vectors are tuples of (vector, label) pairs.
    unique_vectors = []
    for v in vectors: #loop through pairs of vectors
        append = True
        for u in unique_vectors:
            if (abs(np.inner(v[0],u[0]) - 1.0) < epsilon or abs(np.inner (v[0],u[0]) + 1.0) < epsilon): #check if v is in the same direction as u for any u.
                #print u
                #print v
                u[1] = u[1]+ "," + v[1] #add the label of v to the label of u
                append = False
        if append:
            unique_vectors.append([v[0],v[1]])
    return unique_vectors

def createSageGraph(graph): #expects an adjacency list as a dict.
    G = Graph(graph)
    G.allow_multiple_edges(False)

    g = G.plot()
    save(g,'/tmp/dom.png',axes=False,aspect_ratio=True)
    os.system('display /tmp/dom.png')
    print G.vertex_connectivity()
    print G.complement().vertex_connectivity()
    print G.automorphism_group().order()
    print G.is_vertex_transitive()
    print G.is_edge_transitive()
    print G.is_distance_regular()
    #print G.is_distance_transitive()
    print len(G.vertices())
    print len(G.edges())
    print G.lovasz_theta()
    print len(G.independent_set())


yu_oh_vectors = [[0,1,-1], [1,0,-1], [1,-1,0], [0,1,1], [1,0,1], [1,1,0], #ys
                [-1,1,1], [1,-1,1], [1,1,-1], [1,1,1], #hs
                 [1,0,0],[0,1,0],[0,0,1]] #zs #see https://arxiv.org/pdf/1109.4396.pdf

cab_vects =   [[1,0,0,0],
              [0,1,0,0,],
              [0,0,1,1,],
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
cab_cliques = [ [0, 1, 2, 3],
                [3, 4, 5, 6],
                [6, 7, 8, 9],
                [9, 10,11,12],
                [12,13,14,15],
                [15,16,17,0],
                [17,1,8,10],
                [2,4,11,13],
                [5,7,14,16]]
def build_cab_graph(cab_vects, cab_cliques): #builds the ortho_graph for the cab vectors. But there are extra orthogonalities which aren't edges.
    vectors = cab_vects
    labels = map(lambda x: str(x), cab_vects)
    adjacency_list = {}
    for v in range(len(vectors)):
        current_list = []
        for b in range(len(vectors)):
            if v!=b and any(v in q and b in q for q in cab_cliques):
                current_list.append(labels[b])
        adjacency_list[labels[v]]=current_list
    return adjacency_list

#cab_vects = zip(map(lambda x: np.array(x)/np.linalg.norm(x), cab_vects ), map(lambda x: str(x), cab_vects))



yu_oh_vectors = map(lambda x: x/np.linalg.norm(x), yu_oh_vectors)
labels = ['y1m', 'y2m', 'y3m', 'y1p', 'y2p', 'y3p', 'h1', 'h2', 'h3', 'h0', 'z1', 'z2', 'z3']
yu_oh_vectors = list(map(list, zip(yu_oh_vectors, labels)))
#completed_yu_oh_vectors=complete_3d_vectors(yu_oh_vectors)
#clean_yu_oh_labels(completed_yu_oh_vectors)

obs_list = [np.kron(I,Z), np.kron(Z,I), np.kron(Z,Z), np.kron(X,I), np.kron(I,X), np.kron(X,X), np.kron(X,Z), np.kron(Z,X), np.kron(Y,Y)] #Mermin Square
# obs_list = [np.kron(np.kron(I,I), X),#0
#             np.kron(np.kron(X,Z), X),#1
#             np.kron(np.kron(X,I), I),#2
#             np.kron(np.kron(I,Z), I),#3
#             np.kron(np.kron(Z,Z), Z),#4
#             np.kron(np.kron(Z,X), X),#5
#             np.kron(np.kron(X,X), Z),#6
#             np.kron(np.kron(Z,I), I),#7
#             np.kron(np.kron(I,I), Z),#8
#             np.kron(np.kron(I,X), I)] #Mermin Star
ortho_graph=build_cab_graph(cab_vects,cab_cliques)
#ortho_graph = build_ortho_graph(cab_vects)
#print ortho_graph
createSageGraph(ortho_graph)

def example_with_obs_list(obs_list): #This code is kept for reference. It shows how to use an observable list to build a graph.
    DIMENSION = preprocessObservables(obs_list)
    comm_graph = build_graph(obs_list)
    clique_list = bronkerbosch.find_cliques(comm_graph)

    max_clique_size = len(max(clique_list, key=lambda x: len(x)))
    clique_list = [map(lambda x: obs_list[x], list(q))  for q in clique_list]
    if USE_LARGEST_CLIQUES:
            clique_list = [w for w in clique_list if len(w)==max_clique_size]


    vectors = map(lambda q: simultaneouslyDiagonalize(q) , clique_list)

    vectors_copy = [] #vectors get flattened
    for v in vectors:
        for i in range(len(vectors[0])):
            vectors_copy.append(v[i])
    vectors_copy = createLabels(vectors_copy, DIMENSION)
    vectors_copy = removeRedundantVectors(vectors_copy)
    ortho_graph = build_ortho_graph(vectors_copy)
    createSageGraph(ortho_graph)
