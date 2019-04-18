#This file assumes the ramunjan graphs have been created and we can unpickle them.
import pickle
from sage.all import *
import numpy as np
def retrieve_graph(p,q):
    adj_list = pickle.load( open( "/home/owg/Documents/projects/Contextuality/ramanujan_graphs/example_(%d)_(%d).pickle" %(p, q), "rb" ) )
    return adj_list

def test_ramanujan_property(p,q):
    graph = Graph(retrieve_graph(p,q))
    mat = graph.adjacency_matrix()
    eigenvalues = np.linalg.eigh(mat)[0]
    if kronecker(p,q)==1:
        if (-eigenvalues[0] < 2*sqrt(p) and eigenvalues[-2] < 2*sqrt(p)):
            print "test passed!"
    elif kronecker(p,q)==-1:
        if (-eigenvalues[1] < 2*sqrt(p) and eigenvalues[-2] < 2*sqrt(p)):
            print "test passed"
#test_ramanujan_property(5,17)
def get_vertex_connectivity():
    graph = Graph(retrieve_graph(13,17))
    print graph.complement().vertex_connectivity()

get_vertex_connectivity()
#
# adj_list = retrieve_graph(13,17)
# print len(adj_list.keys())
# G = Graph(adj_list)
# print len(G.independent_set())
# G = G.complement()
# print len(G.independent_set())
