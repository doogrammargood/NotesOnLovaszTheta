import numpy as np
import pickle

##This file uses the gram schmidt process to create an orthonormal representation, according to
##Theorem 5.2.2 of http://web.cs.elte.hu/~lovasz/kurzusok/orth16-2.pdf
graph = pickle.load( open( "./ramanujan_graphs/example_(13)_(17).pickle", "rb" ) )


# def keys_to_ints(graph):#This was supposed to have been done in Canonize dict. Bad code smell.
#     #This creates a dictionary where you can look up the key to int conversion.
#     keys = sorted(graph.keys())
#     #Returns a dictionary from keys of the graph to consecutive integers.
#     to_return = {}
#     counter = 0
#     for k in keys:
#         to_return[k]=counter
#         counter +=1
#     return to_return


#The theorem garauntees success when k is the connectivity of the complement of our graph
a=8
vectors = [np.random.rand(2448-a) for i in range(len(graph.keys()))]

epsilon = 10**-2
def straighten(vector, list):
    for l in list:
        vector = vector - np.dot(vector,l)*l
    if np.linalg.norm(vector)<epsilon:
        print "too small"
    vector = vector/np.linalg.norm(vector)
    return vector


for k in graph.keys():
    orthogonal_vectors = graph[k]
    orthogonal_vectors = filter(lambda x, k=k: x<k ,orthogonal_vectors)
    orthogonal_vectors = [ vectors[o] for o in orthogonal_vectors]
    #print orthogonal_vectors
    #print kti_dict[5126]

    vectors[k] = straighten(vectors[k], orthogonal_vectors)
    # vectors[k] = straighten(vectors[k], orthogonal_vectors)
    # vectors[k] = straighten(vectors[k], orthogonal_vectors)
#print vectors

for k in graph.keys():
    neighbors = graph[k]
    for n in neighbors:
        if np.dot(vectors[k],vectors[n]) > epsilon:
            print "some error"
