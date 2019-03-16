import numpy as np
from numpy import linalg as LA
import itertools
import bronkerbosch
import re
from sage.all import *


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

def calculate_separation(G):
    #given a graph G, compute the ratio of a/t.
    a = len(G.independent_set())
    t = G.lovasz_theta()
    return a/t

def print_table_of_separation_values(max):
    best_sep = 1
    best_n = 0
    for n in range(2,max):
        if n in Primes():
            G = construct_chorded_cycles(n)
            sep = calculate_separation(G)
            if sep < best_sep:
                best_sep = sep
                best_n = n
            print n, " . . . ", sep
    print "best values:"
    print best_n, " . . .", best_sep


def rotation_map_from_dict(dict):
    return(lambda x,y: [dict[x][y], dict[dict[x][y]].index(x)])

def dict_from_rotation(r, num_vertices, degree):
    dict = {}
    for v in range(num_vertices):
        neighbors=[]
        for d in range(degree):
            neighbors.append(r(v,d)[0])
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
        print key
        new_key = replace(key, replacement_mapping, counter)
        print counter.tick()
        print new_key
        new_neighbors = [replace(k, replacement_mapping,counter) for k in dict[key] ]
        new_dict[new_key] = new_neighbors
    return new_dict

def zig_zag(G,H):
    lambda v,a,i,j:

    #assumes that G,H


z = rotation_map_from_dict({0: [1,2], 1: [2,0], 2: [1,0]})
#print z(1,0)
#print dict_from_rotation(z, 3,2)

print canonize_dict({'a': ['b','c'], 'b': ['c','a'], 'c': ['b','a']})
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
