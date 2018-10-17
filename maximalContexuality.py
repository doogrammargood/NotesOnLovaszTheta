from sage.all import *
from sage.graphs.modular_decomposition import modular_decomposition
from sage.graphs.independent_sets import IndependentSets
import itertools


def add_vertex_without_increasing_independence(g):
    #expect g to be a dict.
    n = len(g.vertices())
    new_vertex = str(n)
    target_vertices = g.vertices()
    indyNumber = g.independent_set(value_only=True)
    maximal_independent_sets = [v for v in IndependentSets(g) if len(v) == indyNumber]
    target_vertices = []
    max_indep_sets = 0
    g_copy = g.copy()
    while(g_copy.independent_set(value_only=True) == indyNumber):
        for vertex in g_copy.vertices():
            indep_sets_for_vertex = len([indy for indy in maximal_independent_sets if vertex in indy])
            if indep_sets_for_vertex > max_indep_sets:
                indep_sets_for_vertex = max_indep_sets
                max_vertex = vertex
        target_vertices.append(max_vertex)
        g_copy.delete_vertex(max_vertex)
    g.add_vertex(new_vertex)
    g.add_edges([(new_vertex,v) for v in target_vertices])

#g = Graph({'0':['1','2'], '1':['0','2','3'], '2':['0','1'], '3':['1','4']})
g = graphs.CycleGraph(8)
for k in range(10):
    add_vertex_without_increasing_independence(g)
print "the lovasz number is"
print g.lovasz_theta()
print "the independence number"
print g.independent_set(value_only=True)
print "is it perfect?"
print g.is_perfect()

g.allow_multiple_edges(False)
s = g.plot()
save(s,'/tmp/dom.png',axes=False,aspect_ratio=True)
os.system('display /tmp/dom.png')
