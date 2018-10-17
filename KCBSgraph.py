from sage.all import *
from sage.graphs.modular_decomposition import modular_decomposition
import itertools

def compatible_list(vertex, list_of_verticies):
    #Finds the neighbors of vertex in list_of_verticies.
    a = len(vertex)/2 #Same a as above
    neighbors = []
    for i in range(a):
        neighbors += [v for v in list_of_verticies if vertex[i+a] in v[a:] and v[v[a:].index(vertex[i+a])] != vertex[i]]
    return neighbors

experiments = range(5)
experiments = [(str(e),str((e+1)%5)) for e in experiments]

outcomes = itertools.product(['0','1'],repeat = 2)
verticies = [''.join(o+e) for o,e in itertools.product(outcomes, experiments)]
graph = {node: compatible_list(node, verticies) for node in verticies}
G = Graph(graph)
G.allow_multiple_edges(False)

H = G.subgraph(('0112', '0023', '1040', '0001', '0134', '1123', '0140', '1034', '1012', '1101'))

g = H.plot()
save(g,'/tmp/dom.png',axes=False,aspect_ratio=True)
os.system('display /tmp/dom.png')
print len(verticies)

d = graph
counter = 0
length_subgraph = 0
max_discrepancy = 0
found_graphs = []
max_subgraph = []
# for n in range(len(d.keys())):
#     for H in itertools.combinations(d.keys(), n):
#         if len(H) > length_subgraph:
#             length_subgraph = len(H)
#             found_graphs = [] #clear found graphs each time H increases
#         counter += 1
#         if counter %1000 == 0:
#             print counter
#         S = G.subgraph(H)
#         label = S.canonical_label()
#         if not label in found_graphs:
#             found_graphs.append(label)
#             discrepancy = abs(S.lovasz_theta() - len(S.independent_set()))
#             if discrepancy > max_discrepancy:
#                 max_discrepancy = discrepancy
#                 max_subgraph = H
#                 print "Contextuality found.\n"
#                 print discrepancy
#                 print "Subgraph is:\n"
#                 print H

print 'we found a subgraph with discrepancy:'
print max_discrepancy
print 'the subgraph is'
print max_subgraph
