from sage.all import *
from sage.graphs.modular_decomposition import modular_decomposition
import itertools

def compatible_list(vertex, list_of_verticies):
    #Finds the neighbors of vertex in list_of_verticies.
    a = len(vertex)/2 #Same a as above
    neighbors = []
    for i in range(a):
        neighbors += [v for v in list_of_verticies if v[i] != vertex[i] and v[i+a]==vertex[i+a]]
    return neighbors

def create_bell_graph(a,b):
    '''
    There are a qubits distributed among a people, each of which can perform b binary experiments.
    When a=b=2, we have Bell's experiment.
    This function returns the exclusivity graph of the experiment.
    The verticies are experiments and their outcomes represented as a string of 1's & 0's, (for the outcomes) followed by a string of length a with letters [0..b]
    They are connected if the two outcomes could not have occured together.
    More precisely, they are connected if the two verticies share an experiment which yeilds different measurements.
    '''
    experiment_choices = range(b)
    experiment_choices = [str(exp) for exp in experiment_choices] #convert to string
    total_experiments = itertools.product(experiment_choices, repeat=a)
    total_outcomes = itertools.product(['0','1'], repeat=a)
    verticies = [''.join(o+e) for o,e in itertools.product(total_outcomes, total_experiments)]
    graph = {node: compatible_list(node, verticies) for node in verticies}
    return graph

def addZeros(v):
    while(len(v)<4):
        v='0'+v
    return v

def filterMaximalSubgraph(vertex):
    #returns the verticies which are in the subgraph which maximazes discrepancy.
    length = len(vertex)
    a = length/2
    aor3 = min(a,3)
    experimental_parity = reduce(lambda x,y: (int(x)+int(y)) % 2, vertex[0:aor3], 0)
    #if len(list(set(vertex[a:a+aor3]))) != 1 and experimental_parity == 0: #The first 3 tests are different.
    if vertex[a:a+2] != '11' and experimental_parity == 0:
    #    print vertex[a:a+2]
        return True
    #elif len(list(set(vertex[a:a+aor3]))) == 1 and experimental_parity == 1:
    elif vertex[a:a+2] == '11' and experimental_parity == 1:
        return True
    else:
        return False

def bound_theta(S):
    return S.order()*max((-1*S.adjacency_matrix()).eigenvalues())/(max(S.adjacency_matrix().eigenvalues()) + max((-1*S.adjacency_matrix()).eigenvalues()))

#G = Graph(d)
d = create_bell_graph(2,2)
G = Graph(d)
G.allow_multiple_edges(False)
H = [v for v in d.keys() if filterMaximalSubgraph(v)]
S = G.subgraph(H)
S.add_vertex('c')
s.add_edges
Scomp = S.complement()
print 'subrgaph search'
print Scomp.subgraph_search(S, induced = false)

#print 'bound'
#print bound_theta(S)

'''if S.is_vertex_transitive():
    print 'vetex transitive'
else:
    print 'not vt'
if S.canonical_label() == graphs.CubeGraph(3).canonical_label():
    print 'cube'
else:
    print 'not a cube'
'''
#g = S.plot()
#save(g,'/tmp/dom.png',axes=False,aspect_ratio=True)
#os.system('display /tmp/dom.png')
#print verticies
max_discrepancy = 0
length_subgraph = 0
print 'theta'
print S.lovasz_theta()
print 'independent set'
print len(S.independent_set())

counter = 0 #number of subgraphs searched
found_graphs = []

'''
for n in range(len(d.keys())):
    for H in itertools.combinations(d.keys(), n):
        if len(H) > length_subgraph:
            length_subgraph = len(H)
            found_graphs = [] #clear found graphs each time H increases
        counter += 1
        if counter %1000 == 0:
            print counter
        S = G.subgraph(H)
        label = S.canonical_label()
        if not label in found_graphs:
            found_graphs.append(label)
            discrepancy = abs(S.lovasz_theta() - len(S.independent_set()))
            if discrepancy > max_discrepancy:
                max_discrepancy = discrepancy
                print "Contextuality found.\n"
                print discrepancy
                print "Subgraph is:\n"
                print H



'''
#print verticies
#print G.lovasz_theta()
#print 'hello'
