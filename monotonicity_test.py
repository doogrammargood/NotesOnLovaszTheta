
from sage.all import *
import itertools
import multiprocessing
import profile
#vertices = ['0','1','2','3','4','5','6','7','8','9','10']
vertices = [0,1,2,3,4,5,6]



def canonize_from_edgeset(edgeset):
    #turns a list of tuples into a pair [graph, canonical_label]
    d = {node: [] for node in vertices}
    for e in edgeset:
        d[e[0]].append(e[1])
    g = Graph(d)
    l = g.canonical_label()
    return [g,l]
all_edges = powerset(itertools.combinations(vertices, 2))
pool = multiprocessing.Pool()
canonized_graphs = pool.map(canonize_from_edgeset, all_edges)
#all_graphs = [g for g in canonized_graphs]
canonized_graphs = [g[0] for g in canonized_graphs if g[0]==g[1]]
thetas = map(lambda r: r.lovasz_theta(), canonized_graphs)
table = zip(canonized_graphs, thetas)
table = sorted(table, key = lambda t: t[0].size())


lower_critical_graphs = [] #removing any edge decreases theta.
for i in range(len(table)):
    critical = True
    for j in range(0,i):
        if table[j][0].is_subgraph(table[i][0], induced=False) and table[j][1] == table[i][1] and i!=j:
            critical = False
    if critical:
        lower_critical_graphs.append(table[i])

upper_critical_graphs = []
for i in range(len(table)):
    critical = True
    for j in range(i+1,len(table) ):
        if table[i][0].is_subgraph(table[j][0], induced=False) and table[j][1] == table[i][1] and i!=j:
            critical = False
    if critical:
        upper_critical_graphs.append(table[i])
print "lower_critical"
print len(lower_critical_graphs)
print "upper_critical"
print len(upper_critical_graphs)
low=[]
for i in range(len(lower_critical_graphs)):
    low.append( lower_critical_graphs[i][0].plot() )
    save(low[i],'/tmp/critical_graphs/%flower%d.png' % (lower_critical_graphs[i][1], i),axes=False,aspect_ratio=True)
    #os.system('display /tmp/critical_graphs/%flower%d.png' % (lower_critical_graphs[i][1], i))
up=[]
for i in range(len(upper_critical_graphs)):
    up.append( upper_critical_graphs[i][0].plot() )
    save(up[i],'/tmp/critical_graphs/%fupper%d.png' % (upper_critical_graphs[i][1], i),axes=False,aspect_ratio=True)
    #os.system('display /tmp/critical_graphs/%fupper%d.png' % (upper_critical_graphs[i][1], i))

common = [pair for pair in lower_critical_graphs if pair in upper_critical_graphs]
for i in range(len(common)):
    up.append( common[i][0].plot() )
    save(up[i],'/tmp/critical_graphs/%fboth%d.png' % (common[i][1], i),axes=False,aspect_ratio=True)
    #os.system('display /tmp/critical_graphs/%fboth%d.png' % (common[i][1], i))
