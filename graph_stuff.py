from sage.all import *

g = {0: [1, 2, 3, 5, 7, 9, 10, 14, 15],
1: [0, 2, 3, 4, 6, 9, 10, 12, 13],
2: [0, 1, 3, 5, 7, 8, 11, 14, 15],
3: [0, 1, 2, 4, 6, 8, 11, 12, 13],
4: [1, 3, 5, 6, 7, 12, 14, 15, 17, 19, 21, 22],
5: [0, 2, 4, 6, 7, 12, 13, 14, 16, 18, 21, 22],
6: [1, 3, 4, 5, 7, 13, 14, 15, 16, 18, 20, 23],
7: [0, 2, 4, 5, 6, 12, 13, 15, 17, 19, 20, 23],
8: [2, 3, 9, 10, 11, 18, 19, 20, 22],
9: [0, 1, 8, 10, 11, 16, 17, 20, 22],
10: [0, 1, 8, 9, 11, 18, 19, 21, 23],
11: [2, 3, 8, 9, 10, 16, 17, 21, 23],
12: [1, 3, 4, 5, 7, 13, 14, 15, 16, 18, 20, 23],
13: [1, 3, 5, 6, 7, 12, 14, 15, 17, 19, 21, 22],
14: [0, 2, 4, 5, 6, 12, 13, 15, 17, 19, 20, 23],
15: [0, 2, 4, 6, 7, 12, 13, 14, 16, 18, 21, 22],
16: [5, 6, 9, 11, 12, 15, 17, 18, 19],
17: [4, 7, 9, 11, 13, 14, 16, 18, 19],
18: [5, 6, 8, 10, 12, 15, 16, 17, 19],
19: [4, 7, 8, 10, 13, 14, 16, 17, 18],
20: [6, 7, 8, 9, 12, 14, 21, 22, 23],
21: [4, 5, 10, 11, 13, 15, 20, 22, 23],
22: [4, 5, 8, 9, 13, 15, 20, 21, 23],
23: [6, 7, 10, 11, 12, 14, 20, 21, 22]}

G = Graph(g)

G.allow_multiple_edges(False)

g = G.plot()
save(g,'/tmp/dom.png',axes=False,aspect_ratio=True)
os.system('display /tmp/dom.png')
print len(G.vertices())
print len(G.edges())
print G.lovasz_theta()
print G.independent_set()