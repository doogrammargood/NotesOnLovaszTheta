import numpy as np
import itertools

#The point of this script is to determine the structure of unitary matrices
#which transform the cliques in Cabello's 18-vector graph into one another.
EPSILON = 10**-3
base_vects = [[1,0,0,0],
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
base_vects = map(lambda x: np.array(x)/np.linalg.norm(x), base_vects )
print len(base_vects)

list_of_cliques = [ [0, 1, 2, 3],
                    [3, 4, 5, 6],
                    [6, 7, 8, 9],
                    [9, 10,11,12],
                    [12,13,14,15],
                    [15,16,17,0],
                    [17,1,8,10],
                    [2,4,11,13],
                    [5,7,14,16]]

def test_cliques(): #makes sure the cliques are pairwise orthogonal.
    for c in list_of_cliques:
        vectors_in_clique = [ b[1] for b in enumerate(base_vects) if (b[0] in c)]
        for pair in itertools.combinations(vectors_in_clique,2):
            if np.dot( pair[0], pair[1]) != 0:
                print "test failed."
def clean_matrix(M):#The matrix is hard to look at.
    EPSILON = 10**-3
    pretty = M.tolist()
    for i in range(len(pretty)):
        for j in range(len(pretty[0])):
            if (pretty[i][j]>EPSILON):
                pretty[i][j]=1
            elif (pretty[i][j]<-EPSILON):
                pretty[i][j]=-1
            else: pretty[i][j]=0
    return pretty
def pretty_print(M):
    clean = clean_matrix(M)
    for c in clean:
        print c

def get_unitary(i,j):
    #For cliques i,j returns the unitary transformation from i to j.
    c0 = list_of_cliques[i]
    c1 = list_of_cliques[j]
    # print c0
    # print c1
    v0 = [ b[1] for b in enumerate(base_vects) if (b[0] in c0)]
    v1 = [ b[1] for b in enumerate(base_vects) if (b[0] in c1)]
    # print v0
    # print v1
    m0 = np.matrix(v0)
    m1 = np.matrix(v1)
    # print m0
    # print m1
    unitary = m0.T*m1
    return unitary

def print_orbit_zero():
    c0 = list_of_cliques[0]
    v0 = [ b[1] for b in enumerate(base_vects) if (b[0] in c0)]

    for j in range(9):
        for k in range(9):
            for i in range(4):
                pretty_print(v0[i]*get_unitary(j,k))
            print "-----"
        print "next set:",i," ",j

def verify_clique(M):
    #tests whether the rows of M are one of the cliques.
    clique_found = False
    for i in range(9): #loops through the cliques
        c = list_of_cliques[i]
        v = [ b[1] for b in enumerate(base_vects) if (b[0] in c)]
        for j in range(4): #loops through the rows of M.
            vec_found = False
            for k in range(4): #loops through the rows of v
                #print np.linalg.norm(M[j]-v[k])
                if  np.linalg.norm(M[j]-v[k]) < EPSILON:
                    vec_found = True
            if not vec_found:
                break
            elif j==3:
                clique_found = True
    return clique_found


def test_orbit():
    #The claim is that the clique*get_unitary(j,k) is a clique in the graph.
    for i in range(9):
        c = list_of_cliques[i]
        clique = [ b[1] for b in enumerate(base_vects) if (b[0] in c)]
        for k in range(9):
            for j in range(9):
                if verify_clique(np.matrix(clique)*get_unitary(j,k)):
                    # print "out of the matrix!"
                    print i,j,k
                    # print clique
                    # print np.matrix(clique)
                    # print get_unitary(j,k)
                    pretty_print(np.matrix(clique)*get_unitary(j,k))

def test_orbit_zero():
        c = list_of_cliques[0]

        for i in range(9):
            pretty_print(np.matrix(clique)*get_unitary(0,1)**i)
            print "--"
            if not verify_clique(np.matrix(clique)*get_unitary(0,1)**i):
                print i," not there"
def check_for_automorphisms():
    #looks at all unitaries to see if any are automorphisms.
    for i in range(9):
        for j in range(9):
            u = get_unitary(i,j)
            automorphism = True
            for k in range(9):
                v = [ b[1] for b in enumerate(base_vects) if b[0] in list_of_cliques[k]]
                if not verify_clique(np.matrix(v)*u):
                    automorphism = False
            if automorphism:
                print "automorphism found"
                print i, j, get_unitary(i,j)


clique = [ b[1] for b in enumerate(base_vects) if (b[0] in list_of_cliques[0])]
#print np.matrix(clique)*get_unitary(7,8)
#print verify_clique(np.matrix(clique)*get_unitary(7,8))
check_for_automorphisms()
#test_orbit()
#print get_unitary(8,6)
# test_orbit_zero()
# for i in range(10):
#     print i
#     pretty_print(get_unitary(2,5)**i)
#     print "--"
