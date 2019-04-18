import itertools
import numpy as np
from numpy import linalg as LA
import itertools
import math
import profile
from multiprocessing import Process, Manager, Pool
import pickle

class int_counter:
    #produces a integer each time it's called.
    i = 0
    def tick(self):
        self.i+=1
        return self.i -1

def canonize_dict(dict):
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
        replace(key, replacement_mapping, counter)
    for key in dict.keys():
        new_neighbors = []
        print key
        print replacement_mapping[key]
        print dict[key]
        new_dict[replacement_mapping[key]]=[replacement_mapping[x] for x in dict[key]]
    return new_dict
'''
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
        #print key
        new_key = replace(key, replacement_mapping, counter)
        new_neighbors = []
        for k in dict[key]:
            new_neighbors.append(replace(k, replacement_mapping,counter))
        #new_neighbors = [replace(k, replacement_mapping,counter) for k in dict[key] ]
        new_dict[new_key] = new_neighbors
    #print replacement_mapping
    return replacement_mapping
    #return new_dict
    '''
def kronecker(p,q):
    #The legendre symbol. Assume q an odd prime.
    if p**((q-1)/2) %q == 1:
        return 1
    else:
        return -1

#Next two functions come from https://stackoverflow.com/questions/4798654/modular-multiplicative-inverse-function-in-python
def egcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)
def modinv(a, m):
    g, x, y = egcd(a, m)
    if g != 1:
        raise Exception('modular inverse does not exist')
    else:
        return x % m
def add_signs(x):
    to_return = []
    if x == []:
        return [[]]
    if x[0] >0:
        to_return.extend( [[x[0]] + h for h in add_signs(x[1:])])
        to_return.extend( [[-x[0]] + h for h in add_signs(x[1:])])
    else:
        to_return.extend( [[x[0]] + h for h in add_signs(x[1:])])
    return to_return
def sum_of_four_squares(n):
    list = []
    top = int(math.ceil(n**0.5) + 1)
    for i in range(top):
        for j in range(top):
            for k in range(top):
                for h in range(top):
                    if i**2 + j**2 +k**2 +h**2 == n:
                        new_elements = add_signs([i,j,k,h])
                        list.extend(new_elements)
    return list

def sum_divisors(n):
    #sums the divisors of n
    current_sum = 0
    for i in range(1,n+1):
        if n%i==0:
            current_sum += i
    return current_sum
def test_sum_of_four_squares():
    for i in range(1,30):
        if i%2==1:
            sos = sum_of_four_squares(i)
            if not len(sos) == 8*sum_divisors(i):
                print "failed sum of squares test"
                print i
                print sos
                print len(sos)
                print sum_divisors(i)
                return
            for s in sos:
                if s[0]**2+s[1]**2+s[2]**2+s[3]**2 !=i:
                    print "failed sos test"
                    return
    print "passed test"
    return

#test_sum_of_four_squares()
#print sum_of_four_squares(26)
class glement(object):#implements the general linear group on Z
    def __init__(self,q,matrix):
        self.q = q
        self.matrix=[[matrix[0][0]%q, matrix[0][1]%q], [matrix[1][0]%q, matrix[1][1]%q]] #make sure the entries are modulo q.
    def determinant(self):
        return (self.matrix[0][0]*self.matrix[1][1] - self.matrix[1][0]*self.matrix[0][1])%self.q

    def string_name(self):
        return str(self.canonize_matrix())
    def __mul__(self, another_matrix):
        if self.q != another_matrix.q:
            print "can't multiply these matricies!"
        a = self.matrix
        b = another_matrix.matrix
        q = self.q
        e1= (a[0][0]*b[0][0]+a[0][1]*b[1][0]) %q
        e2= (a[0][0]*b[0][1]+a[0][1]*b[1][1]) %q
        e3= (a[1][0]*b[0][0]+a[1][1]*b[1][0]) %q
        e4= (a[1][0]*b[0][1]+a[1][1]*b[1][1]) %q
        new_matrix = [[ e1,e2 ],[ e3,e4 ]]
        return new_matrix


class pslement(glement):
        #pselement is exactly the same as glement, except that we identify multiplication by [[-1,0][0,-1]]
        #We also need to check that the determinant is 1. This is currently a weakness in this code.
    def __init__(self, q, matrix):
        glement.__init__(self, q, matrix)
    def __eq__(self, another_matrix):
        if self.q != another_matrix.q:
            return False
        a = self.matrix
        b = another_matrix.matrix
        q = self.q
        if (a[0][0] - b[0][0]) %q ==0 and (a[0][1]-b[0][1])%q==0 and (a[1][0]-b[1][0])%q==0 and (a[1][1]-b[1][1])%q==0:
            return True
        elif (a[0][0] + b[0][0]) %q ==0 and (a[0][1]+b[0][1])%q==0 and (a[1][0]+b[1][0])%q==0 and (a[1][1]+b[1][1])%q==0:
            return True
        return False
    def inv(self):
        a = self.matrix[0][0]
        b = self.matrix[0][1]
        c = self.matrix[1][0]
        d = self.matrix[1][1]
        q = self.q
        det = (a*d-b*c)%q
        i = modinv(det,q)
        return pslement(q,[[(i*d)%q, (i*-b)%q],[(i*-c)%q,(i*a)%q]])
    def canonize_matrix(self):
        #returns a canonical matrix. We require that self.matrix[0][0] be in [1,floor (q/2) ] if it isn't 0.
        if self.matrix[0][0]!=0:
            if self.matrix[0][0] in range(1,self.q/2+1):
                return self.matrix
            else:
                return [[(-self.matrix[0][0])%self.q, (-self.matrix[0][1])%self.q ],[(-self.matrix[1][0])%self.q, (-self.matrix[1][1])%self.q]]
        else:
            if self.matrix[0][1] in range(1,self.q/2+1):
                return self.matrix
            else:
                return [[(-self.matrix[0][0])%self.q, (-self.matrix[0][1])%self.q ],[(-self.matrix[1][0])%self.q, (-self.matrix[1][1])%self.q]]
    def __mul__(self, another_matrix):
        if self.q != another_matrix.q:
            print "can't multiply these matricies!"
        a = self.matrix
        b = another_matrix.matrix
        q = self.q
        e1= (a[0][0]*b[0][0]+a[0][1]*b[1][0]) %q
        e2= (a[0][0]*b[0][1]+a[0][1]*b[1][1]) %q
        e3= (a[1][0]*b[0][0]+a[1][1]*b[1][0]) %q
        e4= (a[1][0]*b[0][1]+a[1][1]*b[1][1]) %q
        new_matrix = [[ e1,e2 ],[ e3,e4 ]]
        return pslement(self.q,new_matrix)


class pglement(glement):
    def __init(self, q, matrix):
        glement.__init__(self, q, matrix)
    def __eq__(self, another_matrix):
        if self.q != another_matrix.q:
            return False
        a = self.matrix
        b = another_matrix.matrix
        q = self.q
        for r in range(q):
            if (a[0][0] - r*b[0][0]) %q ==0 and (a[0][1] - r*b[0][1])%q==0 and (a[1][0] - r*b[1][0])%q==0 and (a[1][1] - r*b[1][1])%q==0:
                return True
        return False
    def canonize_matrix(self):
        #We make the first entry self.matrix[0][0]==1 unless it is 0.
        #In this case, we make self.matrix[0][1]==1
        if self.matrix[0][0]%self.q!=0:
            r = modinv(self.matrix[0][0],self.q)
        else:
            r = modinv(self.matrix[0][1],self.q)
        return [[(r*self.matrix[0][0])%self.q, (r*self.matrix[0][1])%self.q ],[(r*self.matrix[1][0])%self.q, (r*self.matrix[1][1])%self.q]]
    def __mul__(self, another_matrix):
        if self.q != another_matrix.q:
            print "can't multiply these matricies!"
        a = self.matrix
        b = another_matrix.matrix
        q = self.q
        e1= (a[0][0]*b[0][0]+a[0][1]*b[1][0]) %q
        e2= (a[0][0]*b[0][1]+a[0][1]*b[1][1]) %q
        e3= (a[1][0]*b[0][0]+a[1][1]*b[1][0]) %q
        e4= (a[1][0]*b[0][1]+a[1][1]*b[1][1]) %q
        new_matrix = [[ e1,e2 ],[ e3,e4 ]]
        return pglement(self.q,new_matrix)

def _total_list_psl(q):#used for testing only.
    list = []
    for a,b,c,d in itertools.product(range(q), repeat = 4):
        c = pslement(q,[[a,b],[c,d]])
        if c.determinant() ==1:
            list.append(c)
    return list

def list_psl(q):
    #returns a list of all the elements of psl.
    list = []
    for b,d in itertools.product(range(1,q/2+1), range(q)):
        c = (-modinv(b,q))%q
        mat = pslement(q, [[0,b],[c,d]])
        list.append(mat)
    for a,b,c in itertools.product(range(1,q/2+1), range(q), range(q)):
        d = ((c*b+1)*modinv(a,q) )%q
        mat = pslement(q,[[a,b],[c,d]])
        list.append(mat)
    # for a,b,c,d in itertools.product(range(q), repeat = 4):
    #     mat = pslement(q,[[a,b],[c,d]])
    #     if mat.determinant()==1:# and not mat in list:
    #         if mat.matrix[0][0]==0 and mat.matrix[0][1] in range(1,q/2+1):
    #             list.append(mat)
    #         elif mat.matrix[0][0] in range(1,q/2+1):
    #             list.append(mat)
    return list
# def add_to_list_gsl(q, tuple, list):
#     a,b,c,d = tuple
#     mat = pglement(q,[[a,b],[c,d]])
#     if mat.determinant()!=0:# and not mat in list:
#         if mat.matrix[0][0]==0 and mat.matrix[0][1]==1:
#             list.append(mat)
#         elif mat.matrix[0][0]==1:
#             list.append(mat)
#
# def list_gsl(q):
#     man = Manager()
#     list = man.list([])
#     pool = Pool(processes=8)
#     for tuple in itertools.product(range(q), repeat = 4):
#         pool.apply_async(add_to_list_gsl, args = (q, tuple, list))
#     return list

# #
def list_pgl(q):
    #lists the elements of pgl
    list = []
    for b,c,d in itertools.product(range(q), repeat =3):
        mat = pglement(q,[[1,b],[c,d]])
        if not mat.determinant() == 0: #determinant must not be 0.
            list.append(mat)
    for c,d in itertools.product(range(1,q), range(q)):
        mat = pglement(q, [[0,1],[c,d]])
        list.append(mat)
    # for a,b,c,d in itertools.product(range(q), repeat = 4):
    #     mat = pglement(q,[[a,b],[c,d]])
    #     if mat.determinant()!=0:# and not mat in list:
    #         if mat.matrix[0][0]==0 and mat.matrix[0][1]==1:
    #             list.append(mat)
    #         elif mat.matrix[0][0]==1:
    #             list.append(mat)
    return list

def list_special_generators(p,q):
    #These come from the number of ways to write p as the sum of 4 squares.
    answers = sum_of_four_squares(p)
    i=0
    S=[]
    for x in range(q):
        if x**2 % q == q-1:
            i=x
            break
    if kronecker(p,q)==1:
        for c in range(q):
            if (c**2 -p) %q ==0:
                normalizer = modinv(c,q)
                break
        for a in answers:
            if a[0]>0 and a[0]%2==1 and a[1]%2==0 and a[2]%2==0 and a[3]%2==0:
                el = pslement(q,[ [ normalizer*(a[0]+i*a[1]), normalizer*(a[2]+i*a[3])], [ normalizer*(-a[2]+i*a[3]), normalizer*(a[0]-i*a[1])]] )
                S.append(el)
        return S
    elif kronecker(p,q)==-1:
        for a in answers:
            if a[0]>0 and a[0]%2==1 and a[1]%2==0 and a[2]%2==0 and a[3]%2==0:
                el = pglement(q,[ [ (a[0]+i*a[1]), (a[2]+i*a[3])], [ (-a[2]+i*a[3]), (a[0]-i*a[1])]] )
                S.append(el)
        return S

# def add_nodes(dict, node1, S):
#     neighbors = []
#     for s in S:
#         neighbors.append((node1*s).string_name())
#     dict[node1.string_name()]= neighbors
#     return

def createRamanujanGraph(p,q):
    #assume p and q are primes congruent to 1 mod 4
    if kronecker(p,q)==1:
        total_nodes = list_psl(q)
    elif kronecker(p,q)==-1:
        total_nodes = list_pgl(q)
    S = list_special_generators(p,q) #this set defines the edges in the graph
    adjacency_list = {}
    # manager = Manager()
    # dict = manager.dict()

    # pool = Pool(processes=4)
    for node in total_nodes:
    #     pool.apply_async(add_nodes, args=(dict,node,S))
        #add_nodes(dict, node)
        neighbors = []
        for s in S:
            neighbors.append((node*s).string_name())
        adjacency_list[node.string_name()]= neighbors
    #dict[1]=1
    #print dict
    return adjacency_list

#def test_psl_items(p,q):

    #R = createRamanujanGraph(13,17)


def test_list_psl(q):
    list = list_psl(q)
    if len(list) != q*(q**2-1)/2:
        print "failed to list everything"
    for q1, q2 in itertools.combinations(list,2):
        if q1 == q2:
            print "messed up"
    print "test complete"

def test_canonize_dict():
    R = createRamanujanGraph(13,17)
    q = 17
    #print R
    if not len(R.keys()) == q*(q**2-1)/2:
        print "wrong number of vertices!"
    C = canonize_dict(R)
    print C.keys()
    if not len(C.keys()) == q*(q**2-1)/2:
        print "We dropped some vertices!"
    if C.keys()[-1] != q*(q**2-1)/2 -1:
        print "cannonized wrong."
    #print C.keys()
    print "test canonize dict complete!"

def test_canonize_dict2():
    #have modified canonize_dict to return the replacement mapping
    R = createRamanujanGraph(13,17)
    #print R['[[0, 13], [13, 11]]']
    #print len(R.keys())
    q = 17
    print filter(lambda x: x[2]=='0', R.keys())
    if not len(R.keys()) == q*(q**2-1)/2:
        print "wrong number of vertices!"
    C = canonize_dict(R)
    total_strings = C.keys()
    print len(total_strings)

def test_rmg(p,q):
    R = createRamanujanGraph(p,q)
    #print R
    for r in R.keys():
        if not len(R[r])==p+1:
            print "test failed regularity"
        if not r in R.keys():
            print "test failed"
    if not len(R.keys()) == q*(q**2-1)/2:
        print "wrong number of vertices!"
    print "test complete"
def profile_test():
    R= createRamanujanGraph(17,59)
    R = canonize_dict(R)
def is_prime(num):
    for i in range(2,int(math.sqrt(num)+1)):
        if num%i==0:
            return False
    return True

def write_graphs():
    sum = 18
    while(True):
        for p in range(2,sum/2+1):
            if p%4 == 1 and is_prime(p):
                q = sum - p
                if is_prime(q) and q>2*math.sqrt(p) and p!= q:
                    print p,q
                    R = createRamanujanGraph(p,q)
                    R = canonize_dict(R)
                    print R
                    pickle_out = open("/home/owg/Documents/projects/Contextuality/ramanujan_graphs/example_(%d)_(%d).pickle" %(p, q), "wb")
                    pickle.dump(R,pickle_out)
                    pickle_out.close()

        sum+=4
def test_psl_again():
    list = _total_list_psl(17)
    print len(list)
    for a,b in itertools.combinations(list,2):
        if a == b:
            if not a.string_name() == b.string_name():
                print "string name problem"
    print "test psl complete"


write_graphs()
#profile.run('profile_test()')
#test_list_psl(17)
#test_rmg(13,17)

#test_canonize_dict()
#test_psl_again()
