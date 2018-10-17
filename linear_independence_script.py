import numpy.linalg as LA
def dot(K,L):
    return sum(i[0] * i[1] for i in zip(K, L))
#Orthonormal Representation of the Petersen Graph
a=1
b=-4-15**0.5
c=6+15**0.5
d= (3*a**2 + 6*b**2 + c**2) ** 0.5
a = a/d
b = b/d
c = c/d
#          1  2  3  4  5  6  7  8  9 10
matrix = [[c, b, b, b, b, b, b, a, a, a], #1
          [b, c, b, b, b, a, a, b, b, a], #2
          [b, b, c, b, a, b, a, b, a, b], #3
          [b, b, b, c, a, a, b, a, b, b], #4
          [b, b, a, a, c, b, b, b, b, a], #5
          [b, a, b, a, b, c, b, b, a, b], #6
          [b, a, a, b, b, b, c, a, b, b], #7
          [a, b, b, a, b, b, a, c, b, b], #8
          [a, b, a, b, b, a, b, b, c, b], #9
          [a, a, b, b, a, b, b, b, b, c]] #10
handle = [1/10**0.5]*10
theta = 1/dot(matrix[0],handle)**2 #should be 4.

print "adjancent angle"
print dot(matrix[0],matrix[9])
print LA.matrix_rank(matrix) #should be 6
coloring_from_m = []
#The coloring from this is coloring_from_m.
for row in matrix:
     coloring_from_m.append( [(theta**0.5*r-handle[i])/(theta-1)**0.5 for i,r in enumerate(row)] )
print "coloring from m"
print coloring_from_m
print LA.matrix_rank(coloring_from_m)
coloring_from_m.append(handle)
print "coloring with handle"
print LA.matrix_rank(coloring_from_m)
a=-1/18**0.5
b=-a
c=-1/2**0.5

a = a/d
b = b/d
c = c/d

print LA.matrix_rank(matrix)
for i in range(10):
    if 1/dot(matrix[i],handle)**2 != 4.0:
        print "problem at"
        print i
        print 1/dot(matrix[i],handle)**2

a = -3
b = 2
d = (2*a**2 +3*b**2)**0.5


a = a/d
b = b/d
chi = [[a, a, b, b, b],
       [a, b, a, b, b],
       [a, b, b, a, b],
       [a, b, b, b, a],
       [b, a, a, b, b],
       [b, a, b, a, b],
       [b, a, b, b, a],
       [b, b, a, a, b],
       [b, b, a, b, a],
       [b, b, b, a, a]]
print LA.matrix_rank(chi)
