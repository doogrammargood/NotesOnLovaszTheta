import numpy.linalg as LA
import math
def dot(K,L):
    return sum(i[0] * i[1] for i in zip(K, L))
#Orthonormal Representation of the Petersen Graph
a=-1
b=4+15**0.5
c=-6-15**0.5
d= (3*a**2 + 6*b**2 + c**2) ** 0.5
print d
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
print dot(matrix[0],matrix[0])
handle = [1/10**0.5, 1/10**0.5, 1/10**0.5, 1/10**0.5, 1/10**0.5, 1/10**0.5, 1/10**0.5, 1/10**0.5, 1/10**0.5, 1/10**0.5 ]
print dot(handle,handle)
theta = 1/dot(matrix[9],handle)**2 #should be 4.
print theta
for i in range(10):
    print dot(matrix[i],handle)
# for row in matrix:
#      coloring_from_m.append( [(theta**0.5*r-handle[i])/(theta-1)**0.5 for i,r in enumerate(row)] )
print "dot"
coloring_from_m = [[ ((theta**0.5)*r-1/10**0.5) / ((theta-1)**0.5) for r in row] for row in matrix]
print dot(coloring_from_m[0],coloring_from_m[0])
print LA.matrix_rank(coloring_from_m)
a=-1/(18)**0.5
b=1/(18)**0.5
c=1/2**0.5
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

print "other"
print dot(matrix[0],matrix[1])
