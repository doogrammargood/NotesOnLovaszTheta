from cvxopt import matrix, solvers
import numpy as np
#c = matrix([1.,-1.,1.])
#G = [ matrix([[-7., -11., -11., 3.],
#                  [ 7., -18., -18., 8.],
#                  [-2.,  -8.,  -8., 1.]]) ]
#G += [ matrix([[-21., -11.,   0., -11.,  10.,   8.,   0.,   8., 5.],
#                   [  0.,  10.,  16.,  10., -10., -10.,  16., -10., 3.],
#                   [ -5.,   2., -17.,   2.,  -6.,   8., -17.,  8., 6.]]) ]
#h = [ matrix([[33., -9.], [-9., 26.]]) ]
#h += [ matrix([[14., 9., 40.], [9., 91., 10.], [40., 10., 15.]]) ]
#sol = solvers.sdp(c, Gs=G, hs=h)
#print(sol['x'])
#print(sol['zs'][0])
#print(sol['zs'][1])

A = np.array([[0., 1., 0., 0., 1., 1., 0., 0., 0., 0. ],
              [1., 0., 1., 0., 0., 0., 1., 0., 0., 0. ],
              [0., 1., 0., 1., 0., 0., 0., 1., 0., 0. ],
              [0., 0., 1., 0., 1., 0., 0., 0., 1., 0. ],
              [1., 0., 0., 1., 0., 0., 0., 0., 0., 1. ],
              [1., 0., 0., 0., 0., 0., 0., 1., 1., 0. ],
              [0., 1., 0., 0., 0., 0., 0., 0., 1., 1. ],
              [0., 0., 1., 0., 1., 1., 0., 0., 0., 1. ],
              [0., 0., 0., 1., 0., 1., 1., 0., 0., 0. ],
              [0., 0., 0., 0., 1., 0., 1., 1., 0., 0. ]])

I = np.array([[1., 0., 0., 0., 0., 0., 0., 0., 0., 0. ],
              [0., 1., 0., 0., 0., 0., 0., 0., 0., 0. ],
              [0., 0., 1., 0., 0., 0., 0., 0., 0., 0. ],
              [0., 0., 0., 1., 0., 0., 0., 0., 0., 0. ],
              [0., 0., 0., 0., 1., 0., 0., 0., 0., 0. ],
              [0., 0., 0., 0., 0., 1., 0., 0., 0., 0. ],
              [0., 0., 0., 0., 0., 0., 1., 0., 0., 0. ],
              [0., 0., 0., 0., 0., 0., 0., 1., 0., 0. ],
              [0., 0., 0., 0., 0., 0., 0., 0., 1., 0. ],
              [0., 0., 0., 0., 0., 0., 0., 0., 0., 1. ]])


#print np.ravel(A)
ONES = np.array([[1. for i in range(10)] for j in range(10)])
M = ONES - A - I
tedious_arrays = []
for i in range(10):
    for j in range(i):
        if M[i][j] > 0.1:
            tedious_arrays.append([[ 1. if ( (r==j and k==i) or (r==i and k==j) )else 0. for r in range(10) ] for k in range(10)])
tedious_arrays= map(lambda t: np.array(t), tedious_arrays)
#sum = A + I
#print tedious_arrays[0].type()
#for T in tedious_arrays:
#    sum = sum + T
#print sum
constraintMatrix = []
constraintMatrix.append(np.ravel(-I).tolist())
for T in tedious_arrays:
    constraintMatrix.append(np.ravel(-T).tolist())
G=[ matrix(constraintMatrix)]
H=[ matrix(-A - I)]
d=[1.]
for T in tedious_arrays:
    d.append(0.)
d = matrix(d)
print G
sol = solvers.sdp(matrix(d), Gs=G, hs=H)
print sol['x']
