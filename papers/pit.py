from cvxpy import *
from numpy import *

yy = array([[21, 23.5, 23, 24, 21, 25, 21.5, 22, 19, 23.5,25]]).T

xx = array([[8, 8, 8, 10, 10, 10, 12, 12, 12, 14, 14]]).T

x = Variable(11)
obj = Minimize(Pnorm(yy - x, 2))

b1 = diag([1./3, 1./3, 1./3, 0., 0., 0., 0, 0, 0, 0, 0])
b2 = diag([0, 0, 0, 1./3, 1./3, 1./3, 0., 0., 0., 0, 0])
b3 = diag([0, 0, 0, 0, 0, 0, 1./3, 1./3, 1./3, 0., 0.])
b4 = diag([0, 0, 0, 0., 0., 0., 0, 0, 0, 0.5, 0.5]) 

A = array([[0., 0., 0., 1., 1., 1.]])
x = Variable(6)
obj = Minimize(A * x)

B  = array([[1., 1., 2, 0., 0., 0.]])
c1 = B * x == 1.

C = 
c2 = C * x >= 0.

D = array([[1., 1., 0., 0., 0., 0.]])
E = array([[0., 0., 0., 1., 0., 0.]])
c3 = norm(D * x, 2) <= E * x

F1 = array([[0., 0., 0., 0., 1., 1.]])
F2 = diag([0., 0., 0., 0., 0., 0.])
F2[2, 2] = sqrt(2.0)
F2[4, 4] = 1.
F2[4, 5] = -1.
c4 = norm(F2 * x, 2) <= F1 * x

constraints = [c1, c2, c3, c4]

p = Problem(obj, constraints)

## b MOSEK.solve
p.solve(verbose = True, solver = "MOSEK")
p.value
x.value


