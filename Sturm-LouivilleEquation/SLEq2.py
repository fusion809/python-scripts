#!/usr/bin/env python3
# Second attempt at solving the Sturm-Liouville problem:
#  - d2y/dx2 + kxy = lambda y
# where lambda is the eigenvalue and y is the eigenfunction. 

from numpy import zeros
from numpy import cos
from numpy import sin
from numpy import arccos
from numpy import pi
from numpy import sqrt
from numpy import diag
from scipy import linalg as LA
N  = 1000
a  = 0
b  = 10
ya = 0
yb = 0
k  = 1
x = zeros(N+1)
T = zeros((N+1,N+1), dtype=float)
d2T = zeros((N+1,N+1), dtype=float)

x[0] = -1
x[N] = 1
for m in range(0,N+1):
    T[0,m]   = (-1)**m
    T[N,m]   = 1
    d2T[0,m] = (m**2)*((m**2 - 1)/3)*(-1)**m
    d2T[N,m] = (m**2)*(m**2 - 1)/3
    for n in range(1,N):
        x[n] = -cos(pi*n/N)
        T[n,m] = cos(m*arccos(x[n]))
        d2T[n,m] = (1-x[n]**2)*((m*x[n]*sin(m*arccos(x[n])))/sqrt(1-x[n]**2)-(m**2)*T[n,m])

y = (b-a)/2*x + (a+b)/2
H = -4/((b-a)**2)*d2T + k*diag(y)*T
Hn = LA.inv(T) * H
LAM, Y = LA.eig(Hn)
