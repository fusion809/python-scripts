#!/usr/bin/env python3
# Import modules
import numpy as np
from scipy import linalg as LA
from matplotlib import pyplot as plt
from scipy.special import airy
import time

# Start timer
start = time.time()

# Number of steps
N            = 10000
# The domain we're transforming to is [a,b]
b            = 700
a            = 0
na           = np.arange(N+1)
# Chebyshev extreme grid
xa           = -np.cos(np.pi*na/N)
# Matrix of T_n(x_m) values
T            = np.cos(np.arccos(np.reshape(xa,(N+1,1)))*np.reshape(na,(1,N+1)))
# Matrix of T'_n(x_m) values
dT           = np.zeros((N+1,N+1))
dT[1:N,:]    = np.matmul(np.matmul(np.diag(np.divide(1,np.sqrt(1-np.power(xa[1:N],2)))),np.sin(np.arccos(np.reshape(xa[1:N],(N-1,1)))*np.reshape(na,(1,N+1)))),np.diag(na))
dT[0,:]      = np.multiply(np.power(-1,na+1),np.power(na,2))
dT[N,:]      = np.power(na,2)

# Differentiation matrices
D1           = np.matmul(dT,LA.inv(T))
# Save RAM by clearing dT which is unused after this
del dT
# D2 = D1^2
D2           = np.matmul(D1,D1)
# Save RAM by clearing D1 which is unused from here on
del D1
# Linearly transformed independent variable to the problem
ysub         = (b-a)/2*xa[1:N]+(a+b)/2
# Save even more RAM by clearing na and xa
del na, xa
# Eigenvalues and eigenvectors
values, vecs = LA.eig(-(4/(b-a)**2)*D2[1:N,1:N]+np.diag(ysub))
# Order eigenvalues and eigenvectors
values       = np.abs(values)
idx          = values.argsort()[::1]   
values       = values[idx]
vecs         = vecs[:,idx]

def f(xinput):
    x0=xinput
    xoutput=x0
    Ai=airy(-xoutput)[0]

    while np.abs(Ai)>1e-12:
        ai=airy(-xoutput)
        Ai=ai[0]
        Aip=ai[1]
        xoutput=xoutput+Ai/Aip

    return xoutput

NN=3000
from numpy import zeros
exact_values=zeros([NN+1,1])
error_in_values=zeros([NN+1,1])
for i in range(0,NN+1):
    exact_values[i] = f(values[i])
    error_in_values[i] = exact_values[i]-values[i]

error_in_values = np.reshape(error_in_values,(NN+1,1))
rms = np.sqrt(np.dot(np.reshape(error_in_values,(1,NN+1)),error_in_values)/(NN+1))[0][0]

# Clear idx to save RAM
del idx

# Plot first eigenvector
plt.figure(1)
plt.plot(ysub,vecs[:,0])

# Plot second eigenvector
plt.figure(2)
plt.plot(ysub,vecs[:,1])

# Plot the third eigenvector
plt.figure(3)
plt.plot(ysub,vecs[:,2])

# Plot the fourth eigenvector
plt.figure(4)
plt.plot(ysub,vecs[:,3])

# Plot the fifth eigenvector
plt.figure(5)
plt.plot(ysub,vecs[:,4])

# Plot the sixth eigenvector
plt.figure(6)
plt.plot(ysub,vecs[:,5])

# Plot the seventh eigenvector
plt.figure(7)
plt.plot(ysub,vecs[:,6])

# Plot the eighth eigenvector
plt.figure(8)
plt.plot(ysub,vecs[:,7])

# Plot the nineth eigenvector
plt.figure(9)
plt.plot(ysub,vecs[:,8])

# Plot the tenth eigenvector
plt.figure(10)
plt.plot(ysub,vecs[:,9])

# Plot the eleventh eigenvector
plt.figure(11)
plt.plot(ysub,vecs[:,10])

# Plot errors
n=np.linspace(0,NN,NN+1)
plt.figure(12)
plt.plot(n,np.abs(error_in_values))

# Time taken to run script
print("It took", round(time.time()-start, ndigits=2), "seconds for this script to perform the integration.")