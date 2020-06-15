#!/usr/bin/env python3
# Import modules
# Import all Numpy functions
from numpy import *
# Import all SciPy (NumPy can also be used)
# linear algebra functions
from scipy.linalg import *
# Import PyPlot for plotting
from matplotlib import pyplot as plt
# Import Airy function
from scipy.special import airy
# Timer
import time

# Start timer
start = time.time()

# Number of steps
N            = 10000
# The domain we're transforming to is [a,b]
a            = 0
b            = 870
# na=0,1,2,3,...,N
na           = arange(N+1)
# nad is the transpose of na
nad          = reshape(na,(N+1,1))
# t are values of t (the parameterization for the
# extrema grid mentioned in Boyd, 2000)
t            = pi*(1-na/N)
# td is the transpose of t
td           = reshape(t,(N+1,1))
# Chebyshev extrema grid
xa           = cos(t)
# Matrix of T_n(x_m) values
# m is the row index and n is the column index
T            = cos(td*reshape(na,(1,N+1)))
# Matrix of T'_n(x_m) values
# Initiate it
dT           = zeros((N+1,N+1))
# Define the inner rows
dT[1:N,:]    = matmul(matmul(diag(divide(1,sqrt(1-power(xa[1:N],2)))),sin(td[1:N]*na)),diag(na))
del t, td
# Define the endpoints
dT[0,:]      = multiply(power(-1,na+1),power(na,2))
dT[N,:]      = power(na,2)

# Differentiation matrices
D1           = matmul(dT,inv(T))
# Save RAM by clearing dT which is unused after this
del dT
# D2 = D1^2
D2           = matmul(D1,D1)
del D1
# Linearly transformed independent variable to the
# problem
ysub         = (b-a)/2*xa[1:N]+(a+b)/2
# Save even more RAM by clearing na and xa
del na, xa
# Eigenvalues and eigenvectors (values of the eigenfunctions on
# the transformed Chebyshev extrema grid)
values, vecs = eig(-(4/(b-a)**2)*D2[1:N,1:N]+diag(ysub))
y            = vstack([[a], reshape(ysub,(N-1,1)), [b]])
del ysub, D2
# Order eigenvalues and eigenvectors by absolute value of the
# eigenvalues
values       = abs(values)
idx          = values.argsort()[::1]
values       = values[idx]
vecs         = vecs[:,idx]
del idx

# Uses Newton's method to more precisely estimate eigenvalues
def f(xinput):
    x0      = xinput
    xoutput = x0
    # Initial value of Ai(-x) and Ai'(-x)
    ai      = airy(-xoutput)
    Ai      = ai[0]
    Aip     = ai[1]

    # Keep using Newton's until the error is acceptably small
    while abs(Ai/Aip)>1e-12:
        ai      = airy(-xoutput)
        # Ai(-x)
        Ai      = ai[0]
        # Ai'(-x)
        Aip     = ai[1]
        # x_(n+1) = x_n + Ai(-x_n)/Ai'(-x_n)
        xoutput = xoutput+Ai/Aip

    return xoutput

# Includes zeros of the eigenfunctions at the beginning and end
# of the domain
vecs            = vstack([zeros((1,N-1)), vecs, zeros((1,N-1))])

# How many eigenfunctions and eigenvalues we are going to analyse
NN              = 5360
# Initiate vectors
# Exact eigenvalues vector
exact_values    = zeros([NN+1,1])
# Error in eigenvalues vector
error_in_values = zeros([NN+1,1])
# Exact eigenvectors, to be calculated using our analytical
# solution
exact_vecs      = zeros([N+1,NN+1])
# Error in our Chebyshev-approximated eigenfunctions
error_in_vecs   = zeros([N+1,NN+1])
# Root-mean square error in the
# Chebyshev-approximated eigenfunctions
vecs_rms        = zeros([NN+1,1])
for i in range(0,NN+1):
    # Use Newton's method to compute the eigenvalues more precisely
    exact_values[i]    = f(values[i])
    # Analytical eigenfunctions
    exact_vecs[:,i]    = reshape(airy(y-exact_values[i,0])[0],(N+1))
    # Fixing the sign of our eigenfunctions
    vecs[:,i]          = vecs[:,i]*sign(vecs[1,i])/(sign(exact_vecs[1,i]))
    # Our Chebyshev-approximated eigenfunctions may be off from our
    # analytical ones by a constant multiplier.
    vecs[:,i]          = exact_vecs[1,i]/vecs[1,i]*vecs[:,i]
    # Error in our eigenfunctions
    error_in_vecs[:,i] = abs((vecs[:,i]-exact_vecs[:,i])/max(abs(exact_vecs[:,i])))
    # Root mean square of the eigenfunction error
    vecs_rms[i]        = sqrt(dot(reshape(error_in_vecs[:,i],(1,N+1)),error_in_vecs[:,i])/(N+1))[0]

# Compute the relative error in eigenvalues
error_in_values        = abs(divide(exact_values-reshape(values[0:NN+1],(NN+1,1)),exact_values))

# Eigenvalue root mean square error
values_rms             = sqrt(dot(reshape(error_in_values,(1,NN+1)),error_in_values)/(NN+1))[0]

# Print relevant numbers
print("N is:\n ", N)
print("NN is:\n ", NN)
print("b is:\n ", b)
print("Root mean square relative error in the eigenvalues is:\n ", values_rms[0])

# The following cannot be replaced with a loop, I tried to do so on
# 5 Jun 2020 and I received the error shown in ../logs/Error_from_adding_plotting_loop_to_SLEq_with_fixes.py.log
# Plot first eigenvector
plt.figure(1)
plt.plot(y[0:1000],vecs[0:1000,0])

# Plot second eigenvector
plt.figure(2)
plt.plot(y[0:1000],vecs[0:1000,1])

# Plot the third eigenvector
plt.figure(3)
plt.plot(y[0:1000],vecs[0:1000,2])

# Plot the fourth eigenvector
plt.figure(4)
plt.plot(y[0:1000],vecs[0:1000,3])

# Plot the fifth eigenvector
plt.figure(5)
plt.plot(y[0:1000],vecs[0:1000,4])

# Plot the sixth eigenvector
plt.figure(6)
plt.plot(y[0:1000],vecs[0:1000,5])

# Plot the seventh eigenvector
plt.figure(7)
plt.plot(y[0:1000],vecs[0:1000,6])

# Plot the eighth eigenvector
plt.figure(8)
plt.plot(y[0:1000],vecs[0:1000,7])

# Plot the nineth eigenvector
plt.figure(9)
plt.plot(y[0:1000],vecs[0:1000,8])

# Plot the tenth eigenvector
plt.figure(10)
plt.plot(y[0:1000],vecs[0:1000,9])

# Plot the eleventh eigenvector
plt.figure(11)
plt.plot(y[0:1000],vecs[0:1000,10])

# Plot the twelfth eigenvector
plt.figure(12)
plt.plot(y[0:1000],vecs[0:1000,11])

# Plot NN+1th eigenvector
plt.figure(13)
plt.plot(y,vecs[:,NN+1])

# Plot eigenvalue errors
n=linspace(0,NN,NN+1)
plt.figure(14)
plt.plot(n,error_in_values[0:NN+1])

# Plot eigenfunction errors
plt.figure(15)
plt.plot(n,vecs_rms)

# Time taken to run script
print("It took:\n ", round(time.time()-start, ndigits=2), " seconds for this script to perform the computation.")
