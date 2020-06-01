#!/usr/bin/python3
# This script essentially finds the time taken for a
# simple pendulum that starts with zero velocity at the 
# positive x-axis to reach the negative x-axis (theta=-pi)
# using Simpson's rule
# t is time and y is theta 
# Import the required NumPy functions
from numpy import sin, pi, sqrt, zeros, abs, linspace
# Import time
import time
start = time.time() # start time
# Import SciPy
from scipy.integrate import quad
# Import matplotlib
from matplotlib import pyplot as plt
# Acceleration rate due to gravity in metres per second
# squared.
g=9.8
# Length of the pendulum in metres.
l=1

# Define our integrand
def f(y):
    return 1/sqrt(-2*g/l*sin(y))

# Number of steps
N=1000000
NN=100
tol=linspace(1/(5*N),1/(20*N),NN+1)
error=zeros([NN+1,1])

# Our SciPy approximation to the integral
T=abs(quad(f, 0, -pi))
T=T[0]

for j in range(0,NN+1):
    # Integration interval, [y0, yend]
    y0=-tol[j]
    yend=-pi+tol[j]
    # Step size
    h=(yend-y0)/N
    # Initiate t and y
    t=0
    y=y0

    # Our integration function
    def Simpson(h,y,i,N):
        if i == 0 or i == N:
            return h/3*f(y)
        elif (i % 2) == 0:
            return 2*h/3*f(y)
        else:
            return 4*h/3*f(y)

    # The actual integration
    for i in range(0,N):
        t = t + Simpson(h,y,i,N)
        y = y + h

    # Difference between our SciPy approximation and our 
    # approximation using Simpson's method
    error[j,0] = abs(T + t)

# Time taken to run script
print("It took", round(time.time()-start, ndigits=2), "seconds for this script to perform the integration.")
plt.semilogy(tol,error)