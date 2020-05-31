#!/usr/bin/python3
# This script essentially finds the time taken for a
# simple pendulum that starts with zero velocity at the 
# positive x-axis to reach the negative y-axis (theta=-pi)
# using Simpson's rule
# t is time and y is theta 
# Import the required NumPy functions
from numpy import sin, pi, sqrt, zeros, abs
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
N=100000000
# How close to theta=0 we start our integration and how
# close to theta=-pi we end our intgration.
# At theta=0,-pi our function becomes undefined, hence
# why we cannot start at exact theta=0 or theta=-pi.
tol=1/(18.823*N)
# Integration interval, [y0, yend]
y0=-tol
yend=-pi+tol
# Step size
h=(yend-y0)/N
# Initiate the t and y vectors
t=zeros([N+1,1])
y=zeros([N+1,1])
y[0,0]=y0

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
    t[i+1,0] = t[i,0] + Simpson(h, y[i,0],i,N)
    y[i+1,0] = y[i,0] + h

# Number of steps
print("N is", N)

# Time taken to run script
print("It took", round(time.time()-start, ndigits=2), "seconds for this script to perform the integration.")

# Our SciPy approximation to the integral
T=abs(quad(f, 0, -pi))
T=T[0]

# Difference between our SciPy approximation and our 
# approximation using Simpson's method
error = abs(T + t[N,0])
print("The difference between what SciPy's quad computes and what we computed is", error)

# Plot theta against t
plt.plot(y,t)