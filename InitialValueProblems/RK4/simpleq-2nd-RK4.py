#!/usr/bin/env python3
import time
start = time.time() # start time
# Import the required modules
from numpy import cos, matrix, pi, log10, zeros
from matplotlib import pyplot as plt
from matplotlib import rc
plt.rc('text', usetex=True)

# The ODE we're solving is (note d2 = d^2 in the following equations):
# d2y/dx2 = -g/l cos(y)
# (which is the equation of the Simple Pendulum where y is the angle
# from the x axis; and x is the time in seconds,
# g is the acceleration due to gravity, namely ~9.8 m/s2
# l is the length of the rod)
# this equation can be re-written as:
# d2y/dx2 = f(x, y, dy)
# where f(x, y, dy) is the function defined below
def f(x, y, dy):
    g = 9.785
    l = 1
    return -g/l * cos(y)

# Fourth-order Runge-Kutta (RK4) function
# You can look it up here
# http://mathworld.wolfram.com/Runge-KuttaMethod.html
def RK4(h, x, y, dy):
    K1=h * f(x, y, dy)
    L1=h * dy
    K2=h * f(x + h/2, y + L1/2, dy + K1/2)
    L2=h * (dy + K1/2)
    K3=h * f(x + h/2, y + L2/2, dy + K2/2)
    L3=h * (dy + K2/2)
    K4=h * f(x + h, y + L3, dy + K3)
    L4=h * (dy + K3)
    return matrix([1/6 * (L1 + 2 * L2 + 2 * L3 + L4), 1/6 * (K1 + 2 * K2 + 2 * K3 + K4)])

# integrating on x in [0,3]; x0 is 0; x1 is 3
x0, x1 = 0, 3
# N points the solution is being integrated over.
N = 10000
# h is the step size, that is, the distance between each individual
# integration point
h = (x1-x0)/N
# Relative tolerance in x. Thanks to floating-point arithmetic this
# needs to be used to adjust x1 in the while loop below
rtol = h/2
# y[x=0] = y0
y0 = 0
# dy/dx [x = 0] = dy0
dy0 = 0
# initializing the x variable at x0
x = zeros([N+1, 1])
x[0,0] = x0
# initializing the y variable at y0
y = zeros([N+1, 1])
y[0,0] = y0
# initializing the dy variable at dy0
dy = zeros([N+1, 1])
dy[0,0] = dy0
# miny is the minimum y value. This should be equal to exactly negative pi!
miny = -3.14159
# this is how much the minimum y value differs from its expected value, -pi
err = abs(miny + pi)

# integrate y until x => x1 - rtol
for i in range(1,N+1):
    # RK is a 1x2 vector of the incremental changes to be made to y and dy since
    # the previous x value
    RK = RK4(h, x[i-1,0], y[i-1,0], dy[i-1,0])
    # Add RK correction to y
    y[i,0] = y[i-1,0] + RK[0,0]
    # Add RK correction to dy
    dy[i,0] = dy[i-1,0] + RK[0,1]
    # Add h to x
    x[i,0] = x[i-1,0] + h
    # This is designed to determine miny and the error
    while abs(y[i,0]+pi) < err:
        miny = y[i,0]
        err = abs(y[i,0]+pi)
        logerr = log10(err)

# print x[N], i.e., the final x value calculated in the above while loop
print("x[N] is", x[N,0])
# print y[N], the final y value calculated in the above loop
print("y[N] is", y[N,0])
# print error
print("err is", err)
# error log to base 10
print("log10 of err is", logerr)
print("It took", round(time.time()-start, ndigits=2), "seconds for this script to run.")

# Plot of angle against time
plt.figure(1)
plt.plot(x,y)
plt.xlabel(r'$t$',fontsize=16)
plt.ylabel(r'$\theta \hspace{0.2cm}$',fontsize=16,rotation=0)

# Phase plot of theta dot against theta
plt.figure(2)
plt.plot(y,dy)
plt.xlabel(r'$\theta$',fontsize=16)
plt.ylabel(r'$\frac{d\theta}{dt}\hspace{0.2cm}$',fontsize=20,rotation=0)