#!/usr/bin/env python3
import time
start = time.time() # start time
# Import the required modules
from numpy import matrix, zeros
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d

# The ODE we're solving is (note d2 = d^2 in the following equations):
# dx/dt = sigma (y - x)
# dy/dt = x (rho - z) - y
# dz/dt = xy - beta z
# (which is the equation of the Simple Pendulum where y is the angle
# from the x axis; and x is the time in seconds,
# g is the acceleration due to gravity, namely ~9.8 m/s2
# l is the length of the rod)
# this equation can be re-written as:
# dr/dt = f(t, x, y, z)
def f(t, x, y, z):
    sigma = 10
    beta  = 8/3
    rho   = 28
    return matrix([sigma*(y-x), x*(rho-z) - y, x*y - beta*z])

# Fourth-order Runge-Kutta (RK4) function
# You can look it up here
# http://mathworld.wolfram.com/Runge-KuttaMethod.html
def RK4(h, t, x, y, z):
    F1 = f(t,x,y,z)
    K1 = h * F1[0,0]
    L1 = h * F1[0,1]
    M1 = h * F1[0,2]
    F2 = h * f(t + h/2, x + K1/2, y + L1/2, z + M1/2)
    K2 = F2[0,0]
    L2 = F2[0,1]
    M2 = F2[0,2]
    F3 = h * f(t + h/2, x + K2/2, y + L2/2, z + M2/2)
    K3 = F3[0,0]
    L3 = F3[0,1]
    M3 = F3[0,2]
    F4 = h * f(t + h, x + K3, y + L3, z + M3)
    K4 = F4[0,0]
    L4 = F4[0,1]
    M4 = F4[0,2]
    return matrix([1/6 * (K1 + 2 * K2 + 2 * K3 + K4), 1/6 * (L1 + 2 * L2 + 2 * L3 + L4), 1/6 * (M1 + 2 * M2 + 2 * M3 + M4)])

# integrating on t in [0,3]; x0 is 0; x1 is 3
t0, t1 = 0, 100
# N points the solution is being integrated over.
N = 1000000
# h is the step size, that is, the distance between each individual
# integration point
h = (t1-t0)/N
# Relative tolerance in x. Thanks to floating-point arithmetic this
# needs to be used to adjust x1 in the while loop below
rtol = h/2
# x[t=t0] = x0
x0 = 1
# y[t=t0] = y0
y0 = 1
# z[t=t0] = z0
z0 = 1

# Initializing matrixes
t = zeros([N+1, 1])
x = zeros([N+1, 1])
y = zeros([N+1, 1])
z = zeros([N+1, 1])
# initializing the t matrix at t0
t[0] = t0
# initializing the x variable at x0
x[0] = x0
# initializing the y variable at y0
y[0] = y0
# initializing the dy variable at dy0
z[0] = z0

# integrate y until x => x1 - rtol
for i in range(1,N+1):
    # RK is a 1x2 vector of the incremental changes to be made to y and dy since
    # the previous x value
    RK = RK4(h, t[i-1,0], x[i-1,0], y[i-1,0], z[i-1,0])
    # Add RK correction to y
    x[i,0] = x[i-1,0] + RK[0,0]
    # Add RK correction to y
    y[i,0] = y[i-1,0] + RK[0,1]
    # Add RK correction to z
    z[i,0] = z[i-1,0] + RK[0,2]
    # Add h to x
    t[i,0] = t[i-1,0] + h

# print t[N], i.e., the final t value calculated in the above while loop
print("t[N] is", t[N,0])
# print y[N], the final x value calculated in the above loop
print("x[N] is", x[N,0])
print("It took", round(time.time()-start, ndigits=2), "seconds for this script to run.")

# Call the plot function if you want to plot the data
def plot():
    fig = plt.figure(1)
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_wireframe(x, y, z, rstride=10, cstride=10)
    ax.set_xlabel('x(t)')
    ax.set_ylabel('y(t)')
    ax.set_zlabel('z(t)')
    plt.show()

plot()