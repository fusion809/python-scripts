#!/usr/bin/env python
# Simple Pendulum SciPy solver
# Based on http://www.nathantypanski.com/blog/2014-08-23-ode-solver-py.html
# d^2 theta / dt^2 = -g/l cos(theta)
import numpy as np
from scipy import integrate
from math import cos
from matplotlib import pyplot as plt
a, b = 0, 5
g, l = 9.8, 1
N = 1000000000
h = (b-a) / float(N)

def f(x, t):
    return [x[1], -g/l*cos(x[0])]

t    = np.arange(a, b, h)
sol  = integrate.odeint(f, [0, 0], t, rtol=1e-15)
th   = sol[:,0]
dth  = sol[:,1]
Err=[abs(min(th)+np.pi),max(th)]
print(max(Err))
print("The minimum value of dth is",min(dth))
print("The maximum value of dth is",max(dth))
print("The minimum value of th is",min(th))
print("The maximum value of th is",max(th))
# Call the plot function if you want to plot theta against t
# and a phase plot of d{theta}/dt against theta.
def plot():
    plt.figure(1)
    plt.plot(t,sol)

    plt.figure(2)
    plt.plot(th,dth)
    plt.show()

#plot()
