#!/usr/bin/env python3
# Rossler Equations SciPy solver
import numpy as np
from scipy import integrate
from mpl_toolkits.mplot3d import axes3d
from math import cos
from matplotlib import pyplot as plt
tmin, tmax = 0, 150
x0, y0, z0 = 0, 1, 0
a, b, c = 0.1, 0.1, 14
N = 1000000
h = (tmax - tmin) / float(N)

def f(x, t):
    return [-x[1]-x[2], x[0]+a*x[1], b+x[2]*(x[0]-c)]

t    = np.arange(tmin, tmax, h)
sol  = integrate.odeint(f, [x0, y0, z0], t)
x    = sol[:,0]
y    = sol[:,1]
z    = sol[:,2]

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
