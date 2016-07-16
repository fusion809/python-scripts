#!/usr/bin/env python
# Lorenz Equations SciPy solver
import numpy as np
from scipy import integrate
from mpl_toolkits.mplot3d import axes3d
from math import cos
from matplotlib import pyplot as plt
tmin, tmax = 0, 200
x0, y0, z0 = 0, 1, 0
sigma, rho, beta = 10, 28, 8/3
N = 1000000
h = (tmax - tmin) / float(N)

def f(X, t):
    x, y, z = X[0], X[1], X[2]
    return [sigma*(y-x), x*(rho-z)-y, x*y-beta*z]

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
