# Lorenz Equations SciPy solver
import numpy as np
from scipy import integrate
from mpl_toolkits.mplot3d import axes3d
from math import cos
from matplotlib import pyplot as plt
a, b = 0, 100
sigma, rho, beta = 10, 28, 8/3
N = 1000000
h = (b-a) / float(N)

def f(x, t):
    return [sigma*(x[1]-x[0]), x[0]*(rho-x[2])-x[1], x[0]*x[1]-beta*x[2]]

t    = np.arange(a, b, h)
asol = integrate.odeint(f, [0, 1, 1], t)
x    = asol[:,0]
y    = asol[:,1]
z    = asol[:,2]

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
