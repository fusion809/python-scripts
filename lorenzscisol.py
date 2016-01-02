# Lorenz Equations SciPy solver
import numpy as np
from scipy import integrate
from mpl_toolkits.mplot3d import axes3d
from math import cos
from matplotlib import pyplot as plt
a, b = 0, 200
sigma, rho, beta = 10, 28, 8/3
N = 1000000
h = (b-a) / float(N)

def solvr(Y, t):
    return [sigma*(Y[1]-Y[0]), Y[0]*(rho-Y[2])-Y[1], Y[0]*Y[1]-beta*Y[2]]

t    = np.arange(a, b, h)
asol = integrate.odeint(solvr, [0, 1, 1], t)
x    = asol[:,0]
y    = asol[:,1]
z    = asol[:,2]

# Call the plot function if you want to plot the data
def plot():
    fig = plt.figure(1)
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_wireframe(x, y, z, rstride=10, cstride=10)
    plt.show()

plot()
