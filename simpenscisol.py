# Simple Pendulum SciPy solver
# d^2 theta / dt^2 = -g/l cos(theta)
import numpy as np
from scipy import integrate
from math import cos
from matplotlib import pyplot as plt
a, b = 0, 10
g, l = 9.8, 1
N = 10000000
h = (b-a) / float(N)

def solvr(Y, t):
    return [Y[1], -g/l*cos(Y[0])]

t    = np.arange(a, b, h)
asol = integrate.odeint(solvr, [0, 0], t)
th   = asol[:,0]
dth  = asol[:,1]

# Call the plot function if you want to plot theta against t
# and a phase plot of d{theta}/dt against theta.
def plot():
    plt.figure(1)
    plt.plot(t,th)

    plt.figure(2)
    plt.plot(th,dth)
    plt.show()
