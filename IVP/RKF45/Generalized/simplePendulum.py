from matplotlib import pyplot as plt
import matplotlib
from numpy import *
import time
from RKF45 import RKF45

class paramObj:
    def __init__(self, g, l):
        self.g = g
        self.l = l

def simpPen(params, t, vars, dt):
    g, l = params.g, params.l
    theta, thetaDot = vars[0], vars[1]
    thetaDDot = (-g/l)*cos(theta)
    return [dt * thetaDot, dt * thetaDDot]

# Problem parameters
t0 = 0.0
tf = 10
epsilon = 1e-9
theta0, thetaDot0 = 0, 0
dtInitial = 1e-2       
conds = array([[theta0, thetaDot0]])    
params = paramObj(g=9.81, l=1.0)
t, vars = RKF45(simpPen, t0, tf, dtInitial, params, conds, epsilon)

# Plot data
matplotlib.use('svg')
plt.figure(1)
plt.plot(vars[:,0], vars[:,1])
plt.xlabel(r"$\theta$")
plt.ylabel(r"$\dot{\theta}$")
plt.savefig("simplePendulum-phasePlot.svg")
plt.figure(2)
plt.plot(t, vars[:,0], label=r"$\theta$")
plt.plot(t, vars[:,1], label=r"$\dot{\theta}$")
plt.legend();
plt.savefig("simplePendulum-timePlot.svg")