from matplotlib import pyplot as plt
import matplotlib
from numpy import *
import time
from RKF45 import RKF45

def f(params, t, vars, dt):
    """
    """
    sigma, rho, beta = params[0], params[1], params[2]
    x, y, z = vars[0], vars[1], vars[2]
    dx = dt * (sigma*(y-x))
    dy = dt * (x*(rho-z)-y)
    dz = dt * (x*y - beta*z)
    return array([dx, dy, dz])

# Problem parameters
t0 = 0.0
tf = 60
epsilon = 1e-9
x0 = 1
y0 = 1
z0 = 1
dtInitial = 0.1               
conds = array([[x0, y0, z0]])
params = array([10, 28, 8/3])

# Solve problem and time it
start_time = time.perf_counter()
t, vars = RKF45(f, t0, tf, dtInitial, params, conds, epsilon)
end_time = time.perf_counter()
difference = end_time - start_time
print("It took: ", difference, " seconds for this command to run")

# Plot data
matplotlib.use('svg')
plt.figure(1)
plt.plot(vars[:,0], vars[:,1], vars[:,2])
plt.savefig("Lorenz-3D.svg")
plt.figure(2)
plt.plot(t, vars[:,0], label="x")
plt.plot(t, vars[:,1], label="y")
plt.plot(t, vars[:,2], label="z")
plt.legend()
plt.savefig("Lorenz-time.svg")