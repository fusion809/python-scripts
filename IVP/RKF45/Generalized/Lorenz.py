from matplotlib import pyplot as plt
import matplotlib
import numpy as np
import time
from RKF45 import RKF45

class paramObj:
    def __init__(self, sigma, rho, beta):
        self.sigma = sigma
        self.rho = rho
        self.beta = beta

def Lorenz(params, t, vars, dt):
    """Function that, for the Lorenz system, returns dx, dy, dz.

    params:   An object of parameters created using paramObj.
    t:        Time value.
    vars:     NumPy array of x, y and z values.
    dt:       Time step.

    returns:  Array of dx, dy and dz.
    """
    sigma = params.sigma
    rho = params.rho
    beta = params.beta
    x, y, z = vars[0], vars[1], vars[2]
    dx = dt * (sigma*(y-x))
    dy = dt * (x*(rho-z)-y)
    dz = dt * (x*y - beta*z)
    return np.array([dx, dy, dz])

def main():
    # Problem parameters
    t0 = 0.0
    tf = 60.0
    epsilon = 1e-9
    x0, y0, z0 = 1.0, 1.0, 1.0
    dtInitial = 1e-2       
    conds = np.array([[x0, y0, z0]])
    params = paramObj(sigma = 10, rho = 28, beta = 8/3)

    # Solve problem and time it
    start_time = time.perf_counter()
    t, vars = RKF45(Lorenz, t0, tf, dtInitial, params, conds, epsilon)
    end_time = time.perf_counter()
    difference = end_time - start_time
    print("It took: ", difference, " seconds for this command to run")

    # Plot data
    matplotlib.use('TkAgg')
    import matplotlib as mpl
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure(1)
    ax = fig.gca(projection='3d')
    ax.plot(vars[:,0], vars[:,1], vars[:,2])
    #plt.savefig("Lorenz-3D.svg")
    plt.show()
    plt.figure(2)
    plt.plot(t, vars[:,0], label="x")
    plt.plot(t, vars[:,1], label="y")
    plt.plot(t, vars[:,2], label="z")
    plt.legend()
    #plt.savefig("Lorenz-time.svg")
    plt.show()

if __name__ == "__main__":
    main()