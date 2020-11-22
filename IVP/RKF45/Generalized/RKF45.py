#!/usr/bin/python3
# Import required Python modules
from numpy import array, append, abs, multiply

def RKF45(f, t0, tf, dtInitial, params, conds, epsilon):
    """Solve the specified initial value ordinary differential equation using 
    the Runge-Kutta-Fehlberg 4/5th order method.

    f:         Function that takes arguments params, time value, a row vector 
    of dependent variable values, and a final argument of time step size and 
    returns dt * the time derivatives of the dependent variables of the 
    problem.
    t0:        Start time of the computation.
    tf:        End time of the computation.
    dtInitial: Initial step size.
    params:    An object or vector containing parameter values.
    conds:     Initial conditions of the problem contained within the only row
    of a 2D array.
    epsilon:   Error tolerance.

    returns:   An array of time values and a 2D-array of dependent variable 
    values with each column representing a different dependent variable and 
    each row corresponding to a different time value.
    """
    i = 0
    t = array([t0])
    depVars = conds
    dt = dtInitial

    # Loop over time values
    while (t[i] < tf):
        dt = min([dt, tf-t[i]])

        # Approximators of depVars[i+1, :]
        K1 = f(params, t[i], depVars[i,:], dt)
        K2 = f(params, t[i] + dt/4, depVars[i,:] + multiply(1/4, K1), dt)
        K3 = f(params, t[i] + 3*dt/8, depVars[i,:] + multiply(3/32, K1) + multiply(9/32, K2), dt)
        K4 = f(params, t[i] + 12*dt/13, depVars[i,:] + multiply(1932/2197, K1) + multiply(-7200/2197, K2) + multiply(7296/2197, K3), dt)
        K5 = f(params, t[i] + dt, depVars[i,:] + multiply(439/216, K1) + multiply(-8, K2) + multiply(3680/513, K3) + multiply(-845/4104, K4), dt)
        K6 = f(params, t[i] + dt/2, depVars[i,:] + multiply(-8/27, K1) + multiply(2, K2) + multiply(-3544/2565, K3) + multiply(1859/4104, K4) + multiply(-11/40, K5), dt)

        # 4th order approx for depVars[i+1,:]
        depVars1 = depVars[i,:] + multiply(25/216, K1) + multiply(1408/2565, K3) 
        depVars1 += multiply(2197/4104, K4) + multiply(-1/5, K5)
        # 5th order approx for depVars[i+1,:]
        depVars2 = depVars[i,:] + multiply(16/135, K1) + multiply(6656/12825, K3) 
        depVars2 += multiply(28561/56430, K4) + multiply(-9/50, K5);
        depVars2 += multiply(2/55, K6)
        
        # Rarr is an array of error estimates for depVars[i+1,:]
        Rarr = abs((depVars2-depVars1)/dt)
        R = max(Rarr)
        
        # Correct step size, if error estimate is too big
        s = (epsilon/(2*R))**(1/4)
        if R <= epsilon:
            t = append(t, [t[i]+dt], axis=0)
            depVars = append(depVars, [depVars1], axis=0)
            i += 1
        
        dt *= s

    return t, depVars
