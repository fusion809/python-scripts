#!/usr/bin/python3
# Import required Python modules
from numpy import array, append, abs, multiply

def RKF45(f, t0, tf, dtInitial, params, conds, epsilon):
    i = 0
    t = array([t0])
    vars = conds
    dt = dtInitial

    while (t[i] < tf):
        dt = min([dt, tf-t[i]])
        K1 = f(params, t[i], vars[i,:], dt)
        K2 = f(params, t[i] + dt/4, vars[i,:] + multiply(1/4, K1), dt)
        K3 = f(params, t[i] + 3*dt/8, vars[i,:] + multiply(3/32, K1) + multiply(9/32, K2), dt)
        K4 = f(params, t[i] + 12*dt/13, vars[i,:] + multiply(1932/2197, K1) + multiply(-7200/2197, K2) + multiply(7296/2197, K3), dt)
        K5 = f(params, t[i] + dt, vars[i,:] + multiply(439/216, K1) + multiply(-8, K2) + multiply(3680/513, K3) + multiply(-845/4104, K4), dt)
        K6 = f(params, t[i] + dt/2, vars[i,:] + multiply(-8/27, K1) + multiply(2, K2) + multiply(-3544/2565, K3) + multiply(1859/4104, K4) + multiply(-11/40, K5), dt)

        vars1 = vars[i,:] + multiply(25/216, K1) + multiply(1408/2565, K3) + multiply(2197/4104, K4) + multiply(-1/5, K5)
        vars2 = vars[i,:] + multiply(16/135, K1) + multiply(6656/12825, K3) + multiply(28561/56430, K4) + multiply(-9/50, K5) + multiply(2/55, K6)
        
        Rarr = abs((vars2-vars1)/dt)
        R = max(Rarr)
        s = (epsilon/(2*R))**(1/4)

        if (R <= epsilon):
            t = append(t, [t[i]+dt], axis=0)
            vars = append(vars, [vars1], axis=0)
            i += 1
        
        dt *= s

    return t, vars
