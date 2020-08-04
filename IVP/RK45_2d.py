#!/usr/bin/python3
# Import required Python modules
from numpy import *
from matplotlib import pyplot as plt
from scipy.special import ellipk

# d2y/dx2 = f(x,y,dy)
def f(x,y,dy):
    return [dy, -9.8*cos(y)]

def RK45(f, x0, xf, dxInitial, y0, dy0, epsilon):
    """
    Our 2nd order Runge-Kutta-Fehlberg solver function.

    f:         the RHS of our 1st-order problem.
    x0:        the initial independent variable value for our problem.
    xf:        the final independent variable value for our problem.
    dxInitial: our initial step size.
    y0:        y(x0) - our y initial condition.
    dy0:       dy/dx(x0) - our dy initial condition.
    epsilon:   our error tolerance.

    returns    [x,y,dy] where x is a 1d array of our independent 
    variable values, y is a 1d array of our y solution values and dy is
    a 1d array of our dy/dx solution values.
    """
    i = 0
    x = array([x0])
    y = array([y0])
    dy = array([dy0])
    dx = dxInitial
    while x[i] < xf:
        # dx should be the smallest out of xf-x[i] and the dx determined
        # the last iteration, as otherwise we won't finish at exactly xf
        dx = min([dx, xf-x[i]])
        
        # y, dy approximators
        K1 = multiply(dx,f(x[i],y[i],dy[i]))
        k1 = K1[0]
        l1 = K1[1]
        K2 = multiply(dx,f(x[i]+dx/4,y[i]+k1/4,dy[i]+l1/4))
        k2 = K2[0]
        l2 = K2[1]
        K3 = multiply(dx,f(x[i]+3*dx/8, y[i]+3*k1/32+9*k2/32, dy[i]+3*l1/32+9*l2/32))
        k3 = K3[0]
        l3 = K3[1]
        K4 = multiply(dx,f(x[i]+12*dx/13, y[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, dy[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197))
        k4 = K4[0]
        l4 = K4[1]
        K5 = multiply(dx,f(x[i]+dx, y[i]+439*k1/216 - 8*k2 + 3680*k3/513 - 845*k4/4104, dy[i]+439*l1/216 - 8*l2 + 3680*l3/513 - 845*l4/4104))
        k5 = K5[0]
        l5 = K5[1]
        K6 = multiply(dx,f(x[i]+dx/2, y[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, dy[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40))
        k6 = K6[0]
        l6 = K6[1]

        # Solution approximations
        y1 = y[i] + 25*k1/216 + 1408*k3/2565+2197*k4/4104-k5/5
        y2 = y[i] + 16*k1/135 + 6656*k3/12825 + 28561*k4/56430 - 9*k5/50 + 2*k6/55
        dy1 = dy[i] + 25*l1/216 + 1408*l3/2565+2197*l4/4104-l5/5
        dy2 = dy[i] + 16*l1/135 + 6656*l3/12825 + 28561*l4/56430 - 9*l5/50 + 2*l6/55

        # A representation of the error in our y approximation
        Ry = abs(y1-y2)/dx
        Rdy = abs(dy1-dy2)/dx

        # What our step size should be multiplied by in order to achieve 
        # an error tolerance of epsilon
        sy = 0.84*(epsilon/Ry)**(1/4)
        sdy = 0.84*(epsilon/Rdy)**(1/4)
        
        # If R is less than or equal to epsilon, move onto the next step, 
        # otherwise repeat the iteration with our corrected dx value
        if ( Ry <= epsilon) & (Rdy <= epsilon ):
            x = append(x,[x[i]+dx], axis=0)
            y = append(y,[y1], axis=0)
            dy = append(dy,[dy1], axis=0)
            i += 1
            dx *= min([sy, sdy])
        else:
            dx *= min([sy, sdy])
    
    return [x,y,dy]

# Domain of integration is [x0, xf]
x0 = 0.0
xf = 2*ellipk(1/2)/sqrt(2.45)
# Error tolerance
epsilon = 4e-12
# Initial condition
y0 = 0.0
dy0 = 0.0
[x,y,dy] = RK45(f, x0, xf, (xf-x0)/100, y0, dy0, epsilon)

# Plot our solution
plt.rcParams["figure.figsize"] = (25,14)
plt.figure(1)
plt.plot(x,y)
plt.figure(2)
plt.plot(y,dy)