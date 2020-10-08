#!/usr/bin/python3
# Import required Python modules
from numpy import *
from matplotlib import pyplot as plt

# dy/dx = f(x,y)
def f(x,y):
    return 1+y**2

def RK45(f, x0, xf, dxInitial, y0, epsilon):
    """
    Our Runge-Kutta-Fehlberg solver function.

    f:         the RHS of our 1st-order problem.
    x0:        the initial independent variable value for our problem.
    xf:        the final independent variable value for our problem.
    dxInitial: our initial step size.
    y0:        y(x0) - our initial condition.
    epsilon:   our error tolerance.

    returns    [x,y] where x is a 1d array of our independent variable values,
    and y is a 1d array of our solution values.
    """
    i = 0
    x = array([x0])
    y = array([y0])
    dx = dxInitial
    while x[i] < xf:
        # dx should be the smallest out of xf-x[i] and the dx determined
        # the last iteration, as otherwise we won't finish at exactly xf
        dx = min([dx, xf-x[i]])
        
        # dy approximators
        k1 = dx*f(x[i],y[i])
        k2 = dx*f(x[i]+dx/4,y[i]+k1/4)
        k3 = dx*f(x[i]+3*dx/8, y[i]+3*k1/32+9*k2/32)
        k4 = dx*f(x[i]+12*dx/13, y[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197)
        k5 = dx*f(x[i]+dx, y[i]+439*k1/216 - 8*k2 + 3680*k3/513 - 845*k4/4104)
        k6 = dx*f(x[i]+dx/2, y[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40)
        
        # y[i+1] approximations
        y1 = y[i] + 25*k1/216 + 1408*k3/2565+2197*k4/4104-k5/5
        y2 = y[i] + 16*k1/135 + 6656*k3/12825 + 28561*k4/56430 - 9*k5/50 + 2*k6/55

        # A representation of the error in our y approximation
        R = abs(y1-y2)/dx
        
        # What our step size should be multiplied by in order to achieve 
        # an error tolerance of epsilon
        s = 0.84*(epsilon/R)**(1/4)
        
        # If R is less than or equal to epsilon, move onto the next step, 
        # otherwise repeat the iteration with our corrected dx value
        if R <= epsilon:
            x = append(x,[x[i]+dx], axis=0)
            y = append(y,[y1], axis=0)
            i += 1
            dx *= s
        else:
            dx *= s
    
    return [x,y]

# Domain of integration is [x0, xf]
x0 = 0.0
xf = pi/4
# Error tolerance
epsilon = 5e-12
# Initial condition
y0 = 0.0
[x,y] = RK45(f, x0, xf, (xf-x0)/100, y0, epsilon)

# Our exact solution evaluated at the grid used for our approximation
yexact = tan(x)

# Perform a little error analysis
error = abs(yexact-y)
print("The maximum error is", max(error))
rmsError = sqrt(dot(error,error)/len(error))
print("rmsError is", rmsError)

# Plot our solution
plt.rcParams["figure.figsize"] = (25,14)
plt.plot(x,y)