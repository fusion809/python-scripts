#!/usr/bin/python3
# Import required Python modules
from numpy import array, append, abs, multiply
from matplotlib import pyplot as plt
import matplotlib

# dr/dt = f(t,x,y,z)
# where r = (x,y,z)
def f(t,x,y,z):
    rho = 28
    sigma = 10
    beta = 8.0/3.0
    return [sigma*(y-x), x*(rho-z)-y, x*y-beta*z]

def RK45(f, t0, tf, dtInitial, x0, y0, z0, epsilon):
    """
    Our Runge-Kutta-Fehlberg solver function for a group of three coupled
    ODEs.

    f:         the RHS of our problem.
    t0:        the initial independent variable value for our problem.
    tf:        the final independent variable value for our problem.
    dtInitial: our initial step sixe.
    x0:        x(t0) - our x initial condition.
    y0:        y(t0) - our y initial condition.
    z0:        z(t0) - our z initial condition.
    epsilon:   our error tolerance.

    returns    [t,x,y,z] where t is a 1d array of our independent 
    variable values, x is a 1d array of our x solution values, y is
    a 1d array of our y solution values and z is a 1d array of our z solution values.
    """
    i = 0
    t = array([t0])
    x = array([x0])
    y = array([y0])
    z = array([z0])
    dt = dtInitial
    while t[i] < tf:
        # dt should be the smallest out of tf-t[i] and the dt determined
        # the last iteration, as otherwise we won't finish at exactly tf
        dt = min([dt, tf-t[i]])
        
        # r approximators
        K1 = multiply(dt,f(t[i], x[i], y[i], z[i]))
        k1 = K1[0]
        l1 = K1[1]
        m1 = K1[2]

        K2 = multiply(dt,f(t[i]+dt/4, x[i]+k1/4, y[i]+l1/4, z[i]+m1/4))
        k2 = K2[0]
        l2 = K2[1]
        m2 = K2[2]

        K3 = multiply(dt,f(t[i]+3*dt/8, x[i]+3*k1/32+9*k2/32, y[i]+3*l1/32+9*l2/32, z[i]+3*m1/32+9*m2/32))
        k3 = K3[0]
        l3 = K3[1]
        m3 = K3[2]

        K4 = multiply(dt,f(t[i]+12*dt/13, x[i]+1932*k1/2197-7200*k2/2197+7296*k3/2197, y[i]+1932*l1/2197-7200*l2/2197+7296*l3/2197, z[i]+1932*m1/2197-7200*m2/2197+7296*m3/2197))
        k4 = K4[0]
        l4 = K4[1]
        m4 = K4[2]

        K5 = multiply(dt,f(t[i]+dt, x[i]+439*k1/216-8*k2+3680*k3/513-845*k4/4104, y[i]+439*l1/216-8*l2+3680*l3/513-845*l4/4104, z[i]+439*m1/216-8*m2+3680*m3/513-845*m4/4104))
        k5 = K5[0]
        l5 = K5[1]
        m5 = K5[2]

        K6 = multiply(dt,f(t[i]+dt/2, x[i]-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40, y[i]-8*l1/27+2*l2-3544*l3/2565+1859*l4/4104-11*l5/40, z[i]-8*m1/27+2*m2-3544*m3/2565+1859*m4/4104-11*m5/40))
        k6 = K6[0]
        l6 = K6[1]
        m6 = K6[2]

        x1 = x[i] + 25*k1/216+1408*k3/2565+2197*k4/4104-k5/5
        y1 = y[i] + 25*l1/216+1408*l3/2565+2197*l4/4104-l5/5
        z1 = z[i] + 25*m1/216+1408*m3/2565+2197*m4/4104-m5/5
        x2 = x[i] + 16*k1/135+6656*k3/12825+28561*k4/56430-9*k5/50+2*k6/55
        y2 = y[i] + 16*l1/135+6656*l3/12825+28561*l4/56430-9*l5/50+2*l6/55
        z2 = z[i] + 16*m1/135+6656*m3/12825+28561*m4/56430-9*m5/50+2*m6/55

        # A representation of the error in our x approximation
        Rx = abs(x1-x2)/dt
        Ry = abs(y1-y2)/dt
        Rz = abs(z1-z2)/dt
        
        # What our step sixe should be multiplied by in order to achieve 
        # an error tolerance of epsilon
        sx = 0.84*(epsilon/Rx)**(1/4)
        sy = 0.84*(epsilon/Ry)**(1/4)
        sz = 0.84*(epsilon/Rz)**(1/4)

        # If R is less than or equal to epsilon, move onto the next step, 
        # otherwise repeat the iteration with our corrected dt value
        if ( Rx <= epsilon) & (Ry <= epsilon) & (Rz <= epsilon):
            t = append(t,[t[i]+dt], axis=0)
            x = append(x,[x1], axis=0)
            y = append(y,[y1], axis=0)
            z = append(z,[z1], axis=0)
            i += 1
            dt *= min([sx, sy, sz])
        else:
            dt *= min([sx, sy, sz])

    
    return [t,x,y,z]

# Domain of integration is [t0, tf]
t0 = 0.0
tf = 200
# Error tolerance
epsilon = 1e-8
# Initial condition
x0 = 10
y0 = 10
z0 = 10
[t,x,y,z] = RK45(f, t0, tf, (tf-t0)/1000, x0, y0, z0, epsilon)

# Plot our solution
matplotlib.use('WXAgg')
plt.figure(1)
plt.plot(x,y)
plt.figure(2)
plt.plot(x,z)
plt.figure(3)
plt.plot(y,z)
fig = plt.figure(4)
ax = fig.gca(projection='3d')
plt.rcParams["figure.figsize"] = (25,14)
ax.plot(x,y,z)