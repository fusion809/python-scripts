# Import the required modules
from numpy import cos, matrix, pi

# The ODE we're solving is (note d2 = d^2 in the following equations):
# d2y/dx2 = -g/l cos(y)
# (which is the equation of the Simple Pendulum where y is the angle
# from the x axis; and x is the time in seconds)
# g is the acceleration due to gravity, namely ~9.8 m/s2
# l is the length of the rod
def f(x, y, dy):
    g = 9.8
    l = 1
    return -g/l * cos(y)

# Fourth-order Runge-Kutta (RK4) function
# You can look it up here 
# http://mathworld.wolfram.com/Runge-KuttaMethod.html
def RK4(h, x, y, dy):
    K1=h * f(x, y, dy)
    L1=h * dy
    K2=h * f(x + h/2, y + L1/2, dy + K1/2)
    L2=h * (dy + K1/2)
    K3=h * f(x + h/2, y + L2/2, dy + K2/2)
    L3=h * (dy + K2/2)
    K4=h * f(x + h, y + L3, dy + K3)
    L4=h * (dy + K3)
    return matrix([1/6 * (L1 + 2 * L2 + 2 * L3 + L4), 1/6 * (K1 + 2 * K2 + 2 * K3 + K4)])

# integrating on x in [0,3]; x0 is 0; x1 is 3
x0, x1 = 0, 3 
# N points the solution is being integrated over. 
N = 1000000
# h is the step size, that is, the distance between each individual
# integration point
h = (x1-x0)/N
# Relative tolerance in x. Thanks to floating-point arithmetic this
# needs to be used to adjust x1 in the while loop below
rtol = h/2
# y[x=0] = y0
y0 = 0
# dy/dx [x = 0] = dy0
dy0 = 0
# initializing the x variable at x0
x = x0
# initializing the y variable at y0
y = y0
# initializing the dy variable at dy0
dy = dy0
# miny is the minimum y value. This should be equal to exactly negative pi!
miny = -3.14159265
# this is how much the minimum y value differs from its expected value, -pi
err = abs(miny + pi)

# integrate y until x => x1 - rtol
while x < x1 - rtol:
    # RK is a 1x2 vector
    RK = RK4(h, x, y, dy)
    # Add RK correction to y
    y = y + RK[0,0]
    # Add RK correction to dy
    dy = dy + RK[0,1]
    # Add h to x
    x = x + h
    # This is designed to determine miny and the error
    while abs(y+pi) < err:
        miny = y
        err = abs(y+pi)

# print x[N]
print("x[N] is", x)
# print y[N]
print("y[N] is", y)
# print error
print("err is", err)