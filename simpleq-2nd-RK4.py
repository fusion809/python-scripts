# Import the required modules
from numpy import cos
from numpy import matrix
from numpy import pi

def f(x, y, dy):
    return - 9.8 * cos(y)

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

x0 = 0
x1 = 3
rtol = 1e-4
y0 = 0
dy0 = 0
N = 1000000
h = (x1-x0)/N
x = x0
y = y0
dy = dy0
miny = -3.14159265
err = abs(miny + pi)

while x < x1 - rtol:
    RK = RK4(h, x, y, dy)
    y = y + RK[0,0]
    dy = dy + RK[0,1]
    x = x + h
    while abs(y+pi) < err:
        err = abs(y+pi)

print("x[N] is", x)
print("y[N] is", y)
print("err is", err)