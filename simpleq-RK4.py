from numpy import cos

def f(x, y):
    return x * cos(y)

def RK4(h, x, y):
    K1=h*f(x, y)
    K2=h*f(x+h/2, y+K1/2)
    K3=h*f(x+h/2, y+K2/2)
    K4=h*f(x+h, y+K3)
    return 1/6*(K1+2*K2+2*K3+K4)

x0=0
x1=100
rtol=1e-4
y0=0
N=10000
h=(x1-x0)/N
x=x0
y=y0

while x < x1 - rtol:
    y = y + RK4(h, x, y)
    x = x + h


print(x)
print(y)
