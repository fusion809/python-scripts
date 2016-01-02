import numpy as np
from matplotlib import pyplot as plt
a=0
b=np.pi
g=9.8
l=1
N=1000

def RK4(f):
    return lambda t, y, dt: (
            lambda dy1: (
            lambda dy2: (
            lambda dy3: (
            lambda dy4: (dy1 + 2*dy2 + 2*dy3 + dy4)/6
            )( dt * f( t + dt  , y + dy3   ) )
	    )( dt * f( t + dt/2, y + dy2/2 ) )
	    )( dt * f( t + dt/2, y + dy1/2 ) )
	    )( dt * f( t       , y         ) )

from math import sqrt
dy = RK4(lambda t, y: -y)

t, y, dt = 0., 1., (b-a)/N
i=0
T=np.zeros((N+1,1))
DY=T.copy()
Y=T.copy()
while t < b:
    T[i]=t
    DY[i]=dy(t,y,dt)
    t, y = t + dt, y + dy( t, y, dt )
    Y[i]=y
    i=i+1

plt.figure(1)
plt.plot(T,Y)
plt.show()
