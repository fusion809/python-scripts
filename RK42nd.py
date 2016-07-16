#!/usr/bin/env python
def h(g,l,t, x1, v1):
    from math import cos
    return ( -g/l*cos(x1) )

def f(g,l, t, x1, v1):
    return( [v1, h(g,l,t,x1,v1)] )

def run(data={}):
    xi, vi, ti = data[2], data[3], data[4]
    g=data[0]
    l=data[1]

    while ti[len(ti)-1] <= data[5]:
        xn = xi[-1]
        vn = vi[-1]
        tn = ti[-1]

        K1 = f(g, l, t = tn, x1 = xn, v1 = vn)
        K1 = [dt*K1[i] for i in range(len(K1))]

        K2 = f(g, l, t = tn + 0.5*dt, x1 = xn + 0.5*K1[0], v1 = vn + 0.5*K1[1])
        K2 = [dt*K2[i] for i in range(len(K2))]

        K3 = f(g, l, t = tn + 0.5*dt, x1 = xn + 0.5*K2[0], v1 = vn + 0.5*K2[1])
        K3 = [dt*K3[i] for i in range(len(K3))]

        K4 = f(g, l, t = tn + dt, x1 = xn + K3[0], v1 = vn + K3[1])
        K4 = [dt*K4[i] for i in range(len(K4))]

        xn = xn + (K1[0] + 2*K2[0] + 2*K3[0] + K4[0])/6
        vn = xn + (K1[1] + 2*K2[1] + 2*K3[1] + K4[1])/6

        ti.append(tn+dt)
        xi.append(xn)
        vi.append(vn)

    return(ti, xi, vi)
