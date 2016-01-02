import RK42nd
g=9.8
l=1
N=1000
u1=0
v1=0
t_ini=0
t_end=4
dt=(t_end-t_ini)/N
data=[g, l, u1, v1, t_ini, t_end, dt]
RK42nd.run(data)
