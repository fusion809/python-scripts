#!/usr/bin/env python3
### Not functioning ###
import numpy as np
from scipy import linalg as LA
from matplotlib import pyplot as plt

N=10000
b=10
a=0
n=np.arange(0,N+1)
nsub=np.arange(1,N)
n=np.asmatrix(n)
x=np.reshape(np.cos(np.pi*n/N),(N+1,1))
xsub=x[1:N,:]
T=np.cos(np.arccos(x)*n)
na=np.arange(N+1)
xa=np.cos(np.pi*na/N)
xsuba=xa[1:N]
dT=np.zeros((N+1,N+1))
dTsub=np.zeros((N-1,N-1))
Un=np.diag(na)

Ux = np.diag(np.divide(1,np.sqrt(1-np.power(np.cos(np.pi*nsub/N),2))))

dTsub=np.matmul(np.matmul(Ux,np.sin(np.arccos(xsub)*n)),Un)
dTm1=np.multiply(np.power(-1,n+1),np.power(n,2))
dTp1=np.power(n,2)
dT[0,:] = dTm1
dT[N,:] = dTp1
dT[1:N,:] = dTsub
del dTsub, dTm1, Ux, Un, xsub, x, n
D1=np.matmul(dT,LA.inv(T))
del dT
D2=np.matmul(D1,D1)
del D1
#E1=D1[1:N,1:N]
E2=D2[1:N,1:N]
ysub=(b-a)/2*xsuba+(a+b)/2
E2=(4/(b-a)**2)*E2
values, vecs = LA.eig(-E2+np.diag(ysub))
values = np.abs(values)
idx = values.argsort()[::-1]   
values = values[idx]
vecs = vecs[:,idx]

plt.plot(ysub,vecs[:,N-3])