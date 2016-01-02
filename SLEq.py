import numpy as np
from scipy import linalg as LA

N=100
n=np.arange(0,N+1)
nsub=np.arange(0,N-1)
n=np.asmatrix(n)
x=np.cos(np.pi*n/N)
x=x.T
xsub=x[1:N]
T=np.cos(np.arccos(x) * n)
na=np.asarray(n).T
Un=np.zeros((N+1,N+1))
Ux=np.zeros((N-1,N-1))
for i in na:
    Un[i,i] = i

for i in nsub:
    Ux[i,i] = np.divide(1,np.sqrt(1-np.power(x[i+1],2)))

Usub=Ux*np.sin(np.arccos(xsub)*n)*Un
Um1=np.multiply(np.power(-1,n+1),np.power(n,2))
Up1=np.power(n,2)
U=[Um1, Usub, Up1]
