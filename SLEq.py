import numpy as np
from scipy import linalg as LA

N=100
n=np.arange(0,N+1)
nsub=np.arange(1,N)
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
    print(i)
    Ux[i-1,i-1] = np.divide(1,np.sqrt(1-np.power(np.cos(np.pi*i/N),2)))

dTsub=Ux*np.sin(np.arccos(xsub)*n)*Un
dTm1=np.multiply(np.power(-1,n+1),np.power(n,2))
dTp1=np.power(n,2)
for i in n:
    for j in n:
        if i == 0:
            dT[i,j] = dTm1[j]
        elif i == N:
            dT[i,j] = dTp1[j]
        else:
            dT[i,j]=dTsub[i-1,j]

D1=LA.inv(dT)
