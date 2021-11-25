#!/usr/bin/env python
import numpy as np

# Loop limits
count = 0
itMax = 1e3
diff = 1e-1
tol = 1e-14

# Parameter of problem
a = 3

# Initial guess
x0 = 1.5
x = x0

# Function we're setting to zero
def f(x, params):
    return 2**x - params["a"]

# Function and first derivative
def ffd(x, dx, params):
    fx = f(x, params)
    # 2nd order approximation to derivative
    fdx = (f(x + dx, params) - f(x - dx, params))/(2*dx)
    return fx, fdx

# Apply Newton's until we have a satisfactory approximation
while (np.abs(diff) > tol and count < itMax):
    fx, fdx = ffd(x, 1e-10, {"a": a})
    diff = fx/fdx
    x -= diff
    count += 1

# Error in our approximation, print relevant vars
error = x - np.log(3)/np.log(2)
print("x = {}, error = {}, count = {}".format(x, error, count))