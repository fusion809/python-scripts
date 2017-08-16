#!/usr/bin/env python3
# -*- coding: utf-8 -*-

N = 20       # Number of steps
x = 2        # Initial guess
S = 6        # Number of square root

print("i is 0, x is %s" %x)

for i in range(1,N+1):
    x = 1/2.0 * (x + S/x)
    print("i is %d, x is %s" %(i,x))