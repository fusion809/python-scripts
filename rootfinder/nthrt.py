#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 17:15:26 2017

@author: fusion809
"""

N = 20        # Number of steps
x = 18/17     # Initial guess
n = 12        # Root order
S = 2         # Number of square root

print("i is 0, x is %s" %x)

for i in range(1,N+1):
    x = (1/n) * ((n-1)*x + S/(x**(n-1)))
    print("i is %d, x is %s" %(i,x))