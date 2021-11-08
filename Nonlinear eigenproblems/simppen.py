#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 14 21:43:50 2021

@author: fusion809
"""

import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as la

N = 1000
t = np.linspace(0, 2*np.pi, num=N+1)
x = -np.cos(t);