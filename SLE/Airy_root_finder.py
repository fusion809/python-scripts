#!/usr/bin/env python3
from scipy.special import airy
from numpy import abs

def f(xinput):
    x0=xinput
    xoutput=x0
    Ai=abs(airy(-xoutput)[0])

    while Ai>1e-12:
        ai=abs(airy(-xoutput))
        Ai=ai[0]
        Aip=ai[1]
        xoutput=xoutput+Ai/Aip

    return Ai, xoutput