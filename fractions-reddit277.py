#!/usr/bin/env python3
# https://www.reddit.com/r/dailyprogrammer/comments/4uhqdb/20160725_challenge_277_easy_simplifying_fractions/
def gcd(a,b):
    while b: a, b = b, a%b
    return a

def f(x,y):
    return x // gcd(x,y), y // gcd(x,y)

print(f(105,25))
