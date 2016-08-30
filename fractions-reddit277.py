def gcd(a,b):
    while b: a, b = b, a%b
    return a

def f(x,y):
    return x // gcd(x,y), y // gcd(x,y)

print(f(105,25))
