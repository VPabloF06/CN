import numpy as np
import scipy
import matplotlib.pyplot as plt

def f(x):
    return np.exp(-x**2)

def trapeziCompost(f, a, b, m):
    x = np.linspace(a, b, m+1)
    suma = 0
    
    h = (b-a)/m
    for i in range(1, m):
        suma += 2*f(x[i])

    return (h/2)*(suma + f(a) + f(b))


def SimpsonCompost(f, a, b, m):
    x = np.linspace(a, b, 2*m+1)
    suma = 0
    
    h = (b-a)/(2*m)
    for i in range(1, m+1):
        suma += f(x[2*i-2]) + 4*f(x[2*i-1]) + f(x[2*i])

    return (h/3)*suma


print(trapeziCompost(f, -1, 1, 10))
print(SimpsonCompost(f, -1, 1, 10))