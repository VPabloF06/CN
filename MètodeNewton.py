# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import matplotlib.pyplot as plt


def f(x):
    return x**5 - 4*x**4 + 7*x**3 - 21*x**2 + 6*x + 18

def df(x):
    return 5*x**4 - 16*x**3 + 21*x**2 - 42*x + 6

def grafica_convergencia(errors, niter):
    """
    Parameters
    ----------
    errors : Lista amb residus ordenats per iteració
    niter : Nombre d'iteracions'

    Returns
    -------
    Gràfica de convergència

    """
    x = np.linspace(0, len(errors), niter)
    y = np.log10(errors)
    plt.plot(x,y)
    plt.grid(True)
    plt.show()    

def newton_iter(x, niter, tolx = 1e-4, tolf = 1e-4):
    """
    Parameters
    ----------
    x : TYPE
        DESCRIPTION.
    niter : TYPE
        DESCRIPTION.

    Returns
    -------
    x : TYPE
        DESCRIPTION.
    r : TYPE
        DESCRIPTION.
    errors : TYPE
        DESCRIPTION.

    """
    errors = []
    i = 0
    
    
    x1 = 10
    r = 1
    
    while i < niter and r > tolx and abs(f(x)) > tolf:
        
        x1 = x - f(x)/df(x)
        r = abs(abs(x-x1)/x1)
        errors.append(r)
        x = x1
        i += 1
        
    grafica_convergencia(errors, i)
    return x, r, errors, i
    
niter = 10

x0, r0, errors0, i0 = newton_iter(-1, niter)
x1, r1, errors1, i1 = (newton_iter(2, niter))
x2, r2, errors2, i2 = (newton_iter(3, niter))
x3, r3, errors3, i3 = (newton_iter(2.5, niter))

x0 = np.linspace(0, len(errors0), i0)
y0 = np.log10(errors0)
plt.plot(x0,y0, marker = 'o', label="x0")

x1 = np.linspace(0, len(errors1), i1)
y1 = np.log10(errors1)
plt.plot(x1,y1, marker = 'o', label="x1")

x2 = np.linspace(0, len(errors2), i2)
y2 = np.log10(errors2)
plt.plot(x2,y2, marker = 'o', label="x2")

x3 = np.linspace(0, len(errors3), i3)
y3 = np.log10(errors3)
plt.plot(x3,y3, marker = 'o', label="x3")

plt.grid(True)
plt.legend()
plt.show()    



