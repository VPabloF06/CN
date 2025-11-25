#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 30 12:43:38 2025

@author: victor.pablo
"""

import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt


def f(x):
    x1, x2, x3 = x
    return np.array([
        6*x1 - 2*np.cos(x2*x3) - 1,
        9*x2 + np.sqrt(x1**2 + np.sin(x3) + 1.06) + 0.9,
        60*x3 + 3*np.exp(x1*x2) + 10*np.pi - 3
    ])

def J(x):
    x1, x2, x3 = x
    return np.array([
        [6, 2*x3*np.sin(x2*x3), 2*x2*np.sin(x2*x3)],
        [x1/np.sqrt(x1**2 + np.sin(x3) + 1.06), 9, np.cos(x3)/(2*np.sqrt(x1**2 + np.sin(x3) + 1.06))],
        [3*x2*np.exp(x1*x2), 3*x1*np.exp(x1*x2), 60]
    ])


def grafica_convergencia(errors,niter):
    x = np.linspace(0, len(errors), niter)
    y = np.log10(errors)
    plt.plot(x,y)
    plt.grid(True)
    plt.show()  


def Newton_Raphson(f, J, x0, niter=1e4, tol=1e-6):
    x = np.array(x0, dtype=float)
    fx = f(x)
    Jx = J(x)
    i = 0
    errors = []
   
    while LA.norm(fx) > tol and i < niter:
       
        delta = np.linalg.solve(Jx, fx)
        x1 = x - delta
        r = LA.norm(x1 -x)/LA.norm(x1)
        errors.append(r)
        x = x1
        fx = f(x)
        Jx = J(x)
       
       
        i += 1
   
    grafica_convergencia(errors, i)
   
    return x, i, LA.norm(fx)

# Proves
x0 = [0,0,0]
x1 = [1,1,1]
x2 = [5,5,5]
x3 = [-15,15,-15]

print(Newton_Raphson(f, J, x0))
print(Newton_Raphson(f, J, x1))
print(Newton_Raphson(f, J, x2))
#print(Newton_Raphson(f, J, x3))
