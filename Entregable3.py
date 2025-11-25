import numpy as np
import matplotlib.pyplot as plt


R = 0.00132 
g = np.array([0, -9,8])

def Euler(f, a, b, alpha, nOfSteps):
    h = (b - a) / nOfSteps
    t = np.linspace(a, b, nOfSteps + 1)
    y = np.zeros((nOfSteps + 1, len(alpha)))
    y[0, :] = alpha
    for i in range(nOfSteps):
        val = f(t[i], y[i, :])
        y[i+1, :] = y[i, :] + h * val.flatten()
    return t, y


def odef1(x):
    return v

def odef2(v):
    return -R*np.linalg.norm(v)*v + g

alphax = np.array([0,0])
alphav = np.array([0])

#--- EXERCICI 1--- 

# Resolem EDO per a x