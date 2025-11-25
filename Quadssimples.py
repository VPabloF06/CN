import numpy as np
import scipy
import matplotlib.pyplot as plt
from quadraturaGauss import quadraturaGauss

def f(x):
    return 1/(1+16*x**2)

def g(x):
    return np.exp(-x**2)

def Newton_Cotes(x,a,b):
    n = len(x)
    A = np.zeros((n,n))
    c = np.zeros(n)
    
    for i in range(n):
        A[i] = [y**i for y in x]
        c[i] = 1/(i+1)*(b**(i+1)-a**(i+1))
        
    
    return np.linalg.solve(A,c)

def Quadratura_simple(y,w, f):
    n = len(y)
    suma = 0
    for i in range(n):
        suma += w[i]*y[i]
    r = abs(suma -  scipy.integrate.quad(f, -1, 1)[0])    
    return suma, r
    




a = -1
b = 1

errors1 = []
errors2 = []
for i in range(1, 14):
    "f(x)"
    x1 = np.linspace(a, b, i)
    y1 = f(x1)
    
    w1 = Newton_Cotes(x1, a, b)
    I1, r1 = Quadratura_simple(y1, w1, f)
    errors1.append(r1)
    
    z2, w2 = quadraturaGauss(i)
    y2 = f(z2)
    I2, r2 = Quadratura_simple(y2, w2, f)
    errors2.append(r2)
    



plt.title("Gràfica error")
xf = np.linspace(1, 14, 13)
yf1 = np.log10(errors1)
plt.plot(xf,yf1, marker = "o", label="Newton Cotes")
xf = np.linspace(1, 14, 13)
yf2 = np.log10(errors2)
plt.plot(xf,yf2, marker = "o", label = "Gauss")
plt.grid(True)
plt.show()  


errors1 = []
errors2 = []
for i in range(1, 14):

    "g(x)"
    x1 = np.linspace(a, b, i)
    y1 = g(x1)
    
    w1 = Newton_Cotes(x1, a, b)
    I1, r1 = Quadratura_simple(y1, w1, g)
    errors1.append(r1)
    
    z2, w2 = quadraturaGauss(i)
    y2 = g(z2)
    I2, r2 = Quadratura_simple(y2, w2, g)
    errors2.append(r2)
    
plt.title("Gràfica error")
xf = np.linspace(1, 14, 13)
yf1 = np.log10(errors1)
plt.plot(xf,yf1, marker = "o", label="Newton Cotes")
xf = np.linspace(1, 14, 13)
yf2 = np.log10(errors2)
plt.plot(xf,yf2, marker = "o", label = "Gauss")
plt.grid(True)
plt.show()  