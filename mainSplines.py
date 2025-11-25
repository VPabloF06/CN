import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate

dataX = np.array([0, 1, 3, 4, 5, 7])
dataY = np.array([1, 1.25, 1, 0.5, 0, 0.4])


def calculaDerivada(x,y):
    n = len(dataX) -1
    der = np.zeros(n+1)
    der[0] = (dataY[1] - dataY[0]) / (dataX[1] - dataX[0])
    # for i in range(1,len(dataX)-1):
    #    der[i] = (dataY[i+1] - dataY[i-1]) / (dataX[i+1] - dataX[i-1])
    i = np.arange(1, n)
    der[i] = (dataY[i+1] - dataY[i-1]) / (dataX[i+1] - dataX[i-1])
    der[n] = (dataY[-1] - dataY[-2]) / (dataX[-1] - dataX[-2])
    return der

n = len(dataX)
h = np.zeros(n-1)
t = np.zeros(n-1)

for i in range(len(h)):
    h[i] = dataX[i+1] - dataX[i]
    t[i] = dataY[i+1] - dataY[i]

lamda  = np.zeross(n)
ei = np.zeros(n)
mu = np.zeross(n)

for i in range(1, len(ei)):
    lamda[i] = h[i]/(h[i]+h[i-1])
    mu[i] = h[i-1]/(h[i]+h[i-1])
    ei[i] = 3/(h[i]+h[i-1])*(h[i]*t[i-1]/h[i-1]+h[i-1]*t[i]/h[i])

A = np.zeros(n,n)

#------------------------------------------------------------
# Spline cubic amb aproximació de les derivades
#------------------------------------------------------------
der = calculaDerivada(dataX, dataY)
splineCubicC1 = scipy.interpolate.CubicHermiteSpline(dataX, dataY, der)
print('Coeficients del spline: ', splineCubicC1.c)

x = np.arange(0,7.02,0.05)
y1 = splineCubicC1(x)
dy1 = splineCubicC1(x, nu = 1)   # primera derivada del spline
d2y1 = splineCubicC1(x, nu = 2)  # segona derivada del spline
plt.plot(dataX, dataY,'k*')
plt.plot(x,y1,label = 'S')
plt.plot(x,dy1,label = 'dS')
plt.plot(x,d2y1, label = 'd2S')
plt.title('Spline cubic C1 i les seves derivades')
plt.grid()
plt.legend()
plt.show()

