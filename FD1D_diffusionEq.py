# Solution of the 1D parabolic equation u_t=nu*u_xx for x in (a,b), t in (0,T)
# with constant boundary conditions

import numpy as np
from scipy.sparse import diags
import matplotlib.pyplot as plt
# from fweulerbweuler import explicitDirichlet, implicitDirichlet

# Definicio del problema
a = 0
b = 1
finalT = 0.1
nu = 1
BC = 'DD'
ua = 0
ub = 1
# Discretitzacio
nOfIntervals = 10
nOfTimeSteps = 40



def initialCondition(x):
    return 2*np.sin(np.pi*x)+x
    # return 1-2*np.abs(x-0.5)


def buildSystem(a, b, nOfIntervals, nu, BCtype, BCvalue_a, BCvalue_b):
    Dx = (b - a) / nOfIntervals
    if BCtype == 'DD':
        dm = nOfIntervals - 1
    elif BCtype == 'NN':
        dm = nOfIntervals + 1
    A = (nu/Dx**2) * diags([1, -2, 1], [-1, 0, 1], shape=(dm,dm)).toarray()
    F = np.zeros(dm)
    if BCtype == 'DD':
        F[0] = nu * BCvalue_a / Dx**2
        F[-1] = nu * BCvalue_b / Dx ** 2
    elif BCtype == 'NN':
        A[0, 1] = 2 * nu / Dx**2
        A[-1, -2] = 2 * nu / Dx**2
        F[0] = 2 * nu * BCvalue_a / Dx
        F[-1] = 2 * nu * BCvalue_b / Dx
    return A, F


def fwEuler(A, F, U0, T, nOfTimeSteps, BCtype, BCvalue_a, BCvalue_b):
    dm = len(U0)
    U = np.zeros((dm, nOfTimeSteps+1))
    U[:, 0] = U0
    if BCtype == 'DD':
        U[0, :] = BCvalue_a
        U[-1, :] = BCvalue_b
        ind = np.arange(1, dm-1)
    elif BCtype == 'NN':
        ind = np.arange(dm)

    Dt = T / nOfTimeSteps
    for n in range(nOfTimeSteps):
        Un = U[ind, n]
        U[ind, n+1] = Un + Dt*(A@Un + F)
    return U


A, F = buildSystem(a, b, nOfIntervals, nu, 'DD', ua, ub)
x = np.linspace(a, b, nOfIntervals+1)
U0 = initialCondition(x)
U_fw = fwEuler(A, F, U0, finalT, nOfTimeSteps, 'DD', ua, ub)
print(U_fw[4, :])
# U_bw = bwEuler(A, F, U0, finalT, nOfTimeSteps, 'DD', ua, ub)

plt.plot(x, U_fw)
plt.show()

# Q1: Dona un resultat raonable
# Q2: nu=4 sortim de la regió d'estabilitat i el mètode no és estable i amb nu=0.25 no sortim de la
# regió d'estabilitat
# Q3: En aquest cas, continuem a la regió d'estabilitat (1*(10)/100 = 0.1)/Continua sent estable, en augmentar el temps
# la temperatura s'estabilitza/si incrementem la T final, la \delta t també incrementa
# Q4: Necessitem 40*2² passos
# Q5: Si T=.5, estem fent l'interval inicial 5 cops més grans, per tant, cada step és 5 cops més gran
#  => steps = 5*steps, ídem per T=1, T=2 / La solució és estable
# Q6: 
