import numpy as np
import matplotlib.pyplot as plt


R = 0.00132 
g = -9.8

def Euler(f, a, b, alpha, nOfSteps):
    h = (b - a) / nOfSteps
    t = np.linspace(a, b, nOfSteps + 1)
    y = np.zeros((nOfSteps + 1, len(alpha)))
    y[0, :] = alpha
    for i in range(nOfSteps):
        val = f(t[i], y[i, :])
        y[i+1, :] = y[i, :] + h * val.flatten()
    return t, y

def RK4(f, a, b, alpha, nOfSteps):
    h = (b - a) / nOfSteps
    t = np.linspace(a, b, nOfSteps + 1)
    y = np.zeros((nOfSteps + 1, len(alpha)))
    y[0, :] = alpha
    for i in range(1, nOfSteps + 1):
        k1 = f(t[i-1], y[i-1,:])
        k2 = f(t[i-1] + h/2, y[i-1,:] + h/2*k1)
        k3 = f(t[i-1] + h/2, y[i-1,:] + h/2*k2)
        k4 = f(t[i-1] + h, y[i-1,:] + h*k3)

        y[i, :] = y[i-1, :] + h/6 * (k1+ 2*k2 + 2*k3 + k4)
    return t, y


def odef(t, x):
    v = np.array([x[2], x[3]])
    mod = np.linalg.norm(v)
    return np.array([x[2], x[3], -R*mod*x[2], -R*mod*x[3] + g])

alpha = np.array([0,0,100*np.cos(np.pi/4),100*np.sin(np.pi/4)])
a = 0
b = 10
nOfSteps = 200


#--- EXERCICI 1--- 

t, y = Euler(odef, a, b, alpha, nOfSteps)
t, y2 = RK4(odef, a, b, alpha, nOfSteps)
print("Posició amb EULER a l'instant t = 10:", y[-1,0], y[-1,1] )
print("Posició amb RK4 a l'instant t = 10:", y2[-1,0], y2[-1,1])



plt.figure()
plt.plot(y[:,0], y[:, 1], '*-', label="Euler")
plt.plot(y2[:,0], y2[:, 1], '*-', label="RK4")
plt.legend()
plt.grid()
plt.title('Entregable 3')
plt.show()
# Resolem EDO per a x

#%%
#--- EXERCICI 2 ---
t, y_E = Euler(odef, a, b, alpha, nOfSteps)
t, y2_E = Euler(odef, a, b, alpha, nOfSteps*2)
# Erorr euler
Ea_e = np.linalg.norm(y_E[-1,0:1] - y2_E[-1,0:1])
Er_e = Ea_e/np.linalg.norm(y2_E[-1, 0:1])

t, y_R = RK4(odef, a, b, alpha, nOfSteps)
t, y2_R = RK4(odef, a, b, alpha, nOfSteps*2)
# Error RK4
Ea_R = np.linalg.norm(y_R[-1,0:1] - y2_R[-1,0:1])
Er_R = Ea_e/np.linalg.norm(y2_R[-1, 0:1])
print("Error EULER amb 200 passos (Absolut i relatiu)", Ea_e, Er_e)
print("Error RK amb 200 passos (Absolut i relatiu)", Ea_R, Er_R)