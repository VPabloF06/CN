import numpy as np
import matplotlib.pyplot as plt
import scipy


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
#nofSteps = 200
nOfSteps = 200


#--- EXERCICI 1--- 
print("EXERCICI 1")

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
print()
print("EXERCICI 2")
t, y_E = Euler(odef, a, b, alpha, nOfSteps)
t, y2_E = Euler(odef, a, b, alpha, nOfSteps*2)
# Erorr euler
Ea_e = np.linalg.norm(y_E[-1,0:2] - y2_E[-1,0:2])
Er_e = Ea_e/np.linalg.norm(y2_E[-1, 0:2])

t, y_R = RK4(odef, a, b, alpha, nOfSteps)
t, y2_R = RK4(odef, a, b, alpha, nOfSteps*2)
# Error RK4
Ea_R = np.linalg.norm(y_R[-1,0:2] - y2_R[-1,0:2])
Er_R = Ea_R/np.linalg.norm(y2_R[-1, 0:2])
print("Error EULER amb 200 passos (Absolut i relatiu)", Ea_e, Er_e)
print("Error RK amb 200 passos (Absolut i relatiu)", Ea_R, Er_R)

#%%
#--- EXERCICI 3 ---
print()
print("EXERCICI 3")
print("El mètode d'euler usa una avaluació per pas i RK4 4, és a dir, per a cada pas de R4 afegim 4 per al de Euler")


errors_E = []
errors_R = []

for i in range(10, 200, 10):
    #Euler
    t, y_E = Euler(odef, a, b, alpha, i*4)
    t, y2_E = Euler(odef, a, b, alpha, i*8)
    # Erorr euler
    Ea_e = np.linalg.norm(y_E[-1,0:2] - y2_E[-1,0:2])
    Er_e = Ea_e/np.linalg.norm(y2_E[-1, 0:2])
    #RK4
    t, y_R = RK4(odef, a, b, alpha, i)
    t, y2_R = RK4(odef, a, b, alpha, i*2)
    # Error RK4
    Ea_R = np.linalg.norm(y_R[-1,0:2] - y2_R[-1,0:2])
    Er_R = Ea_R/np.linalg.norm(y2_R[-1, 0:2])
    
    errors_E.append(Er_e)
    errors_R.append(Er_R)

xx = np.linspace(10, 200, 19)
xxx = np.log10(xx)
yyy_E = np.log10(errors_E)
yyy_R = np.log10(errors_R)
plt.title("Gràfica de convergència")
plt.plot(xxx, yyy_E, label="Euler amb 4*N passos")
plt.plot(xxx, yyy_R, label="RK4 amb N passos")
plt.grid()
plt.legend()
plt.show()

#%%
#--- EXERCICI 4 ---

print()
print("EXERCICI 4")

sol = scipy.integrate.solve_ivp(odef, (a, b), alpha, rtol=1e-9)
sol2 = scipy.integrate.solve_ivp(odef, (a,b), alpha, rtol=1e-14)

#Calculem l'error
Ea_sol = np.linalg.norm(sol.y[-1,0:2] - sol2.y[-1,0:2])
Er_sol = Ea_sol/np.linalg.norm(sol2.y[-1, 0:2])
print("Error EULER amb 1e-9 rtol passos (Absolut i relatiu)", Ea_sol, Er_sol)
print("El nombre d'avaluacions necessitat ha sigut:", sol.nfev)

    