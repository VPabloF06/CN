import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def sol(t):
    return 1/3*np.exp(-2*t) + 2/3*np.exp(-t), 2/3*np.exp(-2*t)- 2/3*np.exp(-t)

def odef(t,y):
    return np.array([-1/3*(4*y[0] + y[1]), -1/3*(2*y[0] + 5*y[1])])
#    return np.array([y[0]+y[1], 2*y[0]*y[1]])

#def odef2(t,y):
 #   return np.array([-4*y[0] - 3*y[1], -6*y[0] - 7*y[1]])

A2 = np.array([[-4., -3.], [-6., -7.]])
def odef2(t,y):
    return A2@y


alpha = np.array([1,0]) #initial condition y(0)
finalTime=2

#Euler method(explícit)
nOfSteps=10
h=finalTime/nOfSteps
t=np.arange(0,finalTime+h,h)
y = np.zeros((nOfSteps+1,len(alpha)))
y[0,:] = alpha
for i in np.arange(0,nOfSteps):
    y[i+1,:] = y[i,:] + h*odef(t[i], y[i,:])    
# Plots
plt.plot(t,y[:,0],'*-')
plt.plot(t,y[:,1],'*-')
plt.legend(['y1','y2'])
plt.grid()
plt.title('Euler method')
plt.show()
endError=y[-1,0]
print('Error final time = %0.1e' %(endError))





alpha2 = np.array([1,1]) #initial condition y(0)

finalTime=2

#Euler method(explícit)
nOfSteps=10
h=finalTime/nOfSteps
t=np.arange(0,finalTime+h,h)
y = np.zeros((nOfSteps+1,len(alpha2)))
y[0,:] = alpha2
for i in np.arange(0,nOfSteps):
    y[i+1,:] = y[i,:] + h*odef2(t[i], y[i,:])    
# Plots
plt.plot(t,y[:,0],'*-')
plt.plot(t,y[:,1],'*-')
plt.legend(['y1','y2'])
plt.grid()
plt.title('Euler method')
plt.show()
endError=y[-1,0]
print('Error final time = %0.1e' %(endError))

#Euler enrere
y1 = np.zeros((nOfSteps+1,len(alpha)))
y1[0,:] = alpha2
M = np.eye(len(alpha)) - h*A2
for i in np.arange(0, nOfSteps):
    y1[i+1,:] = np.linalg.solve(M, y1[i,:])

plt.plot(t,y1[:,0],'*-')
plt.plot(t,y1[:,1],'*-')
plt.legend(['y1','y2'])
plt.grid()
plt.title('Euler enrere method')
plt.show()



#Obs: amb odef2 i nOfSteps = 10 no es estable
# IMPORTANT
# Metode euler estable si |lambda*y| < 2 IMPORTANT
# On lambda esta definit per y' = lambda*y

# EULER ENRERE
#Y_i+1 = Y_i + h*f(t_i+1,Y_i+1)



#ODE solution with RK45
# sol = solve_ivp(odef, [0, finalTime], alpha,method='RK45') #,rtol=1.e-8)
# tref=sol.t
# yref=sol.y
# #Comparison plots
# plt.plot(t,y[:,0],'*-')
# plt.plot(tref,yref[0,:],'-o')
# plt.legend(['Euler','RK45'])
# plt.title('Comparison')
# plt.show()