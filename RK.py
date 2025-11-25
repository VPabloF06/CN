import numpy as np
import matplotlib.pyplot as plt

def odef(t, y):
    return np.array([y[1], 2*y[1]-2*y[0]])

def g(x):
    return np.exp(x)*np.sin(x)


alpha = np.array([0., 1.]) 
a = 0
b = 1
yEx_b = np.exp(b)*np.sin(b)
nOfSteps = 20

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
    for i in range(1, nOfSteps):
        k1 = f(t[i-1], y[i,:])
        k2 = f(t[i-1] + h/2, y[i-1,:] + h/2*k1)
        k2 = f(t[i-1] + h/2, y[i-1,:] + h/2*k2)
        k2 = f(t[i-1] + h, y[i-1,:] + h*k3)

        y[i+1, :] = y[i, :] + h/6 * (k1+ 2*k2 + 2*k3 + k4)
    return t, y

# --- Part 1 ---
t, y = Euler(odef, a, b, alpha, nOfSteps)
t, y2 = RK4(odef, a, b, alpha, nOfSteps)
plt.figure()
plt.plot(t, y, '*-')
plt.plot(t, g(t), label="f(x)")
plt.grid()
plt.title('Euler method')
plt.show()

endError = np.linalg.norm(yEx_b - y[-1,0])
print('Error final time = %0.1e' % (endError))

# --- Part 2 ---
errors = []
steps_list = [2**i for i in range(1, 10)]

for n in steps_list:
    t, y = Euler(odef, a, b, alpha, n)
    endError = np.linalg.norm(yEx_b - y[-1,0])
    errors.append(endError)
    
xx = np.log2(np.array(steps_list))
yy = np.log2(errors)

plt.figure()
plt.plot(xx, yy, '-o')
plt.title('Euler method convergencia')
plt.xlabel('log2(Steps)')
plt.ylabel('log2(Error)')
plt.grid()
plt.show()


