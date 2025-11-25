import numpy as np
import matplotlib.pyplot as plt

def odef(t, y):
    return np.array([-y])

A2 = np.array([[-4., -3.], [-6., -7.]])
def odef2(t, y):
    return A2 @ y

alpha = np.array([1.]) 
a = 0
b = 1
yEx_b = np.exp(-b)
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

# --- Part 1 ---
t, y = Euler(odef, a, b, alpha, nOfSteps)

plt.figure()
plt.plot(t, y, '*-')
plt.grid()
plt.title('Euler method')
plt.show()

endError = np.linalg.norm(yEx_b - y[-1])
print('Error final time = %0.1e' % (endError))

# --- Part 2 ---
errors = []
steps_list = [2**i for i in range(1, 10)]

for n in steps_list:
    t, y = Euler(odef, a, b, alpha, n)
    endError = np.linalg.norm(yEx_b - y[-1])
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

# --- Part 3 ---
alpha2 = np.array([1., 1.])
finalTime = 2
nOfSteps = 20 

# Euler Explícito (Sistema)
t, y = Euler(odef2, 0, finalTime, alpha2, nOfSteps)

plt.figure()
plt.plot(t, y[:, 0], '*-', label='y1')
plt.plot(t, y[:, 1], '*-', label='y2')
plt.legend()
plt.grid()
plt.title('Euler method (Sistema)')
plt.show()

# Euler Implícito (Enrere)
h = finalTime / nOfSteps
t = np.linspace(0, finalTime, nOfSteps + 1)
y1 = np.zeros((nOfSteps + 1, len(alpha2)))
y1[0, :] = alpha2
M = np.eye(len(alpha2)) - h * A2

for i in range(nOfSteps):
    y1[i+1, :] = np.linalg.solve(M, y1[i, :])

plt.figure()
plt.plot(t, y1[:, 0], '*-', label='y1')
plt.plot(t, y1[:, 1], '*-', label='y2')
plt.legend()
plt.grid()
plt.title('Euler enrere method')
plt.show()