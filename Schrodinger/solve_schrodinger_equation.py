import matplotlib.pyplot as plt
import matplotlib
import numpy as np

from scipy.integrate import cumtrapz

matplotlib.use('qt5Agg')
plt.rcParams.update({'font.size':25})


def y_r(x):
    return A*np.exp(-0.5*(((x-x_c)/s)**2))*np.cos((2*np.pi*(x-x_c))/wavelength)


def y_i(x):
    return A*np.exp(-0.5*(((x-x_c)/s)**2))*np.sin((2*np.pi*(x-x_c))/wavelength)

def U(x):
    if L/8 > 0 and x < L/2:
        return 500
    else:
        return 0

# h_bar = 1
# e = 1
# m = 1

m = 9.10938291e-31
h_bar = 1.054571726e-34
e = 1.602176565e-19

N = 1000
dt = 1e-12
L = 5.00e-9
A = 1
wavelength = L/25
s = L/40
X = np.linspace(0, L, N)
x_c = X[int(N/2)]
dx = X[1] - X[0]
y_real = np.zeros(N)
y_image = np.zeros(N)
y_real[0] = 0
y_image[0] = 0
C1 = (dt*h_bar)/(2*m*dx*dx)
C2 = (e*dt)/h_bar

for i in range(0, N-1):
    y_real[i] = y_real[i] - C1*(y_i(X[i+1]) - 2*y_i(X[i]) + y_i(X[i-1])) + C2*(U(X[i])*y_i(X[i]))
    y_image[i] = y_i(X[i]) + C1 * (y_r(X[i + 1]) - 2 * y_r(X[i]) + y_r(X[i - 1])) - C2 * (U(X[i])*y_r(X[i]))


Y = y_real**2 + y_image**2
A = np.sum(cumtrapz(Y, X, initial=0))
y_real = y_real/np.sqrt(A)
y_image = y_image/np.sqrt(A)
fig, (ax1, ax2) = plt.subplots(ncols=2)
ax1.plot(X, y_real, 'b-')
ax1.plot(X, y_image, 'r-')
plt.grid(True)

psi = y_real**2 + y_image**2
ax2.plot(X[::1], psi[::1], 'k-')
plt.grid(True)
plt.show()