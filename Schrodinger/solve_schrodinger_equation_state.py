import matplotlib.pyplot as plt
import matplotlib
import numpy as np

matplotlib.use('qt5Agg')
plt.rcParams.update({'font.size':25})


def gauss_function(x, time_shift, sigma):
    return np.exp(-(x-time_shift)**2/(2*sigma**2))

N = 1200
T = 5*N
dx = 1
m = 1
h_bar = 1
sigma = 40.0
time_shift = round(N/2) - 5*sigma
k0 = np.pi/20
E = (h_bar**2/2.0/m)*(k0**2+0.5/sigma**2)
X = dx*np.linspace(0, N, N)

V = np.zeros(N) # free particle

tau = h_bar/(2*h_bar**2/(m*dx**2)+max(V))   # critical time

C1 = (tau*h_bar)/(2*m*dx*dx)
C2 = 2*tau/h_bar
C2V = C2*V

# wave function Three states past, present and future
psi_real = np.zeros((3, N))
psi_imaginary = np.zeros((3, N))
psi_probability = np.zeros((3, N))

x_n = range(1, int(N/2))
x = X[x_n]/dx
gauss = gauss_function(x, time_shift, sigma)

cos_x = np.cos(k0*x)
sin_x = np.sin(k0*x)

psi_real[0, x_n] = cos_x*gauss
psi_imaginary[0, x_n] = sin_x * gauss
psi_real[1, x_n] = cos_x*gauss
psi_imaginary[1, x_n] = sin_x*gauss

psi_probability = psi_real[1]**2 + psi_imaginary[1] ** 2
P = dx*psi_probability.sum()
normalize = np.sqrt(P)
psi_real = psi_real/normalize
psi_imaginary = psi_imaginary / normalize
psi_probability = psi_probability/normalize

fig, ax1 = plt.subplots()
ax1.plot(X, psi_real[1], 'r-', label="real")
ax1.plot(X, psi_imaginary[1], 'b-', label="imaginary")
ax1.plot(X, psi_probability, 'k-', label="probability")
plt.legend()
plt.grid(True)
plt.show()
