import matplotlib.pyplot as plt
import matplotlib
import matplotlib.animation as animation
import numpy as np

matplotlib.use('qt5Agg')
plt.rcParams.update({'font.size':25})


def gauss_function(x, time_shift, sigma):
    return np.exp(-(x-time_shift)**2/(2*sigma**2))


def barrier(n , v0, thickness):
    "Barrier potential"
    v = np.zeros(n)
    v[int(n/2):int(n/2)+thickness] = v0
    return v


N = 300
T = 8*N
dx = 1
m = 1
h_bar = 1
sigma = 40.0
time_shift = round(N/2) - 5*sigma
k0 = np.pi/10
E = (h_bar**2/2.0/m)*(k0**2+0.5/sigma**2)
X = dx*np.linspace(0, N, N)

V = np.zeros(N) # free particle
V = barrier(N, 0.80, 20)
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

# fig, ax1 = plt.subplots()
# ax1.plot(X, psi_real[1], 'r-', label="real")
# ax1.plot(X, psi_imaginary[1], 'b-', label="imaginary")
# ax1.plot(X, psi_probability, 'k-', label="probability")
# ax2 = ax1.twinx()
# ax2.plot(X, V, 'g-')



idex1 = range(1, N-1)
idex2 = range(2, N)
idex3 = range(0, N-2)

fig, ax = plt.subplots()
ax2 = ax.twinx()
ax2.plot(X, V, 'g-')
x = np.arange(0, 2 * np.pi, 0.01)  # defining 'x'

line, = ax.plot(x, np.sin(x))
ax.set_xlim([0, N])
ax.set_ylim([-1, 1])

def run(t):
    print(t)
    psi_real_present = psi_real[1]
    psi_imaginary_present = psi_imaginary[1]

    psi_imaginary[2, idex1] = psi_imaginary[0, idex1] + \
                              C1 * (psi_real_present[idex2] - 2 * psi_real_present[idex1] +
                                    psi_real_present[idex3])
    psi_imaginary[2] -= C2V * psi_real[1]

    psi_real[2, idex1] = psi_real[0, idex1] - \
                         C1 * (psi_imaginary_present[idex2] - 2 * psi_imaginary_present[idex1] +
                               psi_imaginary_present[idex3])
    psi_real[2] += C2V * psi_imaginary[1]


    psi_real[0] = psi_real_present
    psi_real[1] = psi_real[2]
    psi_imaginary[0] = psi_imaginary_present
    psi_imaginary[1] = psi_imaginary[2]
    line.set_xdata(X)
    line.set_ydata(psi_real[1])
    return line,


# psi_probability = psi_real[2]**2 + psi_imaginary[2] ** 2
# P = dx*psi_probability.sum()
# normalize = np.sqrt(P)
# psi_real = psi_real/normalize
# psi_imaginary = psi_imaginary / normalize
# psi_probability = psi_probability/normalize
# ax1.plot(X, psi_real[2], 'r-', label="real")
# ax1.plot(X, psi_imaginary[2], 'b-', label="imaginary")
# ax1.plot(X, 6*psi_probability, 'k-', label="probability")
# plt.legend()
# plt.grid(True)
# plt.show()




# Init function ti initialize variables
def init():
    plt.xlabel("x")
    plt.ylabel("y")
    return line,  # return the variables that will updated in each frame


# def animate(i):  # 'i' is the number of frames
#     line.set_ydata(i*x)  # update the data
#     line.set_xdata(x)
#     print(i)
#     return line,


ani = animation.FuncAnimation(fig, run, T, init_func=init, interval=0.00000001, blit=False)
plt.show()




plt.show()