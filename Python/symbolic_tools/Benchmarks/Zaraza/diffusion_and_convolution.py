from sympy import fourier_transform, exp, sqrt, pi, cos, simplify
from sympy.abc import x, k, t, symbols
from sympy import init_printing
init_printing(use_unicode=False, wrap_line=False)

import numpy as np
from matplotlib import pyplot    #and the useful plotting library
import matplotlib.pyplot as plt

nx = 128
domain_length = 64
dx = domain_length / (nx-1)
xspace = np.linspace(0, domain_length, nx)

nt = 5000                # the number of timesteps we want to calculate
nu = 5                  # the value of viscosity
sigma = .2              # sigma is a parameter, we'll learn more about it later
dt = sigma * dx**2 / nu # dt is defined using sigma ... more later!

u_IC = 0*np.ones(nx)                 # numpy function ones()
u_IC[int((nx-1)/4):int(nx/2 + 1)] = 1 # setting u = 2 between 0.5 and 1 as per our I.C.s

pyplot.plot(xspace, u_IC)
pyplot.show()

def calc_diffusion_FD(IC,nx,nt,nu,dt):
    u = IC.copy()
    un = IC.copy() #our placeholder array, un, to advance the solution in time
    beta = nu * dt / dx**2
    for n in range(nt):  #iterate through time
        un = u.copy() ##copy the existing values of u into un
        for i in range(0, nx):
            if i == nx-1:
                u[i] = beta*un[i-1]+ (1-2*beta)*un[i] + beta*un[0] # periodic BC
            else:
                u[i] = beta*un[i-1]+ (1-2*beta)*un[i] + beta*un[i+1]
    return u

u_FD = calc_diffusion_FD(u_IC,nx,nt,nu,dt)
pyplot.plot(xspace, u_FD)
pyplot.show()


def calc_diffusion_single_convolution(IC, x, nt, nu, dt):
    u = IC.copy()

    def get_gaussian(x, alfa, t):
        g = -(x - domain_length / 2.) ** 2
        g /= (4 * alfa * t)
        g = np.exp(g)
        g /= np.sqrt(4 * np.pi * alfa * t)
        g *= domain_length / nx  # normalize --> sum(g)=1
        return g

    time_spot = dt * nt
    fundamental_solution = get_gaussian(x, nu, time_spot)

    u = np.convolve(fundamental_solution, u, 'same')
    # pyplot.plot(x, fundamental_solution, marker='v', linestyle="", markevery=5)
    return u, fundamental_solution




from sympy.matrices import Matrix
import numpy as np
import sympy as sp

class GaussianHillAnal:
    def __init__(self, C0, X0, Sigma2_0, k, U, D):
        """
        :param C0: initial concentration
        :param X0: initial position of the hill's centre = Matrix([x0, y0])
        :param U:  velocity = Matrix([ux, uy])
        :param Sigma2_0: initial width of the Gaussian Hill
        :param k: conductivity
        :param dimenions: number of dimensions
        """
        self.C0 = C0
        self.X0 = X0
        self.U = U
        self.Sigma2_0 = Sigma2_0
        self.k = k
        self.dim = D

    def get_concentration_ND(self, X, t):
        decay = 2.*self.k*t
        L = X - self.X0 - self.U*t
        C = self.C0
        C *= pow(2. * np.pi * self.Sigma2_0, self.dim / 2.)
        C /= pow(2. * np.pi * (self.Sigma2_0 + decay), self.dim / 2.)
        C *= sp.exp(-(L.dot(L)) / (2.*(self.Sigma2_0 + decay)))
        return C



conductivity = nu

nt = 200  #the number of timesteps we want to calculate
nu = 5   #the value of viscosity
sigma = .2 #sigma is a parameter, we'll learn more about it later
dt = sigma * dx**2 / conductivity #dt is defined using sigma ... more later!

time_0    = dt*nt/2     # initial contidion for FD
time_spot = dt*nt       # time to be simulated (by FD and analytically)

X0 = Matrix([domain_length/2.]) # center of the hill
C0 = 1.                 # concentration
variance = 30           # initial variance
reference_level = 0

T_0 = np.zeros(nx)
T_anal = np.zeros(nx)

gha = GaussianHillAnal(C0, X0, variance, conductivity, Matrix([0]), D=1)

for i in range(nx):
    T_0[i] = reference_level + gha.get_concentration_ND(Matrix([xspace[i]]), time_0)
    T_anal[i] = reference_level + gha.get_concentration_ND(Matrix([xspace[i]]), time_spot)

T_FD = calc_diffusion_FD(T_0,nx,nt,nu,dt)
T_single_conv, fs = calc_diffusion_single_convolution(T_0,xspace,nt,nu,dt)

plt.rcParams.update({'font.size': 16})
figure, axis = plt.subplots(1, 1, figsize=(8, 6))
plt.subplots_adjust(hspace=1)
axis.set_title('Diffusion of a Gaussian Hill')
axis.plot(xspace, T_0, label=r'$T_{0}$')
axis.plot(xspace, T_anal, label=r'$T_{anal}$')
axis.plot(xspace, T_FD, label=r'$T_{FD}$', marker='x', linestyle="", markevery=1)
axis.plot(xspace, T_single_conv, label=r'$T_{conv}$', marker='v', linestyle="", markevery=1)
axis.set_xlabel('x')
axis.set_ylabel('Concentration')
axis.legend(loc="upper right")
pyplot.show()

u_single_conv, fs = calc_diffusion_single_convolution(u_IC, xspace, nt, nu, dt)
pyplot.plot(xspace, u_single_conv)
pyplot.show()