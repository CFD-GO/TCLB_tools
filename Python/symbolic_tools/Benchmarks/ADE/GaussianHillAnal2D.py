"""
See Advection-Diffusion of a Gaussian Hill Section 8.6.1  p322
from 'Lattice Boltzmann Method Principles and Practise'
"""

from sympy.matrices import Matrix, diag
import numpy as np
import sympy as sp


class GaussianHillAnal2D:
    def __init__(self, C0, X0, U, Sigma02, k):
        """
        :param C0: initial concentration
        :param X0: initial position of the hill's centre = Matrix([x0, y0])
        :param U:  velocity = Matrix([ux, uy])
        :param Sigma02: initial width of the Gaussian Hill
        :param k: conductivity
        """
        self.C0 = C0
        self.X0 = X0
        self.U = U
        self.Sigma02 = Sigma02
        self.k = k

    def get_concentration(self, X, t):
        Sigma_D2 = np.sqrt(2*self.k*t)
        L = X - self.X0 - self.U*t

        C = self.C0*self.Sigma02/(self.Sigma02+Sigma_D2)

        C *= sp.exp(-(L.dot(L))/(2*(self.Sigma02 + Sigma_D2)))
        return C
