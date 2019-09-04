"""
See Advection-Diffusion of a Gaussian Hill Section 8.6.1  p322
from 'Lattice Boltzmann Method Principles and Practise'
"""

from sympy.matrices import Matrix, diag
import numpy as np
import sympy as sp
import os

class GaussianHillAnal2D:
    def __init__(self, C0, X0, Sigma2_0, k, U=Matrix([0, 0])):
        """
        :param C0: initial concentration
        :param X0: initial position of the hill's centre = Matrix([x0, y0])
        :param U:  velocity = Matrix([ux, uy])
        :param Sigma2_0: initial width of the Gaussian Hill
        :param k: conductivity
        """
        self.C0 = C0
        self.X0 = X0
        self.U = U
        self.Sigma2_0 = Sigma2_0
        self.k = k

    def get_concentration(self, X, t):
        Sigma2_D = 2*self.k*t
        L = X - self.X0 - self.U*t
        C = self.C0 * self.Sigma2_0 / (self.Sigma2_0 + Sigma2_D)
        C *= sp.exp(-(L.dot(L)) / (2 * (self.Sigma2_0 + Sigma2_D)))
        return C


def prepare_anal_data_ADE_Gaussion_Hill(gha, ux, time_spot, ySIZE, xSIZE, dump_file_path, shall_recalculate_results=False):

    if os.path.isfile(dump_file_path) and not shall_recalculate_results:
        print(f'{dump_file_path} found, loading results from disc')
        (xx, yy, T_anal) = np.load(dump_file_path)
        return xx, yy, T_anal
    else:
        print(f'recalculating results for: {dump_file_path}')

        x_grid = np.linspace(0, xSIZE, xSIZE, endpoint=False) + 0.5
        y_grid = np.linspace(0, ySIZE, ySIZE, endpoint=False) + 0.5
        xx, yy = np.meshgrid(x_grid, y_grid)

        total_time = 1
        T_anal = np.zeros((ySIZE, xSIZE, total_time))

        for i in range(ySIZE):
            print(f"running i/ySIZE = {i}/{ySIZE}...")
            for j in range(xSIZE):
                T_anal[i][j][0] = gha.get_concentration(Matrix([xx[i][j], yy[i][j]]), int(time_spot))  # lets cheat

        T_anal = T_anal[:, :, 0]  # take time slice

        shift = ux * int(time_spot) % xSIZE
        if float(shift).is_integer():
            shift = int(shift)
        else:
            raise Exception('Choose such ux and time_step that the obtained shift will be an int')
        T_anal = np.roll(T_anal, shift, axis=1)

        dump_folder = os.path.dirname(dump_file_path)
        if not os.path.exists(dump_folder):
            os.makedirs(dump_folder)
        np.save(dump_file_path, (xx, yy, T_anal))
        return xx, yy, T_anal