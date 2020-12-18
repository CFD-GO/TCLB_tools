"""
See Advection-Diffusion of a Gaussian Hill Section 8.6.1  p322
from 'Lattice Boltzmann Method Principles and Practise'
"""

from sympy.matrices import Matrix, diag
import numpy as np
import sympy as sp
import os


class GaussianHillAnal:
    def __init__(self, C0, X0, Sigma2_0, k, U, D=2):
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

    def get_concentration_2D(self, X, t):
        Sigma2_D = 2*self.k*t
        L = X - self.X0 - self.U*t
        C = self.C0 * self.Sigma2_0 / (self.Sigma2_0 + Sigma2_D)
        C *= sp.exp(-(L.dot(L)) / (2 * (self.Sigma2_0 + Sigma2_D)))
        return C

    def get_concentration_3D(self, X, t):
        decay = 2*self.k*t
        L = X - self.X0 - self.U*t
        C = self.C0
        C *= pow(2 * np.pi * self.Sigma2_0, self.dim / 2.)
        C /= pow(2 * np.pi * (self.Sigma2_0 + decay), self.dim / 2.)
        C *= sp.exp(-(L.dot(L)) / (2*(self.Sigma2_0 + decay)))
        return C

def prepare_anal_data_ADE_Gaussion_Hill_2D(gha: GaussianHillAnal, ux, time_spot, ySIZE, xSIZE, dump_file_path, shall_recalculate_results=False, reference_level=0):

    if os.path.isfile(dump_file_path) and not shall_recalculate_results:
        print(f'{dump_file_path} found, loading results from disc')
        (xx, yy, T_anal) = np.load(dump_file_path)
        return xx, yy, T_anal
    else:
        print(f'recalculating results for: {dump_file_path}')

        x_grid = np.linspace(0, xSIZE, xSIZE, endpoint=False) + 0.5
        y_grid = np.linspace(0, ySIZE, ySIZE, endpoint=False) + 0.5
        xx, yy = np.meshgrid(x_grid, y_grid)

        T_anal = np.zeros((ySIZE, xSIZE, 1))

        for i in range(ySIZE):
            print(f"running i/ySIZE = {i}/{ySIZE}...")
            for j in range(xSIZE):
                T_anal[i][j][0] = reference_level + gha.get_concentration_2D(Matrix([xx[i][j], yy[i][j]]), int(time_spot))  # lets cheat

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


def prepare_anal_data_ADE_Gaussion_Hill_3D(
        gha: GaussianHillAnal,
        ux, time_spot, ySIZE, xSIZE, zSIZE, dump_file_path,
        shall_recalculate_results=False, reference_level=0):

    if os.path.isfile(dump_file_path) and not shall_recalculate_results:
        print(f'{dump_file_path} found, loading results from disc')
        (xx, yy, zz, T_anal) = np.load(dump_file_path)
        return xx, yy, zz,  T_anal
    else:
        print(f'recalculating results for: {dump_file_path}')

        x_grid = np.linspace(0, xSIZE, xSIZE, endpoint=False) + 0.5
        y_grid = np.linspace(0, ySIZE, ySIZE, endpoint=False) + 0.5
        z_grid = np.linspace(0, zSIZE, zSIZE, endpoint=False) + 0.5
        xx, yy, zz = np.meshgrid(x_grid, y_grid, z_grid)

        T_anal = np.zeros((ySIZE, xSIZE, zSIZE))

        for i in range(ySIZE):
            print(f"running i/ySIZE = {i}/{ySIZE}...")
            for j in range(xSIZE):
                print(f"running j/xSIZE = {j}/{xSIZE}...")
                for k in range(zSIZE):
                    point = Matrix([xx[i][j][k], yy[i][j][k], zz[i][j][k]])
                    value = reference_level + gha.get_concentration_3D(point, int(time_spot))
                    T_anal[i][j][k] = value

        shift = ux * int(time_spot) % xSIZE
        if float(shift).is_integer():
            shift = int(shift)
        else:
            raise Exception('Choose such ux and time_step that the obtained shift will be an int')
        T_anal = np.roll(T_anal, shift, axis=1)

        dump_folder = os.path.dirname(dump_file_path)
        if not os.path.exists(dump_folder):
            os.makedirs(dump_folder)
        np.save(dump_file_path, (xx, yy, zz, T_anal))
        return xx, yy, zz, T_anal