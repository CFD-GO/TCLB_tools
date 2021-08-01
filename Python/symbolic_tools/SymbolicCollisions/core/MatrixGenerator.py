from SymbolicCollisions.core.cm_symbols import ux, uy, uz
from SymbolicCollisions.core.printers import round_and_simplify
from sympy.matrices import Matrix, diag
import numpy as np
import pandas as pd


# def get_m_order_as_in_r(x, y, z):
#     if z is not None:
#         yG, zG, xG = np.meshgrid(x, y, z)  # create the actual grid
#         xG = xG.flatten()  # make the grid 1d
#         yG = yG.flatten()  # same
#         zG = zG.flatten()
#         df = pd.DataFrame({'x': xG, 'y': yG, 'z': zG})  # prepare a dataframe
#     else:
#         xG, yG = np.meshgrid(x, y)  # create the actual grid
#         xG = xG.flatten()  # make the grid 1d
#         yG = yG.flatten()  # same
#         df = pd.DataFrame({'x': xG, 'y': yG})
#     return df.to_numpy()


# def get_e_as_in_r(x, y, z):
#     if z is not None:
#         yG, zG, xG = np.meshgrid(x, y, z)  # create the actual grid
#         xG = xG.flatten()  # make the grid 1d
#         yG = yG.flatten()  # same
#         zG = zG.flatten()
#         ex_D3Q27 = Matrix(xG)
#         ey_D3Q27 = Matrix(yG)
#         ez_D3Q27 = Matrix(zG)
#         e_D3Q27 = ex_D3Q27.col_insert(1, ey_D3Q27)
#         e_D3Q27 = e_D3Q27.col_insert(2, ez_D3Q27)
#         return ex_D3Q27, ey_D3Q27, ez_D3Q27, e_D3Q27
#     else:
#         xG, yG = np.meshgrid(x, y)  # create the actual grid
#         xG = xG.flatten()  # make the grid 1d
#         yG = yG.flatten()  # same
#         ex_D3Q27 = Matrix(xG)
#         ey_D3Q27 = Matrix(yG)
#         e_D3Q27 = ex_D3Q27.col_insert(1, ey_D3Q27)
#         return ex_D3Q27, ey_D3Q27, e_D3Q27

def get_m_order_as_in_r(x, y, z):
    yG, zG, xG = np.meshgrid(x, y, z)  # create the actual grid
    xG = xG.flatten()  # make the grid 1d
    yG = yG.flatten()  # same
    zG = zG.flatten()
    df = pd.DataFrame({'x': xG, 'y': yG, 'z': zG})  # prepare a dataframe
    return df.to_numpy()


def get_e_as_in_r(x, y, z):
    yG, zG, xG = np.meshgrid(x, y, z)  # create the actual grid
    xG = xG.flatten()  # make the grid 1d
    yG = yG.flatten()  # same
    zG = zG.flatten()
    ex_D3Q27 = Matrix(xG)
    ey_D3Q27 = Matrix(yG)
    ez_D3Q27 = Matrix(zG)
    e_D3Q27 = ex_D3Q27.col_insert(1, ey_D3Q27)
    e_D3Q27 = e_D3Q27.col_insert(2, ez_D3Q27)
    return ex_D3Q27, ey_D3Q27, ez_D3Q27, e_D3Q27


def get_reverse_direction_idx(e, index):
    rev_e = -1*e
    q, _ = e.shape
    for i in range(q):
        if rev_e[i, :] == e[index, :]:
            return i


def get_reverse_indices(e):
    q, _ = e.shape
    rev_indices = []
    for i in range(q):
        rev_i = get_reverse_direction_idx(e, i)
        rev_indices.append(rev_i)

    return np.array(rev_indices)


class MatrixGenerator:
    def __init__(self, ex, ey, ez, order_of_moments):
        self.ex = ex
        self.ey = ey
        self.ez = ez
        self.order_of_moments = order_of_moments

    def __matrix_maker(self, row_maker_fun):
        M = [row_maker_fun(*row) for row in self.order_of_moments]
        return M

    def _check_dimensions(self):
        qx = len(self.ex)

        if self.ez is None:
            self.ez = Matrix([0 for i in range(0, qx)])
            d = '2'
        else:
            d = '3'

        if not (qx == len(self.ey) == len(self.ez)):
            raise Exception('ex, ey, ez dimensions mismatch')

        return d, qx

    def get_raw_moments_matrix(self):
        """
        :param ex_: lattice vector
        :param ey_:
        :param ez_:
        :return: transformation matrix from DF to raw moments
        """
        d, q = self._check_dimensions()

        def get_row(m, n, o):
            row = [pow((self.ex[i]), m) * pow((self.ey[i]), n) * pow((self.ez[i]), o) for i in range(0, q)]
            return row

        m_ = self.__matrix_maker(get_row)
        return Matrix(m_)

    def get_shift_matrix(self, K=None):
        """
        See 'Generalized local equilibrium in the cascaded lattice Boltzmann method' by P. Asinari, 2008
        or Incorporating forcing terms in cascaded lattice Boltzmann approach by method of central moments' by Kannan N. Premnath, Sanjoy Banerjee†, 2009
        :param ex_: lattice vector
        :param ey_:
        :param ez_:
        :return: the shift matrix for passing from the frame at rest to the moving frame
        """
        if K is None:
            K = self.get_raw_moments_matrix().inv()  # transformation matrix, from moments to physical DF
        d, q = self._check_dimensions()

        def get_row(m, n, o):
            def get_entry(m, n, o, column):
                coeff = lambda i, m_, n_, o_: pow((self.ex[i] - ux), m_) * pow((self.ey[i] - uy), n_) * pow((self.ez[i] - uz), o_)
                entry = sum([K[i, column] * coeff(i, m, n, o) for i in range(0, q)])
                return round_and_simplify(entry)

            row = [get_entry(m, n, o, i) for i in range(0, q)]
            return row

        cm_ = self.__matrix_maker(get_row)
        return Matrix(cm_)

    def get_raw_x_shift_moments_matrix(self):
        """
        :param ex_: lattice vector
        :param ey_:
        :param ez_:
        :return: transformation matrix from DF to central moments
        """
        d, q = self._check_dimensions()

        def get_row(m, n, o):
            row = [pow((self.ex[i] - ux), m) * pow((self.ey[i] - uy), n) * pow((self.ez[i] - uz), o) for i in range(0, q)]
            row = [round_and_simplify(r) for r in row]
            return row

        m_ = self.__matrix_maker(get_row)
        return Matrix(m_)