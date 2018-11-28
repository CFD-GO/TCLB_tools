
from SymbolicCollisions.core.cm_symbols import ux, uy
from SymbolicCollisions.core.sym_col_fun import round_and_simplify
from sympy.matrices import Matrix, diag

class MatrixGenerator:

    def get_cm_coeff_diag_matrix(self, m, n, ex_, ey_):
        N = len(ex_)
        diagonala = [pow((ex_[i] - ux), m) * pow((ey_[i] - uy), n) for i in range(0, N)]
        return diag(*diagonala)


    def __init__(self):
        pass

    def m_maker(self, d2qN, row_maker_fun):
        matrix_dict = {'d2q5': [row_maker_fun(0, 0),
                                row_maker_fun(1, 0),
                                row_maker_fun(0, 1),
                                row_maker_fun(2, 0),
                                row_maker_fun(0, 2)],

                       'd2q9': [row_maker_fun(0, 0),
                                row_maker_fun(1, 0),
                                row_maker_fun(0, 1),
                                row_maker_fun(2, 0),
                                row_maker_fun(0, 2),
                                row_maker_fun(1, 1),
                                row_maker_fun(2, 1),
                                row_maker_fun(1, 2),
                                row_maker_fun(2, 2)]
                       }
        M = matrix_dict[d2qN]
        return M


    def get_raw_moments_matrix(self, ex_, ey_):

        N = len(ex_)

        def get_row(m, n):
            row = [pow((ex_[i]), m) * pow((ey_[i]), n) for i in range(0, N)]
            return row

        m_ = self.m_maker(f'd2q{N}',get_row)

        return Matrix(m_)


    def get_shift_matrix(self, K, ex_, ey_):
        """
        See 'Generalized local equilibrium in the cascaded lattice Boltzmann method' by P. Asinari, 2008
        or Incorporating forcing terms in cascaded lattice Boltzmann approach by method of central moments' by Kannan N. Premnath, Sanjoy Banerjee†, 2009
        :param K: transformation matrix, from orthogonal moments to physical DF
        :param ex_: lattice vector
        :param ey_:
        :return: the shift matrix for passing from the frame at rest to the moving frame
        """

        N = len(ex_)

        def get_row(m, n):
            def get_entry(m, n, column):
                coeff = lambda i, m_, n_: pow((ex_[i] - ux), m_) * pow((ey_[i] - uy), n_)
                entry = sum([K[i, column] * coeff(i, m, n) for i in range(0, N)])
                return round_and_simplify(entry)

            row = [get_entry(m, n, i) for i in range(0, N)]
            return row

        cm_ = self.m_maker(f'd2q{N}', get_row)
        return Matrix(cm_)
