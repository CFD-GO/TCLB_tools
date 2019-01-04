from SymbolicCollisions.core.cm_symbols import ux, uy, uz, moments_dict
from SymbolicCollisions.core.sym_col_fun import round_and_simplify
from sympy.matrices import Matrix, diag


def get_cm_coeff_diag_matrix(m, n, ex_, ey_):
    N = len(ex_)
    diagonala = [pow((ex_[i] - ux), m) * pow((ey_[i] - uy), n) for i in range(0, N)]
    return diag(*diagonala)


def __matrix_maker(DdQq, row_maker_fun):  # TODO: moments_dict - dependency injection
    M = [row_maker_fun(*row) for row in moments_dict[DdQq]]
    return M


def _check_dimensions(ex, ey, ez):
    qx = len(ex)

    if ez is None:
        ez = Matrix([0 for i in range(0, qx)])
        d = '2'
    else:
        d = '3'

    if not (qx == len(ey) == len(ez)):
        raise Exception('ex, ey, ez dimensions mismatch')

    return d, qx, ez


def get_raw_moments_matrix(ex_, ey_, ez_=None):
    """
    :param ex_: lattice vector
    :param ey_:
    :param ez_:
    :return: transformation matrix from DF to raw moments
    """
    d, q, ez_ = _check_dimensions(ex_, ey_, ez_)

    def get_row(m, n, o):
        row = [pow((ex_[i]), m) * pow((ey_[i]), n) * pow((ez_[i]), o) for i in range(0, q)]
        return row

    m_ = __matrix_maker(f'D{d}Q{q}', get_row)
    return Matrix(m_)


def get_shift_matrix(K, ex_, ey_, ez_=None):
    """
    See 'Generalized local equilibrium in the cascaded lattice Boltzmann method' by P. Asinari, 2008
    or Incorporating forcing terms in cascaded lattice Boltzmann approach by method of central moments' by Kannan N. Premnath, Sanjoy Banerjee†, 2009
    :param K: transformation matrix, from moments to physical DF
    :param ex_: lattice vector
    :param ey_:
    :param ez_:
    :return: the shift matrix for passing from the frame at rest to the moving frame
    """

    d, q, ez_ = _check_dimensions(ex_, ey_, ez_)

    def get_row(m, n, o):
        def get_entry(m, n, o, column):
            coeff = lambda i, m_, n_, o_: pow((ex_[i] - ux), m_) * pow((ey_[i] - uy), n_) * pow((ez_[i] - uz), o_)
            entry = sum([K[i, column] * coeff(i, m, n, o) for i in range(0, q)])
            return round_and_simplify(entry)

        row = [get_entry(m, n, o, i) for i in range(0, q)]
        return row

    cm_ = __matrix_maker(f'D{d}Q{q}', get_row)
    return Matrix(cm_)
