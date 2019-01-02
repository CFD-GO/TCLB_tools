from SymbolicCollisions.core.cm_symbols import ux, uy
from SymbolicCollisions.core.sym_col_fun import round_and_simplify
from sympy.matrices import Matrix, diag


def get_cm_coeff_diag_matrix(m, n, ex_, ey_):
    N = len(ex_)
    diagonala = [pow((ex_[i] - ux), m) * pow((ey_[i] - uy), n) for i in range(0, N)]
    return diag(*diagonala)


def _matrix_maker(DdQq, row_maker_fun):
    matrix_dict = {
        # order of 2D (central) moments as in
        # `Modeling incompressible thermal flows using a central-moments-based lattice Boltzmann method`
        # by Linlin Fei, Kai Hong Luo Chuandong Lin , Qing Li. 2017
        'D2Q5': [(0, 0, 0),
                 (1, 0, 0),
                 (0, 1, 0),
                 (2, 0, 0),
                 (0, 2, 0)],

        'D2Q9': [(0, 0, 0),
                 (1, 0, 0),
                 (0, 1, 0),
                 (2, 0, 0),
                 (0, 2, 0),
                 (1, 1, 0),
                 (2, 1, 0),
                 (1, 2, 0),
                 (2, 2, 0)],

        # order of 3D (central) moments as in
        # `Three-dimensional cascaded lattice Boltzmann method:
        # Improved implementation and consistent forcing scheme`
        # by Linlin Fei, Kai H.  Luo,  Qing Li. 2018
        'D3Q7': [(0, 0, 0),
                 (1, 0, 0),
                 (0, 1, 0),
                 (0, 0, 1),
                 (2, 0, 0),
                 (0, 2, 0),
                 (0, 0, 2)],

        'D3Q15': [(0, 0, 0),
                  (1, 0, 0),
                  (0, 1, 0),
                  (0, 0, 1),
                  (2, 0, 0),
                  (0, 2, 0),
                  (0, 0, 2),
                  (1, 1, 1),  # skipped in      D3Q19, D3Q7
                  (2, 1, 1),  # skipped in      D3Q19, D3Q7
                  (1, 2, 1),  # skipped in      D3Q19, D3Q7
                  (1, 1, 2),  # skipped in      D3Q19, D3Q7
                  (1, 2, 2),  # skipped in      D3Q19, D3Q7
                  (2, 1, 2),  # skipped in      D3Q19, D3Q7
                  (2, 2, 1),  # skipped in      D3Q19, D3Q7
                  (2, 2, 2),  # skipped in      D3Q19, D3Q7
                  ],

        'D3Q19': [(0, 0, 0),
                  (1, 0, 0),
                  (0, 1, 0),
                  (0, 0, 1),
                  (1, 1, 0),  # skipped in D3Q15, D3Q7
                  (1, 0, 1),  # skipped in D3Q15, D3Q7
                  (0, 1, 1),  # skipped in D3Q15, D3Q7
                  (2, 0, 0),
                  (0, 2, 0),
                  (0, 0, 2),
                  (1, 2, 0),  # skipped in D3Q15, D3Q7
                  (1, 0, 2),  # skipped in D3Q15, D3Q7
                  (2, 1, 0),  # skipped in D3Q15, D3Q7
                  (2, 0, 1),  # skipped in D3Q15, D3Q7
                  (0, 1, 2),  # skipped in D3Q15, D3Q7
                  (0, 2, 1),  # skipped in D3Q15, D3Q7
                  (2, 2, 0),  # skipped in D3Q15, D3Q7
                  (2, 0, 2),  # skipped in D3Q15, D3Q7
                  (0, 2, 2),  # skipped in D3Q15, D3Q7
                  ],

        'D3Q27': [(0, 0, 0),
                  (1, 0, 0),
                  (0, 1, 0),
                  (0, 0, 1),
                  (1, 1, 0),  # skipped in D3Q15, D3Q7
                  (1, 0, 1),  # skipped in D3Q15, D3Q7
                  (0, 1, 1),  # skipped in D3Q15, D3Q7
                  (2, 0, 0),
                  (0, 2, 0),
                  (0, 0, 2),
                  (1, 2, 0),  # skipped in D3Q15, D3Q7
                  (1, 0, 2),  # skipped in D3Q15, D3Q7
                  (2, 1, 0),  # skipped in D3Q15, D3Q7
                  (2, 0, 1),  # skipped in D3Q15, D3Q7
                  (0, 1, 2),  # skipped in D3Q15, D3Q7
                  (0, 2, 1),  # skipped in D3Q15, D3Q7
                  (1, 1, 1),  # skipped in      D3Q19, D3Q7
                  (2, 2, 0),  # skipped in D3Q15, D3Q7
                  (2, 0, 2),  # skipped in D3Q15, D3Q7
                  (0, 2, 2),  # skipped in D3Q15, D3Q7
                  (2, 1, 1),  # skipped in      D3Q19, D3Q7
                  (1, 2, 1),  # skipped in      D3Q19, D3Q7
                  (1, 1, 2),  # skipped in      D3Q19, D3Q7
                  (1, 2, 2),  # skipped in      D3Q19, D3Q7
                  (2, 1, 2),  # skipped in      D3Q19, D3Q7
                  (2, 2, 1),  # skipped in      D3Q19, D3Q7
                  (2, 2, 2),  # skipped in      D3Q19, D3Q7
                  ],
    }

    M = [row_maker_fun(*row) for row in matrix_dict[DdQq]]
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

    m_ = _matrix_maker(f'D{d}Q{q}', get_row)
    return Matrix(m_)


def get_shift_matrix(K, ex_, ey_, ez_=None):
    """
    See 'Generalized local equilibrium in the cascaded lattice Boltzmann method' by P. Asinari, 2008
    or Incorporating forcing terms in cascaded lattice Boltzmann approach by method of central moments' by Kannan N. Premnath, Sanjoy Banerjee†, 2009
    :param K: transformation matrix, from orthogonal moments to physical DF
    :param ex_: lattice vector
    :param ey_:
    :param ez_:
    :return: the shift matrix for passing from the frame at rest to the moving frame
    """

    d, q, ez_ = _check_dimensions(ex_, ey_, ez_)

    def get_row(m, n, o):
        def get_entry(m, n, o, column):
            coeff = lambda i, m_, n_: pow((ex_[i] - ux), m_) * pow((ey_[i] - uy), n_)
            entry = sum([K[i, column] * coeff(i, m, n) for i in range(0, q)])
            return round_and_simplify(entry)

        row = [get_entry(m, n, o, i) for i in range(0, q)]
        return row

    cm_ = _matrix_maker(f'D{d}Q{q}', get_row)
    return Matrix(cm_)
