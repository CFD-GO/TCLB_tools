"""
M - distributions to raw moment transformation matrix
N - raw moments to central moments transformation matrix

based on:
'Modeling incompressible thermal flows using a central-moment-based lattice Boltzmann method'
Linlin Fei, Kai Hong Luo, Chuandong Lin, Qing Li
2017
"""


from sympy import Symbol
from sympy.matrices import Matrix
from SymbolicCollisions.core.printers import round_and_simplify

from SymbolicCollisions.core.cm_symbols import \
    rho, w_D2Q9, m00

from joblib import Parallel, delayed
import multiprocessing


class DiscreteCMTransforms:
    def __init__(self, e, u, F, rho):
        """
        :param e: direction (x,y,z)
        :param u: velocity (x,y,z)
        :param u: Force (x,y,z)
        :param rho: density (not necessarily m00, for instance in multiphase flows)
        """
        self.e = e
        self.u = u
        self.F = F
        self.rho = rho

    def get_gamma_first_order(self, i):
        cs2 = 1. / 3.
        # cs2 = Symbol('cs2')

        """ 
         OMG, sympy...
         Matrix([1]) + 1

         Traceback (most recent call last):
           File "/home/grzegorz/Downloads/pycharm-professional-2018.3.1/pycharm-2018.3.1/helpers/pydev/_pydevd_bundle/pydevd_exec2.py", line 3, in Exec
             exec(exp, global_vars, local_vars)
           File "<input>", line 1, in <module>
           File "/home/grzegorz/GITHUB/TCLB_tools/Python/symbolic_tools/venv/lib/python3.6/site-packages/sympy/core/decorators.py", line 132, in binary_op_wrapper
             return func(self, other)
           File "/home/grzegorz/GITHUB/TCLB_tools/Python/symbolic_tools/venv/lib/python3.6/site-packages/sympy/matrices/common.py", line 1976, in __add__
             raise TypeError('cannot add %s and %s' % (type(self), type(other)))
         TypeError: cannot add <class 'sympy.matrices.dense.MutableDenseMatrix'> and <class 'int'>
         """

        eu = self.e[i, :] * self.u
        gamma = w_D2Q9[i] * (Matrix([1]) + eu / cs2)
        return gamma[0]

    def get_gamma(self, i):
        cs2 = 1. / 3.
        # cs2 = Symbol('cs2')

        """ 
         OMG, sympy...
         Matrix([1]) + 1

         Traceback (most recent call last):
           File "/home/grzegorz/Downloads/pycharm-professional-2018.3.1/pycharm-2018.3.1/helpers/pydev/_pydevd_bundle/pydevd_exec2.py", line 3, in Exec
             exec(exp, global_vars, local_vars)
           File "<input>", line 1, in <module>
           File "/home/grzegorz/GITHUB/TCLB_tools/Python/symbolic_tools/venv/lib/python3.6/site-packages/sympy/core/decorators.py", line 132, in binary_op_wrapper
             return func(self, other)
           File "/home/grzegorz/GITHUB/TCLB_tools/Python/symbolic_tools/venv/lib/python3.6/site-packages/sympy/matrices/common.py", line 1976, in __add__
             raise TypeError('cannot add %s and %s' % (type(self), type(other)))
         TypeError: cannot add <class 'sympy.matrices.dense.MutableDenseMatrix'> and <class 'int'>
         """

        eu = self.e[i, :] * self.u
        u2 = Matrix([self.u.dot(self.u)])
        gamma = w_D2Q9[i] * (Matrix([1]) + eu / cs2 + eu * eu / (2 * cs2 * cs2) - u2 / (2 * cs2))
        return gamma[0]

    def get_EDF_hydro(self, i):
        gamma = self.get_gamma(i)
        g = m00 * w_D2Q9[i] + gamma - w_D2Q9[i]
        return g

    def get_m(self, mno, DF, q):
        k = 0
        for i in range(q):
            pop = DF(i)
            for e_i, mno_i in zip(self.e[i, :], mno):
                pop *= pow(e_i, mno_i)
            k += pop
        return round_and_simplify(k)

    def get_cm(self, mno, DF, q):
        k = 0
        for i in range(q):
            pop = DF(i)
            for e_i, u_i, mno_i in zip(self.e[i, :], self.u, mno):
                pop *= pow((e_i - u_i), mno_i)
            k += pop

        return round_and_simplify(k)

    # TODO: add class attribute weights, moments order?
    def get_force_Guo(self, i):
        """
        'Discrete lattice effects on the forcing term in the lattice Boltzmann method',  Guo et al., 2001
        version for 'Improved locality of the phase-field lattice-Boltzmann model for immiscible fluids at high density ratios' A. Fakhari et. al., 2017
        """
        # extended version with second order terms
        cs2 = 1. / 3.
        # cs2 = Symbol('cs2')

        eu_terms = self.e[i, :] - self.u.transpose() + self.e[i, :].dot(self.u)*self.e[i, :]/cs2
        result = w_D2Q9[i] * self.F.dot(eu_terms) / (self.rho * cs2)
        return result

    def get_force_He(self, i):
        """
        'Discrete Boltzmann equation model for the incompressible Navier-Stokes equation', He et al., 1998
        """
        cs2 = 1. / 3.
        # cs2 = Symbol('cs2')

        eu_dot_f = (self.e[i, :] - self.u.transpose()).dot(self.F)
        pop_eq = m00 * self.get_gamma(i)
        result = pop_eq * eu_dot_f / (rho * cs2)
        return result


def get_mom_vector_from_discrete_def(fun, discrete_transform, moments_order, serial_run=False):
    """
    # obviously 2D is faster
    # However 3D works for 2D as well
    :param fun:
    :param discrete_transform:
    :param moments_order:
    :param serial_run:
    python debugger may crash in parallel mode.
    Moreover code coverage doesn't work multiprocessing, since the processes are independent beings,
    :return:
    """
    # for example: discrete_transform=get_discrete_cm
    q = len(moments_order)

    # row = moments_order[3]
    # result = discrete_transform(row, fun, q)
    if serial_run:
        result = [discrete_transform(row, fun, q) for row in moments_order]  # serial run
    else:  # run in parallel
        # if you experience debugger crashing then run a serial version
        # /pycharm-2018.3.1/helpers/pydev/pydevd.py", line 1487, in dispatch
        #     host = setup['client']
        # TypeError: 'NoneType' object is not subscriptable
        # run in parallel:
        num_cores = multiprocessing.cpu_count()
        result = Parallel(n_jobs=num_cores)(delayed(discrete_transform)(row, fun, q) for row in moments_order)

    return Matrix([result])


def get_mom_vector_from_shift_mat(fun, mat):
    pop = Matrix([fun(i) for i in range(mat.cols)])
    # pop = Matrix(9, 1, lambda i,j: i+j)  # column vect
    cm_ = mat * pop  #for example: Mat=Nraw * Mraw)
    cm_ = round_and_simplify(cm_)
    return Matrix([cm_])


def get_DF(q=9, print_symbol='default_symbol2'):
    symbols_ = [Symbol("%s[%d]" % (print_symbol, i)) for i in range(0, q)]
    return Matrix(symbols_)


def get_m00(q=9, print_symbol='default_symbol3'):
    m00_ = Symbol("%s[%d]" % (print_symbol, 0))

    for i in range(1, q):
        m00_ += Symbol("%s[%d]" % (print_symbol, i))

    return m00_


def get_e(ex, ey, ez):
    q = len(ex)
    symbols_ = [Matrix([ex[i], ey[i], ez[i]]) for i in range(q)]
    return Matrix([symbols_])
