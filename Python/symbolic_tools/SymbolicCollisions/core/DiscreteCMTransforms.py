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
    def __init__(self, e, u, F, rho, cs2=1./3., w=w_D2Q9):
        """
        :param e: direction (x,y,z)
        :param u: velocity (x,y,z)
        :param u: Force (x,y,z)
        :param rho: density (not necessarily m00, for instance in multiphase flows)
        :param cs2: (speed of sound)^2, for isothermal LB cs2=1./3;
            otherwise  cs2 = Symbol('RT', positive=True)  # positive, negative, real, nonpositive, integer, prime and commutative.
        :param w: weight coefficients
        """
        self.e = e
        self.u = u
        self.F = F
        self.rho = rho
        self.cs2 = cs2
        self.w = w

    def get_gamma_first_order(self, i):

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
        # Symbol
        eu = self.e[i, :] * self.u
        gamma = self.w[i] * (Matrix([1]) + eu / self.cs2)
        return gamma[0]

    def get_gamma_first_order_cht(self, i):
        """
        checks the moments from
        'Boundary Conditions at two-phase interface in the LBM for the convection diffusion equation'
        by H. Yoshida et al., 2014
        """
        gamma = Symbol("Gamma")
        w_D3Q7 = Matrix([1. - 6. * gamma, gamma, gamma, gamma, gamma, gamma, gamma])

        eu = self.e[i, :] * self.u
        result = self.rho*w_D3Q7[i] * (Matrix([1]) + eu / (2.*gamma))
        return result[0]

    def get_velocity_bc(self, i):
        eu = self.e[i, :] * self.u
        gamma = self.w[i] * (eu / self.cs2)
        return -2*gamma[0]

    def get_pressure_bc(self, i):
        eu = self.e[i, :] * self.u
        u2 = Matrix([self.u.dot(self.u)])
        gamma = self.w[i] * (Matrix([1]) + eu * eu / (2 * self.cs2 * self.cs2) - u2 / (2 * self.cs2))
        return 2*gamma[0]

    def get_heat_flux_bc(self, i):
        # TODO: check it

        # eu_dot_f = (self.e[i, :] - self.u.transpose()).dot(self.F)
        # pop_eq = m00 * self.get_gamma(i)
        # result = pop_eq * eu_dot_f / ( self.cs2)
        # return
        # eu = self.e[i, :] * self.u
        eu = self.e[i, :] * Matrix([Symbol('nx', positive=True), Symbol('ny', positive=True)])
        gamma = self.w[i] * (eu / self.cs2)
        return -2*gamma[0]

    def get_gamma(self, i):
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
        gamma = self.w[i] * (Matrix([1]) + eu / self.cs2 + eu * eu / (2 * self.cs2 * self.cs2) - u2 / (2 * self.cs2))

        return gamma[0]

    def get_gamma_stabilized(self, i):
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

        preconditioner = Symbol('preconditioner', positive=True)
        eu = self.e[i, :] * self.u
        u2 = Matrix([self.u.dot(self.u)])
        gamma = self.w[i] * (Matrix([1]) +
                             eu / self.cs2 + eu * eu / (2 * preconditioner * self.cs2 * self.cs2)
                             - u2 / (2 * preconditioner * self.cs2))
        return gamma[0]

    def get_EDF(self, i):
        gamma = self.get_gamma(i)
        return m00 * gamma

    def get_EDF_incompressible(self, i):
        gamma = self.get_gamma(i)
        g = m00 * self.w[i] + gamma - self.w[i]
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
        eu_terms = self.e[i, :] - self.u.transpose() + self.e[i, :].dot(self.u)*self.e[i, :]/self.cs2
        result = self.w[i] * self.F.dot(eu_terms) / (self.rho * self.cs2)
        return result

    def get_force_He(self, i):
        """
        'Discrete Boltzmann equation model for the incompressible Navier-Stokes equation', He et al., 1998
        """
        eu_dot_f = (self.e[i, :] - self.u.transpose()).dot(self.F)
        pop_eq = m00 * self.get_gamma(i)
        result = pop_eq * eu_dot_f / (rho * self.cs2)
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

#
# def get_DF(q=9, print_symbol='default_symbol2'):
#     symbols_ = [Symbol("%s[%d]" % (print_symbol, i)) for i in range(0, q)]
#     return Matrix(symbols_)


def get_m00(q=9, print_symbol='default_symbol3'):
    m00_ = Symbol("%s[%d]" % (print_symbol, 0))

    for i in range(1, q):
        m00_ += Symbol("%s[%d]" % (print_symbol, i))

    return m00_

