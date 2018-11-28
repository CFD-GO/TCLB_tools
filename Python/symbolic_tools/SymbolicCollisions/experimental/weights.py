
from SymbolicCollisions.core.sym_col_fun import *
from SymbolicCollisions.core.cm_symbols import w
from sympy.matrices import Matrix
from sympy import pretty_print, exp
from SymbolicCollisions.core.cm_symbols import ex, ey
from SymbolicCollisions.core.printers import print_as_vector
import numpy as np
from sympy.utilities.iterables import flatten
from sympy import pretty_print

print("\n\n=== weights ===\n")


def get_w_i(i):
    """
    PhD Thesis: `The lattice Boltzmann method: Fundamentals and acoustics`
    by Erlend Magnus Viggen
    4.1  The discrete-velocity Boltzmann equation, pp75
    :param i: i-th lattice direction
    :return: returns weight in i-th lattice direction
    """
    e2 = ex[i]*ex[i] + ey[i]*ey[i]
    cs2 = 1./3.
    D = 2  # dimension of the space
    w_ = 1./pow((2*np.pi*cs2), D/2.)
    w_ *= exp(-e2/(2*cs2))

    return w_


def check_lattice_isotropy(weigths):
    """
    MSc Thesis: `The Lattice Boltzmann Method with Applications in Acoustics`
    by Erlend Magnus Viggen
    Chapter 4.3, pp26  Lattice isotropy
    Appendix C.2, pp92 - matlab code
    :param weigths:
    :return:
    """
    # To obtain correct weights, one should use Hermite polynomials for integration...

    print("\n\n=== Verification of LBM lattice isotropy ===")

    print("0th order: = %f" % sum(weigths))

    print("1st order: = %f" % sum([weigths[i]*ex[i] for i in range(9)]))
    print("1st order: = %f" % sum([weigths[i]*ey[i] for i in range(9)]))

    print("2nd order: = %f" % sum([weigths[i]*ex[i]*ex[i] for i in range(9)]))
    print("2nd order: = %f" % sum([weigths[i]*ey[i]*ey[i] for i in range(9)]))
    print("2nd order: = %f" % sum([weigths[i]*ex[i]*ey[i] for i in range(9)]))


w_test = [get_w_i(i) for i in range(9)]
print(w_test)

check_lattice_isotropy(w)
check_lattice_isotropy(w_test)

