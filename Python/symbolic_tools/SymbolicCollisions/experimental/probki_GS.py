

from SymbolicCollisions.core.cm_symbols import *
from sympy.matrices import Matrix
from sympy import pretty_print, exp
from SymbolicCollisions.core.cm_symbols import ex_D2Q9, ey_D2Q9
from SymbolicCollisions.core.printers import print_as_vector, get_print_symbols_in_indx_notation
import numpy as np
from sympy.utilities.iterables import flatten


DF = get_print_symbols_in_indx_notation('g')
# pretty_print(M_ortho_GS*DF)


print("\n\n=== is orthogonal and orthonormal? ===\n")
pretty_print(M_ortho_GS*M_ortho_GS.transpose())  # TODO: why not M_ortho_GS.transpose()* M_ortho_GS ?!
pretty_print(K_ortho_Geier.transpose()*K_ortho_Geier)



print("\n\n=== from raw moments to ortho moments ===\n")
T_raw_to_ortho = M_ortho_GS * Mraw_D2Q9.inv()
pretty_print(T_raw_to_ortho)

print("\n\n=== relax raw moments in ortho space and go back to raw moments ===\n")
S_relax_ortho = T_raw_to_ortho.inv() * S_relax_MRT_GS * T_raw_to_ortho
pretty_print(S_relax_ortho)


# print("\n\n=== normalize matrix  ===\n")
# from sklearn.preprocessing import normalize
# Column_Normalized = normalize(M_ortho_GS)#, norm='l1', axis=0)
# """axis = 0 indicates, normalize by column and if you are
# interested in row normalization just give axis = 1"""
# pretty_print(Column_Normalized*Column_Normalized.transpose())

print("\n\n=== CM ===\n")

# N_ortho = T_raw_to_ortho*N
# pretty_print(N_ortho)
XXX_ortho = T_raw_to_ortho * Mraw_D2Q9
pretty_print(XXX_ortho * XXX_ortho.transpose())
print_as_vector(XXX_ortho * XXX_ortho.transpose(), outprint_symbol='xxx')



# print("\n\n=== MRT ===\n")


# momencik = get_discrete_cm(2, 0, lambda i: Symbol('f[%d]' % i))
# momencik = get_discrete_m(2, 0, lambda i: Symbol('f[%d]' % i))
# print_as_vector_raw(Matrix([momencik]), print_symbol='momencik')




L = [Matrix([2,3,5]), Matrix([3,6,2]), Matrix([8,3,2])]


def to_GS_matrix(stuff, shape=(3, 3), orthonormal=False, ):
    m_ = GramSchmidt(stuff, orthonormal)
    m_ = Matrix([i for i in flatten(m_)])
    m_ = Matrix(np.reshape(m_, shape)).transpose()
    return m_

#
# m = to_GS_matrix(L, orthonormal=False)
# pretty_print(m)
#
# print("\n\n\ntest ortho")
# pretty_print(m.transpose()*m)



NN = [
    Matrix([1, 0, 0, 0, 0, 0, 0, 0, 0]),
    Matrix([-ux, 1, 0, 0, 0, 0, 0, 0, 0]),
    Matrix([-uy, 0, 1, 0, 0, 0, 0, 0, 0]),
    Matrix([ux * ux, -2 * ux, 0, 1, 0, 0, 0, 0, 0]),
    Matrix([uy * uy, 0, -2 * uy, 0, 1, 0, 0, 0, 0]),
    Matrix([ux * uy, -uy, -ux, 0, 0, 1, 0, 0, 0]),
    Matrix([-ux * ux * uy, 2 * ux * uy, ux * ux, -uy, 0, -2 * ux, 1, 0, 0]),
    Matrix([-uy * uy * ux, uy * uy, 2 * ux * uy, 0, -ux, -2 * uy, 0, 1, 0]),
    Matrix([ux * ux * uy * uy, -2 * ux * uy * uy, -2 * uy * ux * ux, uy * uy, ux * ux, 4 * ux * uy, -2 * uy, -2 * ux, 1])]  # noqa


# GramSchmidt(NN)

MMraw = [
        Matrix([1, 1, 1, 1, 1, 1, 1, 1, 1]),
        Matrix([0, 1, 0, -1, 0, 1, -1, -1, 1]),
        Matrix([0, 0, 1, 0, -1, 1, 1, -1, -1]),
        Matrix([0, 1, 0, 1, 0, 1, 1, 1, 1]),
        Matrix([0, 0, 1, 0, 1, 1, 1, 1, 1]),
        Matrix([0, 0, 0, 0, 0, 1, -1, 1, -1]),
        Matrix([0, 0, 0, 0, 0, 1, 1, -1, -1]),
        Matrix([0, 0, 0, 0, 0, 1, -1, -1, 1]),
        Matrix([0, 0, 0, 0, 0, 1, 1, 1, 1])]


def to_matrix(stuff, shape=(9, 9)):
    m_ = Matrix([i for i in flatten(stuff)])
    m_ = Matrix(np.reshape(m_, shape)).transpose()
    return m_

#
# MMraw2 =to_matrix(MMraw)
# m_ = GramSchmidt(MMraw2)

#
hmm = to_GS_matrix(MMraw, shape=(9, 9))
pretty_print(hmm)

pretty_print(hmm.transpose()*hmm)


# flat= [i for i in flatten(out1)]
# pretty_print(flat)



#
# pretty_print(mmm.transpose())
print("\n=========================================")
pretty_print(NrawD2Q9.transpose() * NrawD2Q9)
