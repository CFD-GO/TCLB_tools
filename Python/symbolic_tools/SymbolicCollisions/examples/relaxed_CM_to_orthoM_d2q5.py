
from sympy import Symbol
from sympy.matrices import Matrix, eye, zeros, ones, diag
from sympy import pretty_print
from SymbolicCollisions.core.cm_symbols import Shift_ortho_Straka_d2q5
from SymbolicCollisions.core.printers import get_print_symbols_in_indx_notation

"""
  See 
  'New Cascaded Thermal Lattice Boltzmann Method for simulations for advection-diffusion and convective heat transfer'
  by K.V. Sharma, R. Straka, F.W. Tavares, 2017
"""

# Smat = get_shift_matrix(K_ortho_Straka_d2q5, ex_Straka_d2_q5, ey_Straka_d2_q5)
# pretty_print(Smat)

Smat = Shift_ortho_Straka_d2q5
k = get_print_symbols_in_indx_notation(q=4, print_symbol='k')
Relax = diag(Symbol('w2'), Symbol('w3'), Symbol('w4'), Symbol('w5'))   #
cm_neq = get_print_symbols_in_indx_notation(q=4, print_symbol='cm_neq')

k = Smat.inv()*Relax*cm_neq

pretty_print(k)
