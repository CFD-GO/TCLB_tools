
from sympy.matrices import Matrix
from SymbolicCollisions.core.printers import print_as_vector, print_as_vector_latex
from SymbolicCollisions.core.cm_symbols import e_D2Q9, u2D, F2D, rho, moments_dict
from SymbolicCollisions.core.DiscreteCMTransforms import DiscreteCMTransforms

lattice = 'D2Q9'
dcmt = DiscreteCMTransforms(e_D2Q9, u2D, F2D, rho)
discrete_edf = [dcmt.get_EDF(i) for i in range(0, 9)]
discrete_edf = Matrix(discrete_edf)

print_as_vector(discrete_edf, outprint_symbol='f_eq_from_anal_mom', output_order_of_moments=moments_dict[lattice])
print("Now in LaTeX")
print_as_vector_latex(discrete_edf, outprint_symbol='f^{eq}', output_order_of_moments=moments_dict[lattice])
