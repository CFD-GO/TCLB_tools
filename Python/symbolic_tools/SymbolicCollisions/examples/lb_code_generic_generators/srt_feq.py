
from SymbolicCollisions.core.printers import print_as_vector
from sympy.matrices import Matrix

from SymbolicCollisions.core.cm_symbols import e_D2Q9, u2D, F2D, rho, moments_dict
import time

from SymbolicCollisions.core.DiscreteCMTransforms import \
    DiscreteCMTransforms, get_mom_vector_from_discrete_def, get_mom_vector_from_shift_mat


dcmt = DiscreteCMTransforms(e_D2Q9, u2D, F2D, rho)

start = time.process_time()

edf = [dcmt.get_EDF(i) for i in range(0, 9)]
print_as_vector(Matrix([edf]), outprint_symbol=f"feq")

print(f'\n\n Done in {time.process_time() - start} [s].')
