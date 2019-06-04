
from SymbolicCollisions.core.cumulants import get_cumulant
from SymbolicCollisions.core.printers import print_as_vector
from sympy.matrices import Matrix
from sympy import pprint

cumulant, trunc_cumulant = get_cumulant(2, 0, 0)


pprint(cumulant)
# preview(cumulant, output='dvi')  # open preview in new window

# pprint(trunc_cumulant)
# preview(trunc_cumulant, output='dvi')  # open preview in new window
#


mc = Matrix([cumulant])
print_as_vector(mc, 'c')
