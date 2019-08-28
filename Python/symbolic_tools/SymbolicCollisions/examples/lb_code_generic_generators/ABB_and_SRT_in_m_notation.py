from SymbolicCollisions.core.printers import print_as_vector, get_print_symbols_in_indx_notation, get_print_symbols_in_m_notation
from SymbolicCollisions.core.cm_symbols import dynamic_import, moments_dict
from sympy import Symbol
d = 3
q = 27

# DYNAMIC IMPORTS
ex = dynamic_import("SymbolicCollisions.core.cm_symbols", f"ex_D{d}Q{q}")
ey = dynamic_import("SymbolicCollisions.core.cm_symbols", f"ey_D{d}Q{q}")
if d == 3:
    ez = dynamic_import("SymbolicCollisions.core.cm_symbols", f"ez_D{d}Q{q}")
else:
    ez = None

e = dynamic_import("SymbolicCollisions.core.cm_symbols", f"e_D{d}Q{q}")


populations = get_print_symbols_in_m_notation(moments_dict[f'D{d}Q{q}'], 'h')
populations_eq = get_print_symbols_in_m_notation(moments_dict[f'D{d}Q{q}'], 'heq')
for p, p_eq in zip(populations, populations_eq):
    print(f"\t {p} = {-p} + 2 * {p_eq};")

