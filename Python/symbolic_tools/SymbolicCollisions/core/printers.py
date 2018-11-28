import re
from SymbolicCollisions.core.cm_symbols import ux, uy, uxuy, ux2, uy2, ux3, uy3, uxuy3
from sympy.printing import print_ccode
from sympy import simplify, Float, preorder_traversal
from sympy.core.evalf import N as symbol_to_number

# HELPERS:
def print_u2():
    print("real_t %s = %s*%s;" % (uxuy, ux, uy))
    print("real_t %s = %s*%s;" % (ux2, ux, ux))
    print("real_t %s = %s*%s;" % (uy2, uy, uy))
    print("")


def print_u3():
    print("real_t %s = %s*%s;" % (ux3, ux2, ux))
    print("real_t %s = %s*%s;" % (uy3, uy2, uy))
    print("real_t %s = %s*%s*%s;" % (uxuy3, uxuy, uxuy, uxuy))
    print("")


def round_and_simplify(stuff):
    simplified_stuff = simplify(stuff)
    rounded_stuff = simplified_stuff

    for a in preorder_traversal(rounded_stuff):
        if isinstance(a, Float):
            rounded_stuff = rounded_stuff.subs(a, round(a, 14))

    rounded_and_simplified_stuff = simplify(rounded_stuff)
    return rounded_and_simplified_stuff


def print_as_vector(some_matrix, print_symbol='default_symbol1', regex=False):
    rows = some_matrix._mat

    for i in range(len(rows)):
        row = rows[i]  # evaluate symbolic constants, like pi

        if regex:
            row = symbol_to_number(row)  # evaluate symbolic constants, like pi
            row = str(round_and_simplify(row))
            row = re.sub("%s\*\*2" % ux, '%s' % ux2, row)
            row = re.sub("%s\*\*2" % uy, '%s' % uy2, row)
            row = re.sub("%s\*%s" % (ux, uy), '%s' % uxuy, row)

            row = re.sub("%s\*\*3" % ux, '%s' % ux3, row)
            row = re.sub("%s\*\*3" % uy, '%s' % uy3, row)
            row = re.sub("%s\*\*3" % uxuy, '%s' % uxuy3, row)

            row = re.sub("0.3333333333333[1-9]{0,3}", "1./3.", row)
            row = re.sub("0.1111111111111[1-9]{0,3}", "1./9.", row)
            row = re.sub("0.2222222222222[1-9]{0,3}", "2./9.", row)
            row = re.sub("0.1666666666666[1-9]{0,3}", "1./6.", row)
            row = re.sub("0.6666666666666[1-9]{0,3}", "2./3.", row)
            row = re.sub("3.9999999999999[1-9]{0,3}", "4.0", row)
            row = re.sub("1.0\*", "", row)
        else:
            row = str(row)

        print("%s[%d] = %s;" % (print_symbol, i, row))
