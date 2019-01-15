import re
from SymbolicCollisions.core.cm_symbols import ux, uy, uz, \
    uxuy, uxuz, uyuz, \
    ux2, uy2, uz2,\
    ux3, uy3, uz3,\
    uxuy3
from sympy.printing import print_ccode
from sympy import simplify, Float, preorder_traversal
from sympy.core.evalf import N as symbol_to_number
from fractions import Fraction
from decimal import Decimal

# HELPERS:
def print_u2():
    print(f"real_t {uxuy} = {ux}*{uy};")
    print(f"real_t {uxuz} = {ux}*{uz};")
    print(f"real_t {uyuz} = {uy}*{uz};")
    print(f"real_t {ux2} = {ux}*{ux};")
    print(f"real_t {uy2} = {uy}*{uy};")
    print(f"real_t {uz2} = {uz}*{uz};")
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
            rounded_stuff = rounded_stuff.subs(a, round(a, 10))

    rounded_and_simplified_stuff = simplify(rounded_stuff)
    return rounded_and_simplified_stuff


def print_as_vector(some_matrix, print_symbol='default_symbol1', raw_output=False):
    rows = some_matrix._mat

    for i in range(len(rows)):
        row = rows[i]  # evaluate symbolic constants, like pi

        if raw_output:
            row = str(row)
        else:
            row = symbol_to_number(row)  # evaluate symbolic constants, like pi
            row = str(round_and_simplify(row))
            row = re.sub(r"%s\*\*2" % ux, '%s' % ux2, row)
            row = re.sub(r"%s\*\*2" % uy, '%s' % uy2, row)
            row = re.sub(r"%s\*\*2" % uz, '%s' % uz2, row)

            row = re.sub(r"%s\*\*3" % ux, '%s' % ux3, row)
            row = re.sub(r"%s\*\*3" % uy, '%s' % uy3, row)
            row = re.sub(r"%s\*\*3" % uz, '%s' % uz3, row)

            row = re.sub(r"%s\*%s" % (ux, uy), '%s' % uxuy, row)
            row = re.sub(r"%s\*%s" % (ux, uz), '%s' % uxuz, row)
            row = re.sub(r"%s\*%s" % (uy, uz), '%s' % uyuz, row)

            # row = re.sub(r"1.\*", "", row)

            row = re.sub(r"_{", "", row)  # skip curly brackets from latex
            row = re.sub(r"}", "", row)  #

            # get algebraic fractions from decimal ones
            result = re.findall(r"\d\.\d+", row)  # may return an empty list: []
            while result:
                first_number = result[0]

                row = re.sub(r"\d\.\d+",  # digit, dot, one or more digits
                             str(Fraction(first_number).limit_denominator(max_denominator=1000)) + '.',
                             row, count=1)

                row = re.sub(r"1\.\*", "", row)  # dont multiply by 1.*
                result = re.findall(r"\d\.\d+", row)

        print(f"{print_symbol}[{i}] = {row};")


def print_as_vector_niewypal(some_matrix, print_symbol='default_symbol1', regex=False):
    rows = some_matrix._mat

    for i in range(len(rows)):
        row = rows[i]  # evaluate symbolic constants, like pi
        row = symbol_to_number(row)  # evaluate symbolic constants, like pi
        row = round_and_simplify(row)

        if regex:
            # print(f'parsing row: {row}')
            # parsed = None

            if row.is_zero:
                parsed = str(0)
            elif len(row.args) == 0:
                parsed = str(row)
            else:
                parsed = []

                for r in row.args:
                    # print(f'{r} is numeric {r.is_Number}')
                    if r.is_Number:
                        r = r.evalf()
                        r = float(r)
                        r = Fraction(r).limit_denominator(max_denominator=1000)

                    r = str(r)
                    parsed.append(r)


            # print(f'parsed row: {parsed}')
            print(f"{print_symbol}[{i}] = {parsed};")
            # parsed = str(round_and_simplify(parsed))
            # parsed = sum(parsed, [])
            # parsed = ' '.join(str(p) for p in [parsed])

            print(f'parsed again row: {parsed}')

