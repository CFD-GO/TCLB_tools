import re
from SymbolicCollisions.core.cm_symbols import ux, uy, uz, \
    uxuy, uxuz, uyuz, \
    ux2, uy2, uz2,\
    ux3, uy3, uz3,\
    uxuy3

# from decimal import Decimal
from sympy import simplify, Float, preorder_traversal
from sympy.core.evalf import N as symbol_to_number
from fractions import Fraction


def print_u2(d=3):
    print(f"\treal_t {uxuy} = {ux}*{uy};")
    print(f"\treal_t {ux2} = {ux}*{ux};")
    print(f"\treal_t {uy2} = {uy}*{uy};")

    if d == 3:
        print(f"\treal_t {uxuz} = {ux}*{uz};")
        print(f"\treal_t {uyuz} = {uy}*{uz};")
        print(f"\treal_t {uz2} = {uz}*{uz};")

    print("")


def print_u3():
    print("\treal_t %s = %s*%s;" % (ux3, ux2, ux))
    print("\treal_t %s = %s*%s;" % (uy3, uy2, uy))
    print("\treal_t %s = %s*%s*%s;" % (uxuy3, uxuy, uxuy, uxuy))
    print("")


def round_and_simplify(stuff):
    simplified_stuff = simplify(stuff)
    rounded_stuff = simplified_stuff

    for a in preorder_traversal(rounded_stuff):
        if isinstance(a, Float):
            rounded_stuff = rounded_stuff.subs(a, round(a, 10))

    rounded_and_simplified_stuff = simplify(rounded_stuff)
    return rounded_and_simplified_stuff


def print_as_vector(some_matrix, print_symbol='default_symbol1', raw_output=False, withbrackets=True):
    rows = some_matrix._mat

    for i in range(len(rows)):
        row = rows[i]

        if raw_output:
            row = str(row)
        else:
            row = symbol_to_number(round_and_simplify(row))  # evaluate symbolic constants, like pi
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

            row = re.sub(r"_{", "", row)  # skip curly brackets from latex
            row = re.sub(r"}", "", row)  #

            ugly_fractions = re.findall(r"\d\.\d+", row)   # may return an empty list: []
            while ugly_fractions:

                row = re.sub(r"\d\.\d+",  # digit, dot, one or more digits
                             str(Fraction(ugly_fractions[0]).limit_denominator(max_denominator=1000))
                             + '.',  # add '.' to let the C compiler know that it is a float
                             row, count=1)  # get algebraic fractions from decimal ones
                ugly_fractions = re.findall(r"\d\.\d+", row)

            def find_ugly_operations():
                ugly_operations_patterns = [r"(\w+)\*\*1\.0",
                                            r"(\w+)\*\*1\.",  # dont power by 1.0
                                            r"1\.\*",  # dont multiply by 1.*
                                            r"(\w+)\*\*2\.0",
                                            r"(\w+)\*\*2"  # x**2 --> x*x
                                            ]
                ugly_operations = []
                for ugly_operations_pattern in ugly_operations_patterns:
                    ugly_operations += re.findall(ugly_operations_pattern, row)
                return ugly_operations

            while find_ugly_operations():
                row = re.sub(r"\*\*1\.0", r"", row)  # dont power by 1.0
                row = re.sub(r"\*\*1\.", r"", row)  # dont power by 1.0
                row = re.sub(r"1\.\*", "", row)  # dont multiply by 1.*

                square_patterns = [r"\*\*2\.0", r"\*\*2\.", r"\*\*2"]  # order matters

                for square_pattern in square_patterns:
                    to_be_squared = re.findall(r"(\w+)" + square_pattern, row)
                    if len(to_be_squared) > 1:
                        msg = 'There is to much square patterns and I dont know not how to simplify them yet.'
                        # raise NotImplementedError(msg)
                        print(msg)
                        row = re.sub(square_pattern, "*" + to_be_squared[0], row)

                    elif len(to_be_squared) == 1:
                        row = re.sub(square_pattern, "*" + to_be_squared[0], row)

        if withbrackets:
            print(f"\t{print_symbol}[{i}] = {row};")
        else:
            print(f"\t{print_symbol}{i} = {row};")
