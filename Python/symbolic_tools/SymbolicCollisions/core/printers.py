import re
from SymbolicCollisions.core.cm_symbols import ux, uy, uz, \
    uxuy, uxuz, uyuz, \
    ux2, uy2, uz2,\
    ux3, uy3, uz3,\
    uxuy3,\
    Sigma2

# from decimal import Decimal
from sympy import simplify, Float, preorder_traversal, Matrix
from sympy.core.evalf import N as symbol_to_number
from fractions import Fraction
from sympy import Symbol

import numpy as np
import pandas as pd




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


def print_sigma_cht():
    # print(f"{Sigma2}")
    print_as_vector(Matrix([Sigma2]), outprint_symbol="real_t Sigma2")
    # print(f"real_t Sigma2 = (h_stability_enhancement*1./3.)/(cp*rho);")
    # print(f"real_t Sigma2 = {Sigma2};")


def round_and_simplify(stuff):
    simplified_stuff = simplify(stuff)
    rounded_stuff = simplified_stuff

    for a in preorder_traversal(rounded_stuff):
        if isinstance(a, Float):
            rounded_stuff = rounded_stuff.subs(a, round(a, 10))

    rounded_and_simplified_stuff = simplify(rounded_stuff)
    return rounded_and_simplified_stuff


def get_print_symbols_in_m_notation(moments_order, print_symbol='m_', as_list=False):
    if type(moments_order) != Matrix:
        moments_order = Matrix(moments_order)

    q = moments_order.shape[0]

    print_symbols = []
    for i in range(q):
        direction = moments_order[i, :]
        direction = [str(d) for d in direction]
        direction = ''.join(direction)
        direction = re.sub(r'-1', '2', direction)
        print_symbols.append(f"{print_symbol}{direction}")

    if as_list:
        return print_symbols
    else:
        return Matrix(print_symbols)


def get_print_symbols_in_indx_notation(q=9, print_symbol='default_symbol2', withbrackets=True, as_list=False):
    # symbols_ = [Symbol("%s[%d]" % (print_symbol, i)) for i in range(0, q)]
    # return Matrix(symbols_)

    if q == 1:
        print_symbols = [Symbol("%s" % print_symbol)]
    elif withbrackets:
        print_symbols = [Symbol("%s[%d]" % (print_symbol, i)) for i in range(0, q)]
    else:
        print_symbols = [Symbol("%s%d" % (print_symbol, i)) for i in range(0, q)]

    if as_list:
        return print_symbols
    else:
        return Matrix(print_symbols)


def print_as_vector(some_matrix, outprint_symbol='default_symbol1', raw_output=False, withbrackets=True, output_order_of_moments=None):
    rows = some_matrix._mat
    q = len(rows)

    if output_order_of_moments is not None:
        print_symbols = get_print_symbols_in_m_notation(output_order_of_moments, outprint_symbol, as_list=True)
    else:
        print_symbols = get_print_symbols_in_indx_notation(q, outprint_symbol, withbrackets, as_list=True)

    for i in range(q):
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
                                            r"(\w+)\*\*2",  # x**2 --> x*x
                                            r"(\w+)\*\*3",  # x**3 --> x*x*x
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
                        msg = 'There is to much square patterns and I dont know not how to simplify them yet.\n ' \
                              'Please generate and check the raw, unparsed output!'
                        # raise NotImplementedError(msg)
                        print(msg)
                        row = re.sub(square_pattern, "*" + to_be_squared[0], row)

                    elif len(to_be_squared) == 1:
                        row = re.sub(square_pattern, "*" + to_be_squared[0], row)

                cube_patterns = [r"\*\*3\.0", r"\*\*3\.", r"\*\*3"]  # order matters
                for cube_pattern in cube_patterns:
                    to_be_squared = re.findall(r"(\w+)" + cube_pattern, row)
                    if len(to_be_squared) > 1:
                        msg = 'There is to much cubic patterns and I dont know not how to simplify them yet.\n ' \
                              'Please generate and check the raw, unparsed output!'
                        # raise NotImplementedError(msg)
                        print(msg)
                        row = re.sub(cube_pattern, "*" + to_be_squared[0] + "*" + to_be_squared[0], row)

                    elif len(to_be_squared) == 1:
                        row = re.sub(cube_pattern, "*" + to_be_squared[0] + "*" + to_be_squared[0], row)


        print(f"\t{print_symbols[i]} = {row};")
        # print(f"stuff = re.sub(r'{print_symbols[i]}', '{row}', stuff)")
