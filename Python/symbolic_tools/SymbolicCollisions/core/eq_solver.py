from sympy import symbols, Eq, Matrix, solve, lambdify, ccode, cse, Add
import numpy as np
import re


def extract_real_solution(calc_numerical_solution, *x):
    tmp = calc_numerical_solution(*x)
    for i in range(3):
        if (np.imag(tmp[i]) == 0): return np.real(tmp[i])
    raise Exception("Something wrong")


def simpler(tt, level=0):
    if level == 0:
        testme = tt.simplify()
        ops0 = testme.count_ops()
        tested = simpler(testme, 1)
        ops1 = tested.count_ops()
        if ops1 < ops0:
            return tested
        else:
            return testme
    else:
        if tt.is_Add:
            coeffs = dict()

            for pp in tt.args:
                for qq in pp.args:
                    if qq in coeffs.keys():
                        coeffs[qq] = coeffs[qq] + 1
                    else:
                        coeffs[qq] = 1

            if len(coeffs) > 2:
                for pp in np.argsort(coeffs.values())[::-1]:
                    #  print tt
                    tt = tt.collect(list(coeffs.keys())[pp])
                #    print tt
                tmp = list()
                for pp in tt.args:
                    tmp.append(simpler(pp, 1))
                return Add(*tmp)
            else:
                return tt
        elif len(tt.args) > 1:
            tmp = list()
            for pp in tt.args:
                tmp.append(simpler(pp.expand(), 1))
            return tt.func(*tmp)
        else:
            return tt


def block_simpler(symbols, exprs, tprefix='const real_t ', print_number_of_operations=False):
    opers0 = 0
    for ee in exprs:
        opers0 = opers0 + ee.count_ops() + 1

    if print_number_of_operations:
        print("// Opers0 = ", opers0)

    temp = cse(exprs, optimizations='basic')
    opers = 0
    for nn, ee in temp[0]:
        ee = simpler(ee)
        opers = opers + ee.count_ops() + 1
        output_code = re.sub(r"cbrt\(3\)", r"cbrt(3.)", ccode(ee, standard='c11'))
        print(tprefix, ccode(nn, standard='c11'), '=', output_code, '; //', ee.count_ops())

    for nn, ee in zip(symbols, temp[1]):
        ee = simpler(ee)
        opers = opers + ee.count_ops() + 1
        output_code = re.sub(r"cbrt\(3\)", r"cbrt(3.)", ccode(ee, standard='c11'))
        print(ccode(nn, standard='c11'), '=', output_code, '; //', ee.count_ops())

    if print_number_of_operations:
        print("// Opers = ", opers)
