
from sympy import Function, Symbol, symbols, Derivative, preview, simplify, collect,  Poly
from sympy import diff, ln, pprint
import sympy as sym


def get_cumulant(xn, yn, zn):
    """
    :param xn: order of x-nth cumulant
    :param yn:
    :param zn:
    :return: cumulant_xn_yn_zn

    Order cumulant's subterms: xn + yn + zn
    ex 5th order = m012*m020
    """

    x, y, z = symbols('x y z')
    F = Function('F')(x, y, z)

    mgf = ln(F)
    deriv = Derivative(mgf, (x, xn), (y, yn), (z, zn))

    all_terms = deriv.doit()
    all_terms = all_terms * F  # Geier: C = c * m000

    all_terms = simplify(all_terms)

    # result = collect(result, F)  # Thank you sympy:  NotImplementedError: Improve MV Derivative support in collect
    # so lets start from highest (2,2,2) ...

    all_terms = all_terms.subs({
        sym.Derivative(F, (x, 2), (y, 2), (z, 2)): Symbol('m_{222}'),
    })

    all_terms = all_terms.subs({
        sym.Derivative(F, (x, 1), (y, 2), (z, 2)): Symbol('m_{122}'),
        sym.Derivative(F, (x, 2), (y, 1), (z, 2)): Symbol('m_{212}'),
        sym.Derivative(F, (x, 2), (y, 2), (z, 1)): Symbol('m_{221}'),
    })

    all_terms = all_terms.subs({
        sym.Derivative(F, (x, 0), (y, 2), (z, 2)): Symbol('m_{022}'),
        sym.Derivative(F, (x, 2), (y, 0), (z, 2)): Symbol('m_{202}'),
        sym.Derivative(F, (x, 2), (y, 2), (z, 0)): Symbol('m_{220}'),
        sym.Derivative(F, (x, 2), (y, 1), (z, 1)): Symbol('m_{211}'),
        sym.Derivative(F, (x, 1), (y, 2), (z, 1)): Symbol('m_{121}'),
        sym.Derivative(F, (x, 1), (y, 1), (z, 2)): Symbol('m_{112}'),
    })

    pretty_dict = {
        sym.Derivative(F, (x, 2), (y, 1), (z, 0)): Symbol('m_{210}'),
        sym.Derivative(F, (x, 1), (y, 2), (z, 0)): Symbol('m_{120}'),

        sym.Derivative(F, (x, 1), (y, 1), (z, 1)): Symbol('m_{111}'),

        sym.Derivative(F, (x, 0), (y, 1), (z, 2)): Symbol('m_{012}'),
        sym.Derivative(F, (x, 0), (y, 2), (z, 1)): Symbol('m_{021}'),

        sym.Derivative(F, (x, 1), (y, 0), (z, 2)): Symbol('m_{102}'),
        sym.Derivative(F, (x, 2), (y, 0), (z, 1)): Symbol('m_{201}'),
    }
    all_terms = all_terms.subs(pretty_dict)

    all_terms = all_terms.subs({
        sym.Derivative(F, (x, 2)): Symbol('m_{200}'),
        sym.Derivative(F, (y, 2)): Symbol('m_{020}'),
        sym.Derivative(F, (z, 2)): Symbol('m_{002}'),
        sym.Derivative(F, (x, 1), (y, 1), (z, 0)): Symbol('m_{110}'),
        sym.Derivative(F, (x, 1), (y, 0), (z, 1)): Symbol('m_{101}'),
        sym.Derivative(F, (x, 0), (y, 1), (z, 1)): Symbol('m_{011}'),
    })

    all_terms = all_terms.subs({
        sym.Derivative(F, (x, 1)): Symbol('m_{100}'),
        sym.Derivative(F, (y, 1)): Symbol('m_{010}'),
        sym.Derivative(F, (z, 1)): Symbol('m_{001}'),
        F: Symbol('m_{000}'),
    })

    # all_terms = collect(all_terms, Symbol('m_{100}'))  # does nothing
    all_terms = Poly(all_terms, F).all_terms()
    all_terms = sum(F ** n * term for (n,), term in all_terms)

    # order = 4
    # given_order_terms = sum(F ** n * term for (n,), term in all_terms if n <= order)
    # given_order_terms_inv = sum(F ** (-n) * term for (n,), term in all_terms if n <= order)

    all_terms = simplify(all_terms)  # does nothing agan

    lower_m000_terms = []
    # # PYDEVD_USE_FRAME_EVAL = NO
    for term in all_terms.args:
        # may crash in debug session... sympy - thank you again
        ht = Symbol('m_{000}')
        higher_terms = [ht**(-2), ht**(-3), ht**(-4)]
        is_lower_order = not any(elem in higher_terms for elem in term.args)
        if is_lower_order:
            lower_m000_terms.append(term)
            # pprint(term)
        # print("------------------------")

    # Wolfram Alpha  - Derivatives of Abstract Functions
    # d^2/(dy*dx)(log(f(x,y)))
    return all_terms, lower_m000_terms



