from sympy import Function, Symbol, symbols, Derivative, preview, simplify, collect, Wild
from sympy import diff, ln, sin, pprint, sqrt, latex, Integral
import sympy as sym

x, y, z, t = symbols('x y z t')
f = Function('f')
F = Function('F')(x, y, z)
g = Function('g', real=True)(ln(F))

result = diff(sin(x), x)

print(result)
print(f(x * x).diff(x))
print(g.diff(x))

aaa = g.diff(x, 2)
# aaa._args[1] * aaa._args[0]
# init_printing()

pprint(Integral(sqrt(1 / x), x), use_unicode=True)
print(latex(Integral(sqrt(1 / x), x)))

print("\n\n")
pprint(g.diff((x, 1), (y, 0)), use_unicode=True)
# pprint(g.diff((x, 2),(y, 2)), use_unicode=True)

pprint(sym.diff(sym.tan(x), x))
pprint(sym.diff(g, x))

print("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx")
from sympy import Derivative as D, collect, Function
from sympy.abc import a, b, c, x, y, z

f = Function('f')(x)

pprint(a * D(D(f, x), x) ** 2 + b * D(D(f, x), x) ** 2)
pprint(collect(a * D(D(f, x), x) ** 2 + b * D(D(f, x), x) ** 2, D(f, x)))

# print("----------------------------------------------------")
# pprint(deriv.doit().subs({F: Symbol('F')}))



