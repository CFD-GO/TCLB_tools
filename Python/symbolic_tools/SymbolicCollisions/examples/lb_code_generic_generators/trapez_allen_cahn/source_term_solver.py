# from sympy import *
from sympy import core
# import Exc
from sympy import symbols, Eq, Matrix, solve, lambdify
import numpy as np
import matplotlib.pyplot as plt
from SymbolicCollisions.core.eq_solver import block_simpler, extract_real_solution
import os
# phi = symbols('\phi', real=True) # nieprzesuniete
# Lam = symbols('\lambda', real=True, positive=True, nonzero=True)
# DT = symbols('\delta{t}', real=True, positive=True, nonzero=True)
#
# tilde_phi = symbols('\\tilde{\phi}', real=True)

str_phi = 'phi'
str_lambda = 'lambda'
str_dt = 'dt'
str_tilde_phi = 'tilde_phi'
str_Q = 'Q'

phi = symbols(f'\\{str_phi}', real=True)  # nieprzesuniete
Lambda = symbols(f'\\{str_lambda}', real=True, positive=True, nonzero=True)
DT = symbols(f'\\{str_dt}', real=True, positive=True, nonzero=True)
tilde_phi = symbols(f'\\{str_tilde_phi}', real=True)

given = [tilde_phi]
unknown = [phi]
Q = [Lambda * phi * (1 - phi * phi)]


EQs = Eq(Matrix(given), Matrix(unknown) - DT*Matrix(Q)/2)
solutions = solve(EQs, unknown, dict=True)
symbolic_solutions_as_matrix = Matrix([list(s.values()) for s in solutions])
symbolic_solutions_as_matrix

inputs_as_symbols = [tilde_phi, Lambda, DT]
# inputs_as_str = symbols("tmp[0:%d]" % len(inputs))
inputs_as_str = symbols([str_tilde_phi, str_lambda, str_dt])


##
# make plot to ensure that only real solutions are in the range of interest
calc_numerical_solution = lambdify(
    inputs_as_str, symbolic_solutions_as_matrix.subs(dict(zip(inputs_as_symbols, inputs_as_str))), modules="numpy")

calc_numerical_solution(2, 1, 0.5)  # tilde_phi, Lamda, DT
extract_real_solution(calc_numerical_solution, 2, 1, 0.5)

x = np.linspace(-2, 2, 100)
lambda_num = 1.
phi_num = np.vectorize(lambda t_: extract_real_solution(calc_numerical_solution, t_, lambda_num, 1))(x)
q_num = lambda_num * phi_num * (1 - phi_num*phi_num)

if not os.path.exists('plots'):
    os.makedirs('plots')
fig_name = f'plots/phi_and_Q_lambda_num{lambda_num}'

# -------------------- make dummy plot --------------------
plt.rcParams.update({'font.size': 28})
plt.figure(figsize=(12, 8))

axes = plt.gca()
plt.plot(x, phi_num,
         color="black", marker="", markevery=1, markersize=15, linestyle="--", linewidth=3,
         label=r'$\phi$')

plt.plot(x, q_num,
         color="black", marker="", markevery=1, markersize=15, linestyle=":", linewidth=3,
         label=r'$Q$')

plt.xlabel(r'$\tilde{\phi}$', fontsize=28)
plt.legend()
plt.grid()
fig = plt.gcf()  # get current figure
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.pause(1e-9)  # there is a race condition somewhere in the matplotlib code.
fig.savefig(fig_name + '.pdf', bbox_inches='tight', dpi=200)
plt.show()

# for inp, dum in zip(inputs_as_symbols, inputs_as_str):
#     print(f"{inp} = {dum}")

print(f'const real_t {str_dt} = 1.;')
print(f'real_t {str_phi};')
block_simpler([str_phi], [symbolic_solutions_as_matrix.subs(dict(zip(inputs_as_symbols, inputs_as_str)))[2]])

inputs_as_symbols.append(phi)
inputs_as_str.append(str_phi)

# Qsymbol = symbols('Q', real=True)
print(f'\nreal_t {str_Q};')
block_simpler([str_Q], [Q[0].subs(dict(zip(inputs_as_symbols, inputs_as_str)))])

# Q[0].subs(dict(zip(inputs, dummies)))


