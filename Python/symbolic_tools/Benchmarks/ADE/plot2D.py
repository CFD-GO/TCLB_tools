import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

from numpy import pi
# %% Sympy Plotting

# shows two plots together
from sympy import symbols
from sympy.plotting import plot

# x = symbols('x')
# plot(-x*x + x, (x,0,1), -3*x*x + x, (x,0,1))
# plt.plot(x, -x*x + x, 'r--', x, - np.cos(2*pi*x) - np.cos(4*pi*x) + 1/6, 'b-')

b = 100
x = np.arange(0., b, 0.01)

plt.plot(x, -4 * x * (x - b) / (b * b), 'r--',
         # x, - np.cos(2 * pi * x) / (pi * pi) - np.cos(4 * pi * x) / (4 * pi * pi) + 1 / 6, 'b-',
         # x, -np.sin(pi * x) / (pi * pi) - np.sin(2 * pi * x) / (4 * pi * pi), 'g:'
         )
plt.grid()
plt.show()
