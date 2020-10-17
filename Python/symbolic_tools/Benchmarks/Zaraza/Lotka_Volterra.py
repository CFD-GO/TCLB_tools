import numpy as np

import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


def lotkavolterra(t, z, a, b, c, d):
    x, y = z
    return [a*x - b*x*y, -c*y + d*x*y]


sol = solve_ivp(lotkavolterra, [0, 20], [10, 5], args=(1.5, 1, 3, 1),
                method='BDF',
                dense_output=True)

t = np.linspace(0, 20, 1000)
z = sol.sol(t)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
fig.suptitle('Lotka-Volterra System')

ax1.plot(t, z[0].T, label='x')
ax1.plot(t, z[1].T, label='y')
ax1.set(xlabel='time', ylabel='x, y')
ax1.set_title('x(t), y(t)')
ax1.grid(which='both')


ax2.plot(z[0], z[1], 'tab:green')
ax2.set(xlabel='x', ylabel='y')
ax2.set_title('y(x)')
ax2.grid(which='both')

# plt.plot(t, z.T)
# plt.xlabel('t')
# plt.legend(['x', 'y'], shadow=True)
# plt.title('Lotka-Volterra System')
# plt.show()

# plt.plot(z[0], z[1])
# plt.xlabel('x')
# plt.xlabel('y')
# plt.legend(['x(y)'], shadow=True)
plt.title('Lotka-Volterra System')
plt.show()