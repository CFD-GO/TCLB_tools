import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from Benchmarks.HotKarman_Correlations.HT_Nu_Correlations import get_Nu_cylinder_by_Churchill_Bernstein

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np

fig = plt.figure()
ax = fig.gca(projection='3d')

# Make data.
Re = np.arange(10, 10000, 10)
Pr = np.arange(10, 10000, 10)
Re, Pr = np.meshgrid(Re, Pr)
# R = np.sqrt(Re ** 2 + Pr ** 2)
# Z = np.sin(R)

Z = get_Nu_cylinder_by_Churchill_Bernstein(Re, Pr)
# Plot the surface.
surf = ax.plot_surface(Re, Pr, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

# Customize the z axis.
# ax.set_zlim(-1.01, 1.01)

ax.set_xlabel('Re')
ax.set_ylabel('Pr')
ax.set_zlabel('Nu')
plt.title(f'Nu=Nu(Re,Pr)')
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.2f'))

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
