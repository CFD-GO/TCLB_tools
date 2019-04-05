
import numpy as np
import matplotlib.pyplot as plt
from Benchmarks.ADE.steady_two_layer_cylinder_analytical_2D import PipeWithinPipe

r0 = 1  # inner radius
r1 = 2  # interface between layers
r2 = 3  # outer radius

k1 = 15  # inner layer - heat conductivity for r0 < r < r1
k2 = 15  # outer layer - heat conductivity for r0 < r < r1

T0 = 10  # temperature for r = r0
T2 = 100  # temperature for r = r2


pwp = PipeWithinPipe(r0, r1, r2, k1, k2, T0, T2)
step = 0.01

r = np.arange(r0, r2, step)
y = np.array([pwp.get_temperature_r(r_) for r_ in r])


fig_name = f'pipe_within_pipe_grid_convergence.png'

