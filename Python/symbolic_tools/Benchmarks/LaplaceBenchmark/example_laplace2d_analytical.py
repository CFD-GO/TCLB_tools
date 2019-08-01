from sympy.abc import x
from Benchmarks.LaplaceBenchmark.Laplace_2D_analytical import analytical_laplace_2d
from Benchmarks.LaplaceBenchmark.Laplace_2D_analytical import InputForLaplace2DAnalytical
from Benchmarks.LaplaceBenchmark.Laplace_2D_analytical import make_anal_plot

xSIZE = 64
ySIZE = 64
step = 1
# my_fun = x * (1-x)
my_fun = -4 * x * (x - xSIZE) / (xSIZE * xSIZE)
anal_input = InputForLaplace2DAnalytical(xSIZE, ySIZE, step, my_fun)

xx, yy, zz = analytical_laplace_2d(anal_input)
make_anal_plot(xx, yy, zz)

print("bye")
