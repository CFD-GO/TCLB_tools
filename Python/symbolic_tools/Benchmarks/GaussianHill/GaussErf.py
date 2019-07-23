from scipy import special
import numpy as np
import matplotlib.pyplot as plt
x = np.linspace(-3, 3)
plt.plot(x, special.erfc(x))
plt.xlabel('$x$')
plt.ylabel('$erfc(x)$')
plt.show()


class GaussianErf:
    def __init__(self, Cmax, b, k):
        """
        :param Cmax: max initial concentration (from 0 to Cmax)
        :param b: initial location of the step
        :param k: conductivity
        """
        self.Cmax = Cmax
        self.b = b
        self.k =k

    def get_concentration(self, x, t):
        c = special.erfc((x-self.b)/np.sqrt(4*self.k*t))
        return c
