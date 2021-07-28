from unittest import TestCase
import numpy as np
from Benchmarks.Zaraza.SIR_1D_model import SIR_1D_FD, WSIR_1D_FD
import io
from contextlib import redirect_stdout

import sys
import os
sys.path.append(os.path.join('Python', 'symbolic_tools'))  # allow CI bot to see the stuff from the main repo dir
sys.path.append(os.path.join('.'))  # allow CI bot to see the stuff from the main repo dir


class TestSIR_1D(TestCase):
    def setUp(self):

        self.nx = 128
        self.domain_length = 64
        self.dx = self.domain_length / (self.nx - 1)
        self.xspace = np.linspace(0, self.domain_length, self.nx)

        self.r0 = 15.5  # infectious radius
        self.beta_sir = 3.01  # the average number of contacts per person per time
        self.gamma_sir = 1 / 3.2  # 1 over days to recovery

        self.total_time = 1e0
        self.dt = 1e-3
        self.nt = int(self.total_time / self.dt)

        self.I_IC = np.ones(self.nx) * 0.05  # numpy function ones()
        self.I_IC[int((self.nx - 1) / 4):int(self.nx / 2 + 1)] = 0.6  # setting u = 2 between 0.5 and 1 as per our I.C.s
        self.S_IC = np.ones(self.nx) - self.I_IC
        self.R_IC = np.zeros(self.nx)

        self.N = self.S_IC + self.I_IC + self.R_IC


    def test_sir1D(self):
        S, I, R = SIR_1D_FD(self.S_IC, self.I_IC, self.R_IC,
                            self.nx, self.dx, self.r0, self.beta_sir,
                            self.gamma_sir, self.nt, self.dt)


        np.testing.assert_allclose(self.N - (S+I+R), 0, rtol=1e-14, atol=1e-14)
        np.testing.assert_allclose(S[0] - 0.4126346296533406, 0, rtol=1e-14, atol=1e-14)
        np.testing.assert_allclose(I[0] - 0.5126259921142936, 0, rtol=1e-14, atol=1e-14)
        np.testing.assert_allclose(R[0] - 0.07473937823236734, 0, rtol=1e-14, atol=1e-14)