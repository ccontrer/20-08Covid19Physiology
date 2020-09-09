# Virus load function class

import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

class VirusLoadFunction:

    def __init__(self, tdata, vdata):
        self.tdata = tdata
        self.vdata = vdata
        self.par = np.array([])
        self.p0 = np.array([])
        self.pcov = np.array([])

    def __AssertParameters(self, a1, a2, b1, b2, alpha, minv, maxv):
        assert all(np.array([a1, a2, b1, b2, alpha, minv, maxv]) > 0.),"parameters must be positive"
        assert a1 < a2 < b1 < b2 < np.max(self.tdata),"parameter must satisfy a1 < a2 < b1 < b2 < max time"
        assert minv < maxv,"parameter must satisfy minv < maxv"

    def VirusLoad(self, t, a1, a2, b1, b2, alpha, minv, maxv):
        def v1(t, a1, a2, maxv):
            return 1. + (maxv - 1.)*(np.tanh(6.*(t - (a1 + a2)/2)/(a2 - a1)) - np.tanh(-3.*(a2 + a1)/(a2 - a1)))/2.

        def v2(t, a2, alpha):
            return 1. - np.heaviside(t - a2, 0.5) + np.heaviside(t - a2, 0.5)*np.exp(-alpha*(t - a2))

        def v3(t, b1, b2, minv):
            return 1. - (1. - minv)*(np.tanh(6.*(t - (b1 + b2)/2)/(b2 - b1)) - np.tanh(-3.*(b2 + b1)/(b2 - b1)))/2.

        return np.log10(v1(t, a1, a2, maxv)*v2(t, a2, alpha)*v3(t, b1, b2, minv))        

    def InitialGuess(self, a1, a2, b1, b2, alpha, minv, maxv):
        self.__AssertParameters(a1, a2, b1, b2, alpha, minv, maxv)
        self.p0 = [a1, a2, b1, b2, alpha, minv, maxv]

    def Fit(self, **kwargs):
        maxt = np.max(self.tdata)
        maxv = np.power(np.max(self.vdata), 10)*1e1
        # compute a maximum value for alpha based on data
        maxalpha = 10.0
        minb = 0.
        maxb = [maxt, maxt, maxt, maxt, maxalpha, maxv, maxv]
        par, pcov = curve_fit(self.VirusLoad, self.tdata, self.vdata,
                              p0=self.p0, bounds=(minb, maxb), **kwargs)
        self.__AssertParameters(*par)
        self.par, self.pcov = par, pcov
        self.RSS = sum(np.power(np.ones(4), 2))

    def Plot(self, **kwargs):
        plt.plot(self.tdata, self.vdata, 'ro', label='data', **kwargs)
        if self.par.size:
            tdata = np.linspace(np.min(self.tdata), np.max(self.tdata), num=100)
            plt.plot(tdata, self.VirusLoad(tdata, *self.par), 'b-',
                     label='fit: $a_1$=%2.1f, $a_2$=%2.1f\n $b_1$=%2.1f, $b_2$=%2.1f\n $\\alpha$=%2.1f, min=%1.0e\n max=%1.0e' % tuple(self.par))
        plt.xlabel('time')
        plt.ylabel('V(t)')
        plt.legend()
