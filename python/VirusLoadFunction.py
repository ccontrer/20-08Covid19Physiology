# Virus load function class

import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

class VirusLoadFunction:
    # This class provides and easy way to load, fit, and display virus load data.
    # Usage:
    # >>> vl = VirusLoadFunction(tdate, vdata)
    # >>> vl.InitialGuess(*par0)
    # >>> vl.Plot()
    # >>> vl.Fit()
    # >>> vl.Plot() # after fitting
    # >>> vl.Predict(new_tdata) # evaluates using best estimate

    def __init__(self, tdata, vdata, scale='log10'):
        self.tdata = tdata
        self.vdata = vdata
        self.par = np.array([])
        self.par0 = np.array([])
        self.par_se = np.array([])
        self.scale = scale

    def __AssertParameters(self, a1, a2, b1, b2, alpha, minv, maxv):
        assert all(np.array([a1, a2, b1, b2, alpha, minv, maxv]) > 0.),"parameters must be positive"
        assert a1 < a2 < b1 < b2 < np.max(self.tdata),"parameter must satisfy a1 < a2 < b1 < b2 < max time"
        assert minv < maxv,"parameter must satisfy minv < maxv"

    def __VirusLoad(self, t, a1, a2, b1, b2, alpha, minv, maxv):
        def v1(t, a1, a2, maxv):
            return 1. + (maxv - 1.)*(np.tanh(6.*(t - (a1 + a2)/2)/(a2 - a1)) - np.tanh(-3.*(a2 + a1)/(a2 - a1)))/2.

        def v2(t, a2, alpha):
            return 1. - np.heaviside(t - a2, 0.5) + np.heaviside(t - a2, 0.5)*np.exp(-alpha*(t - a2))

        def v3(t, b1, b2, minv):
            return 1. - (1. - minv)*(np.tanh(6.*(t - (b1 + b2)/2)/(b2 - b1)) - np.tanh(-3.*(b2 + b1)/(b2 - b1)))/2.

        out = v1(t, a1, a2, maxv)*v2(t, a2, alpha)*v3(t, b1, b2, minv)
        if self.scale == 'log10':
            out = np.log10(out)
        return out

    def Eval(self, tdata, par):
        return self.__VirusLoad(tdata, *par)

    def InitialGuess(self, a1, a2, b1, b2, alpha, minv, maxv):
        self.__AssertParameters(a1, a2, b1, b2, alpha, minv, maxv)
        self.par0 = [a1, a2, b1, b2, alpha, minv, maxv]

    def Fit(self, **kwargs):
        maxt = np.max(self.tdata)
        maxv = np.power(np.max(self.vdata), 10)*1e1
        # compute a maximum value for alpha based on data
        maxalpha = 10.0
        minb = 0.
        maxb = [maxt, maxt, maxt, maxt, maxalpha, maxv, maxv]
        par, pcov = curve_fit(self.__VirusLoad, self.tdata, self.vdata,
                              p0=self.par0, bounds=(minb, maxb),
                              method='trf', **kwargs)
        self.__AssertParameters(*par)
        self.par, self.par_se = par, np.sqrt(np.diag(pcov))
        self.RSS = sum(np.power(self.Predict(self.tdata) - self.vdata, 2))

    def Predict(self, ttdata):
        return self.Eval(ttdata, self.par)

    def Plot(self, **kwargs):
        plt.plot(self.tdata, self.vdata, 'ro', label='data', **kwargs)
        if self.par.size:
            tdata = np.linspace(np.min(self.tdata), np.max(self.tdata), num=100)
            plt.plot(tdata, self.Eval(tdata, self.par), 'b-',
                     label='fit: $a_1$=%2.1f, $a_2$=%2.1f\n $b_1$=%2.1f, $b_2$=%2.1f\n $\\alpha$=%2.1f, min=%1.0e\n max=%1.0e' % tuple(self.par))
        plt.xlabel('time')
        plt.ylabel('V(t)')
        plt.legend()
