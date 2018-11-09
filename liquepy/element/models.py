import numpy as np
import scipy.integrate
from liquepy import functions


class ShearTest(object):
    _tau = None
    _gamma = None
    _pp = None
    _esig_v0 = None
    _i_liq = None
    _n_points = 0
    _n_cycles = None
    _ru_limit = None

    def __init__(self, gamma, tau, esig_v0=1, sl=None, pp=None, n_cycles=None):
        self._gamma = np.array(gamma)
        self._tau = np.array(tau)
        self.sl = sl
        self._pp = pp
        if esig_v0 is not None:
            self._esig_v0 = esig_v0
        self._n_points = len(tau)
        self._n_cycles = n_cycles

    @property
    def pp(self):
        return self._pp

    @property
    def tau(self):
        return self._tau

    @property
    def gamma(self):
        return self._gamma

    @property
    def esig_v0(self):
        return self._esig_v0

    @property
    def i_liq(self):
        return self._i_liq

    @property
    def n_points(self):
        return self._n_points

    @property
    def n_cycles(self):
        return self._n_cycles

    @esig_v0.setter
    def esig_v0(self, value):
        self._esig_v0 = value

    @property
    def csr(self):
        try:
            return self.tau / self.esig_v0
        except ValueError:
            return None

    @property
    def epp(self):
        try:
            return self.pp - self.pp[0]
        except ValueError:
            return None

    @property
    def ru(self):
        try:
            return self.epp / self.esig_v0
        except ValueError:
            return None

    def set_pp_via_ru(self, ru, hydrostatic):
        epp = np.array(ru) * self.esig_v0
        self._pp = epp + hydrostatic

    def set_i_liq(self, ru_limit=None, esig_v_limit=None):
        if ru_limit is not None:
            self._ru_limit = ru_limit
            self._i_liq = functions.determine_t_liq_index(self.ru, ru_limit)
        elif esig_v_limit is not None:
            ru_limit = 1 - esig_v_limit / self.esig_v0
            self._ru_limit = ru_limit
            self._i_liq = functions.determine_t_liq_index(self.ru, ru_limit)
        else:
            print("No limit set for set_i_liq")

    @property
    def ru_limit(self):
        return self._ru_limit

    @property
    def av_tau(self):
        average_tau = av_tau = (self.tau[1:] + self.tau[:-1]) / 2
        average_tau = np.insert(average_tau, 0, self.tau[0])  # Include first value
        return average_tau

    @property
    def delta_gamma(self):  # TODO: cache this parameter
        delta_gamma = np.diff(self.gamma)
        delta_gamma = np.insert(delta_gamma, 0, 0)
        return delta_gamma

