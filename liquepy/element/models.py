import numpy as np
import scipy.integrate
from liquepy import functions


class ShearTest(object):
    _tau = None
    _gamma = None
    _pp = None
    _init_vert_eff_stress = None
    _i_liq = None
    _n_points = 0

    def __init__(self, gamma, tau, init_vert_eff_stress=1, sl=None, pp=None):
        self._gamma = np.array(gamma)
        self._tau = np.array(tau)
        self.sl = sl
        self._pp = pp
        if init_vert_eff_stress is not None:
            self._init_vert_eff_stress = init_vert_eff_stress
        self._n_points = len(tau)


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
    def init_vert_eff_stress(self):
        return self._init_vert_eff_stress

    @property
    def i_liq(self):
        return self._i_liq

    @property
    def n_points(self):
        return self._n_points

    @init_vert_eff_stress.setter
    def init_vert_eff_stress(self, value):
        self._init_vert_eff_stress = value

    @property
    def csr(self):
        try:
            return self.tau / self.init_vert_eff_stress
        except ValueError:
            return None

    @property
    def epp(self):
        try:
            return self.pp - self.init_vert_eff_stress
        except ValueError:
            return None

    @property
    def ru(self):
        try:
            return self.epp / self.init_vert_eff_stress
        except ValueError:
            return None

    def set_pp_via_ru(self, ru, hydrostatic):
        epp = np.array(ru) * self.init_vert_eff_stress
        self._pp = epp + hydrostatic

    def set_i_liq(self, ru_limit=None, vert_eff_stress_limit=None):
        if ru_limit is not None:
            self._i_liq = functions.determine_t_liq_index(self.ru, ru_limit)
        elif vert_eff_stress_limit is not None:
            ru_limit = 1 - vert_eff_stress_limit / self.init_vert_eff_stress
            self._i_liq = functions.determine_t_liq_index(self.ru, ru_limit)
        else:
            print("No limit set for set_i_liq")

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

