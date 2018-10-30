import numpy as np
import scipy.integrate
from liquepy import functions


class ShearTest(object):
    _tau = None
    _gamma = None
    _pp = None
    _vert_eff_stress = None
    _liq_index = None

    def __init__(self, gamma, tau, init_vert_eff_stress=1, sl=None, pp=None):
        self._gamma = gamma
        self._tau = tau
        self.sl = sl
        self._pp = pp
        if init_vert_eff_stress is not None:
            self._init_vert_eff_stress = init_vert_eff_stress

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
    def liq_index(self):
        return self._liq_index

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
        epp = ru * self.init_vert_eff_stress
        self._pp = epp + hydrostatic

    def set_liq_index(self, ru_limit=None, vert_eff_stress_limit=None):
        if ru_limit is not None:
            self._liq_index = functions.determine_t_liq_index(self.ru, ru_limit)
        elif vert_eff_stress_limit is not None:
            ru_limit = 1 - vert_eff_stress_limit / self.init_vert_eff_stress
            self._liq_index = functions.determine_t_liq_index(self.ru, ru_limit)
        else:
            print("No limit set for set_liq_index")

    def average_tau(self, initial=0):
        return scipy.integrate.cumtrapz(self.tau, dx=1, initial=initial)

    def delta_gamma(self):
        delta_gamma = np.diff(self.gamma)
        delta_gamma = np.insert(delta_gamma, 0, self.gamma[0])
        return delta_gamma

