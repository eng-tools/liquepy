import numpy as np
from liquepy import functions
import eqsig


class ShearTest(object):
    _stress = None
    _strain = None
    _pp = None
    _esig_v0 = None
    _i_liq = None
    _i_liq_strain = None
    _i_liq_pp = None
    _n_points = 0
    _n_cycles = None
    _ru_limit = None
    _da_strain = None

    def __init__(self, stress, strain, esig_v0=1, sl=None, pp=None, n_cycles=None):
        self._strain = np.array(strain)
        self._stress = np.array(stress)
        self.sl = sl
        self._pp = pp
        if esig_v0 is not None:
            self._esig_v0 = esig_v0
        self._n_points = len(stress)
        self._n_cycles = n_cycles

    @property
    def pp(self):
        return self._pp

    @property
    def stress(self):
        return self._stress

    @property
    def strain(self):
        return self._strain

    @property
    def esig_v0(self):
        return self._esig_v0

    @property
    def i_liq(self):
        return self._i_liq

    @property
    def i_liq_strain(self):
        return self._i_liq_strain

    @property
    def i_liq_pp(self):
        return self._i_liq_pp

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
            return self.stress / self.esig_v0
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

    @n_cycles.setter
    def n_cycles(self, values):
        self._n_cycles = values

    def set_pp_via_ru(self, ru, hydrostatic):
        epp = np.array(ru) * self.esig_v0
        self._pp = epp + hydrostatic

    def set_i_liq(self, ru_limit=None, esig_v_limit=None, strain_limit=None, da_strain_limit=None, or_none=True):
        if ru_limit is not None:
            self._ru_limit = ru_limit
            self._i_liq_pp = functions.determine_t_liq_index(self.ru, ru_limit, return_none=or_none)
        elif esig_v_limit is not None:
            ru_limit = 1 - esig_v_limit / self.esig_v0
            self._ru_limit = ru_limit
            self._i_liq_pp = functions.determine_t_liq_index(self.ru, ru_limit, return_none=or_none)
        elif strain_limit is None:
            pass
            # print("No limit set for set_i_liq")
        if strain_limit is not None:
            self._i_liq_strain = functions.determine_t_liq_index(abs(self.strain), strain_limit, return_none=or_none)
        elif da_strain_limit is not None:
            pinds = eqsig.get_switched_peak_array_indices(self.strain)
            da_strains = self.get_da_strain_series()
            dind = np.where(da_strains > 0.05)
            if len(dind[0]):
                self._i_liq_strain = dind[0][0]
        if self._i_liq_pp is None:
            self._i_liq = self._i_liq_strain
        elif self._i_liq_strain is None:
            self._i_liq = self._i_liq_pp
        else:
            self._i_liq = min(self._i_liq_pp, self._i_liq_strain)

    @property
    def ru_limit(self):
        return self._ru_limit

    @property
    def av_stress(self):
        average_stress = (self.stress[1:] + self.stress[:-1]) / 2
        average_stress = np.insert(average_stress, 0, self.stress[0])  # Include first value
        return average_stress

    @property
    def delta_strain(self):  # TODO: cache this parameter
        delta_strain = np.diff(self.strain)
        delta_strain = np.insert(delta_strain, 0, 0)
        return delta_strain

    def get_da_strain_series(self):
        if self._da_strain is not None:
            return np.array(self._da_strain)
        pinds = eqsig.get_switched_peak_array_indices(self.strain)
        if pinds[-1] != len(self.strain) - 1:
            pinds = np.insert(pinds, len(pinds), len(self.strain) - 1)
        da_strains = [0]
        for j in range(len(pinds) - 1):
            curr_p_strain = self.strain[pinds[j]]
            sgn = np.sign(curr_p_strain)
            if sgn == 0 and j == 0:
                sgn = -1 * np.sign(self.strain[pinds[j + 1]])
            da_strains += list(sgn * -1 * (self.strain[pinds[j] + 1: pinds[j + 1] + 1] - curr_p_strain))
        self._da_strain = np.array(da_strains)
        return np.array(self._da_strain)
