from sfsimodels import models


class PM4Sand(models.CriticalSoil):

    def dilation_angle(self, p_mean):
        critical_relative_density = self._calc_critical_relative_density(p_mean)
        xi_r = critical_relative_density - self.relative_density

    def _calc_critical_relative_density(self, p_mean):
        try:
            return (self.e_max - self.e_critical(p_mean)) / (self.e_max - self.e_min)
        except TypeError:
            return None
