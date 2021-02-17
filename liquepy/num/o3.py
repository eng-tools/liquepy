from liquepy.num.models import PM4Sand as PM4SandBase
from liquepy.num.models import StressDensityModel as StressDensityModelBase
from liquepy.num import models


class PM4Sand(PM4SandBase):
    type = "pm4sand"
    o3_type = 'pm4sand'

    def __init__(self, wmd=None, liq_mass_density=None, liq_sg=1.0, g=9.8, p_atm=101000.0, **kwargs):
        PM4SandBase.__init__(self, wmd=wmd, liq_mass_density=liq_mass_density, liq_sg=liq_sg, g=g, p_atm=p_atm, **kwargs)
        self._extra_class_inputs = []
        self.app2mod = {
            'd_r': 'relative_density',
            'g_o': 'g0_mod',
            'den': 'unit_moist_mass',
            'nu': 'poissons_ratio'
        }

    def __repr__(self):
        return "PM4SandO3 Soil model, id=%i, phi=%.1f, Dr=%.2f" % (self.id, self.phi, self.relative_density)

    def __str__(self):
        return "PM4SandO3 Soil model, id=%i, phi=%.1f, Dr=%.2f" % (self.id, self.phi, self.relative_density)


class ManzariDafaliasModel(models.ManzariDafaliasModel):
    o3_type = 'manzaridafalias_model'

    def __init__(self, wmd=None, liq_mass_density=None, liq_sg=1.0, g=9.8, p_atm=101000.0, **kwargs):
        models.ManzariDafaliasModel.__init__(self, wmd=wmd, liq_mass_density=liq_mass_density, liq_sg=liq_sg, g=g, p_atm=p_atm, **kwargs)
        self._extra_class_inputs = []
        self.app2mod = {
            'den': 'unit_moist_mass',
            'nu': 'poissons_ratio'
        }

    def __repr__(self):
        return f"ManzariDafaliasModelO3 Soil model, id={self.id}, m_c={self.m_c:.1f}, e_curr={self.e_curr:.2f}"

    def __str__(self):
        return f"ManzariDafaliasModelO3 Soil model, id={self.id}, m_c={self.m_c:.1f}, e_curr={self.e_curr:.2f}"


class StressDensityModel(StressDensityModelBase):
    type = "stress_density_model"

    def __init__(self, wmd=None, liq_mass_density=None, liq_sg=1.0, g=9.8, p_atm=101000.0, **kwargs):
        super(StressDensityModel, self).__init__(wmd=wmd, liq_mass_density=liq_mass_density, liq_sg=liq_sg, g=g, p_atm=p_atm, **kwargs)
        self._extra_class_inputs = []
        self.app2mod = {
            'e_init': 'e_curr',
            'den': 'unit_moist_mass',
            'nu': 'poissons_ratio',
            'n': 'a'
        }

    def __repr__(self):
        return "PM4SandO3 Soil model, id=%i, phi=%.1f, Dr=%.2f" % (self.id, self.phi, self.relative_density)

    def __str__(self):
        return "PM4SandO3 Soil model, id=%i, phi=%.1f, Dr=%.2f" % (self.id, self.phi, self.relative_density)

