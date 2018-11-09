from sfsimodels import models
import sfsimodels as sm
import numpy as np
from collections import OrderedDict
from sfsimodels.functions import clean_float
from sfsimodels.build_model_descriptions import build_parameter_descriptions


class FlacSoil(sm.Soil):
    _tension = 0.0  # default to zero
    hyst_str = ""

    def __str__(self):
        return "Base Soil model, id=%i, phi=%.1f" % (self.id, self.phi)

    def __init__(self, pw=9800):
        super(FlacSoil, self).__init__(pw=pw)  # run parent class initialiser function
        self._extra_class_inputs = ["tension"]
        self.inputs = self.inputs + self._extra_class_inputs
        self.flac_parameters = OrderedDict([
            ("bulk", "bulk_mod"),
            ("shear", "g_mod"),
            ("friction", "phi"),
            ("cohesion", "cohesion"),
            ("tension", "tension"),
            ("density", "density"),
            ("dilation", "dilation_angle"),
            ("por", "porosity"),
            ("perm", "flac_permeability")
        ])
        if not hasattr(self, "definitions"):
            self.definitions = OrderedDict()
        self.definitions["density"] = ["Soil mass density", "kg"]
        self.definitions["tension"] = ["Soil strength in tension", "Pa"]
        self.definitions["flac_permeability"] = ["Permeability of soil", "??"]
        self.prop_dict = OrderedDict()

    def set_prop_dict(self):
        plist = []
        for item in self.flac_parameters:
            plist.append(self.flac_parameters[item])
        self.prop_dict = build_parameter_descriptions(self, user_p=self.definitions, output="dict", plist=plist)

    def find_units(self, parameter):
        if parameter in self.prop_dict:
            return self.prop_dict[parameter]["units"]
        else:
            return None

    @property
    def density(self):
        try:
            return self.unit_dry_weight / 9.81
        except TypeError:
            return None

    @property
    def tension(self):
        return self._tension

    @tension.setter
    def tension(self, value):
        value = clean_float(value)
        self._tension = value

    @property
    def flac_permeability(self):
        try:
            return self.permeability / self._pw
        except TypeError:
            return None


class PM4Sand(FlacSoil, sm.SoilStressDependent):
    _hp0 = None
    _csr_n20 = None

    _h_o = None  # not required
    n_b = None
    n_d = None
    a_do = None
    z_max = None
    c_z = None
    c_e = None
    g_degr = None
    c_kaf = None
    q_bolt = None
    r_bolt = None
    m_par = None
    mc_ratio = None
    mc_c = None

    # TODO: add non default inputs here
    type = "pm4sand"

    def __init__(self, pw=9800, p_atm=101000.0):
        FlacSoil.__init__(self, pw=pw)
        sm.SoilStressDependent.__init__(self, pw=pw)
        self._extra_class_inputs = [
            "hp0",
            "csr_n20",
            "p_atm"
        ]
        self.p_atm = p_atm
        self.inputs += self._extra_class_inputs
        self.pm4sand_parameters = OrderedDict([
            ("D_r", "relative_density"),
            ("h_po", "hp0"),
            ("G_o", "g0_mod"),
            ("density", "density"),
            ("porosity", "porosity"),
            ("h_o", "h_o"),
            ("e_min", "e_min"),
            ("e_max", "e_max"),
            ("n_b", "n_b"),
            ("n_d", "n_d"),
            ("c_z", "c_z"),
            ("c_e", "c_e"),
            ("n_d", "n_d"),
            ("k11", "flac_permeability"),
            ("k22", "flac_permeability"),
            ("P_atm", "p_atm"),
            ("phi_cv", "phi"),
            ("pois", "poissons_ratio"),
            ("A_do", "a_do"),
            ("G_degr", "g_degr"),
            ("Ckaf", "c_kaf"),
            ("Q_bolt", "q_bolt"),
            ("R_bolt", "r_bolt"),
            ("MC_ratio", "mc_ratio"),
            ("MC_c", "mc_c")
            ])
        if not hasattr(self, "definitions"):
            self.definitions = OrderedDict()
        self.definitions["hp0"] = ["Contraction rate parameter", "-"]
        self.definitions["G_o"] = ["Normalised shear modulus factor", "-"]
        self.definitions["p_atm"] = ["Atmospheric pressure", "Pa"]

    def __repr__(self):
        return "PM4Sand Soil model, id=%i, phi=%.1f, Dr=%.2f" % (self.id, self.phi, self.relative_density)

    def __str__(self):
        return "PM4Sand Soil model, id=%i, phi=%.1f, Dr=%.2f" % (self.id, self.phi, self.relative_density)

    def set_prop_dict(self):
        plist = []
        for item in self.flac_parameters:
            plist.append(self.flac_parameters[item])
        for item in self.pm4sand_parameters:
            plist.append(self.pm4sand_parameters[item])
        self.prop_dict = build_parameter_descriptions(self, user_p=self.definitions, output="dict", plist=plist)

    @property
    def hp0(self):
        return self._hp0

    @hp0.setter
    def hp0(self, value):
        value = clean_float(value)
        self._hp0 = value

    @property
    def csr_n20(self):
        return self._csr_n20

    @csr_n20.setter
    def csr_n20(self, value):
        value = clean_float(value)
        self._csr_n20 = value

    @property
    def h_o(self):
        """
        Copy description from manual page 79
        :return:
        """
        return self._h_o

    @h_o.setter
    def h_o(self, value):
        self._h_o = value

    def g_mod_at_v_eff_stress(self, sigma_v_eff):  # Override base function since k0 is different
        k0 = 1 - np.sin(self.phi_r)
        return self.g0_mod * self.p_atm * (sigma_v_eff * (1 + k0) / 2 / self.p_atm) ** 0.5

    # def e_critical(self, p):
    #     p = float(p)
    #     return self.e_cr0 - self.lamb_crl * np.log(p / self.p_cr0)
    #
    # def dilation_angle(self, p_mean):
    #     critical_relative_density = self._calc_critical_relative_density(p_mean)
    #     xi_r = critical_relative_density - self.relative_density
    #
    # def _calc_critical_relative_density(self, p_mean):
    #     try:
    #         return (self.e_max - self.e_critical(p_mean)) / (self.e_max - self.e_min)
    #     except TypeError:
    #         return None
