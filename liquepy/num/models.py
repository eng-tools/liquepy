from collections import OrderedDict
import sfsimodels as sm
import numpy as np
from sfsimodels.functions import clean_float
from sfsimodels.build_model_descriptions import build_parameter_descriptions
from liquepy.exceptions import deprecation


class PM4Sand(sm.StressDependentSoil):
    _h_po = None
    _crr_n15 = None

    _h_o = None  # not required
    n_b = None
    n_d = None
    a_do = None
    z_max = None
    c_z = None
    c_e = None
    phi_cv = None
    g_degr = None
    c_dr = None
    c_kaf = None
    q_bolt = None
    r_bolt = None
    m_par = None
    f_sed = None
    p_sed = None
    mc_ratio = None
    mc_c = None

    # TODO: add non default inputs here
    type = "pm4sand"

    def __init__(self, pw=9800, liq_mass_density=None, g=9.8, p_atm=101000.0, **kwargs):
        # Note: pw has deprecated
        _gravity = g  # m/s2
        if liq_mass_density:
            _liq_mass_density = liq_mass_density  # kg/m3
        elif pw is not None and _gravity is not None:
            if pw == 9800 and g == 9.8:
                _liq_mass_density = 1.0e3
            else:
                _liq_mass_density = pw / _gravity
        else:
            _liq_mass_density = None

        sm.StressDependentSoil.__init__(self, liq_mass_density=_liq_mass_density, g=_gravity, **kwargs)
        self._extra_class_inputs = [
            "h_po",
            "crr_n15",
            "p_atm",
            "h_o",
            "n_b",
            "n_d",
            "a_do",
            "z_max",
            "c_z",
            "c_e",
            "phi_cv",
            "g_degr",
            "c_dr",
            "c_kaf",
            "q_bolt",
            "r_bolt",
            "m_par",
            "f_sed",
            "p_sed",
            "mc_ratio",
            "mc_c"
        ]
        self.p_atm = p_atm
        self.inputs += self._extra_class_inputs

        if not hasattr(self, "definitions"):
            self.definitions = OrderedDict()
        self.definitions["crr_n15"] = ["Cyclic resistance ratio for 15 cycles", "-"]
        self.definitions["h_po"] = ["Contraction rate parameter", "-"]
        self.definitions["g0_mod"] = ["Normalised shear modulus factor", "-"]
        self.definitions["p_atm"] = ["Atmospheric pressure", "Pa"]

    def __repr__(self):
        return "PM4Sand Soil model, id=%i, phi=%.1f, Dr=%.2f" % (self.id, self.phi, self.relative_density)

    def __str__(self):
        return "PM4Sand Soil model, id=%i, phi=%.1f, Dr=%.2f" % (self.id, self.phi, self.relative_density)

    @property
    def h_po(self):
        return self._h_po

    @h_po.setter
    def h_po(self, value):
        value = clean_float(value)
        self._h_po = value
        if value is not None:
            self._add_to_stack("h_po", value)

    @property
    def hp0(self):
        deprecation('hp0 is deprecated, used h_po')
        return self._h_po

    @hp0.setter
    def hp0(self, value):
        value = clean_float(value)
        self._h_po = value
        if value is not None:
            self._add_to_stack("h_po", value)

    @property
    def csr_n15(self):
        return self._crr_n15

    @property
    def crr_n15(self):
        return self._crr_n15

    @crr_n15.setter
    def crr_n15(self, value):
        value = clean_float(value)
        self._crr_n15 = value
        if value is not None:
            self._add_to_stack("crr_n15", value)

    @csr_n15.setter
    def csr_n15(self, value):
        value = clean_float(value)
        self._crr_n15 = value
        if value is not None:
            self._add_to_stack("crr_n15", value)

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
        if value is not None:
            self._add_to_stack("h_po", value)

    def g_mod_at_v_eff_stress(self, v_eff_stress):  # Override base function since k0 is different
        return self.get_g_mod_at_v_eff_stress(v_eff_stress)

    def get_g_mod_at_v_eff_stress(self, v_eff_stress, k0=None):  # Override base function since k0 is different
        if k0 is None:
            k0 = self.poissons_ratio / (1 - self.poissons_ratio)
        return self.g0_mod * self.p_atm * (v_eff_stress * (1 + k0) / 2 / self.p_atm) ** 0.5

    def set_g0_mod_from_g_mod_at_v_eff_stress(self, g_mod, v_eff_stress, k0=None):
        if k0 is None:
            k0 = self.poissons_ratio / (1 - self.poissons_ratio)
        self.g0_mod = g_mod / self.p_atm / (v_eff_stress * (1 + k0) / 2 / self.p_atm) ** 0.5

    def get_peak_angle(self, p):
        n_b = self.n_b
        if n_b is None:
            n_b = 0.5
        q_bolt = self.q_bolt
        if q_bolt is None:
            q_bolt = 10
        r_bolt = self.r_bolt
        if r_bolt is None:
            r_bolt = 1.5
        return calc_peak_angle_for_pm4sand(self.relative_density, p, p_atm=self.p_atm, phi_cv=self.phi_cv, n_b=n_b, q_bolt=q_bolt, r_bolt=r_bolt)

    def get_dr_cs(self, p):
        return self.r_bolt / (self.q_bolt - np.log(100 * p / self.p_atm))

    def get_ksi_r(self, p):
        return self.get_dr_cs(p) - self.relative_density

    @property
    def phi_r_cs(self):
        if self.phi_cv is not None:
            return np.radians(self.phi_cv)
    
    @property
    def m_cs(self):
        if self.phi_cv is not None:
            return 2 * np.sin(self.phi_r_cs)  # Eq 14

    def get_m_b(self, p):
        return self.m_cs * np.exp(-self.n_b * self.get_ksi_r(p))

    def get_m_d(self, p):
        return self.m_cs * np.exp(self.n_d * self.get_ksi_r(p))



    def set_all_default_pms(self, p, ops=True):
        """
        Set parameters according to the default settings

        Parameters
        ----------
        p: float
            mean effective confining stress
        ops: bool
            If true then set according to the OpenSees version (slightly different to manual)

        Returns
        -------

        """
        self.h_o = max([(0.25 + self.relative_density) / 2, 0.3])
        if self.e_min is None:
            self.e_min = 0.5
        if self.e_max is None:
            self.e_max = 0.8
        if self.n_b is None:
            self.n_b = 0.5
        if self.n_d is None:
            self.n_d = 0.1
        if self.c_z is None:
            self.c_z = 250.0
        if self.c_e is None:
            if ops:
                self.c_e = np.clip(0.2, 0.55, 0.5 - (self.relative_density - 0.55) * 1.5)
            else:
                self.c_e = np.clip(1, 5, 5 - 4 * (self.relative_density - 0.35) / 0.4)
        if self.phi_cv is None:
            self.phi_cv = 33.0
        if self.poissons_ratio is None:
            self.poissons_ratio = 0.3
        if self.g_degr is None:
            self.g_degr = 2.0
        if self.q_bolt is None:
            self.q_bolt = 10.0
        if self.r_bolt is None:
            self.r_bolt = 1.5
        if self.m_par is None:
            self.m_par = 0.01
        ksi_r0 = self.get_ksi_r(p)
        # dr_cs0 = self.get_dr_cs(p)
        if self.c_dr is None:
            self.c_dr = min([5 + 25 * (self.relative_density - 0.35), 10])
        if self.z_max is None:
            self.z_max = min([0.7 * np.exp(-6.1 * ksi_r0), 20])
        m_b = self.get_m_b(p)
        m_d = self.get_m_d(p)
        if self.a_do is None:
            self.a_do = 1. / 0.4 * (np.arcsin(m_b / 2) - np.arcsin(self.m_cs)) / (m_b - m_d)
        if self.c_kaf is None:
            self.c_kaf = np.clip(4.0, 35.0, 5 + 220.0 * (self.relative_density - 0.26) ** 3)
        if self.f_sed is None:
            self.f_sed = min([0.03 * np.exp(2.6 * self.relative_density), 0.99])
        if self.p_sed is None:
            self.p_sed = -self.p_atm / 5
        if self.mc_ratio is None:
            self.mc_ratio = 0.005
        if self.mc_c is None:
            self.mc_c = 0.0


class PM4Silt(sm.StressDependentSoil):
    _h_po = None
    _crr_n15 = None
    s_u = None
    su_rat = None
    _n_g = None
    e_o = None
    h_o = None
    n_bdry = None
    n_bwet = None
    n_d = None
    a_do = None
    ru_max = None
    z_max = None
    c_z = None
    c_e = None
    g_degr = None
    c_kaf = None
    mc_ratio = None
    mc_c = None
    cg_consol = None

    type = "pm4silt"

    def __init__(self, pw=9800, liq_mass_density=None, g=9.8, p_atm=101000.0, **kwargs):
        # Note: pw has deprecated
        _gravity = g  # m/s2
        if liq_mass_density:
            _liq_mass_density = liq_mass_density  # kg/m3
        elif pw is not None and _gravity is not None:
            if pw == 9800 and g == 9.8:
                _liq_mass_density = 1.0e3
            else:
                _liq_mass_density = pw / _gravity
        else:
            _liq_mass_density = None

        sm.StressDependentSoil.__init__(self, liq_mass_density=_liq_mass_density, g=_gravity, **kwargs)
        self._extra_class_inputs = [
            "s_u",
            "su_rat",
            "h_po",
            "crr_n15",
            "p_atm",
            "g_exp",
            "h_o",
            "e_o",
            "n_bwet",
            "n_bdry",
            "n_d",
            "a_do",
            "ru_max",
            "z_max",
            "c_z",
            "c_e",
            "g_degr",
            "c_kaf",
            "mc_ratio",
            "mc_c",
            "cg_consol"
        ]
        self.p_atm = p_atm
        self.inputs += self._extra_class_inputs

        if not hasattr(self, "definitions"):
            self.definitions = OrderedDict()
        self.definitions["crr_n15"] = ["Cyclic resistance ratio for 15 cycles", "-"]
        self.definitions["h_po"] = ["Contraction rate parameter", "-"]
        self.definitions["g0_mod"] = ["Normalised shear modulus factor", "-"]
        self.definitions["p_atm"] = ["Atmospheric pressure", "Pa"]

    def __repr__(self):
        return f"PM4Silt Soil model, id={self.id}, G0={self.g0_mod:.0f}, su={self.s_u:.1}, su_rat={self.su_rat:.1}"

    def __str__(self):
        return f"PM4Silt Soil model, id={self.id}, G0={self.g0_mod:.0f}, su={self.s_u:.1}, su_rat={self.su_rat:.1}"

    @property
    def h_po(self):
        return self._h_po

    @h_po.setter
    def h_po(self, value):
        value = clean_float(value)
        self._h_po = value
        if value is not None:
            self._add_to_stack("h_po", value)

    @property
    def g_exp(self):
        return self.a

    @g_exp.setter
    def g_exp(self, value):
        self.a = value

    @g_exp.setter
    def g_exp(self, value):
        self.a = value

    @property
    def csr_n15(self):
        return self._crr_n15

    @property
    def crr_n15(self):
        return self._crr_n15

    @crr_n15.setter
    def crr_n15(self, value):
        value = clean_float(value)
        self._crr_n15 = value
        if value is not None:
            self._add_to_stack("crr_n15", value)

    @csr_n15.setter
    def csr_n15(self, value):
        value = clean_float(value)
        self._crr_n15 = value
        if value is not None:
            self._add_to_stack("crr_n15", value)

    # def get_n_g(self):
    #     if self.n_g is not None:
    #         return self.n_g
    #     return 0.75  # default
    @property
    def n_g(self):
        if self._n_g is not None:
            return self._n_g
        return None

    @n_g.setter
    def n_g(self, value):
        self._n_g = value
        self._a = value

    @property
    def a(self):
        if self._n_g is not None:
            return self._n_g
        return 0.75  # default

    @a.setter
    def a(self, value):
        self.n_g = value
        self._a = value

    def g_mod_at_v_eff_stress(self, v_eff_stress):  # Override base function since k0 is different
        return self.get_g_mod_at_v_eff_stress(v_eff_stress)

    def get_g_mod_at_v_eff_stress(self, v_eff_stress, k0=None):  # Override base function since k0 is different
        if k0 is None:
            k0 = self.poissons_ratio / (1 - self.poissons_ratio)
        return self.g0_mod * self.p_atm * (v_eff_stress * (1 + k0) / 2 / self.p_atm) ** self.a

    def set_g0_mod_from_g_mod_at_v_eff_stress(self, g_mod, v_eff_stress, k0=None):
        if k0 is None:
            k0 = self.poissons_ratio / (1 - self.poissons_ratio)
        self.g0_mod = g_mod / self.p_atm / (v_eff_stress * (1 + k0) / 2 / self.p_atm) ** self.a
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

    @property
    def phi_r_cs(self):
        if self.phi_cv is not None:
            return np.radians(self.phi_cv)

    @property
    def m_cs(self):
        if self.phi_cv is not None:
            return 2 * np.sin(self.phi_r_cs)  # Eq 14

    def get_m_b(self, p):
        return self.m_cs * np.exp(-self.n_b * self.get_ksi_r(p))

    def get_m_d(self, p):
        return self.m_cs * np.exp(self.n_d * self.get_ksi_r(p))


def calc_peak_angle_for_pm4sand(d_r, p, p_atm=101.0e3, phi_cv=33.0, n_b=0.5, q_bolt=10, r_bolt=1.5):
    dr_cs = r_bolt / (q_bolt - np.log(100 * p / p_atm))  # Eq 11
    ksi_r = dr_cs - d_r  # Eq 10
    phi_r = np.radians(phi_cv)
    big_m = 2 * np.sin(phi_r)  # Eq 14
    if ksi_r < 0:
        m_b = big_m * np.exp(-n_b * ksi_r)  # Eq 13
    else:
        m_b = big_m * np.exp(-n_b / 4 * ksi_r)  # from source code ?
    return np.degrees(np.arcsin(0.5 * m_b))  # Eq 46


class ManzariDafaliasModel(sm.StressDependentSoil):
    crr_n15 = None
    m_c = None
    c_c = None  # c = Me / Mc, the ratio of critical stress ratios in triaxial compression (Mc) and extension (Me)
    lambda_c = None
    e_0 = None
    ksi = None
    m_yield = None
    h_0 = None
    c_h = None
    n_b = None
    n_d = None
    _a_o = None
    z_max = None
    c_z = None
    _g0 = None
    int_scheme = None

    # TODO: add non default inputs here
    type = "manzaridafalias_model"

    def __init__(self, pw=None, wmd=None, liq_mass_density=None, liq_sg=1, g=9.8, p_atm=101000.0, **kwargs):
        # Note: pw has deprecated
        # _gravity = g  # m/s2
        # if liq_mass_density:
        #     _liq_mass_density = liq_mass_density  # kg/m3
        # elif pw is not None and _gravity is not None:
        #     if pw == 9800 and g == 9.8:
        #         _liq_mass_density = 1.0e3
        #     else:
        #         _liq_mass_density = pw / _gravity
        # else:
        #     _liq_mass_density = None

        sm.StressDependentSoil.__init__(self, pw=pw, wmd=wmd, liq_mass_density=liq_mass_density, liq_sg=liq_sg, g=g, **kwargs)
        self._extra_class_inputs = [
            "crr_n15",
            "m_c",
            "c_c",
            "lambda_c",
            "e_0",
            "ksi",
            "m_yield",
            "h_0",
            "c_h",
            "n_b",
            "n_d",
            "a_o",
            "z_max",
            "c_z",
            "g0",
            'int_scheme'
        ]
        self.p_atm = p_atm
        self.inputs += self._extra_class_inputs

        if not hasattr(self, "definitions"):
            self.definitions = OrderedDict()
        self.definitions["crr_n15"] = ["Cyclic resistance ratio for 15 cycles", "-"]
        self.definitions["g0"] = ["Normalised shear modulus factor", "-"]
        self.definitions["p_atm"] = ["Atmospheric pressure", "Pa"]

    def __repr__(self):
        return "ManzariDafaliasModel Soil model, id=%i, phi=%.1f, Dr=%.2f" % (self.id, self.phi, self.relative_density)

    def __str__(self):
        return "ManzariDafaliasModel Soil model, id=%i, phi=%.1f, Dr=%.2f" % (self.id, self.phi, self.relative_density)

    def get_peak_angle(self, p, opt):
        e_cs = self.e_0 - self.lambda_c * (p / self.p_atm) ** self.ksi
        psi = self.e_curr - e_cs
        m_b_txc = self.m_c * np.exp(-self.n_b * psi)  # Eq 13
        if opt == 'txc':
            return np.degrees(np.arcsin(3 * m_b_txc / (m_b_txc + 6)))
        elif opt == 'ps':
            theta = np.pi / 6
            g = 2 * self.c_c / ((1 + self.c_c) - (1 - self.c_c) * np.cos(3 * theta))
            m_b = g * m_b_txc
            return np.degrees(np.arcsin(0.5 * m_b))  # Eq 46
    
    def get_crit_angle(self):
        phi_r = np.arcsin(self.m_c / 2)
        return np.degrees(phi_r)
    
    def set_big_m_from_phi_cv(self, phi_cv):
        phi_r = np.radians(phi_cv)
        self.m_c = 6 * np.sin(np.radians(phi_cv)) / (3 - np.sin(np.radians(phi_cv)))
        self._add_to_stack("m_c", self.m_c)


    @property
    def g0(self):
        return self._g0

    def recompute_g0_mod(self):
        if self.e_curr is not None and self.g0 is not None:
            self._g0_mod = self.g0 * (2.97 - self.e_curr) ** 2 / (1 + self.e_curr)

    @g0.setter
    def g0(self, value):
        # note g0 * (2.97 - e) ** 2 / (1 + e) = g0_mod
        self._g0 = clean_float(value)
        self.recompute_g0_mod()
        if value is not None:
            self._add_to_stack("g0", value)

    @property
    def a_o(self):
        return self._a_o

    @a_o.setter
    def a_o(self, value):
        self._a_o = value
        if value is not None:
            self._add_to_stack("a_o", value)

    @property
    def a_0(self):
        return self.a_o
    
    @a_0.setter
    def a_0(self, value):
        self.a_o = value
        if value is not None:
            self._add_to_stack("a_o", value)
    
    @property
    def g0_mod(self):
        return self._g0_mod

    @g0_mod.setter
    def g0_mod(self, value):
        value = clean_float(value)
        self._g0_mod = value
        if value is not None:
            self._add_to_stack("g0_mod", value)
        # if self.e_curr is None:
        #     raise ValueError('must set e_curr before setting g0_mod')
        if self.e_curr is not None:
            self._g0 = value * (1 + self.e_curr) / (2.97 - self.e_curr) ** 2

    def recompute_all_weights_and_void(self):
        super(sm.StressDependentSoil, self).recompute_all_weights_and_void()
        self.recompute_g0_mod()

    @property
    def e_curr(self):
        """The current void ratio of the soil"""
        return self._e_curr

    @e_curr.setter
    def e_curr(self, value):
        sm.StressDependentSoil.e_curr.fset(self, value)
        # super(sm.StressDependentSoil, self).e_curr = value
        if self.g0 is not None:
            self._g0_mod = self.g0 * (2.97 - self.e_curr) ** 2 / (1 + self.e_curr)
        elif self.g0_mod is not None:
            self._g0 = self.g0_mod * (1 + self.e_curr) / (2.97 - self.e_curr) ** 2


class StressDensityModel(sm.StressDependentSoil):
    _crr_n15 = None
    big_a = None
    a1 = None
    a2 = None
    a3 = None
    b1 = None
    b2 = None
    b3 = None
    fd = None  # degradation constant
    mu_0 = None  # muNot
    mu_cyc = None
    sc = None
    big_m = None
    ssls = None
    hsl = None
    ps = None

    type = "stress_density_model"

    def __init__(self, pw=9800, liq_mass_density=None, g=9.8, p_atm=101000.0, **kwargs):

        super(StressDensityModel, self).__init__(liq_mass_density=liq_mass_density, g=g, **kwargs)
        self._extra_class_inputs = [
            "crr_n15",
            'big_a',
            # 'a',  # pressure dependency exponent for elastic shear modulus
            'a1',
            'a2',
            'a3',
            'b1',
            'b2',
            'b3',
            'fd',  # degradation constant
            'mu_0',  # muNot
            'mu_cyc',
            'sc',
            'big_m',
            'ssls',
            'hsl',
            'ps',
        ]
        self.p_atm = p_atm
        self.inputs += self._extra_class_inputs

        if not hasattr(self, "definitions"):
            self.definitions = OrderedDict()
        self.definitions["crr_n15"] = ["Cyclic resistance ratio for 15 cycles", "-"]
        self.definitions["big_a"] = ["Elastic shear constant", "-"]
        self.definitions["n_e"] = ["Elastic modulus exponent", "-"]
        self.definitions["big_a"] = ["Elastic shear constant", "-"]
        self.definitions["ssls"] = ["QSS-line void ratios", "-"]
        self.definitions["ps"] = ["QSS-line mean effective pressure", "-"]
        self.definitions["a1"] = ["Peak stress ratio coefficient", "-"]
        self.definitions["b1"] = ["Peak stress ratio coefficient", "-"]
        self.definitions["a2"] = ["Max. shear modulus coefficient", "-"]
        self.definitions["b2"] = ["Max. shear modulus coefficient", "-"]
        self.definitions["a3"] = ["Min. shear modulus coefficient", "-"]
        self.definitions["b3"] = ["Min. shear modulus coefficient", "-"]
        self.definitions["mu_0"] = ["Small strain dilantancy coefficient", "-"]
        self.definitions["big_mm"] = ["Critical state stress ratio", "-"]
        self.definitions["sc"] = ["Dilantancy strain", "-"]
        self.definitions["p_atm"] = ["Atmospheric pressure", "Pa"]

    def __repr__(self):
        return "StressDensityModel Soil model, id=%i, phi=%.1f, Dr=%.2f" % (self.id, self.phi, self.relative_density)

    def __str__(self):
        return "StressDensityModel Soil model, id=%i, phi=%.1f, Dr=%.2f" % (self.id, self.phi, self.relative_density)

    @property
    def crr_n15(self):
        return self._crr_n15

    @crr_n15.setter
    def crr_n15(self, value):
        value = clean_float(value)
        self._crr_n15 = value
        if value is not None:
            self._add_to_stack("crr_n15", value)

    @property
    def n_e(self):
        return self._a

    @n_e.setter
    def n_e(self, value):
        value = clean_float(value)
        self._a = value
        if value is not None:
            self._add_to_stack("a", value)

    @property
    def g0_mod(self):
        return self.big_a * (2.17 - self.e_curr) ** 2 / (1 + self.e_curr)

    @g0_mod.setter
    def g0_mod(self, value):
        value = clean_float(value)
        self._g0_mod = None
        self.big_a = value * (1 + self.e_curr) / (2.17 - self.e_curr) ** 2
        if value is not None:
            self._add_to_stack("big_a", value)

    def g_mod_at_v_eff_stress(self, v_eff_stress):  # Override base function since k0 is different
        return self.get_g_mod_at_v_eff_stress(v_eff_stress)

    def get_g_mod_at_v_eff_stress(self, v_eff_stress):  # Override base function since k0 is different
        k0 = 1 - np.sin(self.phi_r)
        p = v_eff_stress * (1 + k0) / 2
        return self.g0_mod * self.p_atm * (p / self.p_atm) ** self.a

    def set_g0_mod_from_g_mod_at_v_eff_stress(self, g_mod, v_eff_stress):
        k0 = 1 - np.sin(self.phi_r)
        p = v_eff_stress * (1 + k0) / 2
        self.g0_mod = g_mod / self.p_atm / (p / self.p_atm) ** self.a
