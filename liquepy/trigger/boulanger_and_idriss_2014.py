import numpy as np
from liquepy.exceptions import deprecation
import sfsimodels as sm
from liquepy.field import CPT


def calc_void_ratio(unit_dry_weight, specific_gravity, pw):
    return specific_gravity * pw / unit_dry_weight - 1


def calc_unit_dry_weight(fs, q_t, p_a, unit_water_wt):
    """
    Estimate the unit weight of the soil.

    Ref: https://www.cpt-robertson.com/PublicationsPDF/Unit%20Weight%20Rob%20%26%20Cabal%20CPT10.pdf

    Parameters
    ----------
    fs: array_like
        CPT skin friction (kPa)
    q_t: array_like
        CPT cone tip resistance (kPa)
    p_a: float
        Atmospheric pressure (kPa)
    unit_water_wt: float
        Unit weght of water

    Returns
    -------

    """
    # eq Robertson pag 37- CPT guide
    # unit_water_wt = 9.81
    np.clip(q_t, 1e-10, None, out=q_t)
    r_f = np.clip((fs / q_t) * 100, 0.1, None)
    min_unit_weight = 1.5 * unit_water_wt  # minimum value obtained in presented results
    max_unit_weight = 4.0 * unit_water_wt  # maximum value obtained in presented results
    soil_unit_wt = np.clip((0.27 * np.log10(r_f) + 0.36 * np.log10(q_t / p_a) + 1.236) * unit_water_wt, min_unit_weight, max_unit_weight)
    return soil_unit_wt


def calc_unit_weight(e_curr, specific_gravity, saturation, pw):
    unit_void_volume = e_curr / (1 + e_curr)
    unit_dry_weight = (specific_gravity * pw) / (1 + e_curr)
    return unit_dry_weight + unit_void_volume * saturation * pw


def calc_sigma_v(depths, gammas, gamma_predrill=17.0):
    """
    Calculates the vertical stress

    Note: properties are forward projecting
    """
    predrill_depth = depths[0]
    depth_incs = depths[1:] - depths[:-1]
    depth_incs = np.insert(depth_incs, 0, depth_incs[0])
    sigma_v_incs = depth_incs * gammas
    sigma_v = np.cumsum(sigma_v_incs) + predrill_depth * gamma_predrill
    return sigma_v


def calc_pore_pressure(depth, gwl, unit_water_wt):
    pore_pressure = np.where(depth > gwl, (depth - gwl) * unit_water_wt, 0.0)
    return pore_pressure


def calc_sigma_veff(sigma_v, pore_pressure):
    sigma_veff = abs(sigma_v - pore_pressure)
    return sigma_veff


def calc_qt(qc, ar, u2):
    """

    :param qc: kPa, cone tip resistance
    :param ar: -, area ratio
    :param u2: kPa, water pressure beneath cone tip
    :return:
    """
    # qt the cone tip resistance corrected for unequal end area effects, eq 2.3
    return qc + ((1 - ar) * u2)


def calc_f_ic_values(fs, qt, sigma_v):
    # qt is in kPa, so it's not necessary measure unit transformation
    return (fs / (qt - sigma_v)) * 100


def calc_big_q_values(qt, sigma_v, sigma_veff, p_a, n_val=0.5):
    """
    Eq. 2.26
    :param c_n: CN
    :param qt:
    :param sigmav:
    :return:
    """
    # return (qt - sigmav) / 100 * c_n   # this is different to eq 2.26
    return (qt - sigma_v) / p_a * (p_a / sigma_veff) ** n_val


def calc_i_c(big_q, big_f):
    """
    Calculates the index parameter of the soil

    Eq. 2.26

    :param big_q: float or array,
    :param big_f: float or array,
    :return:
    """

    if big_f <= 0.1:
        big_f = 0.1
    if big_q <= 1:
        big_q = 1
    return ((3.47 - np.log10(big_q)) ** 2 + (1.22 + np.log10(big_f)) ** 2) ** 0.5


def calc_fc(ic, cfc):
    fc_tent = 80 * (ic + cfc) - 137
    fc1 = np.where(fc_tent > 100, 100, fc_tent)
    fc = np.where(fc1 <= 137 / 80 - cfc, 0, fc1)
    return fc


def calc_rd(depth, magnitude):
    """
    rd from CPT, Eq 2.14a
    """
    alpha = -1.012 - 1.126 * np.sin((depth / 11.73) + 5.133)
    beta = 0.106 + 0.118 * np.sin((depth / 11.28) + 5.142)
    rd = np.exp(alpha + beta * magnitude)
    return rd


def calc_csr(sigma_veff, sigma_v, pga, rd, gwl, depth):
    """
    Cyclic stress ratio from CPT, Eq 2.2,
    """
    return 0.65 * (sigma_v / sigma_veff) * rd * pga


def calc_cn_values(m, sigma_veff):
    """
    CN parameter from CPT, Eq 2.15a
    """
    c_n = (100 / sigma_veff) ** m
    if c_n > 1.7:
        c_n = 1.7
    return c_n


def calc_m(q_c1ncs):
    """
    m parameter from CPT, Eq 2.15b
    """
    m = 1.338 - 0.249 * q_c1ncs ** 0.264
    if q_c1ncs >= 254:
        m = 0.263823991466759
    if q_c1ncs <= 21:
        m = 0.781756126201587
    return m


def calc_q_c1n(q_c, c_n):
    """
    qc1n from CPT, Eq 2.4
    """
    q_c1n = c_n * q_c * 1000 / 100
    return q_c1n


def calc_crr_m7p5_from_qc1ncs(q_c1n_cs, c_0=2.8):
    """
    Calculation of cyclic resistance ratio for Magnitude 7.5 earthquake

    Parameters
    ----------
    q_c1n_cs: float or array_like
        Clean-sand equivalent, normalised cone tip resistance
    c_0: float (default=2.8)
        Empirical fitting parameter.
         - 2.8=16th percentile (commonly used)
         - 2.6=median response

    Returns
    -------

    """
    return np.exp((q_c1n_cs / 113) + ((q_c1n_cs / 1000) ** 2) -
                        ((q_c1n_cs / 140) ** 3) + ((q_c1n_cs / 137) ** 4) - c_0)


def calc_crr_n15_from_qc1ncs(q_c1ncs, c_0=2.8):
    crr_m7p5 = calc_crr_m7p5_from_qc1ncs(q_c1ncs, c_0=c_0)
    msf_m = 1.09 + (q_c1ncs / 180) ** 3
    msf_max = np.clip(msf_m, None, 2.2)
    b = calc_b_from_msf_max_bi2014(msf_max)
    n_m7p5 = calc_n_cycles_at_m7p5_bi2014(b)
    n_cyc_bi2014 = 15
    msf = (n_m7p5 / n_cyc_bi2014) ** b
    return msf * crr_m7p5


def calc_crr_m7p5_from_qc1ncs_capped(q_c1n_cs, gwl, depth, i_c, i_c_limit=2.6, c_0=2.8):
    """
    cyclic resistance from CPT, Eq. 2.24
    it's not possible to have liquefaction above water table
    """
    crr_values = calc_crr_m7p5_from_qc1ncs(q_c1n_cs, c_0)
    crr_tent = np.where(depth < gwl, 4, crr_values)
    return np.where(i_c <= i_c_limit, crr_tent, 4.)


def crr_m(ksigma, msf, crr_m7p5):
    """Cyclic resistance ratio corrected for m_w and confining stress"""
    return crr_m7p5 * ksigma * msf


def calc_delta_q_c1n(q_c1n, fc):
    """
    delta q_c1n from CPT, Eq 2.22
    """
    delta_q_c1n = (11.9 + (q_c1n / 14.6)) * (np.exp(1.63 - (9.7 / (fc + 2)) - ((15.7 / (fc + 2)) ** 2)))
    return delta_q_c1n


def calc_q_c1ncs(q_c1n, delta_q_c1n):
    """
    q_c1ncs from CPT, Eq 2.10
    """
    q_c1ncs = q_c1n + delta_q_c1n
    return q_c1ncs


def calc_msf(magnitude, q_c1ncs):
    """
    Magnitude scaling factor to correct the cyclic resistance ratio

    Eq. 2.19

    :param magnitude: earthquake m_w
    :param q_c1ncs: clean sand-corrected normalised cone tip resistance
    :return:
    """
    if magnitude == 7.5:
        return np.ones_like(q_c1ncs)
    msf_m = 1.09 + (q_c1ncs / 180) ** 3
    msf_max = np.where(msf_m > 2.2, 2.2, msf_m)
    msf = 1. + (msf_max - 1) * (8.64 * np.exp(-magnitude / 4) - 1.325)
    return msf


def calc_k_sigma(sigma_eff, q_c1ncs, pa=100):
    """
    Overburden correction factor, K_sigma

    Equation 2.16a

    :param sigma_eff: vertical effective stress
    :param q_c1ncs: clean sand-corrected normalised cone tip resistance
    :param pa: atmospheric pressure in kPa
    :return:
    """
    c_sigma_unrestricted = 1. / (37.3 - 8.27 * (q_c1ncs ** 0.264))
    c_sigma = np.where(c_sigma_unrestricted <= 0.3, c_sigma_unrestricted, 0.3)
    k_sigma_unrestricted = 1 - c_sigma * np.log(sigma_eff / pa)
    k_sigma = np.where(k_sigma_unrestricted <= 1.1, k_sigma_unrestricted, 1.1)
    return k_sigma


def _calc_dependent_variables(sigma_v, sigma_veff, q_c, f_s, p_a, q_t, cfc):
    """
    Iteratively calc_volumetric_strain parameters as they are interdependent

    Parameters
    ----------
    sigma_v: array_like [kPa]
        Total vertical stress
    sigma_veff: array_like [kPa]
        Effective vertical stress
    q_c: array_like [kPa]
        Cone tip resistance
    f_s: array_like [kPa]
        Skin friction
    p_a: array_like [kPa]
        Atmospheric pressure
    q_t:
    cfc: float
        correction factor
    :return:
    """
    num_depth = len(sigma_v)
    m_values = np.ones(num_depth)  # create an array of ones
    cn_values = np.zeros(num_depth)
    q_c1n = np.zeros(num_depth)
    delta_q_c1n = np.zeros(num_depth)
    q_c1n_cs = np.zeros(num_depth)

    big_q = np.ones(num_depth)
    ft_values = np.ones(num_depth)
    i_c = np.ones(num_depth)
    fines_content = np.ones(num_depth)

    for dd in range(0, num_depth):
        temp_q_c1n = 1e6
        n_val = 1.0
        for j in range(100):
            cn_values[dd] = min((p_a / sigma_veff[dd]) ** m_values[dd], 1.7)  # Eq 2.15a
            q_c1n[dd] = (cn_values[dd] * q_c[dd] / p_a)  # Eq. 2.4
            big_q[dd] = calc_big_q_values(q_t[dd], sigma_v[dd], sigma_veff[dd], p_a, n_val=n_val)
            ft_values[dd] = calc_f_ic_values(f_s[dd], q_t[dd], sigma_v[dd])
            i_c[dd] = calc_i_c(big_q[dd], ft_values[dd])
            if i_c[dd] < 2.6 and n_val == 1.0:  # See second half of pg 449 of Robertson and Wride (1997)
                n_val = 0.5
                n_val_stable = False
            elif i_c[dd] > 2.6 and n_val == 0.5:
                n_val = 0.75
                n_val_stable = False
            else:
                n_val_stable = True
            fines_content[dd] = calc_fc(i_c[dd], cfc)

            delta_q_c1n[dd] = calc_delta_q_c1n(q_c1n=q_c1n[dd], fc=fines_content[dd])  # Eq. 2.22
            q_c1n_cs[dd] = q_c1n[dd] + delta_q_c1n[dd]
            m_values[dd] = calc_m(q_c1n_cs[dd])
            if abs(q_c1n[dd] - temp_q_c1n) < 0.00001 and n_val_stable:
                break
            temp_q_c1n = q_c1n[dd]
    return q_c1n_cs, q_c1n, fines_content, i_c, big_q, ft_values


class BoulangerIdriss2014CPT(object):

    def __init__(self, cpt, gwl=None, pga=0.25, m_w=None, cfc=0.0, **kwargs):
        """
        Performs the Boulanger and Idriss triggering procedure for a CPT profile

        ref: Boulanger:2014id

        Parameters
        ----------

        gwl: float, m,
            ground water level below the surface
        pga: float, g,
            peak ground acceleration
        m_w: float, -,
            Earthquake magnitude
        a_ratio: float, -, default=0.8
            Area ratio
        cfc: float, -, default=0.0
            Fines content correction factor for Eq 2.29
        magnitude: float, -,
            Earthquake magnitude (deprecated)
        i_c_limit: float, -, default=2.6
            Limit of liquefiable material
        s_g: float or array_like, -, default=2.65
            Specific gravity
        s_g_water: float, -, default=1.0
            Specific gravity of water
        p_a: float, -, kPa, default=101
            Atmospheric pressure
        """

        magnitude = kwargs.get("magnitude", None)
        i_c_limit = kwargs.get("i_c_limit", 2.6)
        self.s_g = kwargs.get("s_g", 2.65)
        self.s_g_water = kwargs.get("s_g_water", 1.0)
        self.p_a = kwargs.get("p_a", 101.)  # kPa
        self.c_0 = kwargs.get("c_0", 2.8)
        saturation = kwargs.get("saturation", None)
        unit_wt_method = kwargs.get("unit_wt_method", "robertson2009")
        gamma_predrill = kwargs.get("gamma_predrill", 17.0)
        if gwl is None and cpt.gwl is not None:
            gwl = cpt.gwl

        if m_w is None:
            if magnitude is None:
                self.m_w = 7.5
            else:
                deprecation('Deprecated input "magnitude" in BoulangerIdriss2014(), should use "m_w"')
                self.m_w = magnitude
        else:
            self.m_w = m_w

        unit_water_wt = self.s_g_water * 9.8
        self.npts = len(cpt.depth)
        self.depth = cpt.depth
        self.cpt = cpt
        # self.q_c = cpt.q_c
        # self.f_s = cpt.f_s
        # self.u_2 = cpt.u_2
        self.gwl = gwl
        self.pga = pga
        self.a_ratio = cpt.a_ratio
        if cpt.a_ratio is None:
            self.a_ratio = 0.8
        self.i_c_limit = i_c_limit

        self.cfc = cfc  # parameter of fines content, eq 2.29
        self.q_t = calc_qt(self.cpt.q_c, self.a_ratio, self.cpt.u_2)  # kPa

        if saturation is None:
            self.saturation = np.where(self.depth < self.gwl, 0, 1)
        else:
            self.saturation = saturation
        if unit_wt_method == "robertson2009":
            self.unit_wt = calc_unit_dry_weight(self.cpt.f_s, self.q_t, self.p_a, unit_water_wt)
        elif unit_wt_method == 'void_ratio':
            self.unit_dry_wt = calc_unit_dry_weight(self.cpt.f_s, self.q_t, self.p_a, unit_water_wt)
            self.e_curr = calc_void_ratio(self.unit_dry_wt, self.s_g, pw=unit_water_wt)
            self.unit_wt = calc_unit_weight(self.e_curr, self.s_g, self.saturation, pw=unit_water_wt)
        else:
            raise ValueError("unit_wt_method should be either: 'robertson2009' or 'void_ratio' not: %s" % unit_wt_method)

        self.sigma_v = calc_sigma_v(self.depth, self.unit_wt, gamma_predrill)
        self.pore_pressure = calc_pore_pressure(self.depth, self.gwl, unit_water_wt)
        self.sigma_veff = calc_sigma_veff(self.sigma_v, self.pore_pressure)
        if self.sigma_veff[0] == 0.0:
            self.sigma_veff[0] = 1.0e-10
        self.rd = calc_rd(self.depth, self.m_w)

        self.q_c1n_cs, self.q_c1n, self.fines_content, self.i_c, self.big_q, self.big_f = _calc_dependent_variables(self.sigma_v,
                                                                                                        self.sigma_veff,
                                                                                                        self.cpt.q_c,
                                                                                            self.cpt.f_s, self.p_a,
                                                                                            self.q_t,
                                                                                            self.cfc)

        np.clip(self.q_c1n_cs, None, 211., out=self.q_c1n_cs)  # pg 11 of BI2014 report
        self.k_sigma = calc_k_sigma(self.sigma_veff, self.q_c1n_cs)
        self.msf = calc_msf(self.m_w, self.q_c1n_cs)
        self.csr = calc_csr(self.sigma_veff, self.sigma_v, pga, self.rd, gwl, self.depth)
        self.crr_m7p5 = calc_crr_m7p5_from_qc1ncs_capped(self.q_c1n_cs, gwl, self.depth, self.i_c, self.i_c_limit, self.c_0)
        self.crr = crr_m(self.k_sigma, self.msf, self.crr_m7p5)  # CRR at set magnitude
        fs_unlimited = self.crr / self.csr
        # fs_fines_limited = np.where(self.fines_content > 71, 2.0, fs_unlimited)  # based on I_c=2.6
        fos = np.where(fs_unlimited > 2, 2, fs_unlimited)
        self.factor_of_safety = np.where(self.i_c <= self.i_c_limit, fos, 2.25)

    @property
    def gammas(self):
        deprecation('Deprecated "BoulangerIdriss2014.gammas", should use "BoulangerIdriss2014.unit_wt"')
        return self.unit_wt

    @property
    def magnitude(self):
        deprecation('Deprecated "BoulangerIdriss2014.magnitude", should use "BoulangerIdriss2014.m_w"')
        return self.m_w


class BoulangerIdriss2014(BoulangerIdriss2014CPT):

    def __init__(self, depth, q_c, f_s, u_2, cpt_gwl=None, pga=0.25, m_w=None, gwl=None, a_ratio=0.8, cfc=0.0, **kwargs):
        """
         depth: array_like m,
            depths measured downwards from surface
        q_c: array_like kPa,
            cone tip resistance
        f_s: array_like kPa,
            skin friction
        u_2: array_like kPa,
            water pressure beneath cone tip
        Parameters
        ----------
        cpt
        pga
        m_w
        gwl
        a_ratio
        cfc
        kwargs
        """

        if gwl is None:
            gwl = cpt_gwl
        cpt = CPT(depth, q_c, f_s, u_2, gwl, a_ratio=a_ratio)
        super(BoulangerIdriss2014CPT, self).__init__(cpt, gwl=gwl, pga=pga, m_w=m_w, cfc=cfc, **kwargs)


class BoulangerIdriss2014SoilProfile(object):  # TODO: validate this properly
    def __init__(self, sp, pga=0.25, m_w=None, **kwargs):
        self.sp = sp
        assert isinstance(self.sp, sm.SoilProfile)
        self.inc = 0.01
        self.sp.gen_split(target=self.inc, props=['csr_n15'])
        split_depths = np.cumsum(self.sp.split['thickness'])
        self.depth = np.arange(0, sp.height + self.inc, self.inc)
        self.npts = len(self.depth)

        self.s_g = kwargs.get("s_g", 2.65)
        self.s_g_water = kwargs.get("s_g_water", 1.0)
        saturation = kwargs.get("saturation", None)

        if m_w is None:
            self.m_w = 7.5
        else:
            self.m_w = m_w

        self.gwl = sp.gwl
        self.pga = pga

        if saturation is None:
            self.saturation = np.where(self.depth < self.gwl, 0, 1)
        else:
            self.saturation = saturation

        self.sigma_v = sp.get_v_total_stress_at_depth(self.depth) / 1e3
        self.pore_pressure = sp.get_hydrostatic_pressure_at_depth(self.depth) / 1e3
        self.sigma_veff = self.sigma_v - self.pore_pressure
        if self.sigma_veff[0] == 0.0:
            self.sigma_veff[0] = 1.0e-10
        self.rd = calc_rd(self.depth, self.m_w)
        crr_unlimited = np.interp(self.depth, split_depths, self.sp.split['csr_n15'])
        self.crr_m7p5 = np.where(self.depth <= self.gwl, 4, crr_unlimited)
        self.q_c1n_cs = calc_q_c1n_cs_from_crr_m7p5(self.crr_m7p5)
        self.k_sigma = calc_k_sigma(self.sigma_veff, self.q_c1n_cs)
        self.msf = calc_msf(self.m_w, self.q_c1n_cs)
        self.crr = crr_m(self.k_sigma, self.msf, self.crr_m7p5)  # CRR at set magnitude
        self.csr = calc_csr(self.sigma_veff, self.sigma_v, pga, self.rd, self.gwl, self.depth)
        fs_unlimited = self.crr / self.csr
        self.factor_of_safety = np.where(fs_unlimited > 2, 2, fs_unlimited)


def run_bi2014(cpt, pga, m_w, gwl=None, p_a=101., cfc=0.0, i_c_limit=2.6, gamma_predrill=17.0, c_0=2.8,
               unit_wt_method='robertson2009', s_g=2.65, s_g_water=1.0, saturation=None):
    """
    Runs the Boulanger and Idriss (2014) triggering method.

    Parameters
    ----------
    cpt: liquepy.field.CPT,
        ground water level below the surface
    pga: float, g,
        peak ground acceleration
    m_w: float, -,
        Earthquake magnitude
    gwl: float, m,
        depth to ground water from surface at time of earthquake
    p_a: float, kPa, default=101
        Atmospheric pressure
    cfc: float, -, default=0.0
        Fines content correction factor for Eq 2.29
    i_c_limit: float, -, default=2.6
        Limit of liquefiable material
    gamma_predrill: float, kN/m3, default=17.0
        Unit weight of soil above pre-drill depth
    c_0: float, -, default=2.8
        Factor that adjusts the CRR-vs-qc1ncs relationship
    unit_wt_method: str, -, default='robertson2009'
        Method used to determine unit weight
    s_g: float or array_like, -, default=2.65
        Specific gravity
    s_g_water: float, -, default=1.0
        Specific gravity of water
    saturation: array_like or None
        Saturation ratio for each depth increment

    Returns
    -------
    BoulangerIdriss2014CPT()
    """

    return BoulangerIdriss2014CPT(cpt, gwl=gwl, pga=pga, m_w=m_w, cfc=cfc, i_c_limit=i_c_limit, s_g=s_g,
                                  s_g_water=s_g_water, p_a=p_a,
                                  saturation=saturation, unit_wt_method=unit_wt_method, gamma_predrill=gamma_predrill,
                                  c_0=c_0)


def _invert_crr_formula(a, b, c, d, e):

    p = (8 * a * c - 3 * b ** 2) / (8 * a ** 2)
    q = (b ** 3 - 4 * a * b * c + 8 * d * (a ** 2)) / (8 * a ** 3)
    delta_zero = c ** 2 - 3 * b * d + 12 * a * e
    delta_one = 2 * c ** 3 - 9 * b * c * d + 27 * e * (b ** 2) + 27 * a * (d ** 2) - 72 * a * c * e
    big_q = ((delta_one + (delta_one ** 2 - 4 * delta_zero ** 3) ** 0.5) / 2) ** (1/3)
    big_s = 0.5 * (- 2/3 * p + 1/(3 * a) * (big_q + (delta_zero / big_q))) ** 0.5
    big_a = (- 4 * big_s ** 2 - 2 * p + q/big_s)
    big_b = (- 4 * big_s ** 2 - 2 * p - q/big_s)
    big_c = - b/(4 * a)

    # Solutions
    x1 = big_c - big_s + 0.5 * big_a ** 0.5
    # x2 = C - big_s - 0.5 * big_a ** 0.5
    x2 = -1  # q_c1ncs would be less than zero

    import warnings  # These solutions are complex for negative
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        x3 = big_c + big_s + 0.5 * big_b ** 0.5
        x4 = big_c + big_s - 0.5 * big_b ** 0.5

        return np.where(big_b < 0, np.where(x1 < 0, x2, x1), np.where(x3 < 0, x4, x3))


def calc_q_c1n_cs_from_crr_m7p5(crr_7p5, c_0=2.8):
    """
    Solves the closed form solution to a quartic to invert the CRR_7p5-vs-q_c1n_cs relationship

    Parameters
    ----------
    crr_7p5: float or array_like
        values of cyclic resistance ratio at m_w 7.5
    c_0: float (default=2.8)
        Empirical fitting parameter.
         - 2.8=16th percentile (commonly used)
         - 2.6=median response
    Returns
    -------
    float or array_like
        value of normalised cone tip resistance corrected to clean sand behaviour
    """
    a = (1 / 137) ** 4
    b = - (1 / 140) ** 3
    c = (1 / 1000) ** 2
    d = (1 / 113)
    e = - (np.log(crr_7p5) + c_0)
    return _invert_crr_formula(a, b, c, d, e)


def calc_n1_60cs_from_crr_m7p5(crr_7p5, c_0=2.8):
    """
    Solves the closed form solution to a quartic to invert the CRR_7p5-vs-N1_60_cs relationship

    Parameters
    ----------
    crr_7p5: float or array_like
        values of cyclic resistance ratio at m_w 7.5
    c_0: float (default=2.8)
        Empirical fitting parameter.
         - 2.8=16th percentile (commonly used)
         - 2.6=median response
    Returns
    -------
    float or array_like
        value of normalised blow-count corrected to clean sand behaviour
    """
    a = (1 / 25.4) ** 4
    b = - (1 / 23.6) ** 3
    c = (1 / 126) ** 2
    d = (1 / 14.1)
    e = - (np.log(crr_7p5) + c_0)
    return _invert_crr_formula(a, b, c, d, e)


def calc_qc_1ncs_from_crr_m7p5(crr_7p5, c_0=2.8):
    deprecation("Use calc_q_c1n_cs_from_crr_m7p5")
    return calc_q_c1n_cs_from_crr_m7p5(crr_7p5, c_0)


def calculate_qc_1ncs_from_crr_7p5(crr_7p5):
    deprecation("Use calc_q_c1n_cs_from_crr_m7p5")
    return calc_qc_1ncs_from_crr_m7p5(crr_7p5)


def calc_n_cycles_at_m7p5_bi2014(b):
    """
    Number of equivalent cycles for Magnitude 7.5 from Fig A.15 in BI2014

    :param b:
    :return:
    """
    b_vals = [0.060, 0.063, 0.068, 0.072, 0.078, 0.086, 0.095, 0.110, 0.129, 0.148,
              0.168, 0.189, 0.213, 0.235, 0.253, 0.273, 0.310, 0.333, 0.372, 0.400]
    n_cyc = [950.789, 628.731, 405.756, 278.274, 197.931, 121.658, 77.552, 47.081, 30.742, 22.947,
             18.878, 16.503, 15.519, 14.593, 14.235, 14.054, 14.384, 14.374, 14.710, 15.060]
    return np.interp(b, b_vals, n_cyc)


def calc_b_from_msf_max_bi2014(msf_max):
    """
    Equivalent b value for a given magnitude scaling factor from Fig A.16 in BI2014

    Where b is the power coefficient for the CSR-vs-n_cycles liquefaction resistance relationship

    :param msf_max:
    :return:
    """
    b_vals = [0.080, 0.101, 0.121, 0.147, 0.170, 0.191, 0.211, 0.250, 0.290, 0.312, 0.333, 0.353, 0.374, 0.399]
    msf_max_vals = [1.000, 1.019, 1.045, 1.084, 1.130, 1.184, 1.242, 1.373, 1.540, 1.651, 1.768, 1.891, 2.025, 2.215]
    return np.interp(msf_max, msf_max_vals, b_vals)


def calc_lrc_from_qc1ncs_bi2014(q_c1n_cs, n_cycles):
    crr_m7p5 = calc_crr_m7p5_from_qc1ncs(q_c1n_cs)
    msf_m = 1.09 + (q_c1n_cs / 180) ** 3
    msf_max = np.clip(msf_m, None, 2.2)
    b = calc_b_from_msf_max_bi2014(msf_max)
    n_m7p5 = calc_n_cycles_at_m7p5_bi2014(b)
    return (n_m7p5 / n_cycles) ** b * crr_m7p5


def calc_k_sigma_w_n1_60cs(sigma_eff, n1_60cs, pa=100):
    """
    Overburden correction factor, K_sigma

    Equation 2.16a

    Parameters
    ----------
    sigma_eff: float or array_like
        Vertical effective stress
    n1_60cs: float or array_like
        Clean-sand equivalent, normalised SPT blow count
    pa: float
        Atmospheric pressure in kPa

    """
    c_sigma = 1. / (18.9 - 2.55 * np.sqrt(n1_60cs))
    c_sigma = np.clip(c_sigma, None, 0.3)
    return np.clip(1 - c_sigma * np.log(sigma_eff / pa), None, 1.1)


def calc_crr_m7p5_from_n1_60cs(n1_60cs, c_0=2.8):
    """
    Calculation of cyclic resistance ratio for Magnitude 7.5 earthquake from SPT

    Parameters
    ----------
    n1_60cs: float or array_like
        Clean-sand equivalent, normalised SPT blow count
    c_0: float (default=2.8)
        Empirical fitting parameter.
         - 2.8=16th percentile (commonly used)
         - 2.6=median response

    Returns
    -------

    """
    return np.exp((n1_60cs / 14.1) + ((n1_60cs / 126) ** 2) - ((n1_60cs / 23.6) ** 3) + ((n1_60cs / 25.4) ** 4) - c_0)


def calc_crr_n15_from_n1_60cs(n1_60cs, c_0=2.8):
    crr_m7p5 = calc_crr_m7p5_from_n1_60cs(n1_60cs, c_0=c_0)
    msf_m = 1.09 + (n1_60cs / 31.5) ** 2
    msf_max = np.clip(msf_m, None, 2.2)
    b = calc_b_from_msf_max_bi2014(msf_max)
    n_m7p5 = calc_n_cycles_at_m7p5_bi2014(b)
    n_cyc_bi2014 = 15
    msf = (n_m7p5 / n_cyc_bi2014) ** b
    return msf * crr_m7p5

if __name__ == '__main__':
    n1_60_cs_values = 3.
    crr_values = calc_crr_m7p5_from_n1_60cs(n1_60_cs_values, c_0=2.6)
    n1_60_cs_back = calc_n1_60cs_from_crr_m7p5(crr_values, c_0=2.6)
    print(n1_60_cs_back)