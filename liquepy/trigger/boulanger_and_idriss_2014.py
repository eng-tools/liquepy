import numpy as np


def calculate_unit_weight(fs, q_t, gwl, depth):
    # eq Robertson pag 37- CPT guide

    pa = 100  # kPa
    n = len(fs)

    gamma_soil = np.zeros(n)  # kN/m^3
    for i in range(0, n):
        if q_t[i] <= 0:
            q_t[i] = zeros_search(q_t, q_t[i])
        r_f = (fs[i] / q_t[i]) * 100
        if r_f <= 0:
            r_f = 0.1
        if depth[i] <= gwl:
            gamma_soil[i] = gamma_dry(q_t[i], r_f, pa)
        else:
            gamma_soil[i] = gamma_wet(q_t[i], r_f, pa)

    return gamma_soil


def gamma_dry(q_t, r_f, pa):
    gamma_water = 9.81
    g = (0.27 * np.log10(r_f) + 0.36 * np.log10(q_t / pa) + 1.236) * gamma_water
    if g < 15:  # for void ratio is 12.985
        gamma_dry = 15
    else:
        gamma_dry = g
    return gamma_dry


def gamma_wet(q_t, r_f, pa):
    gamma_water = 9.81
    g = (0.27 * np.log10(r_f) + 0.36 * np.log10(q_t / pa) + 1.236) * gamma_water
    if g < 15:
        gamma_wet = 15
    else:
        gamma_wet = g
    return gamma_wet


def zeros_search(qt, b):  # it needs to search zeros in qt and so it transorm in number minimum different by 0
    n = len(qt)
    a = np.zeros(n)
    p = max(qt) + 10
    for k in range(0, n):
        if qt[k] == b:
            for i in range(0, n):
                a[i] = qt[i]
                if a[i] <= 0:
                    a[i] = p
            for j in range(0, n):
                if a[j] == p:
                    a[j] = min(a)
    return min(a)


def calculate_sigma_v(depths, gammas):
    """
    Calculates the vertical stress
    """
    depth_incs = depths[1:] - depths[:-1]
    depth_incs = np.insert(depth_incs, 0, depth_incs[0])
    sigma_v_incs = depth_incs * gammas
    sigma_v = np.cumsum(sigma_v_incs)
    return sigma_v


def calculate_pore_pressure(depth, gwl):
    gamma_water = 10
    pore_pressure = np.where(depth > gwl, (depth - gwl) * gamma_water, 0.0)
    return pore_pressure


def calculate_sigma_veff(sigma_v, pore_pressure):
    sigma_veff = abs(sigma_v - pore_pressure)
    return sigma_veff


def calculate_qt(qc, ar, u2):
    """

    :param qc: kPa, cone tip resistance
    :param ar: -, area ratio
    :param u2: kPa, water pressure beneath cone tip
    :return:
    """
    # qt the cone tip resistance corrected for unequal end area effects, eq 2.3
    qt = qc + ((1 - ar) * u2)
    return qt


def calculate_f_ic_values(fs, qt, sigmav):
    # qt is in kPa, so it's not necessary measure unit transormation
    F = (fs / (qt - sigmav)) * 100
    return F


def calculate_big_q_values(CN, qt, sigmav):
    """
    Eq. XXXXX
    :param CN:
    :param qt:
    :param sigmav:
    :return:
    """
    Q = (qt - sigmav) / 100 * CN
    return Q


def calculate_ic(big_q, big_f):
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


def calculate_fc(ic, cfc):
    fc_tent = 80 * (ic + cfc) - 137
    fc1 = np.where(fc_tent > 100, 100, fc_tent)
    fc = np.where(fc1 <= 137 / 80 - cfc, 0, fc1)
    return fc


def calculate_rd(depth, magnitude):
    """
    rd from CPT, Eq 2.14a
    """
    alpha = -1.012 - 1.126 * np.sin((depth / 11.73) + 5.133)
    beta = 0.106 + 0.118 * np.sin((depth / 11.28) + 5.142)
    rd = np.exp(alpha + beta * magnitude)
    return rd


def calculate_csr(sigma_veff, sigma_v, pga, rd, gwl, depth):
    """
    Cyclic stress ratio from CPT, Eq 2.2,
    """
    return np.where(depth <= gwl, 2, 0.65 * (sigma_v / sigma_veff) * rd * pga)


def calculate_cn_values(m, sigma_veff):
    """
    CN parameter from CPT, Eq 2.15a
    """
    CN = (100 / sigma_veff) ** m
    if CN > 1.7:
        CN = 1.7
    return CN


def calculate_m(q_c1ncs):
    """
    m parameter from CPT, Eq 2.15b
    """
    m = 1.338 - 0.249 * q_c1ncs ** 0.264
    if q_c1ncs >= 254:
        m = 0.263823991466759
    if q_c1ncs <= 21:
        m = 0.781756126201587
    return m


def calculate_q_c1n(qc, CN):
    """
    qc1n from CPT, Eq 2.4
    """
    q_c1n = CN * qc * 1000 / 100
    return q_c1n


def crr_7p5_from_cpt(q_c1n_cs, gwl, depth, i_c, i_c_limit=2.6):
    """
    cyclic resistance from CPT, Eq. 2.24
    it's not possible to have liquefaction up water table
    """
    crr_values = np.exp((q_c1n_cs / 113) + ((q_c1n_cs / 1000) ** 2) -
                        ((q_c1n_cs / 140) ** 3) + ((q_c1n_cs / 137) ** 4) - 2.8)
    crr_tent = np.where(depth < gwl, 4, crr_values)
    return np.where(i_c <= i_c_limit, crr_tent, 4.)


def crr_m(ksigma, msf, crr_m7p5):
    """Cyclic resistance ratio corrected for magnitude and confining stress"""
    return crr_m7p5 * ksigma * msf


def calculate_delta_q_c1n(q_c1n, fc):
    """
    delta q_c1n from CPT, Eq 2.22
    """
    delta_q_c1n = (11.9 + (q_c1n / 14.6)) * (np.exp(1.63 - (9.7 / (fc + 2)) - ((15.7 / (fc + 2)) ** 2)))
    return delta_q_c1n


def calculate_q_c1ncs(q_c1n, delta_q_c1n):
    """
    q_c1ncs from CPT, Eq 2.10
    """
    q_c1ncs = q_c1n + delta_q_c1n
    return q_c1ncs


def calculate_msf(magnitude, q_c1ncs):
    """
    Magnitude scaling factor to correct the cyclic resistance ratio

    Eq. 2.19

    :param magnitude: earthquake magnitude
    :param q_c1ncs: clean sand-corrected normalised cone tip resistance
    :return:
    """
    if magnitude == 7.5:
        return np.ones_like(q_c1ncs)
    msf_m = 1.09 + (q_c1ncs / 180) ** 3
    msf_max = np.where(msf_m > 2.2, 2.2, msf_m)
    msf = 1. + (msf_max - 1) * (8.64 * np.exp(-magnitude / 4) - 1.325)
    return msf


def calculate_k_sigma(sigma_eff, qc1ncs, pa=100):
    """
    Overburden correction factor, K_sigma

    Equation 2.16a

    :param sigma_eff: vertical effective stress
    :param qc1ncs: clean sand-corrected normalised cone tip resistance
    :param pa: atmospheric pressure in kPa
    :return:
    """
    c_sigma_unrestricted = 1. / (37.3 - 8.27 * (qc1ncs ** 0.264))
    c_sigma = np.where(c_sigma_unrestricted <= 0.3, c_sigma_unrestricted, 0.3)
    k_sigma_unrestricted = 1 - c_sigma * np.log(sigma_eff / pa)
    k_sigma = np.where(k_sigma_unrestricted <= 1.1, k_sigma_unrestricted, 1.1)
    return k_sigma


# def slow_calculate_dependent_variables(sigma_v, sigma_veff, q_c, f_s, p_a, q_t, cfc):
#     """
#     Iteratively calculate_volumetric_strain parameters as they are interdependent
#
#     :param sigma_v: array, kPa, Total vertical stress
#     :param sigma_veff: array, kPa, Effective vertical stress
#     :param q_c: array, kPa, Cone tip resistance
#     :param f_s: array, kPa, Skin friction
#     :param p_a: array, kPa, Atmospheric pressure
#     :param q_t:
#     :param cfc: float, -, correction factor
#     :return:
#     """
#     num_depth = len(sigma_v)
#     m_values = np.ones(num_depth)  # create an array of ones
#     cn_values = np.zeros(num_depth)
#     q_c1n = np.zeros(num_depth)
#     delta_q_c1n = np.zeros(num_depth)
#     q_c1n_cs = np.zeros(num_depth)
#
#     big_q = np.ones(num_depth)
#     ft_values = np.ones(num_depth)
#     i_c = np.ones(num_depth)
#     fines_content = np.ones(num_depth)
#
#     for dd in range(0, num_depth):
#         for j in range(100):
#             cn_values[dd] = min((p_a / sigma_veff[dd]) ** m_values[dd], 1.7)  # Eq 2.15a
#             temp_cn = cn_values[dd]
#             q_c1n[dd] = (cn_values[dd] * q_c[dd] / p_a)  # Eq. 2.4
#             big_q[dd] = calculate_big_q_values(cn_values[dd], q_t[dd], sigma_v[dd])
#             ft_values[dd] = calculate_f_ic_values(f_s[dd], q_t[dd], sigma_v[dd])
#             i_c[dd] = calculate_ic(big_q[dd], ft_values[dd])
#             fines_content[dd] = calculate_fc(i_c[dd], cfc)
#
#             delta_q_c1n[dd] = calculate_delta_q_c1n(q_c1n=q_c1n[dd], fc=fines_content[dd])  # Eq. 2.22
#             q_c1n_cs[dd] = q_c1n[dd] + delta_q_c1n[dd]
#             m_values[dd] = calculate_m(q_c1n_cs[dd])
#             if abs(q_c1n[dd] - temp_cn) < 0.00001:
#                 break
#     return q_c1n_cs, fines_content, i_c


def calculate_dependent_variables(sigma_v, sigma_veff, q_c, f_s, p_a, q_t, cfc):
    """
    Iteratively calculate_volumetric_strain parameters as they are interdependent

    :param sigma_v: array, kPa, Total vertical stress
    :param sigma_veff: array, kPa, Effective vertical stress
    :param q_c: array, kPa, Cone tip resistance
    :param f_s: array, kPa, Skin friction
    :param p_a: array, kPa, Atmospheric pressure
    :param q_t:
    :param cfc: float, -, correction factor
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
        for j in range(100):
            cn_values[dd] = min((p_a / sigma_veff[dd]) ** m_values[dd], 1.7)  # Eq 2.15a
            q_c1n[dd] = (cn_values[dd] * q_c[dd] / p_a)  # Eq. 2.4
            big_q[dd] = calculate_big_q_values(cn_values[dd], q_t[dd], sigma_v[dd])
            ft_values[dd] = calculate_f_ic_values(f_s[dd], q_t[dd], sigma_v[dd])
            i_c[dd] = calculate_ic(big_q[dd], ft_values[dd])
            fines_content[dd] = calculate_fc(i_c[dd], cfc)

            delta_q_c1n[dd] = calculate_delta_q_c1n(q_c1n=q_c1n[dd], fc=fines_content[dd])  # Eq. 2.22
            q_c1n_cs[dd] = q_c1n[dd] + delta_q_c1n[dd]
            m_values[dd] = calculate_m(q_c1n_cs[dd])
            if abs(q_c1n[dd] - temp_q_c1n) < 0.00001:
                break
            temp_q_c1n = q_c1n[dd]
    return q_c1n_cs, fines_content, i_c


class BoulangerIdriss2014(object):
    def __init__(self, depth, q_c, f_s, u_2, gwl=2.3, pga=0.25, magnitude=7.5, a_ratio=0.8, cfc=0.0, i_c_limit=2.6):
        """
        Performs the Boulanger and Idriss triggering procedure for a CPT profile

        ref: Boulanger:2014id

        :param depth: array, m, depths measured downwards from surface
        :param q_c: array, kPa, cone tip resistance
        :param f_s: array, kPa, skin friction
        :param u_2: array, water pressure beneath cone tip
        :param gwl: float, m, ground water level below the surface
        :param pga: float, g, peak ground acceleration
        :param magnitude: float, -, Earthquake magnitude
        :param a_ratio: float, -, Area ratio
        :param cfc: float, -, Fines content correction factor for Eq 2.29
        :param i_c_limit: float, -, limit of liquefiable material
        :return:
        """
        p_a = 100.  # kPa
        self.depth = depth
        self.q_c = q_c
        self.f_s = f_s
        self.u_2 = u_2
        self.gwl = gwl
        self.pga = pga
        self.magnitude = magnitude
        self.a_ratio = a_ratio
        if a_ratio is None:
            self.a_ratio = 0.8
        self.i_c_limit = i_c_limit

        self.cfc = cfc  # parameter of fines content, eq 2.29
        self.q_t = calculate_qt(self.q_c, self.a_ratio, self.u_2)  # kPa
        self.gammas = calculate_unit_weight(self.f_s, self.q_t, gwl, self.depth)
        self.sigma_v = calculate_sigma_v(self.depth, self.gammas)
        self.pore_pressure = calculate_pore_pressure(self.depth, self.gwl)
        self.sigma_veff = calculate_sigma_veff(self.sigma_v, self.pore_pressure)
        self.rd = calculate_rd(depth, magnitude)

        self.q_c1n_cs, self.fines_content, self.i_c = calculate_dependent_variables(self.sigma_v, self.sigma_veff, q_c,
                                                                                    f_s, p_a,
                                                                                    self.q_t,
                                                                                    self.cfc)

        self.k_sigma = calculate_k_sigma(self.sigma_veff, self.q_c1n_cs)
        self.msf = calculate_msf(magnitude, self.q_c1n_cs)
        self.csr = calculate_csr(self.sigma_veff, self.sigma_v, pga, self.rd, gwl, depth)
        self.crr_m7p5 = crr_7p5_from_cpt(self.q_c1n_cs, gwl, depth, self.i_c, self.i_c_limit)
        self.crr = crr_m(self.k_sigma, self.msf, self.crr_m7p5)  # CRR a magnitudo M
        fs_unlimited = self.crr / self.csr
        # fs_fines_limited = np.where(self.fines_content > 71, 2.0, fs_unlimited)  # based on I_c=2.6
        self.factor_of_safety = np.where(fs_unlimited > 2, 2, fs_unlimited)


def run_bi2014(cpt, pga, magnitude, cfc=0.0, i_c_limit=2.6):
    return BoulangerIdriss2014(cpt.depth, cpt.q_c, cpt.f_s, cpt.u_2, gwl=cpt.gwl, pga=pga, magnitude=magnitude,
                               a_ratio=cpt.a_ratio, cfc=cfc, i_c_limit=i_c_limit)


def calculate_qc_1ncs_from_crr_7p5(crr_7p5):
    """
    Solves the closed form solution to a quartic to invert the CRR_7p5-vs-q_c1n_cs relationship

    :param crr_7p5: float or array, values of cyclic resistance ratio at magnitude 7.5
    :return: float or array, value of normalised cone tip resistance corrected to clean sand behaviour
    """
    x = 5
    a = (1 / 137) ** 4
    b = - (1 / 140) ** 3
    c = (1 / 1000) ** 2
    d = (1 / 113)
    e = - (np.log(crr_7p5) + 2.8)

    # discriminant is unused
    discriminant = 256 * (a ** 3) * (e ** 3) - 192 * ((a * e) ** 2) * b * d - 128 * (a * c * e) **2 + 144 *\
                   ((a * d) ** 2) * c * e - 27 * (a * d) ** 4 + 144 * a * c * ((b * e) ** 2) - 6 * a * e * \
                   ((b * d) ** 2) - 80 * a * b * d * e * (c ** 2) + 18 * a * b * c * (d ** 3) + 16 * a * e * \
                   (c ** 4) - 4 * a * (c ** 3) * (d ** 2) - 27 * (b ** 4) * (e ** 2) + 18 * (b ** 3) * c * d * e - 4 * \
                   (b ** 3) * (d ** 3) - 4 * (b ** 2) * (c ** 3) * e + (b * c * d) ** 2

    p = (8 * a * c - 3 * b ** 2) / (8 * a ** 2)
    q = (b ** 3 - 4 * a *b *c + 8 * d * (a ** 2)) / (8 * a ** 3)
    delta_zero = c ** 2 - 3 * b * d + 12 * a * e
    delta_one = 2 * c ** 3 - 9 * b * c * d + 27 * e * (b ** 2) + 27 * a * (d ** 2) - 72 * a * c * e
    big_q = ((delta_one + (delta_one ** 2 - 4 * delta_zero ** 3) ** 0.5) / 2) ** (1/3)
    big_s = 0.5 * (- 2/3 * p + 1/(3 * a) * (big_q + (delta_zero / big_q))) ** 0.5
    big_a = (- 4 * big_s ** 2 - 2 * p + q/big_s)
    big_b = (- 4 * big_s ** 2 - 2 * p - q/big_s)
    C = - b/(4 * a)

    # Solutions
    x1 = C - big_s + 0.5 * big_a ** 0.5
    # x2 = C - big_s - 0.5 * big_a ** 0.5
    x2 = -1  # q_c1ncs would be less than zero

    import warnings  # These solutions are complex for negative
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        x3 = C + big_s + 0.5 * big_b ** 0.5
        x4 = C + big_s - 0.5 * big_b ** 0.5

    return np.where(big_b < 0, np.where(x1 < 0, x2, x1), np.where(x3 < 0, x4, x3))
