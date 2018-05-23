import numpy as np
from liquepy.trigger import load_cpt_file


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


def crr_7p5_from_cpt(q_c1n_cs, gwl, depth, i_c):
    """
    cyclic resistance from CPT, Eq. 2.24
    it's not possible to have liquefaction up water table
    """
    crr_values = np.exp((q_c1n_cs / 113) + ((q_c1n_cs / 1000) ** 2) -
                        ((q_c1n_cs / 140) ** 3) + ((q_c1n_cs / 137) ** 4) - 2.8)
    crr_tent = np.where(depth < gwl, 4, crr_values)
    return np.where(i_c <= 2.6, crr_tent, 4.)


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


def calculate_msf(magnitude, q_c1ns):
    n = len(q_c1ns)
    msf_max = np.ones(n)
    msf_m = np.ones(n)
    msf = np.ones(n)
    if magnitude == 7.5:
        return msf
    else:
        for i in range(0, n):
            msf_m[i] = 1.09 + (q_c1ns[i] / 180) ** 3
            if msf_m[i] > 2.2:
                msf_max[i] = 2.2
            else:
                msf_max[i] = msf_m[i]
            msf[i] = 1 + (msf_max[i] - 1) * (8.64 * np.exp(-magnitude / 4) - 1.325)
    return msf


def calculate_k_sigma(sigma_eff, qc1ncs):
    n = len(qc1ncs)
    q = np.where(qc1ncs > 211, 211, qc1ncs)

    pa = 100  # kPa
    k = np.ones(n)
    k1 = np.ones(n)
    cs = np.ones(n)
    for i in range(0, n):
        cs[i] = (37.3 - 8.27 * (q[i] ** 0.264)) ** -1
        k1[i] = 1 - cs[i] * np.log(sigma_eff[i] / pa)
        if k1[i] > 1.1:
            k[i] = 1.1
        else:
            k[i] = k1[i]
    return k


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
        for j in range(100):
            cn_values[dd] = min((p_a / sigma_veff[dd]) ** m_values[dd], 1.7)  # Eq 2.15a
            temp_cn = cn_values[dd]
            q_c1n[dd] = (cn_values[dd] * q_c[dd] / p_a)  # Eq. 2.4
            big_q[dd] = calculate_big_q_values(cn_values[dd], q_t[dd], sigma_v[dd])
            ft_values[dd] = calculate_f_ic_values(f_s[dd], q_t[dd], sigma_v[dd])
            i_c[dd] = calculate_ic(big_q[dd], ft_values[dd])
            fines_content[dd] = calculate_fc(i_c[dd], cfc)

            delta_q_c1n[dd] = calculate_delta_q_c1n(q_c1n=q_c1n[dd], fc=fines_content[dd])  # Eq. 2.22
            q_c1n_cs[dd] = q_c1n[dd] + delta_q_c1n[dd]
            m_values[dd] = calculate_m(q_c1n_cs[dd])
            if abs(q_c1n[dd] - temp_cn) < 0.00001:
                break
    return q_c1n_cs, fines_content, i_c


class BoulangerIdriss2014(object):
    def __init__(self, depth, q_c, f_s, u_2, gwl=2.3, pga=0.25, magnitude=7.5, ar=0.8):
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
        :param ar: float, -, Area ratio
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
        self.ar = ar

        self.cfc = 0  # parameter of fines content, eq 2.29
        self.q_t = calculate_qt(self.q_c, self.ar, self.u_2)  # kPa
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
        self.crr_m7p5 = crr_7p5_from_cpt(self.q_c1n_cs, gwl, depth, self.i_c)
        self.crr = crr_m(self.k_sigma, self.msf, self.crr_m7p5)  # CRR a magnitudo M
        fs_unlimited = self.crr / self.csr
        fs_fines_limited = np.where(self.fines_content > 71, 2.0, fs_unlimited)
        self.factor_of_safety = np.where(fs_fines_limited > 2, 2, fs_fines_limited)


def run_standard_bi2014(cpt_file_path):
    depths, q_c, f_s, u_2, gwl = load_cpt_file.load_cpt_data(cpt_file_path)
    return BoulangerIdriss2014(depths, q_c, f_s, u_2, gwl=gwl, pga=0.25, magnitude=7.5, ar=0.8)
