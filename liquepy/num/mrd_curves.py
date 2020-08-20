import numpy as np


def calc_mr_elastic(gamma):
    return np.ones(len(gamma))


def calc_mr_mohr_coloumb(gamma, elastic_mod, cu):
    tau = np.clip(elastic_mod * gamma, None, cu)
    secant_mod_reduction = (tau / gamma) / elastic_mod
    return secant_mod_reduction


def calc_mr_flac_default_curve(gamma, l_1, l_2):
    gamma_percent = gamma * 100
    lg = np.log10(gamma_percent)
    s = np.clip((l_2 - lg) / (l_2 - l_1), None, 1)
    return s ** 2 * (3. - 2. * s)


def calc_mr_seed_and_sun_mod(gamma):
    l_1 = -3.156
    l_2 = 1.904
    return calc_mr_flac_default_curve(gamma, l_1, l_2)


def calc_mr_daren_mod(gamma, a=0.92, gamma_ref=0.001):
    return 1.0 / (1.0 + (gamma / gamma_ref) ** a)


def calc_mr_vardanega_2013_mod(gamma, i_p):
    a = 0.943  # Eq 22b
    j = 3.7  # Eq 23
    gamma_ref = j * (i_p / 1000)
    return 1.0 / (1.0 + (gamma / gamma_ref) ** a)  # Eq. 22b


def calc_ss_ratio_from_ip_vardanega_2013(i_p, gamma_target=0.005):
    """
    Calculates an appropriate strength-stiffness ratio for a soil based on the plasticity index.
    Parameters
    ----------
    i_p
    gamma_target

    Returns
    -------

    """
    vardanega_mr = calc_mr_vardanega_2013_mod(gamma_target, i_p=i_p)
    return vardanega_mr * gamma_target


def calc_ip_from_ss_ratio_vardanega_2013(ss_ratio, gamma_target=0.005):
    """
    Calculates the appropriate plasticity index for a soil based on the strength-stiffness ratio (ss_ratio).

    The strength stiffness ratio uses the undrained strength and the initial stiffness.

    Parameters
    ----------
    ss_ratio: float
        Strength-stiffness ratio
    gamma_target: float
        Strain where backbone response should reach undrained strength
    Returns
    -------

    """
    ips = np.arange(0.05, 1.05, 0.05)
    vardanega_mr = calc_mr_vardanega_2013_mod(gamma_target, i_p=ips)
    ss_ratios = vardanega_mr * gamma_target
    return np.interp(ss_ratio, ss_ratios, ips)


def set_hyp_params_from_op_pimy_or_pdmy_model(sl, p_ref=100.0e3, hyp=True):
    # Octahedral shear stress
    tau_f = (2 * np.sqrt(2.) * np.sin(sl.phi_r)) / (3 - np.sin(sl.phi_r)) * p_ref + 2 * np.sqrt(2.) / 3 * sl.cohesion
    if hasattr(sl, 'get_g_mod_at_m_eff_stress'):
        g_mod_r = sl.get_g_mod_at_m_eff_stress(p_ref)
        if hasattr(sl, 'g_mod_p0'):
            assert sl.g_mod_p0 == 0.0
        d = sl.a
    else:
        g_mod_r = sl.g_mod
        d = 0.0
    print('tau_f: ', tau_f)
    print('cohesion: ', sl.cohesion)
    strain_r = sl.peak_strain * tau_f / (g_mod_r * sl.peak_strain - tau_f)
    sdf = (p_ref / p_ref) ** d
    if hyp:  # hyperbolic model parameters
        sl.strain_curvature = 1.0
        sl.xi_min = 0.01
        dss_eq = 1.  # np.sqrt(3. / 2)  # correct to direct simple shear equivalent
        sl.strain_ref = strain_r / sdf / dss_eq
        sl.sra_type = "hyperbolic"
        sl.inputs += ['strain_curvature', 'xi_min', 'sra_type', 'strain_ref']


def calc_gamma_ref_via_darendeli_2001(i_p, ocr, p_eff, p_atm=101.0e3):
    """
    Calculate the reference strain for a modulus reduction curve from :cite:`Darendeli:2001thesis`

    :param i_p: float
        Plasticity index (as a ratio not percentage)
    :param ocr:
    :param p_eff:
    :return:
    """
    phi_1 = 0.0352
    phi_2 = 0.0010
    phi_3 = 0.3246
    phi_4 = 0.3483
    return (phi_1 + phi_2 * i_p * 100 * ocr ** phi_3) * (p_eff / p_atm) ** phi_4 / 100


def calc_gamma_ref_via_menq_2003(c_u, p_eff, p_atm=101.0e3):
    """
    Calculate the reference strain for a modulus reduction curve from Menq

    :param c_u: float
        Uniformity coefficient (Grain size ratio d_60 / d_10)
    :param p_eff:
    :return:
    """
    big_a_gamma = 0.12 * c_u ** -0.6
    n_gamma = 0.5 * c_u ** -0.15

    return big_a_gamma * (p_eff / p_atm) ** n_gamma / 100


def calc_min_damping_via_menq_2003(c_u, d_50, p_eff, p_atm=101.0e3):
    """
    Minimum damping for hyperbolic model according to Menq 2003
    Parameters
    ----------
    c_u: float
        Uniformity coefficient (Grain size ratio d_60 / d_10)
    d_50: float
        Median grain diameter in mm
    p_eff
    p_atm

    Returns
    -------

    """

    return (0.55 * c_u ** 0.1 * d_mean ** -0.3 * (p_eff / p_atm) ** -0.08) / 100


def calc_alpha_via_menq_2003(p_eff, p_atm=101.0e3):
    """
    Calculate the reference strain for a modulus reduction curve from Menq

    :param p_eff:
    :return:
    """
    return 0.86 + 0.1 * np.log(p_eff / p_atm)
