import numpy as np


def calc_shear_vel_mcgann_2015_cpt(cpt):
    """
    Computes the shear wave velocity profile according to :cite:`McGann:2015fd`

    Parameters
    ----------
    cpt: liquepy.field.CPT object

    Returns
    -------
    array_like
        Shear wave velocity profile corresponding to the depths in the CPT.
    """
    return 18.4 * cpt.q_c ** 0.144 * cpt.f_s ** 0.0832 * cpt.depth ** 0.278


def calc_shear_vel_andrus_and_stokoe_2000_spt_values(n_1_60):
    """
    Eq 98 in Pm4Sand v3.1 manual

    :param n_1_60:
    :return: float or array_like
    """
    return 85.8 * (n_1_60 + 2.5) ** 0.25


def calc_relative_density_idriss_and_boulanger_2008_spt_values(n_1_60, c_d=46):
    """
    Calculate the relative density from SPT normalised blow count (Eq. 35 :cite:`Idriss:2008ua`)

    Parameters
    ----------
    n_1_60: array_like
        Corrected normalised SPT values
    c_d: float
        Correlation factor, default=46 (from :cite:`Idriss:2008ua`)
        other proposals:

         * Meyerhof (1957): 41
         * Skempton (1986): 55 fine and 65 coarse natural norm-consolidated sand, 35 lab, 40 recent fills
         * Cubrinovski and Ishihara (1999): 51 clean sand, 26 silty sand, 39 all samples

    Returns
    -------
    array_like
    """
    return np.sqrt(n_1_60 / c_d)


def calc_relative_density_salgado_et_al_1997_cpt_values(q_c1n, c_dq=0.9):
    """
    Eq. 95 in PM4Sand v3.1 manual

    Parameters
    ----------
    q_c1n: array_like
        Normalised cone penetration resistance
    c_dq: float, default=0.9 (from :cite:`Idriss:2008ua`)
        Correlation factor, (range 0.64-155 from Salgado (1997)

    Returns
    -------
    array_like
    """
    return 0.465 * np.sqrt(q_c1n / c_dq) - 1.063


def calc_relative_density_boulanger_et_al_2014_cpt_values(q_c1n, c_dq=0.9):
    """
    Table 4.1 in PM4Sand v3.1 manual

    Parameters
    ----------
    q_c1n: array_like
        Normalised cone penetration resistance
    c_dq: float, default=0.9 (from :cite:`Idriss:2008ua`)
        Correlation factor, (range 0.64-1.55 from Salgado (1997)

    Returns
    -------
    array_like
    """
    return 0.465 * (q_c1n / c_dq) ** 0.264 - 1.063


def calc_q_c1n_via_inverse_d_r_boulanger_et_al_2014(d_r, c_dq=0.9):
    """
    Table 4.1 in PM4Sand v3.1 manual

    Parameters
    ----------
    d_r: array_like
        Relative density
    c_dq: float, default=0.9 (from :cite:`Idriss:2008ua`)
        Correlation factor, (range 0.64-155 from Salgado (1997)

    Returns
    -------
    array_like
    """
    return c_dq * ((d_r + 1.063) / 0.465) ** (1. / 0.264)


def calc_g0_mod_boulanger_and_ziotopoulou_2015_spt_values(n_1_60):
    """
    Calculate the normalised shear modulus :cite:`Boulanger:2017pm4_v3p1`

    Parameters
    ----------
    n_1_60: array_like
        Corrected normalised SPT values

    Returns
    -------
    array_like
    """
    return 167. * np.sqrt(n_1_60 + 2.5)


def est_permeability_robertson_and_cabal_2012(i_c):
    """
    Estimates the soil permeability based on the soil behaviour index :cite: `Robertson:2012cpt`

    Parameters
    ----------
    i_c: float or array_like

    Returns
    -------

    """
    if np.min(i_c) < 1:
        raise ValueError
    return np.where(i_c < 3.27, 10.0 ** (0.952 - 3.04 * i_c), 10.0 ** (-4.52 - 1.37 * i_c))


def est_shear_vel_hegazy_and_mayne_2006(q_c1n, i_c, esig_v0, p_a):
    """
    Estimates the soil shear wave velocity from CPT :cite: `Hegazy:2012bs`

    Eq 6.

    Parameters
    ----------
    q_c1n
    i_c
    esig_v0
    p_a

    Returns
    -------

    """
    return 0.0831 * q_c1n * np.exp(1.7861 * i_c) * (esig_v0 / p_a) ** 0.25


def est_g_mod_robertson_2009(i_c, big_q, unit_weight):
    """Set normalised shear modulus using Robertson (2009)."""
    alpha_vs = 10 * (0.55 * i_c + 1.68)  # Eq 11 [m/s]
    big_q = np.clip(big_q, 0.001, None)
    vs1 = (alpha_vs * big_q) ** 0.5  # Eq 9 [m/s]
    return unit_weight * vs1 ** 2


def est_g0_mod_robertson_2009(i_c, big_q, unit_weight, esig_v, pa=101000, n=0.5):
    """Set normalised shear modulus using Robertson (2009)."""
    g_mod = est_g_mod_robertson_2009(i_c, big_q, unit_weight)
    return g_mod / pa / (esig_v / pa) ** n


def est_undrained_strength_ratio_robertson_2009(big_q, n_kt=14):
    """determine normalised undrained strength using Robertson (2009)."""
    return big_q / n_kt  # =(su / esig_v)  Eq 33


def set_strength_props(sl, vert_eff_stress, i_c, big_q, n_kt=14):
    if i_c > 2.6:
        sl.cohesion = est_undrained_strength_ratio_robertson_2009(big_q, n_kt=n_kt) * vert_eff_stress
        sl.phi = 0.0
    else:
        sl.phi = np.arctan(est_undrained_strength_ratio_robertson_2009(big_q, n_kt=n_kt))
        sl.cohesion = 0.0