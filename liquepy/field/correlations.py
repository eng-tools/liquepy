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



