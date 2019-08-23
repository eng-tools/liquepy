import numbers
import numpy as np
from liquepy.exceptions import deprecation


def calc_shear_strain_zhang_2004(fs, d_r):
    """
    Calculates the shear strain (not in percentage)

    Parameters
    ----------
    fs
    d_r: float or array_like, [dec]
        Relative density of soil

    Returns
    -------

    """
    if not hasattr(fs, '__len__'):
        if not hasattr(d_r, '__len__'):
            return calc_single_shear_strain(fs, d_r)
        else:
            out_values = []
            for i in range(len(d_r)):
                out_values.append(calc_single_shear_strain(fs, d_r[i]))
            return np.array(out_values)
    else:
        out_values = []
        for i in range(len(fs)):
            if not hasattr(d_r, '__len__'):
                out_values.append(calc_single_shear_strain(fs[i], d_r))
            else:
                out_values.append(calc_single_shear_strain(fs[i], d_r[i]))
        return np.array(out_values)


def calc_shear_strain(fs, d_r):
    deprecation("Use calc_shear_strain_zhang_2004, note that new function returns strain not in percentage!")
    return None


def calc_relative_density_tasuoka_1990(q_c, esig_v0, dr_min=0.0, dr_max=1.0):
    """
    Calculates the relative density according to Tasuoka et al. (1990)

    Parameters
    ----------
    q_c1n
    esig_v0
    dr_min
    dr_max

    Returns
    -------

    """
    return np.clip((-85. + 76. * np.log10(q_c / np.sqrt(esig_v0))), dr_min, dr_max) / 1e2


def calc_relative_density_zhang_2002(q_c1n, dr_min=0.0, dr_max=1.0):
    """
    Calculates the relative density (Eq. 2) :cite:`Zhang:2004el`

    Parameters
    ----------
    q_c1n: Normalised cone tip resistance

    Returns
    -------

    """
    return np.clip((-85. + 76. * np.log10(q_c1n)) / 100, dr_min, dr_max)


def calculate_shear_strain(fos, relative_density):
    deprecation("Use calc_shear_strain")
    return calculate_shear_strain(fos, relative_density)


def calc_single_shear_strain(fs, d_r):
    if d_r == -1:
        return 0
    if d_r > 2.:
        raise ValueError('d_r should be a decimal not a percentage')
    dr_values = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
    for i in range(len(dr_values)):
        if d_r < dr_values[i]:
            if i == 0:
                low_dr = dr_values[0]
                d_r_cur = dr_values[0]
            else:
                low_dr = dr_values[i - 1]
                d_r_cur = d_r
            high_dr = dr_values[i]
            es_low = calc_fixed_dr_gamma_max(fs, low_dr)
            es_high = calc_fixed_dr_gamma_max(fs, high_dr)

            es_actual = np.interp(d_r_cur, [low_dr, high_dr], [es_low, es_high])
            return es_actual
    return 0.0


# Maximum cyclic shear strains
def calc_fixed_dr_gamma_max(fs, relative_density):
    if fs > 2.0:
        return 0.0
    elif relative_density == 1.0:
        return 0.0
    elif relative_density == 0.9:
        if fs >= 0.7:
            gamma_max = 3.26 * fs ** (-1.80)
        else:
            gamma_max = 6.2
    elif relative_density == 0.8:
        if fs >= 0.56:
            gamma_max = 3.22 * fs ** (-2.08)
        else:
            gamma_max = 10
    elif relative_density == 0.7:
        if fs >= 0.59:
            gamma_max = 3.2 * fs ** (-2.89)
        else:
            gamma_max = 14.5
    elif relative_density == 0.6:
        if fs >= 0.66:
            gamma_max = 3.58 * fs ** (-4.42)
        else:
            gamma_max = 22.7
    elif relative_density == 0.5:
        if fs >= 0.72:
            gamma_max = 4.22 * fs ** (-6.39)
        else:
            gamma_max = 34.1
    elif relative_density == 0.4:
        if fs >= 1:
            gamma_max = 3.31 * fs ** (-7.97)
        elif fs >= 0.81:
            gamma_max = 250 * (1 - fs) + 3.5
        else:
            gamma_max = 51.2
    elif relative_density < 0.4:
        gamma_max = 51.2
    else:
        raise ValueError("Relative density (%.4f) not set to standard value" % relative_density)
    return gamma_max / 100
