import numbers
import numpy as np


def calculate_volumetric_strain(factor_of_safety, q_c1n_cs):
    """
    Calculates the Volumetric strain according to Zhang et al. (2002)

    doi: 10.1139/t02-047

    :param factor_of_safety: float or array, factor of safety against liquefaction triggering
    :param q_c1n_cs: float, or array, corrected normalised clean sand cone tip resistance
    :return:
    """
    if isinstance(factor_of_safety, numbers.Real) and isinstance(q_c1n_cs, numbers.Real):
        return _calculate_single_volumetric_strain(factor_of_safety, q_c1n_cs)
    elif not isinstance(q_c1n_cs, numbers.Real) and not isinstance(factor_of_safety, numbers.Real):
        assert len(factor_of_safety) == len(q_c1n_cs)
        out_values = []
        for i in range(len(factor_of_safety)):
            fs_value = factor_of_safety[i]
            ev_value = _calculate_single_volumetric_strain(factor_of_safety[i], q_c1n_cs[i])
            out_values.append(ev_value)
        return np.array(out_values)
    else:
        raise ValueError("Factor of safety and q_c1n_cs must be the same length")


def _calculate_single_volumetric_strain(factor_of_safety, q_c1n_cs):
    """
    Determines the volumetric strain for a single value by interpolation of
    the equations by Zhang et al. (2002).

    :param factor_of_safety:
    :param q_c1n_cs:
    :return:
    """

    fs_values = [0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 2, 300.]
    for i in range(len(fs_values)):
        if factor_of_safety < fs_values[i]:
            if i == 0:
                low_fs = fs_values[0]
            else:
                low_fs = fs_values[i - 1]
            high_fs = fs_values[i]
            ev_low = _calculate_fixed_factor_safety_volumetric_strain(low_fs, q_c1n_cs)
            ev_high = _calculate_fixed_factor_safety_volumetric_strain(high_fs, q_c1n_cs)

            ev_actual = np.interp(factor_of_safety, [low_fs, high_fs], [ev_low, ev_high])
            return ev_actual


def _calculate_fixed_factor_safety_volumetric_strain(factor_of_safety, q_c1n_cs):
    """
    Implementation of the equations from Zhang et al. (2002)

    Note
    ====
    Must have the exact factor of safety
    :param factor_of_safety:
    :param q_c1n_cs:
    :return:
    """

    if q_c1n_cs < 33:
        q_c1n_cs = 33.

    if q_c1n_cs > 200:
        q_c1n_cs = 200.0
    if factor_of_safety == 0.5:
        e_v = 102 * q_c1n_cs ** -0.82
    elif factor_of_safety == 0.6:
        if q_c1n_cs <= 147:
            e_v = 102 * q_c1n_cs ** -0.82
        else:
            e_v = 2411 * q_c1n_cs ** -1.45
    elif factor_of_safety == 0.7:
        if q_c1n_cs <= 110:
            e_v = 102 * q_c1n_cs ** -0.82
        else:
            e_v = 1701 * q_c1n_cs ** -1.42
    elif factor_of_safety == 0.8:
        if q_c1n_cs <= 80:
            e_v = 102 * q_c1n_cs ** -0.82
        else:
            e_v = 1609 * q_c1n_cs ** -1.46
    elif factor_of_safety == 0.9:
        if q_c1n_cs <= 60:
            e_v = 102 * q_c1n_cs ** -0.82
        else:
            e_v = 1403 * q_c1n_cs ** -1.48
    elif factor_of_safety == 1.0:
        e_v = 64 * q_c1n_cs ** -0.93
    elif factor_of_safety == 1.1:
        e_v = 11 * q_c1n_cs ** -0.65
    elif factor_of_safety == 1.2:
        e_v = 9.7 * q_c1n_cs ** -0.69
    elif factor_of_safety == 1.3:
        e_v = 7.6 * q_c1n_cs ** -0.71
    elif factor_of_safety == 2.0:
        e_v1p3 = 7.6 * q_c1n_cs ** -0.71
        e_v = e_v1p3 - (e_v1p3 / 0.7) * (factor_of_safety - 1.3)  # linear interpolate
    elif factor_of_safety > 2.0:
        e_v = 0
    else:
        raise ValueError
    return e_v
