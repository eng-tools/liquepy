import numbers
import numpy as np


def calculate_shear_strain(fs, d_r):
    if isinstance(fs, numbers.Real):
        if isinstance(d_r, numbers.Real):
            return calculate_single_shear_strain(fs, d_r)
        else:
            out_values = []
            for i in range(len(d_r)):
                out_values.append(calculate_single_shear_strain(fs, d_r[i]))
            return np.array(out_values)
    else:
        out_values = []
        for i in range(len(fs)):
            if isinstance(d_r, numbers.Real):
                out_values.append(calculate_single_shear_strain(fs[i], d_r))
            else:
                out_values.append(calculate_single_shear_strain(fs[i], d_r[i]))
        return np.array(out_values)


def calculate_array_shear_strain(fs, d_r):
    if isinstance(fs, numbers.Real) and isinstance(d_r, numbers.Real):
        return calculate_single_shear_strain(fs, d_r)
    elif not isinstance(d_r, numbers.Real) and not isinstance(fs, numbers.Real):
        assert len(fs) == len(d_r)
        out_values = []
        for i in range(len(fs)):
            fs_value = fs[i]
            ev_value = calculate_single_shear_strain(fs[i], d_r[i])
            out_values.append(ev_value)
        return np.array(out_values)
    else:
        raise ValueError("Factor of safety and q_c1n_cs must be the same length")


def calculate_single_shear_strain(fs, d_r):
    if d_r == -1:
        return 0
    dr_values = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
    for i in range(len(dr_values)):
        if d_r < dr_values[i]:
            if i == 0:
                low_dr = dr_values[0]
            else:
                low_dr = dr_values[i - 1]
            high_dr = dr_values[i]
            # print(low_dr, high_dr)
            ev_low = calculate_fixed_dr_gamma_max(fs, low_dr)
            ev_high = calculate_fixed_dr_gamma_max(fs, high_dr)

            ev_actual = np.interp(fs, [low_dr, high_dr], [ev_low, ev_high])
            # print(ev_low, ev_high, ev_actual)
            return ev_actual


# Maximum cyclic shear strains
def calculate_fixed_dr_gamma_max(fs, relative_density):
    if fs > 2.0:
        return 0.0
    elif relative_density == 0.9:
        if fs >= 0.7:
            gamma_max = 3.26 * fs ** (-1.80)
        else:
            gamma_max=6.2
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
    return gamma_max
