__author__ = 'maximmillen'

import numpy as np
import matplotlib.pyplot as plt


def calc_vs(g_mod, density):
    """
    Calculates the soil shear wave velocity

    g_mod: Soil shear modulus
    density: Soil mass density

    """
    return np.sqrt(g_mod / density)  # m/s


def calc_g_mod(vs, density):
    return vs ** 2 * density


def calc_damped_vs_kramer_1996(vs, xi):
    """
    Calculates the damped shear wave velocity
    Ref: Eq 7.9 from Kramer (1996)

    :param vs:
    :param xi:
    :return:
    """
    return vs * (1.0 + 1j * xi)


def calc_damped_vs_dormieux_1990(vs, xi):
    """
    Calculates the damped shear wave velocity

    :param vs:
    :param xi:
    :return:
    """
    return vs * np.sqrt(np.sqrt(1 - 4 * xi ** 2) + 2j * xi)


def calc_alternative_damped_vs(vs, xi):
    """
    Calculates the shear wave velocity with material damping (No citation)

    :param vs:
    :param xi:
    :return:
    """
    sqrt_d = np.sqrt(1.0 + 4 * xi ** 2)
    return vs / sqrt_d * ((1.0 + sqrt_d) / 2 + 1j * xi)


def calc_impedance(density_soil, vs_soil, density_rock, vs_rock):
    """
    Calculates the impedance contrast between the soil and rock
    Ref: Eq 7.23 from Kramer (1996)

    :param density_soil: Soil mass density
    :param g_soil: Soil shear modulus
    :param density_rock: Rock mass density
    :param g_rock: Rock shear modulus
    :return:
    """
    impedance = density_soil * vs_soil / (density_rock * vs_rock)
    return impedance


def calc_tf_rigid_br(h_deposit, vs, omega, absolute=False):
    """
    Calculates the bedrock-to-surface transfer, assuming the bedrock is rigid.
    Ref: Eq 7.8 from Kramer (1996)

    :param h_deposit: Height of the soil deposit
    :param vs: Soil shear wave velocity
    :param omega: [array] Angular frequency
    :return: [array] transfer function
    """
    kh = omega * h_deposit / vs
    if absolute:
        h_surface = np.abs(1.0 / np.cos(kh))
    else:
        h_surface = 1.0 / np.cos(kh)
    return h_surface


def calc_tf_elastic_br(h_deposit, vsi_soil, omega, impedance, absolute=False):
    """
    Surface transfer function for visco-elastic soil and an elastic bedrock
    Ref: Eq 7.27 of Kramer (1996), consistent with Faccioli, 2005

    :param h_deposit: Height of the soil deposit
    :param vs: Soil shear wave velocity
    :param omega: [array] Angular frequency
    :param impedance: Soil-to-bedrock impedance contrast
    :return:
    """

    h_surface = 1.0 / (np.cos(omega * h_deposit / vsi_soil) + 1j * impedance * np.sin(omega * h_deposit / vsi_soil))
    if absolute:
        return np.abs(h_surface)
    return h_surface


def view_all_theoretical(sub_fig):
    h_deposit = 20.0  # m
    xi = 0.15
    g_soil = 80.0e6  # Pa
    density_soil = 1800.0  # kg/m3

    f_low = 0.1  # Hz
    f_high = 15.0  # Hz

    vs_soil = calc_vs(g_soil, density_soil)  # m/s
    vsi_soil = calc_damped_vs(vs_soil, xi)

    omega_soil = vs_soil / (2 * h_deposit)
    f_1 = vs_soil / (4 * h_deposit)
    f_2 = 3 * f_1
    f_3 = 5 * f_1

    frequencies = omega_soil / 10 * np.linspace(1, 100, 100)
    omega = 2 * np.pi * frequencies

    density_rock = 2200.0  # kg/m3
    g_rock = 400.0e6  # Pa
    vs_rock = calc_vs(g_rock, density_rock)

    impedance = calc_impedance(density_soil, vs_soil, density_rock, vs_rock)

    a_ud_rb = calc_tf_rigid_br(h_deposit, vsi_soil, omega)
    sub_fig.plot(frequencies, a_ud_rb, label="theoretical - damped rigid br")

    w_alt = 0
    if w_alt:
        vsi_soil_alt = calc_alternative_damped_vs(vs_soil, xi)
        a_d_rb_alt = calc_tf_rigid_br(h_deposit, vsi_soil_alt, omega)
        sub_fig.plot(frequencies, a_d_rb_alt, label="theoretical - damp rigid (Alt)")
    a_ud_eb = calc_tf_elastic_br(h_deposit, vs_soil, omega, impedance)
    sub_fig.plot(frequencies, a_ud_eb, ls="--", label="theoretical - undamped elastic br")
    a_d_eb = calc_tf_elastic_br(h_deposit, vsi_soil, omega, impedance)
    sub_fig.plot(frequencies, a_d_eb, ls="--", label="theoretical - damped elastic br")


def view_theoretical(h_deposit, g_soil, density_soil=1800, xi=0.1, impedance=0.5, sub_fig=None, freqs=(0.1, 30.0)):
    if sub_fig is None:
        fig, sub_fig = plt.subplots()

    vs_soil = calc_vs(g_soil, density_soil)  # m/s
    vsi_soil = calc_damped_vs(vs_soil, xi)

    frequencies = np.logspace(np.log10(freqs[0]), np.log10(freqs[1]), 200, base=10)
    omega = 2 * np.pi * frequencies

    a_ud_rb = calc_tf_rigid_br(h_deposit, vsi_soil, omega)
    sub_fig.plot(frequencies, a_ud_rb, label="theoretical - damped rigid br")

    a_d_eb = calc_tf_elastic_br(h_deposit, vsi_soil, omega, impedance)
    sub_fig.plot(frequencies, a_d_eb, ls="--", label="theoretical - damped elastic br")


if __name__ == '__main__':
    fig, sub_fig = plt.subplots()
    # theoretical_transfer_function(sub_fig)
    # theoretical_elastic_bedrock_transfer_function(sub_fig)
    view_all_theoretical(sub_fig)
    plt.legend()
    plt.show()
