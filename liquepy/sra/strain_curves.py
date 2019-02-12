import numpy as np


def elastic(gamma):
    return np.ones(len(gamma))


def mohr_coloumb(gamma, elastic_mod, cu):
    tau = []
    for i in range(len(gamma)):
        tau.append(min(elastic_mod * gamma[i], cu))
    tau = np.array(tau)
    secant_mod_reduction = (tau / gamma) / elastic_mod
    return secant_mod_reduction


def flac_default_curve(gamma, l_1, l_2):
    gamma_percent = gamma * 100
    l = np.log10(gamma_percent)
    s = []
    for i in range(len(l)):
        s.append(min((l_2 - l[i]) / (l_2 - l_1), 1))
    s = np.array(s)
    return s ** 2 * (3. - 2. * s)


def default_mod(gamma):
    l_1 = -3.
    l_2 = 1.
    return flac_default_curve(gamma, l_1, l_2)


def ishi_mod(gamma):
    l_1 = -3.7
    l_2 = 0.3
    return flac_default_curve(gamma, l_1, l_2)


def seed_and_sun_mod(gamma):
    l_1 = -3.156
    l_2 = 1.904
    return flac_default_curve(gamma, l_1, l_2)


def daren_mod(gamma):
    a = 0.92
    gamma_ref = 0.001
    return 1.0 / (1.0 + (gamma / gamma_ref) ** a)


def vardanega_2013_mod(gamma, i_p):
    a = 0.943  # Eq 22b
    j = 3.7  # Eq 23
    gamma_ref = j * (i_p / 1000)
    return 1.0 / (1.0 + (gamma / gamma_ref) ** a)  # Eq. 22b


def best1_mod(gamma):
    l_1 = -2.3
    l_2 = 0.63
    return flac_default_curve(gamma, l_1, l_2)


def vardanega_w_default(gamma):
    l_1 = -2.4
    l_2 = 0.48
    return flac_default_curve(gamma, l_1, l_2)


def default_l1_and_l2(gamma, l_1, l_2):
    return flac_default_curve(gamma, l_1, l_2)