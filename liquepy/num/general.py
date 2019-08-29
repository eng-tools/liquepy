import numpy as np


def calc_rayleigh_damping_parameters(f1, f2, d):
    w1 = 2 * np.pi * f1
    w2 = 2 * np.pi * f2

    beta = 2 * (((w1 * d) - (w2 * d)) / ((w1 ** 2) - (w2 ** 2)))
    alpha = beta * w1 * w2

    wmin = np.sqrt(alpha / beta)

    emin = np.sqrt(alpha * beta)
    fmin = wmin / (2 * np.pi)
    return emin, fmin


