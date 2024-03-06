import numpy as np

def calc_k_alpha_for_qc1ncs(alpha, qc1ncs, esig_v0, big_q=10, p_atm=101.0e3, k0=0.45):
    alpha = abs(alpha)
    xi_r = 1 / (big_q - np.log(100 * (1 + 2 * k0) * esig_v0) / (3 * p_atm)) - (0.478 * qc1ncs ** 0.264 - 1.063)
    a = 1267 + 636 * alpha ** 2 - 634 * np.exp(alpha) - 632 * np.exp(-alpha)
    b = np.exp(-1.11 + 12.3 * alpha ** 2 + 1.31 * np.log(alpha + 0.0001))
    c = 0.138 + 0.126 * alpha + 2.52 * alpha ** 3
    k_alpha = a + b * np.exp(-xi_r / c)
    return k_alpha
