import numpy as np


def determine_t_liq_index(ru, ru_limit):
    """
    Finds the index where the pore pressure ratio (ru) exceeds a limit

    :param ru: array_like, pore pressure ratio series
    :param ru_limit: float, limit for liquefaction triggering
    :return:
    """

    ind2 = np.where(ru > ru_limit)
    if len(ind2[0]):
        t_liq_index = ind2[0][0]
    else:
        t_liq_index = -1
    return t_liq_index
