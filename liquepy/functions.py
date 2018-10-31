import numpy as np


def determine_t_liq_index(ru, ru_limit, return_none=False):
    """
    Finds the index where the pore pressure ratio (ru) exceeds a limit

    :param ru: array_like, pore pressure ratio series
    :param ru_limit: float, limit for liquefaction triggering
    :param return_none: bool, if True then returns none if liquefaction does not occur
    :return:
    """

    ind2 = np.where(ru > ru_limit)
    if len(ind2[0]):
        t_liq_index = ind2[0][0]
    else:
        if return_none:
            return None
        t_liq_index = len(ru)
    return t_liq_index
