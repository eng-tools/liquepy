import numpy as np


def calculate_lpi_weighting_factor(depth):
    w_unlimited = (10 - depth / 2) / 100
    w = np.where(depth < 20, w_unlimited, 0)
    return w


def calculate_f_values(factor_safety):
    f_tent = 0
    f = np.where(factor_safety < 1, 1 - factor_safety, f_tent)
    return f


def calculate_lpi_increments(depths, liq_factor_of_safety):
    """
    formulation from 'Soil Dynamics and earthquake engineering', eq. page 317
    """
    w = calculate_lpi_weighting_factor(depths)
    f = calculate_f_values(liq_factor_of_safety)
    return w * f


def calculate_lpi(depths, liq_factor_of_safety):
    """
    formulation from 'Soil Dynamics and earthquake engineering', eq. page 317
    """
    return sum(calculate_lpi_increments(depths, liq_factor_of_safety))


def calculate_lsn_increments(e, depth):
    """
    Calculates the liquefaction severity number (LSN)

    doi: 10.1016/j.soildyn.2015.09.016

    :param e: array, volumetric strain
    :param depth: array, depth from surface
    :return: array, lsn increment at depth
    """
    ds = depth[1:] - depth[:-1]
    ds = np.insert(ds, 0, depth[0])
    lsn = (e * ds) / depth
    return lsn * 10


def calculate_lsn(e, depth):
    """
    Calculates the liquefaction severity number (LSN)

    doi: 10.1016/j.soildyn.2015.09.016

    :param e: array, volumetric strain
    :param depth: array, depth from surface
    :return: float, profile LSN
    """
    return sum(calculate_lsn_increments(e, depth))
