import numpy as np


def calculate_lpi_increments(depth, liq_factor_of_safety):
    """
    formulation from 'Soil Dynamics and earthquake engineering', eq. page 317
    """
    w_unlimited = (10 - depth / 2) / 100
    w = np.where(depth < 20, w_unlimited, 0)
    f = np.where(liq_factor_of_safety < 1, 1 - liq_factor_of_safety, 0)
    return w * f


def calculate_lpi(depths, liq_factor_of_safety):
    """
    Formulation from 'Soil Dynamics and earthquake engineering', eq. page 317
    """
    return np.sum(calculate_lpi_increments(depths, liq_factor_of_safety))


def calculate_lsn_increments(e, depth):
    """
    Calculates the liquefaction severity number (LSN)

    doi: 10.1016/j.soildyn.2015.09.016

    :param e: array, volumetric strain
    :param depth: array, depth from surface
    :return: array, lsn increment at depth
    """
    ds = depth[1:] - depth[:-1]
    depth_av = (depth[1:] + depth[:-1]) / 2
    av_e = (e[1:] + e[:-1]) / 2
    # depth_av = np.insert(depth_av, len(depth_av), depth[-1])
    lsn = (av_e * ds) / depth_av
    lsn = np.insert(lsn, len(lsn), 0)
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
