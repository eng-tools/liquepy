import numpy as np


def calc_lpi_increments(liq_factor_of_safety, depth):
    """
    formulation from 'Soil Dynamics and earthquake engineering', eq. page 317
    """
    ds = depth[1:] - depth[:-1]
    depth_av = (depth[1:] + depth[:-1]) / 2
    w = np.where(depth_av < 20, (10 - 0.5 * depth_av), 0)
    fos_av = (liq_factor_of_safety[1:] + liq_factor_of_safety[:-1]) / 2
    f = np.where(fos_av < 1, 1 - fos_av, 0)
    lpis = np.zeros_like(depth)
    lpis[1:] = w * f * ds
    return lpis  # TODO: Is this correct? or should it be * dz


def calc_lpi(liq_factor_of_safety, depths):
    """
    Formulation from 'Soil Dynamics and earthquake engineering', eq. page 317
    """
    return np.sum(calc_lpi_increments(liq_factor_of_safety, depths))


def calc_lsn_increments(e_v, depth):
    """
    Calculates the liquefaction severity number (LSN)

    doi: 10.1016/j.soildyn.2015.09.016

    :param e_v: array, volumetric strain
    :param depth: array, depth from surface
    :return: array, lsn increment at depth
    """
    ds = np.zeros_like(depth)
    ds[:-1] = depth[1:] - depth[:-1]
    depth_av = np.array(depth)
    depth_av[:-1] = (depth[1:] + depth[:-1]) / 2
    lsn = (e_v * ds) / depth_av
    return lsn * 10


def calculate_lsn_increments(e_v, depth):
    """
    Deprecated: see calc_lsn_increments
    """
    return calc_lsn_increments(e_v, depth)


def calc_lsn(e_v, depth):
    """
    Calculates the liquefaction severity number (LSN)

    doi: 10.1016/j.soildyn.2015.09.016

    :param e_v: array,
        volumetric strain in percentage
    :param depth: array,
        depth from surface
    :return: float,
        LSN for profile
    """
    return sum(calc_lsn_increments(e_v, depth))


def calculate_lsn(e_v, depth):
    return calc_lsn(e_v, depth)


def calc_ldi_increments(e_s, depth):
    """
    Calculates the Lateral Displacement Index :cite:`Zhang:2004el`

    Parameters
    ----------
    e_s
    depth

    Returns
    -------

    """

    ds = depth[1:] - depth[:-1]
    depth_av = (depth[1:] + depth[:-1]) / 2
    av_e = (e_s[1:] + e_s[:-1]) / 2
    # depth_av = np.insert(depth_av, len(depth_av), depth[-1])
    ldi = np.zeros_like(depth)
    ldi[1:] = av_e * ds
    return ldi


def calc_ldi(e_s, depth, z_max=None):
    """
    Calculates the Lateral Displacement Index :cite:`Zhang:2004el`
    
    Parameters
    ----------
    e_s
    depth
    z_max

    Returns
    -------

    """
    if z_max is not None and depth[-1] > z_max:
        indy = np.argmin(abs(depth - z_max))
        e_s = e_s[:indy]
        depth = depth[:indy]
    return np.trapz(y=e_s, x=depth)

