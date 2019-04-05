import numpy as np


def calc_diss_energy_fd(force, disp):
    """
    Calculates the area inside the hysteresis loops

    Parameters
    ----------
    force: array_like
        force values
    disp: array_like
        displacement values
    Returns
    -------
    : array_like
        dissipated energy series
    """
    average_force = (force[1:] + force[:-1]) / 2  # TODO: speed up by pre-allocating array of len(disp), then remove insert statements
    average_force = np.insert(average_force, 0, force[0])  # Include first value
    delta_disp = np.diff(disp)
    delta_disp = np.insert(delta_disp, 0, 0)
    return np.cumsum(average_force * delta_disp)


def calc_diss_energy_et(element_test, to_liq=False, norm=False):
    """
    Calculates the area inside the hysteresis loops of an element test

    Parameters
    ----------
    element_test: ElementTest Object
    to_liq: bool (default=False)
        if True then only go to point of liquefaction
    norm: bool (default=False)
        if True then divide by vertical effective stress

    Returns
    -------
    : array_like
        dissipated energy series
    """
    indy = element_test.n_points
    if to_liq:
        indy = element_test.i_liq
    denom = 1
    if norm:
        denom = element_test.esig_v0
    return np.cumsum(element_test.av_stress * element_test.delta_strain)[:indy] / denom


def calc_abs_delta_tau_fd(force):
    """
    Calculates the absolute change in shear stress

    Parameters
    ----------
    force: array_like
        force values

    Returns
    -------
    : array_like
        change in shear stress series
    """

    delta_force = np.diff(force)
    delta_force = np.insert(delta_force, 0, force[0])
    return np.cumsum(abs(delta_force))


def calc_abs_delta_tau_et(element_test, to_liq=False, norm=False):
    """
    Calculates the absolute change in shear stress of an element test

    Parameters
    ----------
    element_test: ElementTest Object
    to_liq: bool (default=False)
        if True then only go to point of liquefaction
    norm: bool (default=False)
        if True then divide by vertical effective stress

    Returns
    -------
    : array_like
        change in shear stress series
    """
    indy = element_test.n_points
    if to_liq:
        indy = element_test.i_liq
    denom = 1
    if norm:
        denom = element_test.esig_v0
    return calc_abs_delta_tau_fd(element_test.stress)[:indy] / denom


def average_of_absolute_via_trapz(values):
    """
    Calculates the average absolute value of the y-axis of a trapezium

    Parameters
    ----------
    values: array_like
        y-axis values
    Returns
    -------
    : array_like
        average absolute values series
    """
    values_i = values[:-1]
    values_ip1 = values[1:]
    # build trapezoids but be careful of sign changes.
    expected = np.where(values_ip1 * values_i >= 0,
                        (values_i + values_ip1) / 2,
                        (values_i ** 2 + values_ip1 ** 2) / (2 * abs(values_ip1 - values_i)))
    return abs(expected)
