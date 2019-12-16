import numpy as np
import eqsig


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


def calc_stored_energy_abs_incs_fd_peaks_and_indices(forces, disps, peaks_from='disp'):
    """
    Calculates the absolute change stored energy for an oscillating system

    >>> disp = np.array([0, 4, 0, -2, 0])
    >>> force = np.array([0, 4, 0, -4, 0])
    >>> calc_stored_energy_abs_incs_fd_peaks_and_indices(force, disp)
    (array([ 0.,  8., 12.,  4.]), array([0, 1, 3, 4]))

    Parameters
    ----------
    disps: array_like
        displacement time series
    forces: array_like
        force time series
    peaks_from: str or array_like
        Method to determine peaks or array of peak indices
    :return: array_like
    """
    forces = forces.astype(float)
    disps = disps.astype(float)
    if peaks_from == 'disp':
        peak_indices = eqsig.get_peak_array_indices(disps)
    elif peaks_from == 'cyclic_energy':
        peak_indices = get_energy_peaks_for_cyclic_loading(forces, disps)  # slow
    elif not isinstance(peaks_from, str):
        peak_indices = np.array(peaks_from)
    else:
        raise ValueError("'peaks_from' must be either 'disp', 'cyclic_energy' or "
                         "array indices, not {0}".format(peaks_from))
    peak_forces = np.take(forces, peak_indices)
    peak_disps = np.take(disps, peak_indices)

    # Take the average of the absolute average and the signed average
    # accounts for sign change through a step
    peak_abs_av_force = average_of_absolute_via_trapz(peak_forces)
    peak_abs_delta_disps = abs(np.diff(peak_disps))
    peak_abs_delta_work = peak_abs_delta_disps * peak_abs_av_force  # trapezoid
    peak_abs_delta_work = np.insert(peak_abs_delta_work, 0, 0)

    return peak_abs_delta_work, peak_indices


def calc_case_peaks_and_indices_fd(forces, disps):
    """
    Calculates the cumulative change stored energy for an oscillating system at the peaks

    Note: This quantity is double the area of the triangles in a full cycle, since positive and negative
    triangles are counted.

    >>> disp = np.array([0, 4, 0, -2, 0])
    >>> force = np.array([0, 4, 0, -4, 0])
    >>> calc_case_peaks_and_indices_fd(disp, force)
    (array([ 0.,  8., 20.,  24.]), array([0, 1, 3, 4]))

    Parameters
    ----------
    disps: array_like
        displacement time series
    forces: array_like
        force time series
    :return: array_like
    """
    peak_abs_delta_work, peak_indices = calc_stored_energy_abs_incs_fd_peaks_and_indices(forces, disps)
    cum_work = np.cumsum(peak_abs_delta_work)
    return cum_work, peak_indices


def calc_case_fd(forces, disps, stepped=False, peaks_from='disp'):
    """
    Calculates the cumulative change in stored energy for an oscillating system.

    Note: This quantity is double the area of the triangles in a full cycle, since positive and negative
    triangles are counted.

    >>> disp = np.array([0, 4, 0, -2, 0])
    >>> force = np.array([0, 4, 0, -4, 0])
    >>> calc_case_fd(force, disp, stepped=True)
    array([ 0.,  8., 8., 20.,  24.])

    Parameters
    ----------
    disps: array_like,
        displacement time series
    forces: array_like,
        force time series
    stepped: bool,
        if false then values are interpolated between peaks
    :return: array_like
    """

    peak_abs_delta_work, peak_indices = calc_stored_energy_abs_incs_fd_peaks_and_indices(forces, disps,
                                                                                         peaks_from=peaks_from)
    if stepped:
        peak_abs_delta_work_full_series = np.zeros_like(disps)
        np.put(peak_abs_delta_work_full_series, peak_indices, peak_abs_delta_work)
        cum_abs_delta_work_full_series = np.cumsum(peak_abs_delta_work_full_series)
        return cum_abs_delta_work_full_series
    else:
        incs = np.arange(len(disps))
        incs_peaks = np.take(incs, peak_indices)
        cum_abs_delta_work = np.cumsum(peak_abs_delta_work)
        cum_abs_delta_work_full_series = np.interp(incs, incs_peaks, cum_abs_delta_work)
        return cum_abs_delta_work_full_series


def calc_case_et(element_test, stepped=False, to_liq=False, norm=False, peaks_from='disp'):
    """
    Calculates the absolute elastic work (case), cumulative absolute change in stored energy for an element test

    >>> gamma = np.array([0, 4, 0, -2, 0])
    >>> tau = np.array([0, 4, 0, -4, 0])
    >>> element_test = liquepy.element.models.ElementTest(tau, gamma, 1)
    array([ 0.,  8., 8., 20.,  24.])

    Parameters
    ----------
    element_test: ElementTest object
    stepped: bool
        if false then values are interpolated between peaks
    to_liq: bool
        if true then compute only to the point where liquefaction is reached
    norm: bool
        if true the divide by initial vertical effective stress
    :return: array_like
    """
    indy = element_test.n_points
    if to_liq:
        indy = element_test.i_liq
    denom = 1
    if norm:
        denom = element_test.esig_v0
    return calc_case_fd(element_test.stress, element_test.strain, stepped=stepped, peaks_from=peaks_from)[:indy] / denom


def calc_damping_et(element_test, to_liq=False, cumulative=False):
    """
    Calculates the damping ratio during the element test

    Parameters
    ----------
    element_test
    to_liq

    Returns
    -------

    """
    diss_e = calc_diss_energy_et(element_test, to_liq=False, norm=False)
    case = calc_case_et(element_test, to_liq=False, norm=False)
    indy = None
    if to_liq:
        indy = element_test.i_liq
    if cumulative:
        return diss_e / (np.pi * case)[:indy]
    else:
        ci = eqsig.get_zero_crossings_array_indices(element_test.stress)
        loc_damp = np.diff(diss_e[ci]) / (np.pi * np.diff(case[ci]))
        return np.interp(np.arange(element_test.n_points, ci, loc_damp))[:indy]


def get_energy_peaks_for_cyclic_loading(forces, disps):
    zi = eqsig.get_zero_crossings_array_indices(forces)
    if zi[-1] != len(forces):
        zi = np.append(zi, len(forces))
    inds = [0]
    for i in range(len(zi) - 1):
        e = (forces[zi[i]: zi[i + 1]] - forces[inds[i]]) * (disps[zi[i]: zi[i + 1]] - disps[inds[i]])
        inds.append(np.argmax(e) + zi[i])
    if inds[1] == 0:  # remove due to non zero start
        inds = inds[1:]
    return np.array(inds)



if __name__ == '__main__':
    fs = np.array([0, 1., 2., 3., 4., 5., 5.5, 5.5, 4., 3., 2.5, 2.0, 1., 0., -1, -2, -5, 1, 3, 3.5,
                   2.5, 3.5, 2.5, -1, -3])
    ds = np.array([0, 0.5, 1., 1.5, 2.5, 3., 4.25, 5.5, 5.5, 5.25, 5.5, 5.25, 4., 3., 1.5, 0.5, -3, -2, -1, -0.5,
                   -0.75, 1.5, 1., -1.5, -5])
    inds = get_energy_peaks_for_cyclic_loading(-fs, -ds)
    print(inds)
    import matplotlib.pyplot as plt
    plt.plot(ds, fs)
    plt.plot(ds[inds], fs[inds], 'o')
    plt.show()