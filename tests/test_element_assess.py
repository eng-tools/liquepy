import numpy as np
import eqsig
from liquepy.element.models import ShearTest


from liquepy.element import assess


def test_with_one_cycle_no_dissipation():
    strs = np.array([0, -1, -2, -3, -4, -3, -2, -1, 0, 1, 2, 3, 4, 3, 2, 1, 0])
    tau = np.array([0, -2, -4, -6, -8, -6, -4, -2, 0, 2, 4, 6, 8, 6, 4, 2, 0])
    expected_energy = 0
    assert np.isclose(expected_energy, assess.calc_diss_energy_fd(tau, strs)[-1])


def test_with_one_cycle_no_dissipation_with_offset():
    strs = np.array([0, -1, -2, -3, -4, -3, -2, -1, 0, 1, 2, 3, 4, 3, 2, 1, 0]) + 4
    tau = np.array([0, -2, -4, -6, -8, -6, -4, -2, 0, 2, 4, 6, 8, 6, 4, 2, 0])
    expected_energy = 0
    assert np.isclose(expected_energy, assess.calc_diss_energy_fd(tau, strs)[-1])


def test_with_one_cycle_circle():
    angle = np.linspace(0, 2 * np.pi, 3600)
    strs = 4 * np.sin(angle)
    tau = 4 * np.cos(angle)
    expected_energy = 4 ** 2 * np.pi
    assert np.isclose(expected_energy, assess.calc_diss_energy_fd(tau, strs)[-1])


def test_with_one_cycle_circle_with_offset():
    angle = np.linspace(0, 2 * np.pi, 3600)
    strs = 4 * np.sin(angle) + 4
    tau = 4 * np.cos(angle) + 10
    expected_energy = 4 ** 2 * np.pi
    assert np.isclose(expected_energy, assess.calc_diss_energy_fd(tau, strs)[-1])


def test_with_one_cycle_triangles():
    strs = np.array([0, -1, -2, -3, -4, -4, -3, -2, -1, 0, 1, 2, 3, 4, 4, 3, 2, 1, 0])
    tau = np.array([0, -2, -4, -6, -8, 0, 0, 0, 0, 0, 2, 4, 6, 8, 0, 0, 0, 0, 0])
    expected_energy = 8 * 4.
    assert np.isclose(expected_energy, assess.calc_diss_energy_fd(tau, strs)[-1])


def test_average_of_absolute_simple():
    values = np.array([4, -3])
    expected = 12.5 / 7
    av_abs = assess.average_of_absolute_via_trapz(values)
    assert np.isclose(av_abs, expected), (av_abs, expected)


def test_average_of_absolute_matching_neg():
    values = np.array([3, -3, 3])
    expected = 1.5
    av_abs = assess.average_of_absolute_via_trapz(values)
    assert np.isclose(av_abs[0], expected), (av_abs[0], expected)
    assert np.isclose(av_abs[1], expected), (av_abs[1], expected)


def test_determine_cum_stored_energy_series_simple():

    gamma = np.array([0, 4, 0, -3, 0])
    tau = np.array([0, 4, 0, -3, 0])
    two_times_triangle_1 = 2 * (4 * 4 / 2)
    two_times_triangle_2 = 2 * (3 * 3 / 2)
    expected_energy = two_times_triangle_1 + two_times_triangle_2
    et = ShearTest(tau, gamma, 1)
    energy = assess.calc_case_et(et)
    assert energy[-1] == expected_energy, (energy[-1], expected_energy)


def test_small_cycle_behaviour_increases_case():

    gamma_1 = np.array([0, 4, -2, 0])
    tau_1 = np.array([0, 4, -4, 0])
    et_1 = ShearTest(tau_1, gamma_1, 1)
    energy_1 = assess.calc_case_et(et_1)
    gamma_2 = np.array([0, 4, 3, 4, -2, 0])
    tau_2 = np.array([0, 4, 1, 1, -4, 0])
    et_2 = ShearTest(tau_2, gamma_2, 1)
    energy_2 = assess.calc_case_et(et_2)
    assert energy_2[-1] > energy_1[-1]


def skip_test_strain_bulge_behaviour_increases_case():

    gamma_1 = np.array([0, 4, -2, 0])
    tau_1 = np.array([0, 4, -4, 0])
    et_1 = ShearTest(tau_1, gamma_1, 1)
    energy_1 = assess.calc_case_et(et_1)
    gamma_2 = np.array([0, 4, 4.1, -2, 0])
    tau_2 = np.array([0, 4, 1, -4, 0])
    et_2 = ShearTest(tau_2, gamma_2, 1)
    energy_2 = assess.calc_case_et(et_2)
    assert energy_2[-1] > energy_1[-1]


def test_determine_cum_stored_energy_series_simple_up_down():
    """
    /\
    :return:
    """
    gamma = np.array([0., 1., 0.5])
    tau = np.array([0., 1., 0])
    expected_delta_e = 0.75  # two triangles (1x1x0.5 + 1x0.5x0.5)
    energy = assess.calc_case_fd(tau, gamma)
    assert energy[-1] == expected_delta_e, energy
    et = ShearTest(tau, gamma)
    energy = assess.calc_case_et(et)
    assert energy[-1] == expected_delta_e, energy


def test_determine_cum_stored_energy_series_simple_up_down_neg():
    gamma = np.array([0., 1., -1])
    tau = np.array([0., 1., -1])
    expected_delta_e = 1.5
    et = ShearTest(tau, gamma)
    energy = assess.calc_case_et(et)
    assert energy[-1] == expected_delta_e, energy


def test_determine_cum_stored_energy_series_simple_close_loop():
    gamma = np.array([1., -1, 1])
    tau = np.array([1., -1, 1])

    expected_delta_e = 2
    et = ShearTest(tau, gamma)
    energy = assess.calc_case_et(et)
    assert energy[-1] == expected_delta_e, energy


def test_determine_cum_stored_energy_series_simple_4points():

    gamma = np.array([0, 1, -1, 2])
    tau = np.array([0, 1, -1, 1])
    step_1 = (0 + 1) / 2 * (1 - 0)
    step_2 = (0 + 2) / 2 * (1 - 0)
    expected_delta_e = step_1 * 4 + step_2
    et = ShearTest(tau, gamma)
    energy = assess.calc_case_et(et)
    assert energy[-1] == expected_delta_e, (energy, expected_delta_e)


def test_determine_cum_stored_energy_series_simple_trapz_zero():

    gamma = np.array([0, 2, 1])
    tau = np.array([0, 2, 1])
    step_1 = (0 + 2) / 2 * (2 - 0)
    step_2 = (2 + 1) / 2 * abs(2 - 1)
    expected_delta_e = step_1 + step_2
    et = ShearTest(tau, gamma)
    energy = assess.calc_case_et(et)
    assert energy[-1] == expected_delta_e, (energy, expected_delta_e)


def test_determine_cum_stored_energy_series_simple_trapz():

    gamma = np.array([1, 3, 2])
    tau = np.array([1, 2, 0])
    step_1 = (0 + 1) / 2 * (2 - 0)
    step_2 = (2 + 1) / 2 * abs(2 - 0)
    expected_delta_e = step_1 + step_2
    et = ShearTest(tau, gamma)
    energy = assess.calc_case_et(et)
    assert energy[-1] == expected_delta_e, (energy, expected_delta_e)


def test_determine_cum_stored_energy_series_simple_5points():

    gamma = np.array([0, 2, 1, 3, 2])
    tau = np.array([0, 2, 1, 2, 0])
    step_1 = (0 + 2) / 2 * (2 - 0)
    step_2 = (2 + 1) / 2 * abs(2 - 1)
    step_3 = (0 + 1) / 2 * (2 - 0)
    step_4 = (2 + 1) / 2 * abs(2 - 0)
    expected_delta_e = step_1 + step_2 + step_3 + step_4
    et = ShearTest(tau, gamma)
    energy = assess.calc_case_et(et)
    assert energy[-1] == expected_delta_e, (energy, expected_delta_e)


def test_case_et_simple_6points():

    gamma = np.array([0, 1, 0.5, 1.5, -1, 2])
    tau = np.array([0, 1, 0.5, 1, -1, 1])
    expected_delta_e = 4.375
    et = ShearTest(tau, gamma)
    energy = assess.calc_case_et(et)
    assert energy[-1] == expected_delta_e, (energy, expected_delta_e)


def test_get_energy_peaks_for_cyclic_loading():
    fs = np.array([0, 1., 2., 3., 4., 5., 5.5, 5.5, 4., 3., 2.5, 2.0, 1., 0., -1, -2, -5, 1, 3, 3.5,
                   2.5, 3.5, 2.5, -1, -3])
    ds = np.array([0, 0.5, 1., 1.5, 2.5, 3., 4.25, 5.5, 5.5, 5.25, 5.5, 5.25, 4., 3., 1.5, 0.5, -3, -2, -1, -0.5,
                   -0.75, 1.5, 1., -1.5, -5])
    inds = assess.get_energy_peaks_for_cyclic_loading(-fs, -ds)
    expected = np.array([0,  7, 16, 21, 24])
    assert np.sum(abs(inds - expected)) == 0
