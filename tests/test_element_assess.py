import numpy as np


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
