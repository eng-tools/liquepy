import numpy as np

from liquepy.element import models


def test_av_stress():
    tau = np.array([0, -2, -4, -6, -8, -6, -4, -2, 0, 2, 4, 6, 8, 6, 4, 2, 0])
    av_tau = (tau[1:] + tau[:-1]) / 2
    et = models.ShearTest(tau, np.ones_like(tau), 1)
    diff = sum(abs(av_tau - et.av_stress[1:]))
    assert np.isclose(diff, 0)


def test_da_strain():
    strain = np.array(
        [0, -0.01, -0.02, -0.01, 0.015, 0.02, 0.005, -0.017, -0.033, -0.04, -0.02]
    )
    et = models.ShearTest(np.ones_like(strain), strain, 1)
    da_strain = et.get_da_strain_series()
    expected = np.array(
        [0, 0.01, 0.02, 0.01, 0.035, 0.04, 0.015, 0.037, 0.053, 0.06, 0.02]
    )
    abs_diff = np.sum(abs(da_strain - expected))
    if abs_diff > 0.0001:
        for i in range(len(da_strain)):
            assert np.isclose(da_strain[i], expected[i]), (i, da_strain[i], expected[i])
    et.set_i_liq(da_strain_limit=0.05)
    assert et.i_liq == 8, et.i_liq


if __name__ == "__main__":
    test_da_strain()
