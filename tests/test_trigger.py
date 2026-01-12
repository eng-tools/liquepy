import numpy as np
import pytest

import liquepy as lq


@pytest.mark.parametrize(
    "q_c1n, fos, gamma_max",
    [
        (263.04, 2.0, 0.101e-2),
        (40.434, 0.686, 51.2e-2),
        (60.545, 0.59, 33.6007e-2),
        (25.851, 0.45, 51.2e-2),
        (147.94, 1.413, 1.566e-2),
    ],
)
def calc_shear_strain_zhang_2004(q_c1n, fos, gamma_max):
    d_r = lq.trigger.calc_relative_density_zhang_2002(q_c1n) / 1e2
    gamma_inc = lq.trigger.calc_shear_strain_zhang_2004(fos, d_r)
    assert np.isclose(gamma_inc, gamma_max, rtol=0.001), gamma_inc


def test_single_calc_shear_strain():
    fos = 2.0
    q_c1n = 263.04
    d_r = lq.trigger.calc_relative_density_zhang_2002(q_c1n)
    gamma_inc = lq.trigger.calc_shear_strain_zhang_2004(fos, d_r)
    assert np.isclose(gamma_inc, 0.101e-2, rtol=0.001), gamma_inc

    fos = 0.686
    q_c1n = 40.434
    d_r = lq.trigger.calc_relative_density_zhang_2002(q_c1n)
    gamma_inc = lq.trigger.calc_shear_strain_zhang_2004(fos, d_r)
    assert np.isclose(gamma_inc, 51.2e-2), gamma_inc

    fos = 1.413
    q_c1n = 147.94  # not q_c1ncs
    d_r = lq.trigger.calc_relative_density_zhang_2002(q_c1n)
    gamma_inc = lq.trigger.calc_shear_strain_zhang_2004(fos, d_r)
    assert np.isclose(gamma_inc, 1.566e-2, rtol=0.001), gamma_inc


def calc_volumetric_strain_zhang_2004():

    fos = 0.686
    q_c1ncs = 89.667
    e_inc = lq.trigger.calc_volumetric_strain_zhang_2004(fos, q_c1ncs)
    assert np.isclose(e_inc, 2.5553, rtol=0.001), e_inc


if __name__ == "__main__":
    # test_single_calc_shear_strain()
    test_calc_volumetric_strain()
