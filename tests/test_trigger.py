import liquepy as lq
import numpy as np


def test_calc_shear_strain():
    # fos = 0.53
    # q_c1n = 34.06
    # d_r = lq.trigger.calc_relative_density_zhang_2002(q_c1n)
    # gamma_inc = lq.trigger.calc_shear_strain(fos, d_r)
    # assert np.isclose(gamma_inc, 38.48), gamma_inc

    fos = 0.86
    q_c1n = 113.96
    d_r = lq.trigger.calc_relative_density_zhang_2002(q_c1n)
    gamma_inc = lq.trigger.calc_shear_strain(fos, d_r)
    assert np.isclose(gamma_inc, 5.4), gamma_inc


def test_calc_volumetric_strain():
    # fos = 0.53
    # q_c1n = 34.06
    # d_r = lq.trigger.calc_relative_density_zhang_2002(q_c1n)
    # gamma_inc = lq.trigger.calc_shear_strain(fos, d_r)
    # assert np.isclose(gamma_inc, 38.48), gamma_inc

    fos = 0.86
    q_c1n = 113.96
    e_inc = lq.trigger.calc_volumetric_strain(fos, q_c1n)
    assert np.isclose(e_inc, 1.62), e_inc


if __name__ == '__main__':
    # test_calc_shear_strain()
    test_calc_volumetric_strain()
