import liquepy as lq
import numpy as np


def test_calc_q_c1n_via_inverse_d_r_boulanger_et_al_2014():
    q = 70
    d_r = lq.field.correlations.calc_relative_density_boulanger_et_al_2014_cpt_values(q)
    q_new = lq.field.correlations.calc_q_c1n_via_inverse_d_r_boulanger_et_al_2014(d_r)
    assert np.isclose(q, q_new), (q, q_new)


class FakeCPT(object):
    def __init__(self):
        self.f_s = None


def test_calc_shear_vel_mcgann_2015_cpt():
    cpt = FakeCPT()
    cpt.f_s = 15.0  # kPa
    cpt.q_c = 600.0  # kPa
    cpt.depth = 2.9  # m
    vs = lq.field.correlations.calc_shear_vel_mcgann_2015_cpt(cpt)
    assert np.isclose(vs, 77.85187, rtol=0.01), vs

if __name__ == '__main__':
    test_calc_shear_vel_mcgann_2015_cpt()
