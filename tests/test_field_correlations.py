import liquepy as lq
import numpy as np


def test_calc_q_c1n_via_inverse_d_r_boulanger_et_al_2014():
    q = 70
    d_r = lq.field.correlations.calc_relative_density_boulanger_et_al_2014_cpt_values(q)
    q_new = lq.field.correlations.calc_q_c1n_via_inverse_d_r_boulanger_et_al_2014(d_r)
    assert np.isclose(q, q_new), (q, q_new)


if __name__ == '__main__':
    test_calc_q_c1n_via_inverse_d_r_boulanger_et_al_2014()
