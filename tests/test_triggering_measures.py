import liquepy as lq
from tests.conftest import TEST_DATA_DIR
import numpy as np


def test_calc_lsn():
    cpt = lq.field.load_cpt_from_file(TEST_DATA_DIR + "standard_1.csv")
    bi2014 = lq.trigger.run_bi2014(cpt, pga=0.25, m_w=7.5, gwl=cpt.gwl)
    epsilon = lq.trigger.calc_volumetric_strain(bi2014.factor_of_safety, bi2014.q_c1n_cs)
    lsn_direct = lq.trigger.calc_lsn(epsilon, cpt.depth)
    assert np.isclose(lsn_direct, 29.498444636615105, rtol=0.01)


def test_calc_lpi_increments():
    fos = np.array([0.53])
    depth = np.array([8.2, 8.21])
    lpi_inc = lq.trigger.triggering_measures.calc_lpi_increments(depth, fos)
    assert np.isclose(lpi_inc[0], 0.000277, rtol=0.01), lpi_inc[0]  # unvalidated test value


if __name__ == '__main__':
    test_calc_lpi_increments()
