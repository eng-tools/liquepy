import liquepy as lq
from tests.conftest import TEST_DATA_DIR
import numpy as np


def test_can_calculate_lsn():
    cpt = lq.field.load_cpt_from_file(TEST_DATA_DIR + "standard_1.csv")
    bi2014 = lq.trigger.run_bi2014(cpt, pga=0.25, m_w=7.5)
    epsilon = lq.trigger.calculate_volumetric_strain(bi2014.factor_of_safety, bi2014.q_c1n_cs)
    lsn_direct = lq.trigger.calculate_lsn(epsilon, cpt.depth)
    assert np.isclose(lsn_direct, 29.2636, rtol=0.01)
