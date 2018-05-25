import numpy as np
from liquepy import trigger
from liquepy import checking_tools as ct

from tests.conftest import TEST_DATA_DIR


def calculate_fos():
    depths, q_c, f_s, u_2, gwl = trigger.load_cpt_data(TEST_DATA_DIR + "standard_1.csv")
    bi2014 = trigger.BoulangerIdriss2014(depths, q_c, f_s, u_2, gwl=gwl, pga=0.25, magnitude=7.5, ar=0.8)
    factor_safety_values = bi2014.factor_of_safety

    expected_fos_at_40 = 2.0
    expected_fos_at_500 = 0.541205215
    assert factor_safety_values[40] == expected_fos_at_40
    assert ct.isclose(factor_safety_values[500], expected_fos_at_500, rel_tol=0.0001)
    # np.savetxt("standard_1_fos.csv", factor_safety_values)


def compare_fos():
    depths, q_c, f_s, u_2, gwl = trigger.load_cpt_data(TEST_DATA_DIR + "standard_1.csv")
    bi2014 = trigger.BoulangerIdriss2014(depths, q_c, f_s, u_2, gwl=gwl, pga=0.25, magnitude=7.5, ar=0.8)
    factor_safety_values = bi2014.factor_of_safety

    fos_expected = np.loadtxt("standard_1_fos.csv")
    assert np.isclose(fos_expected, factor_safety_values).all()


if __name__ == '__main__':
    compare_fos()
