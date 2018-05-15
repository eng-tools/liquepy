from liquepy import settlements

import numpy as np
from liquepy import checking_tools as ct
from tests.conftest import TEST_DATA_DIR


def test_cavdp():

    fpath = TEST_DATA_DIR + "input_acc.his"
    acc_file = np.loadtxt(fpath, skiprows=4)
    acc = acc_file[:, 1]
    acc = acc / 9.81
    time = acc_file[:, 0]
    dt = time[1] - time[0]

    cav_dp = settlements.calculate_cav_dp(acc, dt)

    # 1.4598176 from liquepy 0.1.0  tested with several motions
    assert ct.isclose(cav_dp, 1.4598176, rel_tol=0.001), cav_dp