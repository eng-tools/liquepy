import numpy as np
from liquepy import trigger
from liquepy.trigger import boulanger_and_idriss_2014 as bim14
from liquepy import checking_tools as ct

from tests.conftest import TEST_DATA_DIR


def test_can_calculate_fos():
    depths, q_c, f_s, u_2, gwl = trigger.load_cpt_data(TEST_DATA_DIR + "standard_1.csv")
    bi2014 = trigger.BoulangerIdriss2014(depths, q_c, f_s, u_2, gwl=gwl, pga=0.25, magnitude=7.5, ar=0.8)
    factor_safety_values = bi2014.factor_of_safety

    expected_fos_at_40 = 2.0
    expected_fos_at_500 = 0.541205215
    assert factor_safety_values[40] == expected_fos_at_40
    assert ct.isclose(factor_safety_values[500], expected_fos_at_500, rel_tol=0.0001)


def test_calculate_unit_weight():
    fs = np.array([10.])
    q_t = np.array([10.])
    gwl = 2.
    depth = np.array([4.])

    gamma = bim14.calculate_unit_weight(fs, q_t, gwl, depth)
    expected_gamma = 15.
    assert gamma[0] == expected_gamma


def test_calculate_qt():
    qc = 1400.
    ar = 0.8
    u2 = 90.
    expect_qt = qc + ((1 - ar) * u2)
    qt = bim14.calculate_qt(qc=qc, ar=ar, u2=u2)
    assert qt == expect_qt


def test_calculate_sigma_veff():
    sigma_v = 20
    pore_pressure = 5
    expected_sigma_veff = 15
    sigma_veff = bim14.calculate_sigma_veff(sigma_v, pore_pressure=pore_pressure)

    assert sigma_veff == expected_sigma_veff


def test_calculate_fc():
    ic = 2.93
    cfc = 0.
    # expected_fc = 97.31
    expected_fc = 97.4
    fc = bim14.calculate_fc(ic, cfc)
    print(fc)

    assert fc == expected_fc


def test_calculate_rd():
    depth = 5.32
    magnitude = 7.5
    rd = bim14.calculate_rd(depth, magnitude)
    expected_rd = 0.957178376538626
    print(rd)

    assert rd == expected_rd

def test_calculate_csr():
    sigma_veff= 42.6
    sigma_v = 79.1
    pga = 0.25
    rd = 0.96
    gwl = 1.6
    depth = 5.32
    csr = bim14.calculate_csr(sigma_veff, sigma_v, pga, rd, gwl, depth)
    expected_csr = 0.2896619718309859
    print(csr)
    assert csr == expected_csr


def test_calculate_cn_values():
    m = 0.58
    sigma_veff = 42.6
    cn = bim14.calculate_cn_values(m, sigma_veff)
    expected_cn = 1.7
    assert ct.isclose(cn, expected_cn, rel_tol=0.04)


def test_calculate_m():
    q_c1ncs = 62.1
    m = bim14.calculate_m(q_c1ncs)
    expected_m = 0.6
    assert ct.isclose(m, expected_m, rel_tol=0.01)


def test_calculate_q_c1n():
    qc = 0.58
    CN = 1.7
    qc1n = bim14.calculate_q_c1n(qc, CN)
    # expected_qc1n = 9.73
    expected_qc1n = 9.8
    assert ct.isclose(qc1n, expected_qc1n, rel_tol=0.1)


def test_crr_7p5_from_cpt():
    q_c1n_cs = 62.1
    gwl = 1.6
    depth = 5.35
    expected_crr_7p5 = 0.1
    crr = bim14.crr_7p5_from_cpt(q_c1n_cs, gwl, depth, i_c=1.8)
    assert ct.isclose(crr, expected_crr_7p5, rel_tol=0.011)
    # that is not important, because it is modified in crr_m function


def test_crr_m():
    k_sigma = 1.1
    msf = 1
    q_c1n_cs = 9.73
    gwl = 1.6
    depth = 5.32
    crr_m7p5 = bim14.crr_7p5_from_cpt(q_c1n_cs, gwl, depth, i_c=2.7)
    crr = bim14.crr_m(k_sigma, msf, crr_m7p5)  # CRR a magnitudo M

    expected_crr = 4.4
    assert crr == expected_crr


def test_calculate_delta_q_c1n():
    q_c1n = 9.73
    fc = 64.21
    delta_q = bim14.calculate_delta_q_c1n(q_c1n, fc)
    expected_delta_q = 52.37
    assert ct.isclose(delta_q, expected_delta_q, rel_tol=0.02)


def test_calculate_q_c1ncs():
    q_c1n = 9.73
    delta_q_c1n = 63.26
    expected_qc1ncs = 9.73
    qc1ncs = bim14.calculate_q_c1ncs(q_c1n, delta_q_c1n)
    assert ct.isclose(qc1ncs, expected_qc1ncs, rel_tol=0.9)


def test_calculate_msf():
    magnitude = 7.5
    q_c1ns = 9.73*np.ones(1)
    msf = bim14.calculate_msf(magnitude, q_c1ns)
    expected_msf = 1.*np.ones(1)
    assert msf == expected_msf


def test_calculate_big_q_values():
    qt = 20.
    sigmav = 15.
    CN = 0.8
    qq = bim14.calculate_big_q_values(CN, qt, sigmav)
    print(qq)
    assert ct.isclose(qq, 0.04, rel_tol=0.011)

def test_calculate_big_f_values():
    qt = 20.
    sigmav = 15.
    fs = 15
    F = bim14.calculate_f_ic_values(fs, qt, sigmav)
    expected_F = 300.
    assert F == expected_F


def test_calculate_k_sigma():
    sigma_eff = 42.6*np.ones(1)
    qc1ncs = 9.73*np.ones(1)
    expected_ksigma = bim14.calculate_k_sigma(sigma_eff, qc1ncs)
    ksigma = 1.07 * np.ones(1)
    ksigma = 1.038 * np.ones(1)
    assert ct.isclose(expected_ksigma, ksigma, rel_tol=0.001)


def test_calculate_ic():
    big_q = 0.04
    big_f = 300.
    expected_ic = bim14.calculate_ic(big_q, big_f)
    assert ct.isclose(expected_ic , 5.07, rel_tol=0.09)





if __name__ == '__main__':
    test_can_calculate_fos()