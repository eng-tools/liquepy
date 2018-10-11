import numpy as np
from liquepy import trigger
from liquepy.trigger import boulanger_and_idriss_2014 as bim14

from tests.conftest import TEST_DATA_DIR


def test_can_calculate_fos():
    cpt = trigger.load_cpt_from_file(TEST_DATA_DIR + "standard_1.csv")
    bi2014 = trigger.run_bi2014(cpt, pga=0.25, magnitude=7.5)
    factor_safety_values = bi2014.factor_of_safety

    expected_fos_at_40 = 2.0
    expected_fos_at_500 = 0.541205215
    assert factor_safety_values[40] == expected_fos_at_40
    assert np.isclose(factor_safety_values[500], expected_fos_at_500, rtol=0.0001)


def test_compare_fos_to_previous_version():
    depths, q_c, f_s, u_2, gwl = trigger.load_cpt_data(TEST_DATA_DIR + "standard_1.csv")
    bi2014 = trigger.BoulangerIdriss2014(depths, q_c, f_s, u_2, gwl=gwl, pga=0.25, magnitude=7.5, a_ratio=0.8)
    factor_safety_values = bi2014.factor_of_safety

    fos_expected = np.loadtxt(TEST_DATA_DIR + "standard_1_fos.csv")  # from liquepy=0.1.6
    assert np.isclose(fos_expected, factor_safety_values).all()


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
    assert fc == expected_fc


def test_calculate_rd():
    depth = 5.32
    magnitude = 7.5
    rd = bim14.calculate_rd(depth, magnitude)
    expected_rd = 0.95717837
    assert np.isclose(rd, expected_rd)


def test_calculate_csr():
    sigma_veff = 42.6
    sigma_v = 79.1
    pga = 0.25
    rd = 0.96
    gwl = 1.6
    depth = 5.32
    csr = bim14.calculate_csr(sigma_veff, sigma_v, pga, rd, gwl, depth)
    expected_csr = 0.28966
    assert np.isclose(csr, expected_csr)


def test_calculate_cn_values():
    m = 0.58
    sigma_veff = 42.6
    cn = bim14.calculate_cn_values(m, sigma_veff)
    expected_cn = 1.7
    assert np.isclose(cn, expected_cn, rtol=0.04)


def test_calculate_m():
    q_c1ncs = 62.1
    m = bim14.calculate_m(q_c1ncs)
    expected_m = 0.6
    assert np.isclose(m, expected_m, rtol=0.01)


def test_calculate_q_c1n():
    qc = 0.58
    c_n = 1.7
    qc1n = bim14.calculate_q_c1n(qc, c_n)
    # expected_qc1n = 9.73
    expected_qc1n = 9.8
    assert np.isclose(qc1n, expected_qc1n, rtol=0.1)


def test_crr_7p5_from_cpt():
    q_c1n_cs = 62.1
    gwl = 1.6
    depth = 5.35
    expected_crr_7p5 = 0.101
    crr = bim14.crr_7p5_from_cpt(q_c1n_cs, gwl, depth, i_c=1.8)
    assert np.isclose(crr.astype(float), expected_crr_7p5, rtol=0.011)


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
    q_c1ns = 84. * np.ones(1)
    msf = bim14.calculate_msf(magnitude, q_c1ns)
    expected_msf = 1.
    assert np.isclose(expected_msf, msf[0], rtol=0.001)

    magnitude = 6.5
    q_c1ns = 84. * np.ones(1)
    msf = bim14.calculate_msf(magnitude, q_c1ns)
    expected_msf = 1.072
    assert np.isclose(expected_msf, msf[0], rtol=0.001)

    magnitude = 8.5
    q_c1ns = 84. * np.ones(1)
    msf = bim14.calculate_msf(magnitude, q_c1ns)
    expected_msf = 0.944
    assert np.isclose(expected_msf, msf[0], rtol=0.001)

    magnitude = 7.5
    q_c1ns = 175. * np.ones(1)
    msf = bim14.calculate_msf(magnitude, q_c1ns)
    expected_msf = 1.
    assert np.isclose(expected_msf, msf[0])

    magnitude = 6.2
    q_c1ns = 175. * np.ones(1)
    msf = bim14.calculate_msf(magnitude, q_c1ns)
    expected_msf = 1.514
    assert np.isclose(expected_msf, msf[0], rtol=0.001)

    magnitude = 5.1
    q_c1ns = 175. * np.ones(1)
    msf = bim14.calculate_msf(magnitude, q_c1ns)
    expected_msf = 2.1
    assert np.isclose(expected_msf, msf[0], rtol=0.01)


def test_calculate_big_q_values():
    qt = 20.
    sigma_v = 15.
    c_n = 0.8
    qq = bim14.calculate_big_q_values(c_n, qt, sigma_v)
    assert ct.isclose(qq, 0.04, rel_tol=0.011)


def test_calculate_big_f_values():
    qt = 20.
    sigmav = 15.
    fs = 15
    f_value = bim14.calculate_f_ic_values(fs, qt, sigmav)
    expected_f = 300.
    assert f_value == expected_f


def test_calculate_k_sigma():
    sigma_eff = 42.6 * np.ones(1)
    qc1ncs = 9.73 * np.ones(1)
    k_sigma = bim14.calculate_k_sigma(sigma_eff, qc1ncs)
    expected_k_sigma = 1.038
    assert np.isclose(expected_k_sigma, k_sigma[0], rtol=0.001)

    sigma_eff = 42.6 * np.ones(1)
    qc1ncs = 84. * np.ones(1)
    k_sigma = bim14.calculate_k_sigma(sigma_eff, qc1ncs)
    expected_k_sigma = 1.080
    assert np.isclose(expected_k_sigma, k_sigma[0], rtol=0.001)

    sigma_eff = 142.6 * np.ones(1)
    qc1ncs = 84. * np.ones(1)
    k_sigma = bim14.calculate_k_sigma(sigma_eff, qc1ncs)
    expected_k_sigma = 0.9667
    assert np.isclose(expected_k_sigma, k_sigma[0], rtol=0.001)

    sigma_eff = 142.6 * np.ones(1)
    qc1ncs = 130 * np.ones(1)
    k_sigma = bim14.calculate_k_sigma(sigma_eff, qc1ncs)
    expected_k_sigma = 0.9521
    assert np.isclose(expected_k_sigma, k_sigma[0], rtol=0.001)

    sigma_eff = 142.6 * np.ones(1)
    qc1ncs = 212 * np.ones(1)  # Should hit maximum limit
    k_sigma = bim14.calculate_k_sigma(sigma_eff, qc1ncs)
    expected_k_sigma = 0.8934
    assert np.isclose(expected_k_sigma, k_sigma[0], rtol=0.001)

    sigma_eff = 142.6 * np.ones(1)
    qc1ncs = 220 * np.ones(1)  # Should hit maximum limit
    k_sigma = bim14.calculate_k_sigma(sigma_eff, qc1ncs)
    expected_k_sigma = 0.8934
    assert np.isclose(expected_k_sigma, k_sigma[0], rtol=0.001)


def test_calculate_ic():
    big_q = 0.04
    big_f = 300.
    expected_ic = bim14.calculate_ic(big_q, big_f)
    assert ct.isclose(expected_ic, 5.07, rel_tol=0.09)


def test_calculate_qc_1ncs_from_crr_7p5():
    q_c1n_cs_values = np.linspace(50, 220, 100)
    crr_values = bim14.crr_7p5_from_cpt(q_c1n_cs_values, depth=10, gwl=0, i_c=1.8)
    q_c1n_cs_back = bim14.calculate_qc_1ncs_from_crr_7p5(crr_values)
    error = np.sum(abs(q_c1n_cs_values - q_c1n_cs_back))
    assert error < 0.01, error
