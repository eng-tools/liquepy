import liquepy as lq
import numpy as np
from tests.conftest import TEST_DATA_DIR


def test_fit_single_layer():
    depths = np.arange(0, 10.1, 0.1)
    crr = 0.6 * np.ones_like(depths)
    crr[10:20] = 0.15
    crr_opts = np.array([0.6, 0.35, 0.2])

    d_liqs, d_nonliqs, csr_n15s, normed_diff = lq.esp.millen_2020.fit_n_layer_profile(depths, crr, n=3, crr_n15_opts=crr_opts)
    assert np.isclose(d_liqs[0], 1.0), d_liqs[0]
    assert np.isclose(d_nonliqs[0], 2.0), d_nonliqs[0]
    assert np.isclose(csr_n15s[0], 0.15), csr_n15s[0]

    d_liqs, d_nonliqs, csr_n15s, normed_diff = lq.esp.millen_2020.fit_n_layer_profile(depths, crr, n=5, crr_n15_opts=crr_opts)

    assert np.isclose(d_liqs[0], 1.0), d_liqs[0]
    assert np.isclose(d_nonliqs[0], 2.0), d_nonliqs[0]
    assert np.isclose(csr_n15s[0], 0.15), csr_n15s[0]

    assert np.isclose(d_liqs[1], 10), d_liqs[1]
    assert np.isclose(d_nonliqs[1], 10), d_nonliqs[1]
    assert np.isclose(csr_n15s[1], 0.6, rtol=1e-3), csr_n15s[1]


def test_fit_many():
    depths = np.arange(0, 10.1, 0.1)
    crr = 0.6 * np.ones_like(depths)
    crr[10:20] = 0.54
    crr[38: 50] = 0.3
    crr[45: 50] = 0.4
    d_liqs, d_nonliqs, csr_n15s, normed_diff = lq.esp.millen_2020.fit_n_layer_profile(depths, crr, n=3)
    assert np.isclose(d_liqs[0], 3.8), d_liqs[0]
    assert np.isclose(d_nonliqs[0], 5.0), d_nonliqs[0]
    assert np.isclose(csr_n15s[0], 0.3), csr_n15s[0]

    d_liqs, d_nonliqs, csr_n15s, normed_diff = lq.esp.millen_2020.fit_n_layer_profile(depths, crr, n=5)
    assert np.isclose(d_liqs[0], 1.0), d_liqs[0]
    assert np.isclose(d_nonliqs[0], 2.0), d_nonliqs[0]
    assert np.isclose(csr_n15s[0], 0.54), csr_n15s[0]

    assert np.isclose(d_liqs[1], 3.8), d_liqs[1]
    assert np.isclose(d_nonliqs[1], 5.0), d_nonliqs[1]
    assert np.isclose(csr_n15s[1], 0.3), csr_n15s[1]


def test_fit_adjacent_layers():
    depths = np.arange(0, 10.1, 0.1)
    crr = 0.6 * np.ones_like(depths)
    crr[38: 45] = 0.3
    crr[46: 60] = 0.4
    crr_opts = np.array([0.6, 0.35, 0.2])
    d_liqs, d_nonliqs, csr_n15s, normed_diff = lq.esp.millen_2020.fit_n_layer_profile(depths, crr, n=3, crr_n15_opts=crr_opts)
    assert np.isclose(d_liqs[0], 3.8), d_liqs[0]
    assert np.isclose(d_nonliqs[0], 6.0), d_nonliqs[0]
    assert np.isclose(csr_n15s[0], 0.4), csr_n15s[0]

    d_liqs, d_nonliqs, csr_n15s, normed_diff = lq.esp.millen_2020.fit_n_layer_profile(depths, crr, n=5, crr_n15_opts=crr_opts)

    assert np.isclose(d_liqs[0], 3.8), d_liqs[0]
    assert np.isclose(d_nonliqs[0], 4.5), d_nonliqs[0]
    assert np.isclose(csr_n15s[0], 0.3), csr_n15s[0]

    assert np.isclose(d_liqs[1], 4.6), d_liqs[1]
    assert np.isclose(d_nonliqs[1], 6.0), d_nonliqs[1]
    assert np.isclose(csr_n15s[1], 0.4), csr_n15s[1]


def test_fit():
    depths = np.arange(0, 10, 0.1)
    crr = 0.6 * np.ones_like(depths)
    crr[10:20] = 0.15
    crr[38: 50] = 0.3
    crr[45: 50] = 0.4
    d_liqs, d_nonliqs, csr_n15s, normed_diff = lq.esp.millen_2020.fit_n_layer_profile(depths, crr, n=3)
    assert np.isclose(d_liqs[0], 1.0), d_liqs[0]
    assert np.isclose(d_nonliqs[0], 2.0), d_nonliqs[0]
    assert np.isclose(csr_n15s[0], 0.15), csr_n15s[0]

    crr = 0.6 * np.ones_like(depths)
    crr[10:30] = 0.15
    crr[43: 50] = 0.3
    crr[45: 50] = 0.4
    d_liqs, d_nonliqs, csr_n15s, normed_diff = lq.esp.millen_2020.fit_n_layer_profile(depths, crr, n=3)
    # h_crusth_liq, p_value, normed_diff
    assert np.isclose(d_liqs[0], 1.0), d_liqs[0]
    assert np.isclose(d_nonliqs[0], 3.0), d_nonliqs[0]
    assert np.isclose(csr_n15s[0], 0.15), csr_n15s[0]


def disable_test_fit_3layer_from_cpt():
    cpt = lq.field.load_mpa_cpt_file(TEST_DATA_DIR + "standard_1.csv")
    bi2014 = lq.trigger.run_bi2014(cpt, pga=0.25, m_w=7.5, gwl=cpt.gwl, unit_wt_method='void_ratio')

    n_layers = 3
    d_liqs, d_nonliqs, csr_n15s, normed_diff = lq.esp.millen_2020.fit_n_layer_profile(bi2014.depth, bi2014.crr_m7p5, n=n_layers)

    print(d_liqs, d_nonliqs, csr_n15s, normed_diff)
    assert np.isclose(d_liqs[0], 14.15), d_liqs[0]
    assert np.isclose(d_nonliqs[0], 16.33), d_nonliqs[0]
    assert np.isclose(csr_n15s[0], 0.15403698, rtol=0.001), csr_n15s[0]


def disabled_test_fit_5layer_from_cpt():
    cpt = lq.field.load_cpt_from_file(TEST_DATA_DIR + "standard_1.csv")
    bi2014 = lq.trigger.run_bi2014(cpt, pga=0.25, m_w=7.5, gwl=cpt.gwl, unit_wt_method='void_ratio')

    n_layers = 5
    d_liqs, d_nonliqs, csr_n15s, normed_diff = lq.esp.millen_2020.fit_n_layer_profile(bi2014.depth, bi2014.crr_m7p5, n=n_layers)
    show = 0
    if show:
        import matplotlib.pyplot as plt
        bf, sps = plt.subplots(nrows=2)
        crr_m7p5 = np.clip(bi2014.crr_m7p5, None, 0.6)
        sps[0].plot(bi2014.depth, crr_m7p5, drawstyle='steps-post')
        sps[0].plot([d_liqs[0], d_nonliqs[0]], [csr_n15s[0], csr_n15s[0]], c="k", ls="--")
        if n_layers == 5:
            sps[0].plot([d_liqs[1], d_nonliqs[1]], [csr_n15s[1], csr_n15s[1]], c="k", ls="--")
        # sps[0].plot([d_nonliq, d_liq2], [crr2, crr2], c="k", ls="--")
        print(d_liqs, d_nonliqs, csr_n15s)
        plt.show()
    assert np.isclose(d_liqs[0], 4.9), d_liqs[0]
    assert np.isclose(d_nonliqs[0], 8.17), d_nonliqs[0]
    assert np.isclose(csr_n15s[0], 0.247002, rtol=0.001), csr_n15s[0]
    assert np.isclose(d_liqs[1], 14.15), d_liqs[1]
    assert np.isclose(d_nonliqs[1], 16.33), d_nonliqs[1]
    assert np.isclose(csr_n15s[1], 0.154037, rtol=0.001), csr_n15s[1]


def run():
    cpt = lq.field.load_mpa_cpt_file(TEST_DATA_DIR + "standard_1.csv")
    bi2014 = lq.trigger.run_bi2014(cpt, pga=0.25, m_w=7.5, gwl=cpt.gwl, unit_wt_method='void_ratio')
    import matplotlib.pyplot as plt
    n_layers = 3
    nsteps = [2, 6, 10, 15]
    for nstep in nsteps:
        crr_opts = np.linspace(0.5, 0.6, nstep)[::-1]
        d_liqs, d_nonliqs, csr_n15s, normed_diff = lq.esp.millen_2020.fit_n_layer_profile(bi2014.depth, bi2014.crr_m7p5, n=n_layers, crr_n15_opts=crr_opts)

        plt.plot(nstep, d_liqs[0], 'o', c='r')
        plt.plot(nstep, d_nonliqs[0], 'o', c='b')
        plt.plot(nstep, csr_n15s[0], 'o', c='g')
        print(normed_diff)
    plt.show()

if __name__ == '__main__':
    run()

