import numpy as np
import pytest

import liquepy as lq
from tests.conftest import TEST_DATA_DIR


def test_calc_lsn():
    cpt = lq.field.load_mpa_cpt_file(TEST_DATA_DIR + "standard_1.csv")
    bi2014 = lq.trigger.run_bi2014(cpt, pga=0.25, m_w=7.5, gwl=cpt.gwl)
    epsilon = lq.trigger.calc_volumetric_strain_zhang_2004(
        bi2014.factor_of_safety, bi2014.q_c1n_cs
    )
    lsn_direct = lq.trigger.calc_lsn(epsilon * 100, cpt.depth)
    assert np.isclose(lsn_direct, 36.0919293469645, rtol=0.01)  # v0.5.5


def test_single_calc_lpi_increments():
    depth = 6.61
    fos = 0.45
    lpi_inc = lq.trigger.triggering_measures.calc_lpi_increments(
        np.ones(2) * fos, np.array([depth, depth + 0.01])
    )
    assert np.isclose(lpi_inc[1], 0.0368, rtol=0.01), lpi_inc[
        1
    ]  # unvalidated test value


@pytest.mark.parametrize(
    "depth, fos, lpi",
    [
        (0.98, 2.0, 0.0),
        (2.18, 0.686, 0.0280),
        (3.2, 0.59, 0.0344),
        (6.61, 0.45, 0.0368),
        (12.390, 1.413, 0.0),
    ],
)
def test_calc_lpi_increments(depth, fos, lpi):
    lpi_inc = lq.trigger.triggering_measures.calc_lpi_increments(
        np.ones(2) * fos, np.array([depth, depth + 0.01])
    )
    assert np.isclose(lpi_inc[1], lpi, rtol=0.01), lpi_inc[1]  # unvalidated test value


if __name__ == "__main__":
    test_calc_lpi_increments()
