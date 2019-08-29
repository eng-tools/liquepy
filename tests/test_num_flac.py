import liquepy as lq
import numpy as np


def test_flac_set_g0_mod():
    sl = lq.num.flac.PM4Sand()
    sl.phi = 30.
    sl.set_g0_mod_from_g_mod_at_v_eff_stress(20.0e6, 50.0e3)
    assert np.isclose(sl.get_g_mod_at_v_eff_stress(50.0e3), 20.0e6)


if __name__ == '__main__':
    test_flac_set_g0_mod()
