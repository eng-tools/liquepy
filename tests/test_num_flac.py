import liquepy as lq
import numpy as np


def test_flac_set_g0_mod():
    sl = lq.num.flac.PM4Sand()
    sl.phi = 30.
    sl.set_g0_mod_from_g_mod_at_v_eff_stress(20.0e6, 50.0e3)
    assert np.isclose(sl.get_g_mod_at_v_eff_stress(50.0e3), 20.0e6)
    assert np.isclose(sl.get_g_mod_at_v_eff_stress(100.0e3), 28284271.24746)
    assert np.isclose(sl.get_g_mod_at_m_eff_stress(50.0e3), 23094010.767585)
    k0 = 1 - np.sin(sl.phi_r)
    m_eff = (1 + k0) / 2 * 50.0e3
    assert np.isclose(sl.get_g_mod_at_m_eff_stress(m_eff), 20.0e6)



if __name__ == '__main__':
    test_flac_set_g0_mod()
