import numpy as np
import pytest
import sfsimodels as sm

import liquepy as lq


def test_flac_set_g0_mod():
    sl = lq.num.flac.PM4Sand()
    # sl.phi = 30.
    sl.poissons_ratio = 0.3
    sl.set_g0_mod_from_g_mod_at_v_eff_stress(20.0e6, 50.0e3)
    assert np.isclose(sl.get_g_mod_at_v_eff_stress(50.0e3), 20.0e6)
    assert np.isclose(sl.get_g_mod_at_v_eff_stress(100.0e3), 28284271.24746)
    assert np.isclose(sl.get_g_mod_at_m_eff_stress(50.0e3), 23664319.132398)


def test_write_obj_to_fis_str():
    sl = lq.num.flac.PM4Sand()
    sl.name = "layer_1"
    sl.phi = 30.0
    sl.g_mod = 50.0e3
    sl.specific_gravity = 2.65
    sl.e_curr = 0.55
    sl.poissons_ratio = 0.3
    sl.dilation_angle = 0.0
    sl.cohesion = 0.0
    sl.permeability = 1.0e-4
    sl.set_g0_mod_from_g_mod_at_v_eff_stress(20.0e6, 50.0e3)
    with pytest.raises(sm.ModelError):
        lq.num.flac.write_obj_to_fis_str(sl, sl.inputs, sl.required_parameters)
    sl.h_po = 0.5
    sl.relative_density = 0.5
    p_str = "".join(
        lq.num.flac.write_obj_to_fis_str(sl, sl.inputs, sl.required_parameters)
    )
    assert "layer_1_relative_density=0.5" in p_str

    sp = sm.SoilProfile()
    sp.add_layer(0, sl)
    sp.height = 10.0
    p_str = "".join(lq.num.flac.write_soil_profile_obj_to_fis_str(sp))
    assert "layer_1_relative_density=0.5" in p_str


def test_write_parameters_to_fis_models():
    sl = lq.num.flac.PM4Sand()
    sl.name = "layer_1"
    sl.phi = 30.0
    sl.g_mod = 50.0e3
    sl.specific_gravity = 2.65
    sl.e_curr = 0.55
    sl.poissons_ratio = 0.3
    sl.dilation_angle = 0.0
    sl.cohesion = 0.0
    sl.permeability = 1.0e-4
    sl.set_g0_mod_from_g_mod_at_v_eff_stress(20.0e6, 50.0e3)
    with pytest.raises(sm.ModelError):
        lq.num.flac.write_obj_to_fis_str(
            sl, sl.all_flac_parameters, sl.required_parameters
        )
    sl.h_po = 0.5
    sl.relative_density = 0.5
    p_str = "".join(
        lq.num.flac.write_parameters_to_fis_models(sl, sl.all_flac_parameters)
    )
    assert "D_r=layer_1_relative_density" in p_str


if __name__ == "__main__":
    test_write_parameters_to_fis_models()
