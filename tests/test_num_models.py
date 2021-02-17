import sfsimodels as sm
import liquepy as lq
import json


def test_save_load_w_diff_wmd():
    sl = lq.num.models.ManzariDafaliasModel(liq_mass_density=1., p_atm=101.)
    sl.e_curr = 0.8
    sl.g0 = 125
    sl.poissons_ratio = 0.3

    sl.m_c = 1.214  # From Cubrinovski 1998
    sl.c_c = 0.712
    sl.lambda_c = 0.019
    sl.e_0 = 0.934
    sl.ksi = 0.7
    sl.p_atm = 101
    sl.m_yield = 0.01
    sl.h_0 = 7.05
    sl.c_h = 0.968
    sl.n_b = 1.1
    sl.a_0 = 0.704
    sl.n_d = 3.5
    sl.z_max = 4
    sl.c_z = 600
    sl.specific_gravity = 2.65
    print(sl.unit_dry_weight)
    sl.cohesion = 0

    eo = sm.Output()
    eo.add_to_dict(sl)
    p_str = json.dumps(eo.to_dict(), indent=4)
    mods = sm.loads_json(p_str, {'soil-manzaridafalias_model': lq.num.models.ManzariDafaliasModel})


if __name__ == '__main__':
    test_save_load_w_diff_wmd()