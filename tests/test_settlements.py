import numpy as np
import sfsimodels.models as sm
import sfsimodels.files as sfsi_files
import geofound as gf
import eqsig

from liquepy import settlements as lqs

from tests.conftest import TEST_DATA_DIR


def test_karamitros():

    # Load ground motion
    fpath = TEST_DATA_DIR + "input_acc.his"
    acc_file = np.loadtxt(fpath, skiprows=4)
    acc = acc_file[:, 1]

    # define the soils and soil profile
    sl_0 = sm.Soil()
    sl_0.phi = 0
    sl_0.cohesion = 50000
    sl_0.unit_dry_weight = 19600
    sl_0.unit_sat_weight = 21000

    sl_1 = sm.Soil()
    sl_1.phi = 33.
    sl_1.cohesion = 0
    sl_1.unit_dry_weight = 19700
    sl_1.unit_sat_weight = 21000

    soil_profile = sm.SoilProfile()
    soil_profile.add_layer(0, sl_0)
    soil_profile.add_layer(4.0, sl_1)
    soil_profile.gwl = 2.

    # Define a foundation
    length = 1000000.
    fd = sm.Foundation()
    fd.length = length
    fd.width = 10.0
    fd.depth = 0.0

    q = 80000  # applied load

    z_c = lqs.cal_z_c(fd, z_liq=4, h0=2)
    vertical_effective_stress = soil_profile.vertical_effective_stress(z_c)
    phi_deg = lqs.calc_degraded_phi(sl_1.phi, vertical_effective_stress, q=q)
    assert np.isclose(phi_deg, 9.9275389), phi_deg  # Not validated
    sl_1.phi = phi_deg
    q_ult = gf.capacity_meyerhof_and_hanna_1978(sl_0, sl_1, h0=2, fd=fd, verbose=0)
    assert np.isclose(q_ult, 107350.07398), q_ult  # Not validated
    dt = 0.005

    sett_dyn = lqs.karamitros_settlement(fd, z_liq=4, q=80000, q_ult=q_ult, acc=acc, dt=dt)
    assert np.isclose(sett_dyn, 0.03242937, rtol=0.001), sett_dyn  # 0.03242937 Not validated, liquepy 0.1.11


def test_calc_degraded_phi():
    phi_deg = lqs.calc_degraded_phi(33., 800, q=800)
    assert np.isclose(phi_deg, 10.12554, rtol=0.001)


def test_bray_and_macedo_settlement():
    # Load ground motion
    fpath = TEST_DATA_DIR + "input_acc.his"
    acc_file = np.loadtxt(fpath, skiprows=4)
    acc = acc_file[:, 1]
    time = acc_file[:, 0]
    dt = time[1] - time[0]
    asig = eqsig.AccSignal(acc, dt)

    models = sfsi_files.load_json(TEST_DATA_DIR + "test_ecp_models.json")
    soil_profile = models["soil_profiles"][0]
    building = models["buildings"][0]
    building.mass_ratio = 1.
    fd = models["foundations"][0]
    fd.vertical_load = building.mass_eff * 9.8
    q_c1ncs = 106
    soil_profile.layer(2).q_c1ncs = q_c1ncs
    asig.magnitude = 6.6

    liq_layers = [2]
    sett_dyn_bray = lqs.bray_and_macedo_settlement_time_series(soil_profile, fd, asig, liq_layers)[-1]
    assert np.isclose(sett_dyn_bray, 0.0828104, rtol=0.001), sett_dyn_bray  # 0.0828104 Not validated, liquepy 0.3.1+


# class FakeSignal(object):
#     fake_sa = None
#
#     def __init__(self, values, dt):
#         self.values = values
#         self.dt = dt
#
#     def generate_response_spectrum(self, **kwargs):
#         pass
#

def test_bray_and_macedo_settlement_paper_values():
    case = [1, 4, 11]
    widths = [29.0, 15.0, 38.0]
    qfs = [100.0e3, 60.0e3, 210.0e3]
    sa1 = [0.9, 0.9, 0.9]
    cavdp = [1.0, 1.0, 1.0]
    hl = [12.0, 10.0, 6.0]
    int_lbs = [71.0, 62.0, 14.0]

    # asig = FakeSignal()
    epsilon = -0.5
    expected = [0.1, 0.07, 0.04]
    for i in range(len(widths)):
        pred_sett = lqs.bray_and_macedo_eq(widths[i], qfs[i], hl[i], sa1[i], cavdp[i], int_lbs[i], epsilon)
        assert np.isclose(pred_sett, expected[i], rtol=0.05), pred_sett
    pass


def test_lu_settlement():

    # Load ground motion
    fpath = TEST_DATA_DIR + "input_acc.his"
    acc_file = np.loadtxt(fpath, skiprows=4)
    acc = acc_file[:, 1]
    acc = acc / 9.81  # TODO: confirm should be in m/s2

    # Define a foundation
    length = 1000000.
    fd = sm.Foundation()
    fd.length = length
    fd.width = 10.0
    fd.depth = 0.0

    sett_dyn = lqs.lu_settlements(q=80, fd=fd, Dr=55, acc=acc)
    assert np.isclose(sett_dyn, 0.366247, rtol=0.001)  # 0.366247 Not validated, liquepy 0.1.0


if __name__ == '__main__':
    test_bray_and_macedo_settlement_paper_values()

# if __name__ == '__main__':
#     test_karamitros()