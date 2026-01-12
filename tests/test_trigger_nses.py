import eqsig
import numpy as np
import sfsimodels as sm

from liquepy.trigger import nses
from tests import conftest


def test_est_nses_millen_et_al_2019():
    # Define soil profile properties
    damp = 0.03
    sp = sm.SoilProfile()  # Soil profile object
    # Top layer
    sl1 = sm.Soil(g_mod=40.0e6, unit_dry_weight=20.0e3)
    sl1.xi = damp
    sp.add_layer(depth=0, soil=sl1)
    # Middle layer - with lower shear modulus
    sl2 = sm.Soil(g_mod=30.0e6, unit_dry_weight=20.0e3)
    sl2.xi = damp
    sp.add_layer(depth=10, soil=sl2)
    # Bottom layer
    sl3 = sm.Soil(g_mod=40.0e6, unit_dry_weight=20.0e3)
    sl3.xi = damp
    sp.add_layer(depth=20, soil=sl3)
    sp.height = 30  # m

    # Load ground motion
    acc = np.loadtxt(conftest.TEST_DATA_DIR + "test_motion_dt0p01.txt", skiprows=2)
    dt = 0.01

    in_signal = eqsig.AccSignal(acc / 2, dt)
    rho = sp.layer(3).unit_dry_weight / 9.8
    in_uke = eqsig.im.calc_unit_kinetic_energy(in_signal)[-1]
    in_cake = in_uke * rho

    # Estimate CASE and CAKE using the input motion with method from Millen et al. (2019)
    odepths = np.array([4.0])
    pred_case = nses.est_case_1d_millen_et_al_2019(sp, in_signal, odepths, xi=damp)[
        0, -1
    ]
    pred_cake = nses.est_case_1d_millen_et_al_2019(
        sp, in_signal, odepths, xi=damp, nodal=False
    )[0, -1]
    assert np.isclose(pred_cake, 729.00213882), pred_cake  # v0.6.3
    assert np.isclose(pred_case, 40.094265976), pred_case  # v0.6.3


if __name__ == "__main__":
    test_est_nses_millen_et_al_2019()
