import liquepy as lq
import numpy as np
import pysra
import eqsig
import sfsimodels as sm
import matplotlib.pyplot as plt
from bwplot import cbox
import engformat as ef


def run_sra(sp, m, odepths, analysis="linear", d_inc=None):
    if d_inc is None:
        d_inc = np.ones(sp.n_layers)

    profile = lq.sra.sm_profile_to_pysra(sp, d_inc=d_inc)
    if analysis == "linear":
        calc = pysra.propagation.LinearElasticCalculator()
    else:
        calc = pysra.propagation.EquivalentLinearCalculator()
    od = {}
    outs = []
    for i, depth in enumerate(odepths):
        od["upACCX_d%i" % i] = len(outs)
        outs.append(pysra.output.AccelerationTSOutput(pysra.output.OutputLocation('incoming_only', depth=depth)))
        od["ACCX_d%i" % i] = len(outs)
        outs.append(pysra.output.AccelerationTSOutput(pysra.output.OutputLocation('within', depth=depth)))
        od["STRS_d%i" % i] = len(outs)
        outs.append(pysra.output.StrainTSOutput(pysra.output.OutputLocation('within', depth=depth), in_percent=False))
        od["TAU_d%i" % i] = len(outs)
        outs.append(pysra.output.StressTSOutput(pysra.output.OutputLocation('within', depth=depth),
                                                normalized=False))

    outputs = pysra.output.OutputCollection(outs)

    # Perform the calculation
    calc(m, profile, profile.location('outcrop', depth=sp.height))
    outputs(calc)

    out_series = {}
    for item in od:
        if "TAU" in item:
            out_series[item] = outputs[od[item]].values
        else:
            out_series[item] = outputs[od[item]].values

    return out_series, profile, calc


def analyse(sp, m):
    odepths = np.arange(0.5, int(sp.height - 1), 0.5)

    oseries, profile, calc = run_sra(sp, m, odepths)

    aew = []
    acc_uke = []
    up_uke = []
    down_uke = []

    for i, depth in enumerate(odepths):
        acc_signal = eqsig.AccSignal(oseries["ACCX_d%i" % i] * 9.8, m.time_step)
        g_mod = profile.location('outcrop', depth=depth).layer.initial_shear_mod
        rho = profile.location('outcrop', depth=depth).layer.unit_wt / 9.8

        tau = oseries["TAU_d%i" % i]
        strain_energy = 0.5 * tau ** 2 / g_mod

        delta_st_energy = np.diff(strain_energy)
        delta_st_energy = np.insert(delta_st_energy, 0, 0)
        cum_delta_st_energy = np.cumsum(abs(delta_st_energy))

        up_acc_signal = eqsig.AccSignal(oseries["upACCX_d%i" % i] * 9.8, m.time_step)
        down_acc_signal = eqsig.AccSignal((oseries["upACCX_d%i" % i] - oseries["ACCX_d%i" % i]) * 9.8, m.time_step)

        acc_uke.append(eqsig.im.calc_unit_kinetic_energy(acc_signal)[-1] * rho)
        up_uke.append(eqsig.im.calc_unit_kinetic_energy(up_acc_signal)[-1] * rho)
        down_uke.append(eqsig.im.calc_unit_kinetic_energy(down_acc_signal)[-1] * rho)

        aew.append(cum_delta_st_energy[-1])

    cake = np.array(aew)
    case = np.array(acc_uke)

    return odepths, case, cake


def build_soil_profile(thicknesses, g_mods, xis, unit_weights, add_rock=0):
    depths = np.cumsum(thicknesses)
    depths = np.insert(depths, 0, 0)

    sp = sm.SoilProfile()
    sp.gwl = 10000
    sp.height = depths[-1]
    sp.hydrostatic = True
    for i in range(len(unit_weights)):
        sl = sm.Soil()
        sl.unit_dry_weight = unit_weights[i]
        sl.g_mod = g_mods[i]
        sl.specific_gravity = 2.65
        sl.xi = xis[i]
        sp.add_layer(depths[i], sl)

    if add_rock:
        rock = sm.Soil()
        rock.unit_dry_weight = 17500.
        rock.specific_gravity = 2.65
        rock.g_mod = 100.0e6
        rock.xi = 0.01

        sp.height = depths[-1] + 1.0
        sp.add_layer(depths[-1], rock)
    return sp


def create(save=0, show=0):

    bf, sps = plt.subplots(ncols=2, nrows=1, figsize=(6.5, 5), squeeze=False)
    sps = sps.flatten()

    # Define soil profile
    damp = 0.03
    thicknesses = np.array([10, 10, 10])
    g_mods = np.array([40.0e6, 30.0e6, 40.0e6])
    unit_weights = [20000., 20000., 20000.]
    xis = np.ones(3) * damp
    sp = build_soil_profile(thicknesses, g_mods, xis, unit_weights)

    sps[-1].axvline(0, c="k", ls="--")

    asig = eqsig.load_asig('test_motion_dt0p01.txt')

    # Should be input as g
    m = pysra.motion.TimeSeriesMotion(filename='test motion', description=None, time_step=asig.dt,
                                      accels=asig.values / 9.8)

    in_signal = eqsig.AccSignal((m.accels / 2) * 9.8, m.time_step)  # Should be input as m/s2
    rho = sp.layer(3).unit_dry_weight / 9.8
    in_uke = eqsig.im.calc_unit_kinetic_energy(in_signal)[-1]
    in_cake = in_uke * rho
    odepths, cake, case = analyse(sp, m)

    sps[1].plot(cake / in_cake, odepths, ls="-", ms=2, label="Linear $\\xi={0:.0f}$%".format(damp * 100), c=cbox(0))
    sps[0].plot(case / in_cake, odepths, ls="-", label="Linear $\\xi={0:.0f}$%".format(damp * 100), c=cbox(0), lw=1)

    pred_case = lq.trigger.nses.est_case_1d_millen_et_al_2019(sp, in_signal, odepths, xi=sp.layer(1).xi, g_scale_limit=4)[:, -1]
    sps[0].plot(pred_case / in_cake, odepths, 's', label="NSES (base)", c=cbox(1))
    pred_cake = lq.trigger.nses.est_case_1d_millen_et_al_2019(sp, in_signal, odepths, xi=sp.layer(1).xi,
                                                              g_scale_limit=4, nodal=False)[:, -1]
    sps[1].plot(pred_cake / in_cake, odepths, 's', label="ASES (base)", c=cbox(1))

    # get surface motion
    oseries, profile, calc = run_sra(sp, m, [0])
    acc_signal = eqsig.AccSignal(oseries["upACCX_d%i" % 0] * 9.8, m.time_step)
    pred_case = lq.trigger.nses.est_case_1d_millen_et_al_2019(sp, acc_signal, odepths, xi=sp.layer(1).xi, in_loc=0, start=False)[:, -1]
    sps[0].plot(pred_case / in_cake, odepths, "x", ms=3, label="NSES (surf)", c=cbox(3))
    pred_cake = lq.trigger.nses.est_case_1d_millen_et_al_2019(sp, acc_signal, odepths, xi=sp.layer(1).xi, in_loc=0,
                                                              start=False, nodal=False)[:, -1]
    sps[1].plot(pred_cake / in_cake, odepths, "x", ms=3, label="ASES (surf)", c=cbox(3))

    sps[0].legend().set_zorder(1012)
    ef.revamp_legend(sps[0], prop={"size": 8}, loc="lower right")

    sps[0].set_ylabel("Depth [m]")
    sps[-1].set_xlabel("$\\frac{CASE - CASE_{pred}}{CASE_{pred}}$")
    ef.revamp_legend(sps[-1], prop={"size": 8}, loc="lower left")
    ef.xy(sps[-1], y_origin=True)
    sps[0].set_xlabel("$\\frac{CASE}{CAKE_{in}}$")
    sps[1].set_xlabel("$\\frac{CAKE}{CAKE_{in}}$")
    for i in range(len(sps)):

        ef.xy(sps[i], x_origin=True, y_origin=True)

        sps[i].set_ylim([0, 30])
        sps[i].invert_yaxis()
        if i != 0:
            sps[i].tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

    plt.tight_layout()
    bf.subplots_adjust(wspace=0.05, hspace=0)

    plt.show()


if __name__ == '__main__':
    create()
