import numpy as np

from scipy.integrate import trapz
import eqsig

import liquepy as lq


def calculate_factor_safety(q_c1ncs, p_a, magnitude, pga, depth, soil_profile):
    """
    Calculate the liquefaction factor of safety at a given depth.

    :param q_c1ncs: float, normalised cone tip resistance corrected to equivalent clean sand
    :param p_a: float, atmospheric pressure
    :param magnitude: float, earthquake m_w
    :param pga: float, peak ground acceleration
    :param depth: float, depth from surface
    :param soil_profile: SoilProfile, A soil profile object
    :return:
    """

    crr_m7p5 = np.exp(q_c1ncs / 113 + (q_c1ncs / 1000) ** 2 - (q_c1ncs / 140) ** 3 + (q_c1ncs / 137) ** 4 - 2.8)

    c_sigma = 1.0 / (37.3 - (8.27 * (q_c1ncs ** 0.264)))

    sigma_v = soil_profile.get_v_total_stress_at_depth(depth)
    sigma_veff = soil_profile.get_v_eff_stress_at_depth(depth)

    k_sigma = np.clip(1.0 - c_sigma * np.log(sigma_veff / p_a), -1000, 1.1)

    msf_max = 1.09 + (q_c1ncs / 180) ** 3
    msf = 1 + ((msf_max - 1) * ((8.64 * np.exp(-magnitude / 4)) - 1.325))

    alpha = -1.012 - 1.126 * np.sin(depth / 11.73 + 5.133)
    beta = 0.106 + 0.118 * np.sin(depth / 11.28 + 5.142)
    r_d = np.exp(alpha + (beta * magnitude))
    csr = 0.65 * pga * sigma_v / sigma_veff * r_d

    fs_liq = crr_m7p5 * k_sigma * msf / csr

    return fs_liq


def calc_degraded_phi(phi, sigma_v_eff, q, a=0.9, ru_ff=1.):
    """
    Equivalent degraded friction angle of liquefied soil under a foundation

    Ref: Cascone and Bouckovalas (1998)

    :param phi: float, friction angle
    :param sigma_v_eff: float, vertical effective stress
    :param q: float, bearing pressure of foundation
    :param a: float, adjustment parameter
    :param ru_ff: float, pore pressure ratio in the free-field
    :return:
    """
    u_foot = a / (1 + (q / sigma_v_eff))
    big_u = (u_foot + ru_ff) / 2  # for strip foundations
    degraded_phi = np.degrees(np.arctan((1 - big_u) * np.tan(np.deg2rad(phi))))
    return degraded_phi


def cal_z_c(fd, z_liq, h0):
    """
    Calculation of characteristic depth from Karamitros et al. (2013)
    :param fd:
    :param z_liq:
    :param h0:
    :return:
    """
    if fd.width > z_liq:
        z_c = h0 + z_liq
    else:
        z_c = h0 + fd.width
    return z_c


def karamitros_settlement(fd, z_liq, q, q_ult, acc, dt):
    """
    Calculate the settlement using the method proposed by Karamitros et al. 2013 - sett

    :param sss:
    :return:
    """
    sett_dyn_ts = karamitros_settlement_time_series(fd, z_liq, q, q_ult, acc, dt)
    return sett_dyn_ts[-1]


def karamitros_settlement_time_series(fd, z_liq, q, q_ult, acc, dt, c_dash=0.003):  # units: m, Pa, s
    """
    Calculate the settlement using the method proposed by :cite:`Karamitros:2013gi`

    Parameters
    ----------
    fd: sfsimodels.Foundation
    z_liq
    q
    q_ult
    acc
    dt
    c_dash

    Returns
    -------

    """
    asig = eqsig.AccSignal(acc, dt)
    fd_q_ult = q_ult
    fd_q_demand = q
    return calc_settlement_karamitros_et_al_2013(fd, asig, fd_q_ult, fd_q_demand, z_liq, c_dash=c_dash)


def calc_settlement_karamitros_et_al_2013(fd, asig, fd_q_ult, fd_q_demand, y_liq, c_dash=0.003):
    """
    Calculate the settlement using the method proposed by :cite:`Karamitros:2013gi`

    Parameters
    ----------
    fd: sfsimodels.Foundation
    asig: eqsig.AccSignal
    fd_q_ult: float or array_like
        Bearing capacity of foundation
    fd_q_demand: float or array_like
        Bearing load on foundation
    y_liq: float
        Depth to liquefaction
    c_dash: float (default 0.003)

    Returns
    -------
    array_like
    """

    c_factor = min(c_dash * (1.0 + 1.65 * fd.length / fd.width), 11.65 * c_dash)  # Karamitros 2013 sett

    int_vel = eqsig.im.calc_integral_of_abs_velocity(asig)
    amax_t2_n = (np.pi ** 2) * int_vel
    fs_deg = fd_q_ult / fd_q_demand
    sett_dyn_ts = c_factor * amax_t2_n * (y_liq / fd.width) ** 1.5 * (1.0 / fs_deg) ** 3
    return sett_dyn_ts


def bray_and_macedo_settlement(soil_profile, fd, asig, liq_layers):
    """
    Calculates foundation settlement using Bray and Macedo (2017)

    :param acc: array, acceleration time series
    :param dt: float, time step of acceleration time series
    :param z_liq:
    :param q: float, foundation bearing pressure
    :param fd: Foundation, foundation object
    :param soil_profile: SoilProfile, soil profile object
    :return:
    """

    sett_dyn_ts = bray_and_macedo_settlement_time_series(soil_profile, fd, asig, liq_layers)
    return sett_dyn_ts[-1]


def bray_and_macedo_settlement_time_series(soil_profile, fd, asig, liq_layers):
    """
    Calculates foundation settlement using Bray and Macedo (2017)

    :param acc: array, acceleration time series
    :param dt: float, time step of acceleration time series
    :param z_liq:
    :param q: float, foundation bearing pressure
    :param fd: Foundation, foundation object
    :param soil_profile: SoilProfile, soil profile object
    """
    gravity = 9.8
    # calculation of CAVdp
    cavdp_time_series = eqsig.im.calc_cav_dp(asig)
    pga_max = max(abs(asig.values)) / gravity
    hl = 0
    q_c1ncs_values = []
    for layer_id in liq_layers:
        hl += soil_profile.layer_height(layer_id)
        q_c1ncs_values.append(soil_profile.layer(layer_id).q_c1ncs)
    q_c1ncs = np.mean(q_c1ncs_values)
    # calculation of LBS

    # calculation of Maximum Cyclic Shear Strains
    z = np.arange((soil_profile.layer_depth(2)) + 0.5, (soil_profile.layer_depth(3) + 0.5), 0.5)
    xmax = len(z)-1
    lbs = []

    for depth in z:

        fs = calculate_factor_safety(q_c1ncs=q_c1ncs, p_a=101000, magnitude=asig.magnitude, pga=pga_max, depth=depth, soil_profile=soil_profile)
        d_r = soil_profile.layer(2).relative_density
        e_shear = lq.trigger.calc_shear_strain(fs=fs, d_r=d_r)
        if depth < fd.depth:
            w = 0
        else:
            w = 1
        lbs1 = w * e_shear / depth
        lbs.append(lbs1)

    x_lower = z[0]  # the lower limit of x
    x_upper = z[xmax]  # the upper limit of x
    x_int = z[np.where((x_lower <= z) * (z <= x_upper))]
    y_int = np.abs(np.array(lbs)[np.where((x_lower <= z) * (z <= x_upper))])
    int_lbs = trapz(y_int, x_int)  # lbs value

    asig.generate_response_spectrum(response_times=np.array([1.]), xi=0.05)
    sa1 = asig.s_a[0] / gravity
    qf = fd.vertical_load / fd.width / fd.length
    return bray_and_macedo_eq(fd.width, qf, hl, sa1, cavdp_time_series, int_lbs)


def bray_and_macedo_eq(width, qf, hl, sa1, cavdp_time_series, int_lbs, epsilon=0.0):
    # calculation of c_1 and c_2
    if int_lbs <= 16:
        c_1 = -8.35
        c_2 = 0.072
    else:
        c_1 = -7.48
        c_2 = 0.014

    q = qf / 1000

    sett_dyn_ts = np.exp(c_1 + (4.59 * np.log(q)) - (0.42 * ((np.log(q)) ** 2)) + (c_2 * int_lbs) + (0.58 * np.log(np.tanh(hl / 6))) - (0.02 * width) + (0.84 * np.log(cavdp_time_series)) + (0.41 * np.log(sa1)) + epsilon)

    sett_dyn_ts = sett_dyn_ts/1000

    return sett_dyn_ts  # TODO: Should return metres not millimetres


def lu_settlements(q, fd, Dr, acc):

    # TODO: q should be in Pa not kPa
    # TODO: DR should be a ratio

    Dr_1b=[30.057, 32.004, 34.065, 35.954, 37.958, 39.962, 41.966, 43.969, 45.973, 47.920, 49.866, 51.985, 53.817, 55.935, 57.939, 59.943]
    N_lr_1b = [204.739, 208.531, 212.322, 216.114, 219.905, 231.280, 242.654, 257.820, 276.777, 299.526, 329.858, 382.938, 428.436, 496.682, 561.137, 636.967]

    Dr_2a = [30.000, 32.004, 33.950, 35.954, 38.187, 40.019, 41.966, 43.969, 46.031, 48.034, 49.981, 51.927, 53.989, 55.992, 57.882, 59.943]
    N_lr_2a = [37.915, 37.915, 45.498, 53.081, 56.872, 60.664, 75.829, 79.621, 90.995, 98.578, 98.578, 109.953, 117.536, 128.910, 140.284, 155.450]

    Dr_2b = [30.115, 31.947, 34.008, 36.298, 37.901, 40.076, 42.137, 43.511, 45.744, 47.748, 49.122, 51.126, 53.130, 55.649, 57.424, 59.943]
    N_lr_2b = [299.526, 299.526, 318.483, 329.858, 337.441, 367.773, 401.896, 424.645, 492.891, 553.555, 587.678, 659.716, 724.171, 800.000, 845.498, 936.493]

    Dr_3a = [30.115, 31.947, 34.179, 35.954, 38.073, 40.019, 42.023, 43.969, 46.088, 48.092, 50.095, 51.985, 53.989, 55.763, 57.767, 59.828]
    N_lr_3a = [60.664, 60.664, 68.246, 75.829, 83.412, 90.995, 98.578, 106.161, 121.327, 128.910, 151.659, 170.616, 189.573, 216.114, 250.237, 288.152]

    Dr_3b = [30.057, 31.897, 33.678, 35.230, 37.356, 39.483, 41.034, 42.414, 43.563, 46.264, 49.080, 50.862, 52.644, 54.713, 57.011, 60.000]
    N_lr_3b = [778.202, 793.461, 820.164, 831.608, 865.940, 900.272, 926.975, 946.049, 972.752, 1022.343, 1079.564, 1121.526, 1163.488, 1205.450, 1251.226, 1319.891]

    x_int = Dr

    abs_acc = abs(acc)
    pga = (max(abs_acc))
    x_pga = [0.1, 0.4]

    if 10 <= q < 30:

        N_lr_1_b = np.interp(x_int, Dr_1b, N_lr_1b)

        y_nlr = [N_lr_1_b, 0]

        logz = np.log10(pga)
        logx = np.log10(x_pga)
        logy = np.log10(y_nlr)
        N_lr = np.power(10.0, np.interp(logz, logx, logy))

    elif 30 <= q < 80:
        N_lr_2_a = np.interp(x_int, Dr_2a, N_lr_2a)
        N_lr_2_b = np.interp(x_int, Dr_2b, N_lr_2b)
        y_nlr = [N_lr_2_b, N_lr_2_a]

        logz = np.log10(pga)
        logx = np.log10(x_pga)
        logy = np.log10(y_nlr)
        N_lr = np.power(10.0, np.interp(logz, logx, logy))

    elif 80 <= q <= 120:
        N_lr_3_a = np.interp(x_int, Dr_3a, N_lr_3a)
        N_lr_3_b = np.interp(x_int, Dr_3b, N_lr_3b)
        y_nlr = [N_lr_3_b, N_lr_3_a]

        logz = np.log10(pga)
        logx = np.log10(x_pga)
        logy = np.log10(y_nlr)
        N_lr = np.power(10.0, np.interp(logz, logx, logy))
    else:
        raise ValueError("q value ({0}) out of range (10-120)".format(q))

    c_d = 1-(fd.depth/(4*fd.width))

    sett_dyn_lu = c_d * (q / N_lr) * ((fd.width / (fd.width + 0.33)) ** 2)

    return sett_dyn_lu



