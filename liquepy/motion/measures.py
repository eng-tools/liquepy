import eqsig
import numpy as np
from numpy import trapz
from liquepy.exceptions import deprecation


def dep_calculate_cav_dp(acc, dt):
    start = 0
    pga_max = 0
    CAVdp = 0
    num_points = (int(1 / dt))
    total_time = int(dt * (len(acc) - 1))

    for i in range(0, total_time):

        end = start + num_points

        interval_total_time = (start * dt) + 1
        interval_time = np.arange(start * dt, interval_total_time, dt)

        acc_interval = []
        for j in range(start, end + 1):
            acc_interval.append(acc[j])

        acc_interval = np.array(acc_interval)
        abs_acc_interval = abs(acc_interval)

        x_lower = start * dt  # the lower limit of x
        x_upper = end * dt  # the upper limit of x
        x_int = interval_time[np.where((x_lower <= interval_time) * (interval_time <= x_upper))]
        y_int = np.abs(np.array(abs_acc_interval)[np.where((x_lower <= interval_time) * (interval_time <= x_upper))])
        int_acc = trapz(y_int, x_int)
        # print (x_lower, x_upper)

        # calculation of pga (g)
        pga = (max(abs_acc_interval))
        if pga > pga_max:
            pga_max = pga

        if (pga - 0.025) < 0:
            h = 0
        elif (pga - 0.025) >= 0:
            h = 1
        else:
            raise ValueError("cannot evaluate pga: {0}".format(pga))

        CAVdp = CAVdp + (h * int_acc)
        start = end

    return CAVdp


def calculate_cav_dp_time_series(acc, dt):
    asig = eqsig.AccSignal(acc, dt)
    return calc_cav_dp_series(asig)


def calculate_cav_dp(acc, dt):
    asig = eqsig.AccSignal(acc, dt)
    return calc_cav_dp_series(asig)[-1]


def calc_cav_dp_series(asig):
    return eqsig.measures.calc_cav_dp(asig)


def calculate_cav_dp_series(asig):
    deprecation("calculate_cav_dp_series is deprecated - use calc_cav_dp_series")
    return calc_cav_dp_series(asig)
