import eqsig
import numpy as np
from numpy import trapz


def integral_of_velocity(acc, dt):
    delta_vel = acc * dt
    vel = np.cumsum(delta_vel)
    abs_vel = abs(vel)
    vel_int = np.cumsum(abs_vel * dt)
    return vel_int


def integral_of_acceleration(acc, dt):

    abs_acc = abs(acc)
    acc_int = np.cumsum(abs_acc * dt)
    return acc_int


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
            H=0
        if (pga - 0.025) >= 0:
            H=1
        #H = 1  # what is x??
        #CAVdp = CAVdp + (H * (pga - 0.025) * int_acc)
        CAVdp = CAVdp + (H * int_acc)
        start = end

    return CAVdp


def calculate_cav_dp_time_series(acc, dt):
    asig = eqsig.AccSignal(acc, dt)
    return calculate_cav_dp_series(asig)


def calculate_cav_dp(acc, dt):
    asig = eqsig.AccSignal(acc, dt)
    return calculate_cav_dp_series(asig)[-1]


def calculate_cav_dp_series(asig):
    return eqsig.measures.calc_cav_dp(asig)