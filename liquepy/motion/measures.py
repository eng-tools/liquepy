import eqsig
import numpy as np
from numpy import trapz
from liquepy.exceptions import deprecation


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
