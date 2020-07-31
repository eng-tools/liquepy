import numpy as np


def est_strain_inc_per_cycle_tasiopoulou_et_al_2020(relative_density, shear_stress, b=28., n_std=0):
    """
    Estimation of the increase in post-liq shear strain per cycle of equal amplitude load
    From Tasiopoulou (2020).

    Parameters
    ----------
    relative_density: float or array_like
        Soil relative density as a fraction (not percentage)
    shear_stress: float or array_like
        Stress amplitude (not normalised), in Pa
    b: float or array_like
        Slope intercept
    n_std: float or array_like
        Multiplier of the standard deviation

    Returns
    -------
    float or array_like
        Increment of increase in shear strain as a fraction (not percentage)
    """
    a = -0.1
    std_dev = 0.5

    return np.exp(a * relative_density * 100 + np.log(b) + n_std * std_dev) / 100 * shear_stress / 1e3

