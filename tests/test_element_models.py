from liquepy.element import models
import numpy as np


def test_av_stress():
    tau = np.array([0, -2, -4, -6, -8, -6, -4, -2, 0, 2, 4, 6, 8, 6, 4, 2, 0])
    av_tau = (tau[1:] + tau[:-1]) / 2
    et = models.ShearTest(tau, np.ones_like(tau), 1)
    diff = sum(abs(av_tau - et.av_stress[1:]))
    assert np.isclose(diff, 0)


if __name__ == '__main__':
    test_av_stress()
