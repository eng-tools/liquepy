import numpy as np
from liquepy.exceptions import deprecation


def load_cpt_data(fname):
    deprecation('Deprecated (load_cpt_data), should use load_cpt_from_file')

    # import data from csv file
    data = np.loadtxt(fname, skiprows=24, delimiter=";")
    depth = data[:, 0]
    q_c = data[:, 1] * 1e3  # should be in kPa
    f_s = data[:, 2]
    u_2 = data[:, 3]
    gwl = None
    infile = open(fname)
    lines = infile.readlines()
    for line in lines:
        if "Assumed GWL:" in line:
            gwl = float(line.split(";")[1])

    return depth, q_c, f_s, u_2, gwl


def load_cpt_from_file(fname):
    # import data from csv file
    data = np.loadtxt(fname, skiprows=24, delimiter=";")
    depth = data[:, 0]
    q_c = data[:, 1] * 1e3  # should be in kPa
    f_s = data[:, 2]
    u_2 = data[:, 3]
    gwl = None
    a_ratio = None
    infile = open(fname)
    lines = infile.readlines()
    for line in lines:
        if "Assumed GWL:" in line:
            gwl = float(line.split(";")[1])
        if "aratio" in line:
            try:
                a_ratio = float(line.split(";")[1])
            except ValueError:
                pass
    return CPT(depth, q_c, f_s, u_2, gwl, a_ratio)


class CPT(object):
    def __init__(self, depth, q_c, f_s, u_2, gwl, a_ratio=None):
        """
        A cone penetration resistance test

        :param depth: array
        :param q_c: array, kPa,
        :param f_s: array, kPa,
        :param u_2: array, kPa,
        :param gwl: float, m, ground water level
        :param a_ratio: float, -, area ratio
        """
        self.depth = depth
        self.q_c = q_c
        self.f_s = f_s
        self.u_2 = u_2
        self.gwl = gwl
        self.a_ratio = a_ratio

