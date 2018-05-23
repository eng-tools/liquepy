import numpy as np


def load_cpt_data(fname):
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
