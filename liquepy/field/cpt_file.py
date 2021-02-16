import numpy as np
from liquepy.exceptions import deprecation
import ntpath


# def load_cpt_data(fname):
#     deprecation('Deprecated (load_cpt_data), should use load_cpt_from_file')
#
#     # import data from csv file
#     data = np.loadtxt(fname, skiprows=24, delimiter=";")
#     depth = data[:, 0]
#     q_c = data[:, 1] * 1e3  # should be in kPa
#     f_s = data[:, 2]
#     u_2 = data[:, 3]
#     gwl = None
#     infile = open(fname)
#     lines = infile.readlines()
#     for line in lines:
#         if "Assumed GWL:" in line:
#             gwl = float(line.split(";")[1])
#
#     return depth, q_c, f_s, u_2, gwl

def load_mpa_cpt_file(ffp, delimiter=",", a_ratio_override=None):
    # import data from csv file
    folder_path, file_name = ntpath.split(ffp)
    ncols = 4
    try:
        data = np.loadtxt(ffp, skiprows=24, delimiter=delimiter, usecols=(0, 1, 2, 3))
    except:
        ncols = 3
        data = np.loadtxt(ffp, skiprows=24, delimiter=delimiter, usecols=(0, 1, 2))
    depth = data[:, 0]
    q_c = data[:, 1] * 1e3  # convert to kPa
    f_s = data[:, 2] * 1e3  # convert to kPa
    if ncols == 4:
        u_2 = data[:, 3] * 1e3  # convert to kPa
    else:
        u_2 = np.zeros_like(depth)
    gwl = None
    a_ratio = 1.0
    pre_drill = None
    infile = open(ffp)
    lines = infile.readlines()
    for line in lines:
        if "Assumed GWL:" in line:
            gwl = line.split(delimiter)[1]
            if gwl == '-':
                gwl = None
            else:
                gwl = float(line.split(delimiter)[1])
        if "aratio" in line:
            try:
                a_ratio = float(line.split(delimiter)[1])
            except ValueError:
                pass
        if "Pre-Drill:" in line:
            val = line.split(delimiter)[1]
            if val != '':
                pre_drill = float(val)
    if a_ratio_override:
        a_ratio = a_ratio_override
    if pre_drill is not None:
        if depth[0] < pre_drill:
            indy = np.argmin(abs(depth - pre_drill))
            depth = depth[indy:]
            q_c = q_c[indy:]
            f_s = f_s[indy:]
            u_2 = u_2[indy:]
    return CPT(depth, q_c, f_s, u_2, gwl, a_ratio, folder_path=folder_path, file_name=file_name, delimiter=delimiter)


def load_cpt_from_file(ffp, delimiter=";"):
    deprecation('Use load_mpa_cpt_file() where file is all in MPa')
    # import data from csv file
    folder_path, file_name = ntpath.split(ffp)
    ncols = 4
    try:
        data = np.loadtxt(ffp, skiprows=24, delimiter=delimiter, usecols=(0, 1, 2, 3))
    except:
        ncols = 3
        data = np.loadtxt(ffp, skiprows=24, delimiter=delimiter, usecols=(0, 1, 2))
    depth = data[:, 0]
    q_c = data[:, 1] * 1e3  # should be in kPa
    f_s = data[:, 2]
    if ncols == 4:
        u_2 = data[:, 3]
    else:
        u_2 = np.zeros_like(depth)
    gwl = None
    a_ratio = None
    infile = open(ffp)
    lines = infile.readlines()
    for line in lines:
        if "Assumed GWL:" in line:
            gwl = float(line.split(delimiter)[1])
        if "aratio" in line:
            try:
                a_ratio = float(line.split(delimiter)[1])
            except ValueError:
                pass
    return CPT(depth, q_c, f_s, u_2, gwl, a_ratio, folder_path=folder_path, file_name=file_name, delimiter=delimiter)


class CPT(object):
    def __init__(self, depth, q_c, f_s, u_2, gwl, a_ratio=None, folder_path="<path-not-set>", file_name="<name-not-set>",
                 delimiter=";"):
        """
        A cone penetration resistance test

        Parameters
        ----------
        depth: array_like
            depths from surface, properties are forward projecting (i.e. start at 0.0 for surface)
        q_c: array_like, [kPa]
        f_s: array_like, [kPa]
        u_2: array_like, [kPa]
        gwl: float, [m]
            ground water level
        a_ratio: float,
            Area ratio
        """
        self.depth = depth
        self.q_c = q_c
        self.f_s = f_s
        self.u_2 = u_2
        self.gwl = gwl
        self.a_ratio = a_ratio
        self.folder_path = folder_path
        self.file_name = file_name
        self.delimiter = delimiter


    @property
    def q_t(self):
        """
        Pore pressure corrected cone tip resistance

        """
        # qt the cone tip resistance corrected for unequal end area effects, eq 2.3
        return self.q_c + ((1 - self.a_ratio) * self.u_2)

    @q_t.setter
    def q_t(self, q_t):
        self.q_c = q_t - ((1 - self.a_ratio) * self.u_2)



def correct_qt_via_boulanger_and_dejong_2018(q_t, depth, d_c):

    z_dash = (depth[:, np.newaxis] - depth[np.newaxis, :]) / d_c
    c_1 = np.where(z_dash >= 0, 1, np.where(z_dash >= -4, 1 + z_dash / 8, 0.5))
    # pg 30 'The parameter C2 is equal to unity for points below the cone tip,
    # and less than unity (0.8 is used herein) for points above the cone tip.
    c_2 = np.where(z_dash >= 0, 1, 0.8)
    z50_dash_ref = 4.0  # pg 30
    m50 = 0.5  # pg 30
    mz = 3.0  # pg 30
    mq = 2  # pg 30
    qt_dash_o_qt = q_t[:, np.newaxis] / q_t[np.newaxis, :]
    qt_dash_o_qt = 100
    w_2 = np.sqrt(2. / (1 + (1 / qt_dash_o_qt) ** mq))
    z50_dash = 1 + 2 * (c_2 * z50_dash_ref - 1) * (1 - 1. / (1 + (qt_dash_o_qt) ** m50))
    w_1 = c_1 / (1 + (z_dash / z50_dash) ** mz)
    w_c = w_1 * w_2 / np.sum(w_1 * w_2, axis=0)
    # import matplotlib.pyplot as plt
    # bf, ax = plt.subplots(ncols=2, sharey='row')
    # ax[1].plot(w_c[400], z_dash[400])
    # ax[0].invert_yaxis()
    # plt.show()
    import matplotlib.pyplot as plt
    bf, ax = plt.subplots(ncols=2, sharey='row')
    ax[0].plot(q_t, depth)
    ax[1].plot(w_c, depth)
    ax[0].invert_yaxis()
    plt.show()
    return w_c


def view_correction():

    z = np.arange(0, 10, 0.01)
    zm = [0, 2.5, 2.9, 3.0, 3.3, 5.0, 5.2, 5.4, 10]
    qm = [3, 3.0, 3.2, 4.8, 5.0, 4.6, 3.2, 3.0, 3.0]
    q_m = np.interp(z, zm, qm)
    qt = np.ones_like(z) * 3
    qt[300:500] = 7
    import matplotlib.pyplot as plt
    bf, ax = plt.subplots(ncols=2, sharey='row')
    ax[0].plot(qt, z)
    ax[0].plot(q_m, z)
    ax[1].plot(q_m / qt, z)
    ax[0].invert_yaxis()
    plt.show()
    # correct_qt_via_boulanger_and_dejong_2018(q_m, z, 0.08)


if __name__ == '__main__':
    view_correction()

