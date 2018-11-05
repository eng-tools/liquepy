import numpy as np
from liquepy.element.models import ShearTest

def load_flac_file_and_dt(fname):
    num_data_k = np.loadtxt(fname, skiprows=4)
    time = num_data_k[:, 0]  # This get the first column
    dt = time[1] - time[0]
    values = num_data_k[:, 1]
    return values, dt


def load_flac_file_and_time(fname):
    num_data_k = np.loadtxt(fname, skiprows=4)
    time = num_data_k[:, 0]  # This get the first column
    values = num_data_k[:, 1]
    return values, time


def load_flac_input_motion(ffp):
    data = np.genfromtxt(ffp, skip_header=1, delimiter=",", names=True, usecols=0)
    # print(ffp)
    # print(data.dtype.names)
    dt = data.dtype.names[0].split("_")[-1]
    dt = "." + dt[1:]
    dt = float(dt)
    acc = data.astype(np.float)
    return acc, dt


def load_flac_element_test(ffp, esig_v0, hydrostatic=0):
    ele_data = np.loadtxt(ffp, delimiter="  ", skiprows=1, usecols=(0, 1, 2, 4))
    n_count = ele_data[:, 0]
    csr_vals = ele_data[:, 1]
    tau = csr_vals * esig_v0
    strs = ele_data[:, 2] / 100
    ru_flac = ele_data[:, 3]
    stest = ShearTest(strs, tau, esig_v0=esig_v0, n_cycles=n_count)
    stest.set_pp_via_ru(ru_flac, hydrostatic=hydrostatic)
    stest.set_i_liq(vert_eff_stress_limit=5000)
    return stest