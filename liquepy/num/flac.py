from sfsimodels import models
import sfsimodels as sm
import numpy as np
from collections import OrderedDict
from sfsimodels.functions import clean_float
from sfsimodels.build_model_descriptions import build_parameter_descriptions
from liquepy.element.models import ShearTest
from liquepy.exceptions import deprecation
from liquepy.num.models import PM4Sand as PM4SandBase


class FlacSoil(sm.Soil):
    _tension = 0.0  # default to zero
    hyst_str = ""

    def __str__(self):
        return "Base Soil model, id=%i, phi=%.1f" % (self.id, self.phi)

    def __init__(self, wmd=None, pw=None, liq_mass_density=None, liq_sg=1.0, g=9.8):
        # _gravity = g  # m/s2
        # if liq_mass_density:
        #     _liq_mass_density = liq_mass_density  # kg/m3
        # elif pw is not None and _gravity is not None:
        #     if pw == 9800 and g == 9.8:
        #         _liq_mass_density = 1.0e3
        #     else:
        #         _liq_mass_density = pw / _gravity
        # else:
        #     _liq_mass_density = None
        # run parent class initialiser function
        super(FlacSoil, self).__init__(wmd=wmd, pw=pw, liq_mass_density=liq_mass_density, liq_sg=liq_sg, g=g)
        self._extra_class_inputs = ["tension"]
        self.inputs = self.inputs + self._extra_class_inputs
        self.app2mod = OrderedDict([
            ("bulk", "bulk_mod"),
            ("shear", "g_mod"),
            ("friction", "phi"),
            ("cohesion", "cohesion"),
            ("tension", "tension"),
            ("density", "density"),
            ("dilation", "dilation_angle"),
            ("por", "porosity"),
            ("perm", "flac_permeability")
        ])
        self.required_parameters = []
        # for item in self.flac_parameters:
        #     param = self.flac_parameters[item]
        self.required_parameters = list(self.app2mod)
        self.optional_parameters = []
        if not hasattr(self, "definitions"):
            self.definitions = OrderedDict()
        self.definitions["density"] = ["Soil mass density", "kg"]
        self.definitions["tension"] = ["Soil strength in tension", "Pa"]
        self.definitions["flac_permeability"] = ["Permeability of soil", "m^2/Pa.sec"]
        self.prop_dict = OrderedDict()

    @property
    def all_flac_parameters(self):
        return self.required_parameters + self.optional_parameters
    
    def sm_to_fis_name(self, sm_name):
        for item in self.app2mod:
            if self.app2mod[item] == sm_name:
                return item
        return None

    def to_fis_mohr_coulomb(self, group_name=None, as_values=False):
        mc_params = OrderedDict([
            ("bulk", "bulk_mod"),
            ("shear", "g_mod"),
            ("friction", "phi"),
            ("cohesion", "cohesion"),
            ("tension", "tension"),
            ("density", "density"),
            ("dilation", "dilation_angle"),
            ("por", "porosity"),
            ("perm", "flac_permeability")
        ])
        if group_name is None:
            group_name = "'{0}'".format(self.name)
        para = ["model mohr notnull group %s" % group_name]
        para.append(write_parameters_to_fis_models(self, mc_params, ncols=1, not_null=True, group_name=group_name,
                                                   as_values=as_values))
        return "\n".join(para)

    def set_prop_dict(self):
        plist = []
        for item in self.app2mod:
            plist.append(self.app2mod[item])
        self.prop_dict = build_parameter_descriptions(self, user_p=self.definitions, output="dict", plist=plist)

    def find_units(self, parameter):
        if parameter in self.prop_dict:
            return self.prop_dict[parameter]["units"]
        else:
            return None

    @property
    def density(self):
        try:
            return self.unit_dry_weight / 9.81
        except TypeError:
            return None

    @property
    def tension(self):
        return self._tension

    @tension.setter
    def tension(self, value):
        value = clean_float(value)
        self._tension = value

    @property
    def flac_permeability(self):
        try:
            return self.permeability / self._uww
        except TypeError:
            return None


class PM4Sand(FlacSoil, PM4SandBase):

    type = "pm4sand"

    def __init__(self, wmd=None, pw=None, liq_mass_density=None, liq_sg=1.0, g=9.8, p_atm=101000.0, **kwargs):
        # Note: pw has deprecated
        _gravity = g  # m/s2
        if liq_mass_density:
            _liq_mass_density = liq_mass_density  # kg/m3
        elif pw is not None and _gravity is not None:
            if pw == 9800 and g == 9.8:
                _liq_mass_density = 1.0e3
            else:
                _liq_mass_density = pw / _gravity
        else:
            _liq_mass_density = None

        FlacSoil.__init__(self, wmd=wmd, pw=pw, liq_mass_density=liq_mass_density, liq_sg=liq_sg, g=g)
        PM4SandBase.__init__(self, wmd=wmd, pw=pw, liq_mass_density=liq_mass_density, liq_sg=liq_sg, g=g, p_atm=p_atm, **kwargs)
        self._extra_class_inputs = []

        additional_dict = OrderedDict([
            ("D_r", "relative_density"),
            ("h_po", "h_po"),
            ("G_o", "g0_mod"),
            ("density", "density"),
            ("porosity", "porosity"),
            ("h_o", "h_o"),
            ("e_min", "e_min"),
            ("e_max", "e_max"),
            ("n_b", "n_b"),
            ("n_d", "n_d"),
            ("c_z", "c_z"),
            ("c_e", "c_e"),
            ("n_d", "n_d"),
            ("k11", "flac_permeability"),
            ("k22", "flac_permeability"),
            ("P_atm", "p_atm"),
            ("phi_cv", "phi_cv"),
            ("pois", "poissons_ratio"),
            ("A_do", "a_do"),
            ("G_degr", "g_degr"),
            ("Ckaf", "c_kaf"),
            ("Q_bolt", "q_bolt"),
            ("R_bolt", "r_bolt"),
            ("MC_ratio", "mc_ratio"),
            ("MC_c", "mc_c"),
            ("z_max", "z_max"),
            ("c_dr", "c_dr"),
            ("m_par", "m_par"),
            ("f_sed", "f_sed"),
            ("p_sed", "p_sed")
        ])
        self.app2mod.update(additional_dict)
        self.pm4sand_parameters = self.app2mod  # deprecated
        self.required_parameters = ['h_po', 'D_r', 'G_o', 'P_atm', 'density']
        self.optional_parameters = [
            "k11",
            "k22",
            "pois",
            "h_o",
            "n_b",
            "n_d",
            "A_do",
            "z_max",
            "c_z",
            "c_e",
            "phi_cv",
            "G_degr",
            "c_dr",
            "Ckaf",
            "Q_bolt",
            "R_bolt",
            "m_par",
            "f_sed",
            "p_sed",
            "MC_ratio",
            "MC_c"
        ]

    def __repr__(self):
        return "PM4SandFLAC Soil model, id=%i, phi=%.1f, Dr=%.2f" % (self.id, self.phi, self.relative_density)

    def __str__(self):
        return "PM4SandFLAC Soil model, id=%i, phi=%.1f, Dr=%.2f" % (self.id, self.phi, self.relative_density)

    def to_fis(self, group_name=None, as_values=False):
        params = self.all_flac_parameters
        if group_name is None:
            group_name = "'{0}'".format(self.name)
        para = [f"model pm4sand notnull group {group_name}"]
        para.append(write_parameters_to_fis_models(self, params, ncols=1, not_null=True, group_name=group_name,
                                                   as_values=as_values))
        return para


def load_element_test(ffp, esig_v0, hydrostatic=0):
    ele_data = np.loadtxt(ffp, delimiter="  ", skiprows=1, usecols=(0, 1, 2, 4))
    n_count = ele_data[:, 0]
    csr_vals = ele_data[:, 1]
    tau = csr_vals * esig_v0
    strs = ele_data[:, 2] / 100
    ru_flac = ele_data[:, 3]
    stest = ShearTest(tau, strs, esig_v0=esig_v0, n_cycles=n_count)
    stest.set_pp_via_ru(ru_flac, hydrostatic=hydrostatic)
    return stest


def load_file_and_dt(fname):
    num_data_k = np.loadtxt(fname, skiprows=4)
    time = num_data_k[:, 0]  # This get the first column
    dt = time[1] - time[0]
    values = num_data_k[:, 1]
    return values, dt


def load_file_and_time(fname):
    num_data_k = np.loadtxt(fname, skiprows=4)
    time = num_data_k[:, 0]  # This get the first column
    values = num_data_k[:, 1]
    return values, time


def load_input_motion_and_dt(ffp):
    """
    Loads acceleration values and time step that were saved in FLAC input format.

    Parameters
    ----------
    ffp: str
        Full file path to output file

    Returns
    -------
    values: array_like
        An array of values
    dt: float
        Time step

    """
    data = np.genfromtxt(ffp, skip_header=1, delimiter=",", names=True, usecols=0)
    dt = data.dtype.names[0].split("_")[-1]
    dt = "." + dt[1:]
    dt = float(dt)
    acc = data.astype(float)
    return acc, dt


def save_input_motion_and_dt(ffp, values, dt, label="unlabelled"):
    """
    Exports acceleration values to the FLAC input format.

    Parameters
    ----------
    ffp: str
        Full file path to output file
    values: array_like
        An array of values
    dt: float
        Time step
    label: str
        A label of the data

    Returns
    -------

    """
    para = [label, "%i %.4f" % (len(values), dt)]
    for i in range(len(values)):
        para.append("%.6f" % values[i])
    ofile = open(ffp, "w")
    ofile.write("\n".join(para))
    ofile.close()


def save_input_motion(ffp, name, values, dt):
    """
    Exports acceleration values to the FLAC input format.

    :param ffp: str, full file path to output file
    :param name: str, name of records
    :param values: array, acceleration values
    :param dt: float, time step
    :return: None
    """
    deprecation("liquepy.num.flac.save_input_motion is deprecated, use liquepy.num.flac.save_input_motion_and_dt")
    para = [name, "%i %.4f" % (len(values), dt)]
    for i in range(len(values)):
        para.append("%.6f" % values[i])
    ofile = open(ffp, "w")
    ofile.write("\n".join(para))
    ofile.close()


def calc_hp0_from_crr_n15_and_relative_density_millen_et_al_2019(crr_n15, d_r):
    return crr_n15 * (2.05 - 2.4 * d_r) / (1. - crr_n15 * (12.0 - (12.5 * d_r)))


def calc_pm4sand_h_po_from_crr_n15_and_relative_density_millen_et_al_2019(crr_n15, d_r):
    return crr_n15 * (2.05 - 2.4 * d_r) / (1. - crr_n15 * (12.0 - (12.5 * d_r)))


def write_parameters_to_fis_models(obj, parameters, ncols=3, not_null=False, group_name=None, as_values=False):
    if group_name is None:
        group_name = obj.name
        
    count = 0
    para = []
    pline = ["prop"]
    for flac_name in parameters:
        if flac_name in obj.app2mod:
            ecp_name = obj.app2mod[flac_name]
        else:
            ecp_name = flac_name
        value = getattr(obj, ecp_name)
        if value is None:
            continue
        if as_values:
            value = f'{value:.6g}'
        else:
            value = f'{obj.name}_{ecp_name}'
        pline.append(f"{flac_name}={value}")
        count += 1
        if count == ncols:
            if not_null:
                pline.append(f"notnull group {group_name}")
            else:
                pline.append(f"group {group_name}")
            para.append(" ".join(pline))
            pline = ["prop"]
            count = 0
    if count:
        para.append(" ".join(pline))

    return "\n".join(para)


def write_obj_to_fis_str(obj, parameters, required=''):
    para = []
    added = []
    for item in parameters:
        if hasattr(obj, "find_units"):
            units = obj.find_units(item)
            if units is None or units == "":
                units_str = ""
            else:
                units_str = "  ; %s" % obj.find_units(item)
        else:
            units_str = ""
        if item in obj.app2mod:
            ecp_name = obj.app2mod[item]
        else:
            ecp_name = item
        val = getattr(obj, ecp_name)
        fis_name = f'fis_name = {obj.name}_{ecp_name}'
        if fis_name not in added:
            added.append(fis_name)
        else:
            continue
        if val is None:
            if item in required:
                raise sm.ModelError(f"{fis_name} is None")
        elif isinstance(val, str):
            para.append(f"  ={val}{units_str}")
        else:
            para.append(f"  {obj.name}_{ecp_name}={val:.5g}{units_str}")

    return para


def write_soil_profile_obj_to_fis_str(soil_profile, auto_name=True):
    """
    This file defines the contents of the soil profile object as parameters in fis language

    :param soil_profile:
    :return:
    """

    num_layers = soil_profile.n_layers

    para = ["", "def soil_profile_parameters"]

    for i in range(1, num_layers + 1):
        sl = soil_profile.layer(i)
        if auto_name:
            sl.name = "layer_%i" % i
        sl.set_prop_dict()
        para.append(" ; LAYER %i" % i)
        para.append(f'; name: {sl.name}')
        para.append(" layer_%i_thickness=%.4g  ; m" % (i, soil_profile.layer_height(i)))
        para += write_obj_to_fis_str(sl, sl.all_flac_parameters, sl.required_parameters)
        para.append("  ;")

    para.append(" ; GROUND WATER TABLE (negative)")
    para.append(" gwl= -%.4g  ; m" % soil_profile.gwl)
    para.append(" gravity=%.4g  ; m/s2" % 9.8)
    para.append(" density_water=%.4g  ; kg/m3" % (soil_profile.layer(1)._uww / 9.8))
    para.append("end")
    para.append("soil_profile_parameters")
    return para


def calc_rayleigh_damping_params(f1, f2, d):
    w1 = 2 * np.pi * f1
    w2 = 2 * np.pi * f2

    beta = 2 * (((w1 * d) - (w2 * d)) / ((w1 ** 2) - (w2 ** 2)))
    alpha = beta * w1 * w2

    wmin = np.sqrt(alpha / beta)

    emin = np.sqrt(alpha * beta)
    fmin = wmin / (2 * np.pi)

    return emin, fmin


def check_max_char_limit(fnames, run_loc):
    """
    Fis files have a maximum character limit of 200. This function checks that the fis files do not exceed it

    :param fnames:
    :param run_loc:
    :return:
    """
    MAX_CHAR_LIMIT = 200  # according to pg 2 of command reference manual this should be 2000
    fname_char_limit = 30
    for fname in fnames:
        if len(fname) > fname_char_limit:
            raise ValueError('File name: %s - exceeds limit of %i' % (fname, fname_char_limit))
        print(fname)
        a = open(run_loc + fname)
        lines = a.read().splitlines()
        for i, line in enumerate(lines):
            if len(line) > MAX_CHAR_LIMIT:
                print(line)
                raise ValueError('%s - Line %i has %i chars and exceeds limit of %i' % (fname, i + 1, len(line),
                                                                                        MAX_CHAR_LIMIT))


def pip_freeze_to_run_loc(run_loc):
    """
    Saves current list of python packages to the run location

    :param run_loc:
    :return:
    """
    import pip
    import datetime
    from pip._internal.operations import freeze
    reqs = freeze.freeze()
    dt = datetime.datetime.today().strftime("%Y%m%d")
    ofile = open(run_loc + f'requirements_{dt}.txt', 'w')
    ofile.write('\n'.join(reqs))
    ofile.close()


# see paper-lateral-spread/nonliq-soil/figure_calibrate_flac_hysteretic.py
def calc_flac_default_l1_from_ip_millen_2020(i_p):
    return 1.7 * i_p ** 0.5 - 3.6


def calc_flac_default_l2_from_ip_millen_2020(i_p):
    return -1.2 * i_p ** -0.25 + 2.25
