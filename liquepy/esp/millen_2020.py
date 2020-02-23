import itertools
import os
import numpy as np
import liquepy as lq
from liquepy.trigger import boulanger_and_idriss_2014 as bi_functions
import sfsimodels as sm
import eqsig
from liquepy.field import correlations as lq_field_correlations
from collections import OrderedDict

# ESP_MODULE_DIR = os.path.dirname(os.path.abspath(__file__)) + "/"


def fit_3_layer_profile(depths, crr_n15s, fitting_values, crr_non_liq=0.6, max_depth=20):
    """
    Determines the best 3-layer profile fit for a given parameter

    Method uses brut force different CRR values. For each CRR is determines the change points and then evaluates which
    change points are best.

    Sets the caps the other layer crr values to the non liquefying value best to use 0.6 as this is approx the upper
    limit of liquefiable soil. Then the error for over-fitting (defining a non-liq layer as liquefiable) incurs an
    error of 0.6-crr_liq, while under-fitting layers less than or equal to crr_liq incurs the same error
    while layers in between crr_liq and the non-liq limit incur an error proportional to the difference from the defined
    and calculated crr.

    Parameters
    ----------
    depths: array,
        distance from surface
    crr_n15s: array,
        actual values to fit to (note CRR_max=4)
    fitting_values: array,
        possible values of layer 2
    crr_non_liq: float,
        value for layers 1 and 3
    :return:
    """

    if depths[-1] > max_depth:
        indy = np.where(depths > max_depth)[0][0]
    else:
        indy = len(depths)
    std_depths = depths[:indy]
    std_crr_n15s = crr_n15s[:indy]
    n_depths = len(std_depths)
    capped_values = np.clip(std_crr_n15s, None, crr_non_liq)
    diffs = capped_values[:, np.newaxis] - fitting_values[np.newaxis, :]
    cdiffs = np.cumsum(np.abs(diffs), axis=0)

    diffs = []
    h_crusts = []
    h_liqs = []
    for ii in range(len(fitting_values)):
        eline = cdiffs[:, 0] - cdiffs[:, ii]
        peak_ids = eqsig.get_peak_array_indices(eline)
        if eline[peak_ids[1]] > 0:
            peak_ids = peak_ids[:-1]
        else:
            peak_ids = peak_ids[1:-1]
        if len(peak_ids):
            loc_min_is = np.take(peak_ids, np.arange(0, len(peak_ids), 2))
            loc_max_is = np.take(peak_ids, np.arange(1, len(peak_ids), 2))
            loc_mins = np.take(eline, loc_min_is)
            loc_maxs = np.take(eline, loc_max_is)
            opts = loc_mins[np.newaxis, :] - loc_maxs[:, np.newaxis]
            opts *= np.tri(*opts.shape)  # remove cases where liq layer is higher than crust

            min_i, min_j = np.unravel_index(opts.argmin(), opts.shape)
            i_crust = loc_min_is[min_j]
            i_liq = loc_max_is[min_i]
            h_crusts.append(depths[i_crust + 1])
            h_liqs.append(depths[i_liq + 1] - depths[i_crust + 1])

            refined_errors = cdiffs[i_crust][0] + (cdiffs[i_liq][ii] - cdiffs[i_crust][ii]) + (cdiffs[-1][0] - cdiffs[i_liq][0])
            diffs.append(refined_errors)
        else:
            h_crusts.append(depths[-1])
            h_liqs.append(0)
            diffs.append(1e6)

    i_best = np.argmin(diffs)

    h_crust = h_crusts[i_best]
    h_liq = h_liqs[i_best]
    p_value = fitting_values[i_best]
    diff = diffs[i_best]
    normed_diff = diff / (n_depths * crr_non_liq)

    return h_crust, h_liq, p_value, normed_diff


def fit_n_layer_profile(depths, crr_n15s, n=3, crr_n15_opts=None, cap=0.6, max_depth=20):
    """
    Determines the best n-layer profile fit for a given parameter

    Method:
    1) Computes the cumulative absolute difference between the upper bound and the actual values
    2) Cycles through a range of guess values (between upper and lower bounds), computes the cumulative absolute
       difference between a guessed value and the actual values
    3) For each guess it determines the change points
    4) For each guess, Evaluates which change points pairs are best (either 1 or 2 depending on n).
    5) For each guess, for each change point pair, compute mean of values between pair (least error)
    6) For each guess, compute the total error as the sum of the error from all of the non-fitted and fitted layers
    7) Select guess which has the least error

    Parameters
    ----------
    depths: array
        Distance from surface
    crr_n15s: array
        Actual values to fit to
    n: int (3 or 5)
        Number of layers to fit
    crr_n15_opts: array_like
        Possibly options for CRR_n15
    max_depth: float
        Maximum distance to consider in fitting
    cap: float
        Clips all actual values to this maximum prior to fitting

    Returns
    -------
    array_like: depths to top of liq layers
    array_like: depths to top of non-liq layers
    array_like: fitted values in liq layers
    float: normalised error
    """
    if not (n == 3 or n == 5):
        raise ValueError("n must be either 3 or 5")
    n_liqs = int((n - 1) / 2)

    # prepare search array
    if crr_n15_opts is None:
        if n_liqs == 1:
            crr_n15_opts = np.array([0.6, 0.5, 0.2, 0.06])
            # q_c1n_cs = np.arange(0, 180., 5.)
            # crr_n15_opts = bi_functions.calc_crr_m7p5_from_qc1ncs(q_c1n_cs)[::-1]
        else:
            q_c1n_cs = np.arange(0, 180., 5.)
            crr_n15_opts = bi_functions.calc_crr_m7p5_from_qc1ncs(q_c1n_cs)[::-1]

    # standardise depth
    if depths[-1] > max_depth:
        indy = np.where(depths > max_depth)[0][0]
    else:
        indy = len(depths)
    std_depths = depths[:indy]
    std_crr_n15s = crr_n15s[:indy]
    n_depths = len(std_depths)

    # enforce cap
    std_crr_n15s = np.clip(std_crr_n15s, None, cap)

    # compute difference between options and actual values
    init_diffs = std_crr_n15s[:, np.newaxis] - crr_n15_opts[np.newaxis, :]
    init_cdiffs = np.cumsum(np.abs(init_diffs), axis=0)

    # prepare output arrays
    normed_diffs = []
    d_liqs = [[] for i in range(n_liqs)]
    d_nonliqs = [[] for i in range(n_liqs)]
    crrs = [[] for i in range(n_liqs)]

    # evaluate each option
    for ii in range(len(crr_n15_opts)):
        opts_lay = []
        eline = init_cdiffs[:, 0] - init_cdiffs[:, ii]
        peak_ids = eqsig.get_peak_array_indices(eline) + 1

        if len(peak_ids) > 4 or (len(peak_ids) > 2 and n_liqs == 1):
            # continue
            if eline[peak_ids[1]] > 0:
                peak_ids = peak_ids[:-1]
            else:
                peak_ids = peak_ids[1:-1]

            poss_i_tops = np.take(peak_ids, np.arange(0, len(peak_ids), 2))
            poss_i_bots = np.take(peak_ids, np.arange(1, len(peak_ids), 2))
            poss_i_tops = np.insert(poss_i_tops, len(poss_i_tops), len(std_crr_n15s) - 1)  # add end
            poss_i_bots = np.insert(poss_i_bots, len(poss_i_bots), len(std_crr_n15s) - 1)
            poss_err_tops = np.take(eline, poss_i_tops)
            poss_err_bots = np.take(eline, poss_i_bots)
            for ll in range(n_liqs):
                if len(poss_i_tops) < n_liqs or len(poss_i_bots) < n_liqs:
                    opts_lay.append([0])
                else:
                    # error occurred for each change point pair
                    opts_lay.append(poss_err_tops[np.newaxis, ll:] - poss_err_bots[ll:, np.newaxis])
                    # remove cases where i_bot is higher than i_top
                    opts_lay[ll] *= np.tri(*opts_lay[ll].shape)
            if n_liqs == 1:
                opts_l1f = opts_lay[0]
                j, i = np.unravel_index(opts_l1f.argmin(), opts_l1f.shape)
                i_tops = [int(poss_i_tops[i])]
                i_bots = [int(poss_i_bots[j])]

            elif n_liqs == 2:
                opts_l1f = opts_lay[0].flatten()
                top1_is = poss_i_tops[np.newaxis, :] * np.ones_like(opts_lay[0], dtype=int)
                bot1_is = poss_i_bots[:, np.newaxis] * np.ones_like(opts_lay[0], dtype=int)

                opts_l2f = opts_lay[1].flatten()
                top2_is = poss_i_tops[np.newaxis, 1:] * np.ones_like(opts_lay[1], dtype=int)
                bot2_is = poss_i_bots[1:, np.newaxis] * np.ones_like(opts_lay[1], dtype=int)

                opts = opts_l1f[np.newaxis, :] + opts_l2f[:, np.newaxis]
                top1_is = top1_is.flatten()[np.newaxis, :] * np.ones_like(opts, dtype=int)
                bot1_is = bot1_is.flatten()[np.newaxis, :] * np.ones_like(opts, dtype=int)

                top2_is = top2_is.flatten()[:, np.newaxis] * np.ones_like(opts, dtype=int)
                bot2_is = bot2_is.flatten()[:, np.newaxis] * np.ones_like(opts, dtype=int)

                # remove cases where i_crust_l2 < i_liq_l1
                opts = np.where(top2_is < bot1_is, 1e10, opts)
                # remove cases where i_crust_l2 > i_liq_l2
                opts = np.where(top2_is > bot2_is, 1e10, opts)
                opts = np.where(top2_is == len(std_depths) - 1, 1e10, opts)

                lay1_i, lay2_i = np.unravel_index(opts.argmin(), opts.shape)
                i_tops = [int(top1_is[lay1_i][lay2_i]), int(top2_is[lay1_i][lay2_i])]
                i_bots = [int(bot1_is[lay1_i][lay2_i]), int(bot2_is[lay1_i][lay2_i])]

            else:
                raise ValueError("n_liqs must be either 1 or 2")

            # total_err = init_cdiffs[i_tops[0]][0]
            crr_profile = np.ones_like(std_crr_n15s) * cap
            for ll in range(n_liqs):
                d_liqs[ll].append(depths[i_tops[ll]])
                d_nonliqs[ll].append(depths[i_bots[ll]])
                if i_tops[ll] == i_bots[ll]:
                    crrs[ll].append(cap)
                    # err = 0
                else:
                    # Median equal to min of absolute deviation
                    crr_median = np.median(std_crr_n15s[i_tops[ll]:i_bots[ll]])
                    crrs[ll].append(crr_median)
                    crr_profile[i_tops[ll]:i_bots[ll]] = crr_median
                # err = np.sum(np.abs(capped_values[i_tops[ll]:i_bots[ll]] - crrs[ll][ii]))
                # total_err += err
            total_err = np.sum(abs(std_crr_n15s - crr_profile))

            normed_diffs.append(total_err / (n_depths * cap))

        else:
            if len(peak_ids) <= 2:
                for ll in range(n_liqs):
                    d_liqs[ll].append(0)
                    d_nonliqs[ll].append(depths[-1])
                    crrs[ll].append(cap)
                normed_diffs.append(1e6)
            else:
                d_liqs[0].append(depths[peak_ids[1]])
                d_nonliqs[0].append(depths[peak_ids[2]])
                crr_median = np.median(std_crr_n15s[peak_ids[1]:peak_ids[2]])
                crrs[0].append(crr_median)
                d_liqs[1].append(depths[-1])
                d_nonliqs[1].append(depths[-1])
                crrs[1].append(cap)
                crr_profile = np.ones_like(std_crr_n15s) * cap
                crr_profile[peak_ids[1]:peak_ids[2]] = crr_median
                total_err = np.sum(abs(std_crr_n15s - crr_profile))
                normed_diffs.append(total_err / (n_depths * cap))

    normed_diffs = np.array(normed_diffs)
    d_liqs = np.array(d_liqs)
    d_nonliqs = np.array(d_nonliqs)
    crrs = np.array(crrs)

    i_best = np.argmin(normed_diffs)
    normed_diff = normed_diffs[i_best]

    # change to ECP logic of depths to top of layers
    d_liqs = d_liqs[:, i_best]
    d_nonliqs = d_nonliqs[:, i_best]

    return d_liqs, d_nonliqs, crrs[:, i_best], normed_diff


def create_profile_array(d_nonliqs, d_liqs, depth, layer_values, other_layers=2.):
    """
    Given equivalent soil profile values it builds an array based on the depths

    :param h_crust: float, height of crust
    :param h_liq: float, height of liquefied layer
    :param depth: array, distance from surface
    :param main_layer_value: float, value for liquefiable layer
    :param other_layers: float, value for non-liquefiable layers
    :return: array, values
    """
    profile = other_layers * np.ones_like(depth)
    inds = np.arange(len(depth))
    for ll in range(len(d_nonliqs)):
        crust_i = int(np.interp(d_liqs[ll], depth, inds))
        lower_i = int(np.interp(d_nonliqs[ll], depth, inds))
        profile[crust_i:lower_i] = layer_values[ll]
    return profile


class EquivalentProfiler(object):
    _e3_d_nonliqs = None
    _e3_d_liqs = None
    _e3_csr_n15 = None
    _e3_norm_diff = None
    _e3_q_c1n_cs = None
    _e5_d_nonliqs = None
    _e5_d_liqs = None
    _e5_csr_n15 = None
    _e5_norm_diff = None
    _e5_q_c1n_cs = None

    def __init__(self, bi2014):
        self._bi2014 = bi2014
        self.depth = bi2014.depth

    @property
    def bi2014(self):
        return self._bi2014

    def compute_e3profile(self, max_depth=20.):
        self._e3_d_liqs, self._e3_d_nonliqs, self._e3_csr_n15, self._e3_norm_diff = fit_n_layer_profile(self.bi2014.depth, self.bi2014.crr_m7p5, n=3, max_depth=max_depth)

    @property
    def e3_h_crust(self):
        if self._e3_d_nonliqs is None:
            self.compute_e3profile()
        return self._e3_d_liqs[0]

    @property
    def e3_h_liq(self):
        if self._e3_d_nonliqs is None:
            self.compute_e3profile()
        return self._e3_d_nonliqs[0] - self._e3_d_liqs[0]

    @property
    def e3_csr_n15(self):
        if self._e3_csr_n15 is None:
            self.compute_e3profile()
        return self._e3_csr_n15[0]

    @property
    def e3_qc1ncs(self):
        qc1ncs_vals = bi_functions.calculate_qc_1ncs_from_crr_7p5(self.e3_csr_n15)
        if self._e3_q_c1n_cs is None:
            self._e3_q_c1n_cs = create_profile_array(self._e3_d_nonliqs, self._e3_d_liqs, self.depth, qc1ncs_vals)
        return self._e3_q_c1n_cs

    @property
    def e3_crr_m7p5(self):
        return create_profile_array(self._e3_d_nonliqs, self._e3_d_liqs, self.depth, self.e3_csr_n15s)

    @property
    def e3_norm_diff(self):
        if self._e3_norm_diff is None:
            self.compute_e3profile()
        return self.e3_norm_diff

    def compute_e5profile(self):
        self._e5_d_nonliqs, self._e5_d_liqs, self._e5_csr_n15, self._e5_norm_diff = fit_n_layer_profile(
            self.bi2014.depth, self.bi2014.crr_m7p5, n=5)

    @property
    def e5_h_nonliqs(self):
        if self._e5_d_nonliqs is None:
            self.compute_e5profile()
        return self._e5_d_nonliqs[0]

    @property
    def e5_h_liqs(self):
        if self._e5_d_nonliqs is None:
            self.compute_e5profile()
        return self._e5_d_liqs[0] - self._e5_d_nonliqs[0]

    @property
    def e5_csr_n15s(self):
        if self._e5_csr_n15 is None:
            self.compute_e5profile()
        return self._e5_csr_n15[0]

    @property
    def e5_qc1ncs(self):
        qc1ncs_vals = bi_functions.calculate_qc_1ncs_from_crr_7p5(self.e5_csr_n15s)
        if self._e5_q_c1n_cs is None:
            self._e5_q_c1n_cs = create_profile_array(self._e5_d_nonliqs, self._e5_d_liqs, self.depth, qc1ncs_vals)
        return self._e5_q_c1n_cs

    @property
    def e5_crr_m7p5(self):
        return create_profile_array(self._e5_d_nonliqs, self._e5_d_liqs, self.depth, self.e5_csr_n15s)

    @property
    def e5_norm_diff(self):
        if self._e5_norm_diff is None:
            self.compute_e5profile()
        return self.e5_norm_diff

    @property
    def total_depth(self):
        return self.depth[-1]


def compute_e5profile(bi2014):
    d_liqs, d_nonliqs, csr_n15s, norm_diff = fit_n_layer_profile(bi2014.depth, bi2014.crr_m7p5, n=5)
    return Equiv5LayerProfile(d_nonliqs, d_liqs, csr_n15s, bi2014.gwl, norm_diff, bi2014.depth, surrogate=True)


def compute_e3profile(bi2014):
    d_liqs, d_nonliqs, csr_n15, norm_diff = fit_n_layer_profile(bi2014.depth, bi2014.crr_m7p5, n=3)
    return Equiv3LayerProfile(d_liqs[0], d_nonliqs[0], csr_n15[0], bi2014.gwl, norm_diff, bi2014.depth, surrogate=True)


def crr_error_variance(eprofile, crr_values):
    return 1 / len(eprofile.depth) * np.sum((crr_values - eprofile.crr_m7p5) ** 2)


def calc_bic(eprofile, bi2014):
    """
    Calculates the Bayesian Information Criterion

    See: https://www.immagic.com/eLibrary/ARCHIVES/GENERAL/WIKIPEDI/W120607B.pdf

    Parameters
    ----------
    eprofile
    bi2014

    Returns
    -------

    """
    err_var = crr_error_variance(eprofile, bi2014.crr_m7p5)
    if eprofile.n_layers == 3:
        k = 3
    else:
        k = 6
    n = len(eprofile.depth)
    bic = n * np.log(err_var) + k * np.log(n)
    return bic


def compute_new_lsn_from_equivalent_profile(eq_profiler, new_pga, new_m_w, n=3):
    """
    Computes the Liquefaction severity number (LSN) for an equivalent soil profile

    :param eq_profile: EquivalentSoilProfile object
    :param new_pga: float, peak ground acceleration to consider LSN at
    :param new_m_w: float, earthquake magnitude to consider LSN at
    :return: float, LSN
    """
    depth = eq_profiler.depth
    gwl = eq_profiler.bi2014.gwl

    if n == 3:
        eq_qc1ncs = eq_profiler.e3_q_c1n_cs
        crr_7p5 = eq_profiler.e3_crr_7p5
    elif n == 5:
        eq_qc1ncs = eq_profiler.e5_q_c1n_cs
        crr_7p5 = eq_profiler.e5_crr_7p5
    else:
        raise ValueError

    # compute msf
    sigma_v = eq_profiler.bi2014.sigma_v
    sigma_v_eff = eq_profiler.bi2014.sigma_veff
    k_sigma = bi_functions.calc_k_sigma(sigma_v_eff, eq_qc1ncs)
    msf = bi_functions.calc_msf(new_m_w, eq_qc1ncs)
    rd = bi_functions.calc_rd(depth, new_m_w)

    # compute CSR
    csr = bi_functions.calc_csr(sigma_v_eff, sigma_v, new_pga, rd, gwl, depth)

    # using the CRR and q_c1n_cs values compute the new FoS
    eq_factor_of_safety = crr_7p5 / csr * k_sigma * msf

    # using the new FoS and q_c1n_cs values compute the volumetric strain
    eq_epsilon = lq.trigger.calculate_volumetric_strain(eq_factor_of_safety, eq_qc1ncs)
    lsn = lq.trigger.calculate_lsn(eq_epsilon, depth)
    # using the volumetric strain and new FoS compute new lsn

    return lsn


class EquivalentProfile(sm.PhysicalObject):  # TODO: create SurrogateObject, always attached
    q_c = None
    f_s = None
    u_2 = None
    pga = None
    magnitude = None
    ar = None

    cfc = 0  # parameter of fines content, eq 2.29
    q_t = None

    q_c1n_cs = None
    fines_content = None
    i_c = None

    k_sigma = None
    msf = None
    csr = None
    crr = None  # Array with depth
    crr_m7p5 = None  # Array with depth
    factor_of_safety = None

    crr_m7p5_layer = None  # This is not an array

    def __init__(self, depth, unit_wt, gwl, h_crust, h_liq):
        # super(EquivalentProfile, self).__init__()
        self.depth = depth
        self.h_crust = h_crust
        self.h_liq = h_liq
        self.unit_wt = unit_wt
        self.gwl = gwl
        self._extra_class_inputs = ["h_crust", "h_liq", "gwl", "csr_n15", "total_depth"]
        self.inputs = [] + self._extra_class_inputs

        self.sigma_v = bi_functions.calc_sigma_v(self.depth, self.unit_wt)
        self.pore_pressure = bi_functions.calc_pore_pressure(self.depth, self.gwl)
        self.sigma_veff = bi_functions.calc_sigma_veff(self.sigma_v, self.pore_pressure)
        depth_increment = depth[2] - depth[1]
        self.crust_i = int(h_crust / depth_increment)
        self.lower_i = min(int((h_crust + h_liq) / depth_increment), len(depth))

    @property
    def gammas(self):
        return self.unit_wt

    @property
    def csr_n15(self):
        return self.crr_m7p5_layer

    @property
    def total_depth(self):
        return self.depth[-1]


def get_esp_names():
    return ['WLS'
            'WLM',
            'WLD',
            'WMS',
            'WMM',
            'WMD',
            'WTS',
            'WTM',
            'WTD',
            'MLS',
            'MLM',
            'MLD',
            'MMS',
            'MMM',
            'MMD',
            'MTS',
            'MTM',
            'MTD',
            'SLX',
            'SMX',
            'STX',
            'RXX']


def get_esp_criteria():
    """
    Builds a dictionary of esp criteria

    :return:
    """
    cd = {'crr_n15':
           {'W': [0.0, 0.15],
            'M': [0.15, 0.25],
            'S': [0.25, 0.45],
            'R': [0.45, 0.60]},
          'h_liq':
                {'T': [0.5, 3.0],
                'M': [3.0, 7.0],
            'L': [7.0, 10.0]},
          'h_crust':
              {'S': [0.0, 2.0],
            'M': [3.0, 7.0],
            'D': [7.0, 10.0]
            }
          }

    return cd


def classify_esp(esp_values: dict) -> str:
    """
    Determine the class of the equivalent soil profile

    :param esp_values: dictionary storing the values for the lq_esp parameters
    :return: str, class of ESP
    """

    esp_class = ""
    clist = ["crr_n15", "h_liq", "h_crust"]
    cd = get_esp_criteria()
    for i, crit in enumerate(clist):
        for key, val in cd[crit]:
            if val[0] <= esp_values[crit] < val[1]:
                esp_class += key
                break
        if len(esp_class) == i:
            esp_class += cd[crit].keys()

    # Adjust for aggregated classes
    if esp_class[0] == "S":
        esp_class = esp_class[:-1] + "X"
    if esp_class[0] == "R":
        esp_class = esp_class[0] + "XX"

    return esp_class


def set_strength_props(sl, vert_eff_stress, i_c, big_q, n_kt=14):
    if i_c > 2.6:
        sl.cohesion = lq.field.correlations.est_undrained_strength_ratio_robertson_2009(big_q, n_kt=n_kt) * vert_eff_stress
        sl.phi = 0.0
    else:
        sl.phi = np.arctan(lq.field.correlations.est_undrained_strength_ratio_robertson_2009(big_q, n_kt=n_kt))
        sl.cohesion = 0.0


def set_soil_props(sl, bi2014, i_top, i_bot, gwl):
    sl.specific_gravity = bi2014.s_g
    val = np.mean(bi2014.unit_dry_wt[i_top:i_bot] * 1e3)
    sl.unit_dry_weight = val
    base_depth = bi2014.depth[i_bot]
    if base_depth > gwl:
        unit_wt = sl.unit_dry_weight
    else:
        unit_wt = sl.unit_sat_weight
    if np.median(bi2014.i_c) > 2.6:  # clay-like
        undrained_strength_ratio_i = lq.field.correlations.est_undrained_strength_ratio_robertson_2009(bi2014.big_q[i_top:i_bot])
        undrained_strength_ratio = np.mean(undrained_strength_ratio_i)
        rel_density_i = lq_field_correlations.calc_relative_density_salgado_et_al_1997_cpt_values(bi2014.q_c1n[i_top:i_bot])
        sl.relative_density = np.mean(rel_density_i)
        av_esig_v = np.mean(bi2014.sigma_veff)
        sl.cohesion = undrained_strength_ratio * av_esig_v
        sl.phi = 0.0
    else:
        sl.phi = 30.0  # TODO: better definition for this
        sl.cohesion = 0.0
    sl.g_mod = np.mean(lq.field.correlations.est_g_mod_robertson_2009(bi2014.i_c, bi2014.big_q, unit_wt))


def build_soil_profile_from_bi2014(esp, bi2014):
    """
    Defined soil profile properties based on an equivalent soil profile define

    Parameters
    ----------
    bi2014: liquepy.trigger.BoulangerIdriss2014CPT
    esp: equivalent_soil object
        Equivalent soil profile of the triggering analysis

    Returns
    -------
    sm.SoilProfile
    """
    assert isinstance(bi2014, lq.trigger.BoulangerIdriss2014CPT)
    sp = sm.SoilProfile()
    sp.gwl = bi2014.gwl
    sp.height = bi2014.depth[-1]
    s_depths = [0, esp.h_crust]
    csr_n15 = [0]
    if esp.n_layers == 3:
        csr_n15 += [esp.csr_n15, 0]
        s_depths.append(esp.d_nonliq)
    else:
        s_depths.append(esp.d_nonliqs[0])
        s_depths.append(esp.d_liqs[1])
        s_depths.append(esp.d_nonliqs[1])
        csr_n15 += [esp.csr_n15s[0], 0, esp.csr_n15s[1], 0]
    s_depths.append(sp.height)

    for i in range(len(s_depths) - 1):
        if s_depths[i + 1] - s_depths[i] > 0:
            sl = sm.Soil()
            set_soil_props(sl, bi2014, i_top=esp.ith(s_depths[i]), i_bot=esp.ith(s_depths[i + 1]),  gwl=sp.gwl)
            if csr_n15[i]:
                sl.csr_n15 = csr_n15[i]
                sl.inputs.append("csr_n15")
            sp.add_layer(s_depths[i], sl)

    sp.set_soil_ids_to_layers()
    return sp


class Equiv3LayerProfile(sm.CustomObject):
    base_type = 'soil_profile'
    type = 'equivalent_3layer'
    n_layers = 3
    _q_c1ncs = None
    _e3_q_c1n_cs = None

    def __init__(self, d_liq, d_nonliq, crr_n15, gwl, norm_diff=None, depth=None, d_inc=None, d_start=None, d_end=None, surrogate=False, **kwargs):
        super(Equiv3LayerProfile, self).__init__()
        if depth is None:
            if d_inc is None or d_start is None or d_end is None:
                raise ValueError
            depth = np.arange(d_start, d_end, d_inc)
        self.depth = depth
        self.d_nonliq = d_nonliq
        self.d_liq = d_liq
        self.crr_n15 = crr_n15
        self.gwl = gwl
        self.norm_diff = norm_diff
        self._extra_class_inputs = ['d_liq', 'd_nonliq', 'crr_n15', 'gwl', 'norm_diff', 'd_inc', 'd_start', 'd_end']
        if surrogate:
            self.inputs = self._extra_class_inputs
        else:
            self.inputs += self._extra_class_inputs

    @property
    def d_end(self):
        return self.depth[-1]

    @property
    def d_inc(self):
        return self.depth[1] - self.depth[0]

    @property
    def d_start(self):
        return self.depth[0]

    @property
    def h_crust(self):
        return self.d_liq

    @property
    def h_liq(self):
        return self.d_nonliq - self.d_liq

    @property
    def h_nonliqs(self):
        return np.array([self.d_liq, self.d_end - self.d_nonliq])

    @property
    def q_c1n_cs(self):
        qc1ncs_vals = bi_functions.calc_q_c1n_cs_from_crr_m7p5(self.crr_n15)
        if self._e3_q_c1n_cs is None:
            self._e3_q_c1n_cs = create_profile_array([self.d_nonliq], [self.d_liq], self.depth, [qc1ncs_vals], 220.0)
        return self._e3_q_c1n_cs

    @property
    def crr_m7p5(self):
        return create_profile_array([self.d_nonliq], [self.d_liq], self.depth, [self.crr_n15], 1.0)

    def ith(self, y_depth):
        return np.where(self.depth >= y_depth)[0][0]


class Equiv5LayerProfile(sm.CustomObject):
    base_type = 'soil_profile'
    type = 'equivalent_5layer'
    n_layers = 5
    _q_c1n_cs = None

    def __init__(self, d_liqs, d_nonliqs, crr_n15s, gwl, norm_diff=None, depth=None, d_inc=None, d_start=None,
                 d_end=None, surrogate=False, **kwargs):

        super(Equiv5LayerProfile, self).__init__()
        if depth is None:
            if d_inc is None or d_start is None or d_end is None:
                raise ValueError
            depth = np.arange(d_start, d_end, d_inc)
        self.depth = depth
        self.d_nonliqs = d_nonliqs
        self.d_liqs = d_liqs
        self.crr_n15s = crr_n15s
        self.gwl = gwl
        self.norm_diff = norm_diff
        self._extra_class_inputs = ['d_liqs', 'd_nonliqs', 'crr_n15s', 'gwl', 'norm_diff', 'd_inc', 'd_start', 'd_end']
        if surrogate:
            self.inputs = self._extra_class_inputs
        else:
            self.inputs += self._extra_class_inputs

    @property
    def d_end(self):
        return self.depth[-1]

    @property
    def d_inc(self):
        return self.depth[1] - self.depth[0]

    @property
    def h_crust(self):
        return self.d_liqs[0]

    @property
    def d_start(self):
        return self.depth[0]

    @property
    def h_nonliqs(self):
        return np.array([self.d_liqs[0], self.d_liqs[1] - self.d_nonliqs[0], self.d_end - self.d_nonliqs[1]])

    @property
    def h_liqs(self):
        return np.array([self.d_nonliqs[0] - self.d_liqs[0], self.d_nonliqs[1] - self.d_liqs[1]])

    @property
    def q_c1n_cs(self):
        q_c1ncs_vals = bi_functions.calculate_qc_1ncs_from_crr_7p5(self.crr_n15s)
        if self._q_c1n_cs is None:
            self._q_c1n_cs = create_profile_array(self.d_nonliqs, self.d_liqs, self.depth, q_c1ncs_vals)
        return self._q_c1n_cs

    @property
    def crr_m7p5(self):
        return create_profile_array(self.d_nonliqs, self.d_liqs, self.depth, self.crr_n15s)

    def ith(self, y_depth):
        return np.where(self.depth >= y_depth)[0][0]


def esp_dict_to_obj(esp_dict):
    return Equiv3LayerProfile(**esp_dict)


def build_soil_profile_from_bi2014_w_layer_dict(bi2014, layers_dict):
    """
    Defined soil profile properties based on an equivalent soil profile define

    Parameters
    ----------
    bi2014: liquepy.trigger.BoulangerIdriss2014CPT
    esp: equivalent_soil object
        Equivalent soil profile of the triggering analysis

    Returns
    -------
    sm.SoilProfile
    """
    assert isinstance(bi2014, lq.trigger.BoulangerIdriss2014CPT)
    d_inc = bi2014.depth[3] - bi2014.depth[2]
    sp = sm.SoilProfile()
    sp.gwl = bi2014.gwl
    sp.height = bi2014.depth[-1]
    s_depths = []
    i_ths = []
    csr_n15 = []
    for crr in layers_dict:
        for depth in layers_dict[crr]:
            s_depths.append(depth)
            csr_n15.append(crr)
            i_ths.append((np.abs(bi2014.depth - depth - d_inc)).argmin())
    i_ths.append(len(bi2014.depth) - 1)
    s_depths.append(sp.height)
    csr_n15.append(-1)  # make same length

    i_ths = np.array(i_ths)
    s_depths = np.array(s_depths)
    csr_n15 = np.array(csr_n15)
    arr_inds = s_depths.argsort()
    i_ths = i_ths[arr_inds]
    s_depths = s_depths[arr_inds]
    csr_n15 = csr_n15[arr_inds]
    s_depths[0] = 0

    crr_values = lq.trigger.boulanger_and_idriss_2014.calc_crr_m7p5_from_qc1ncs(bi2014.q_c1n_cs)
    crr_values = np.clip(crr_values, None, 0.6)

    for i in range(len(s_depths) - 1):
        sl = sm.Soil()
        set_soil_props(sl, bi2014, i_top=i_ths[i] - 1, i_bot=i_ths[i + 1],  gwl=sp.gwl)
        # if csr_n15[i]:
        #     sl.csr_n15 = csr_n15[i]
        #     sl.inputs.append("csr_n15")
        sl.csr_n15 = np.mean(crr_values[i_ths[i] - 1:i_ths[i + 1]])
        sl.inputs.append("csr_n15")
        sp.add_layer(s_depths[i], sl)

    sp.set_soil_ids_to_layers()
    return sp, i_ths


def get_col_dict():
    """Define a dictionary of colors of different ESP (strength-position pairs)"""
    col_dict = OrderedDict([
        ("WT", [(1.0, 0.5, 0.4), "Weak thin"]),
        ("WM", [(1.0, 0.3, 0.25), "Weak mid-size"]),
        ("WL", [(0.8, 0.0, 0.15), "Weak large"]),
        ("MT", [(1.0, 0.9, 0.3), "Mid-strength thin"]),
        ("MM", [(1.0, 0.75, 0.15), "Mid-strength mid-size"]),
        ("ML", [(0.95, 0.55, 0.0), "Mid-strength large"]),
        ("ST", [(0.3, 0.8, 0.3), "Strong thin"]),
        ("SM", [(0.0, 0.6, 0.0), "Strong mid-size"]),
        ("SL", [(0.1, 0.4, 0.1), "Strong large"]),
        ("RX", [(0.6, 0.6, 0.6), "Resistant"]),
    ])
    return col_dict


def get_hatch_dict():
    """Define a dictionary of hatches for different ESP sizes"""
    hatch_dict = OrderedDict([
        ("S", [None, "Shallow"]),
        ("M", ["...", "Mid-height"]),
        ("D", ["o", "Deep"]),
        ("X", ["*", "N/A"])
        ])
    return hatch_dict

def get_esp_classes_legend_elements(include="all"):
    from matplotlib.patches import Patch
    hatch_dict = get_hatch_dict()
    col_dict = get_col_dict()
    legend_elements = []
    for name in col_dict:
        legend_elements.append(Patch(facecolor=col_dict[name][0],
                                     label=col_dict[name][1]))
    for name in hatch_dict:
        legend_elements.append(Patch(facecolor="w", edgecolor='k',
                                     label=hatch_dict[name][1], hatch=hatch_dict[name][0]))

    return legend_elements
