import liquepy as lq
import numpy as np
import eqsig
import pysra
import sfsimodels as sm


class EqlinStockwellAnalysis(object):

    def __init__(self, soil_profile, in_sig, rus=None, wave_field='outcrop', store='surface', gibbs=0, t_inc=1.0, t_win=3.0, strain_at_incs=True, strain_ratio=0.9):
        """
        Equivalent linear Stockwell Analysis

        This method performs the eight step procedure outlined in Millen et al. (2020) [TODO: update cite]
        to obtain the surface acceleration time series from an input motion at the base of a 1D soil profile.

        Note: a soil layer is a layer as defined in the soil_profile object, a slice is a layer from the pysra_profile.

        Parameters
        ----------
        soil_profile: sm.SoilProfile object
        in_sig: eqsig.AccSignal object
            Input motion at base
        rus: array_like or None
            A 2D array of pore pressure ratios of shape `(soil_profile.n_layers, in_sig.npts)`,
            if none the 'total stress' conditions are assumed
        wave_field: str
            If input motion should be used as an `outcrop` or '`within` motion.
        store: str
            if 'surface' (default), it only stores the surface acceleration time series,
            if 'all' then stores Stockwell transforms.
        gibbs: int or None
            If integer then zero-pad input motion to next power of 2, to reduce Gibb's effect
        t_inc: float (default=1s)
            Time increment of interval for determining transfer functions
        t_win: float (default=3s)
            Time window for determining maximum strain value
        strain_at_incs: bool (default=True)
            If true then compute effective strain at time intervals, else use same value for full time series
        strain_ratio: float (default=0.9)
            Ratio between effective strain and peak (maximum) strain
        """

        assert isinstance(soil_profile, sm.SoilProfile)
        org_npts = in_sig.npts
        # Determine number of zeros required for zero padding
        if gibbs is not None:  # If it is an integer then add to the exponent of 2 to remove the Gibbs effect
            nindex = int(np.ceil(np.log2(org_npts))) + gibbs
            new_len = 2 ** nindex
            diff_len = new_len - org_npts
            front = 0  # int(diff_len / 2)
            back = diff_len - front
        else:  # record length must be a factor of 4
            back = int(4 * np.ceil(org_npts / 4) - org_npts)
            front = 0
        # pad the input signal with zeros to make length a factor of 4
        in_sig = eqsig.AccSignal(np.pad(in_sig.values, (front, back), mode='constant'), in_sig.dt)

        self.t_inds = np.arange(0, in_sig.npts - 1, int(t_inc / in_sig.dt), dtype=int)  # indices of time intervals
        self.t_inds = np.insert(self.t_inds, len(self.t_inds), in_sig.npts)  # make sure last value is in list
        ics = np.array((self.t_inds[1:] + self.t_inds[:-1]) / 2, dtype=int)  # halfway between indices of time intervals

        points = int(in_sig.npts / 2)
        freqs = np.arange(0, points) / (points * in_sig.dt * 2)  # All the frequencies in the Stockwell transform
        # freqs_d2 = freqs[:int(points / 2)]  # Frequencies needed to compute the transfer function

        # Steps 1 & 2) Conduct and equivalent linear analysis and obtain strain time series
        pysra_profile, strains = compute_pysra_strain_time_series(soil_profile, in_sig, target_height=0.5,
                                                                  wave_field=wave_field)

        # 3a) Calculate the effective strain in each time interval for each slice of the soil profile
        iside = int(t_win / in_sig.dt)  # width of window in increments
        eff_strains = []
        for i, depth in enumerate(pysra_profile.depth):
            eff_strains.append([])
            for tt in range(len(ics)):
                if strain_at_incs:
                    si = max([ics[tt] - iside, 0])
                    ei = min([ics[tt] + iside, len(strains[i])])
                    max_strain = max(abs(strains[i][si: ei]))
                else:  # Note this is different to the pysra eff. strain -which estimates the peak from the strain tf
                    max_strain = max(abs(strains[i]))
                eff_strains[i].append(strain_ratio * max_strain)

        # 4a) Obtain the reduction in secant stiffness and increase in damping from pore pressure at each time interval

        # Parameters defined in Millen et al. (2020)
        dxi_ld_liq = 0.3  # (Delta-xi-low-density-at-liquefaction)
        dxi_hd_liq = 0.1  # (Delta-xi-high-density-at-liquefaction)
        gr_ld_liq = 0.03  # (secant-shear-modulus-ratio-low-density-at-liquefaction)
        gr_hd_liq = 0.15  # (secant-shear-modulus-ratio-high-density-at-liquefaction)
        x_ld = 0.45  # Low density threshold
        x_hd = 0.8  # high density threshold
        ru_gr_low = 0.3  # Low pore pressure ratio threshold for secant shear modulus change
        ru_gr_high = 0.8  # High pore pressure ratio threshold for secant shear modulus change
        ru_dx_low = 0.5  # Low pore pressure ratio threshold for damping increment change
        ru_dx_high = 1.0  # High pore pressure ratio threshold for damping increment change
        min_g_liq_vs_g0 = 0.001  # Limiting ratio between the shear modulus at liquefaction divided by initial shear mod
        max_xi_liq = 0.3  # Maximum value for damping

        # The arrays for the secant stiffness ratio and damping increase at each time interval from pore pressure
        gred_is = np.ones((soil_profile.n_layers, len(ics)))
        dxi_is = np.zeros((soil_profile.n_layers, len(ics)))
        if rus is not None:  # if pore pressure time series is defined then calculate the pore pressure corrections
            assert len(rus[0]) == org_npts, (len(rus[0]), org_npts)
            gr_liqs = np.ones(soil_profile.n_layers)  # The secant stiffness ratio at liquefaction for each soil layer
            dxi_liqs = np.zeros(soil_profile.n_layers)  # The damping increase at liquefaction for each soil layer
            for i in range(soil_profile.n_layers):
                dr = soil_profile.layer(i + 1).relative_density
                if dr is None:
                    if max(rus[i]):
                        raise ValueError('Relative density must be set for layer: ', i + 1)
                    else:
                        continue
                # Calculate the secant stiffness ratio at liquefaction based on relative density
                gr_liq = np.where(dr < x_ld, gr_ld_liq, gr_ld_liq + (gr_hd_liq - gr_ld_liq) / (x_hd - x_ld) * (dr - x_ld))
                np.clip(gr_liq, None, gr_hd_liq, out=gr_liq)
                gr_liqs[i] = gr_liq
                # Calculate the damping increase at liquefaction based on relative density
                dx_max = np.where(dr < x_ld, dxi_ld_liq, dxi_ld_liq + (dxi_hd_liq - dxi_ld_liq) / (x_hd - x_ld) * (dr - x_ld))
                np.clip(dx_max, dxi_hd_liq, None, out=dx_max)
                dxi_liqs[i] = dx_max

            # zero pad pore pressure time series to be consistent with acceleration time series
            rus = np.pad(rus, [(0, 0), (front, back)], mode='constant')
            # Calculate the secant stiffness ratio at each time step based on pore pressure ratio (ru)
            greds = np.where(rus < ru_gr_low, 1, 1 - (1 - gr_liqs[:, np.newaxis]) / (ru_gr_high - ru_gr_low) * (rus - ru_gr_low))
            np.clip(greds, gr_liqs[:, np.newaxis], None, out=greds)
            # Calculate the damping increase at each time step based on pore pressure ratio (ru)
            dxs = np.where(rus < ru_dx_low, 0, dxi_liqs[:, np.newaxis] / (ru_dx_high - ru_dx_low) * (rus - ru_dx_low))
            np.clip(dxs, None, dxi_liqs[:, np.newaxis], out=dxs)

            # Calculate the secant stiffness ratio and damping increase at each time interval
            for tt in range(len(ics)):
                gred_is[:, tt] = np.mean(greds[:, self.t_inds[tt]: self.t_inds[tt + 1]], axis=1)
                dxi_is[:, tt] = np.mean(dxs[:, self.t_inds[tt]: self.t_inds[tt + 1]], axis=1)

        # 5) Develop input-to-surface transfer functions
        self.tfs = []  # A list to store the transfer functions at each increment

        for tt in range(len(self.t_inds[1:])):
            layers = []
            for i, depth in enumerate(pysra_profile.depth):
                org_layer = pysra_profile.location('outcrop', depth=depth).layer
                org_layer.strain = eff_strains[i][tt]  # Apply effective strain (Step 3b)
                shear_vel0 = org_layer.initial_shear_vel
                shear_vel = org_layer.shear_vel
                damping = org_layer.damping
                slice_thickness = org_layer.thickness
                # get pore pressure effects
                ind = soil_profile.get_layer_index_by_depth(depth) - 1
                dx = dxi_is[ind][tt]
                gred = gred_is[ind][tt]
                # 4b) determine the new shear modulus and damping accounting for strain and pore pressure
                xi_liq = min([damping + dx, max_xi_liq])
                vs_liq = max([np.sqrt(min_g_liq_vs_g0) * shear_vel0, shear_vel * np.sqrt(gred)])
                pysra_sl = pysra.site.SoilType("soil", org_layer.unit_wt, None, xi_liq)
                lay = pysra.site.Layer(pysra_sl, slice_thickness, vs_liq)
                layers.append(lay)

            # rebuild the pysra_profile with the new properties
            strain_comp_profile = pysra.site.Profile(layers, wt_depth=soil_profile.gwl)
            # determine the new transfer function for this interval
            freq1, tf_values = lq.sra.calc_pysra_tf(strain_comp_profile, freqs, wave_field=wave_field)
            # refactor transfer function to be applied to Stockwell transform
            tf_values = np.flipud(np.conj(tf_values))
            # tf_values = np.concatenate((tf_values, np.flipud(np.conj(tf_values))))
            tf_values = tf_values.reshape(len(tf_values), 1)
            self.tfs.append(tf_values)

        # 6) Obtain the Stockwell transform of the input motion
        in_sig.swtf = eqsig.stockwell.transform(in_sig.values)
        # 7) Obtain the surface Stockwell transform by multiplying input Stockwell transform by transfer functions
        ps = []
        for ss in range(len(self.tfs)):
            p1 = self.tfs[ss] * in_sig.swtf[:, self.t_inds[ss]:self.t_inds[ss + 1]]
            ps.append(p1)
        surf_st = np.concatenate(ps, axis=1)

        # 8) Perform the inverse Stockwell transform to obtain the surface acceleration time series
        iacc = eqsig.stockwell.itransform(surf_st)

        # save the surface acceleration series as a parameter
        self.surf_sig = eqsig.AccSignal(iacc, in_sig.dt)

        if store == 'all':  # Store the Stockwell transforms of the input motion if needed
            # self.tfs = tfs
            self.freqs = freqs
            self.in_sig = in_sig
            self.surf_sig.stockwell = surf_st
            self.in_sig.smooth_freqs = np.linspace(0.2, 1 / (4 * in_sig.dt), 30)
            self.surf_sig.smooth_freqs = np.linspace(0.2, 1 / (4 * in_sig.dt), 30)


def compute_pysra_strain_time_series(soil_profile, in_sig, d_inc=None, target_height=1.0, wave_field='outcrop'):
    """
    Perform an equivalent linear analysis and obtain the strain time series at many depths

    Parameters
    ----------
    soil_profile: sm.SoilProfile object
        in_sig: eqsig.AccSignal object
            Input motion at base
    d_inc: float
        Target depth increment for each layer in soil_profile
    target_height: float
        Target depth increment for whole soil profile
    wave_field: str
            If input motion should be used as an `outcrop` or '`within` motion.

    Returns
    -------

    """
    import pysra
    m = pysra.motion.TimeSeriesMotion(filename=in_sig.label, description=None, time_step=in_sig.dt,
                                      accels=in_sig.values / 9.8)
    if d_inc is None:
        d_inc = 1.0 * np.ones(soil_profile.n_layers)
    profile = lq.sra.sm_profile_to_pysra(soil_profile, target_height=target_height, d_inc=d_inc)

    calc = pysra.propagation.EquivalentLinearCalculator()
    calc(m, profile, profile.location(wave_field, depth=soil_profile.height))

    outs = []
    for i, depth in enumerate(profile.depth):
        outs.append(pysra.output.StrainTSOutput(pysra.output.OutputLocation('within', depth=depth),
                                                in_percent=False))

    outputs = pysra.output.OutputCollection(outs)
    outputs(calc)
    strains = []
    for i, depth in enumerate(profile.depth):
        strains.append(outputs[i].values)

    return profile, strains
