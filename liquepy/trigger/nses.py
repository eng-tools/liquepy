import numpy as np
import eqsig


def calc_energy_ratio_w_time(xi, total_time, time, av_period):
    up_reduction = np.exp(-xi * (total_time - time) * 2 * np.pi / av_period) ** 2
    surf_reduction = np.exp(-xi * total_time * 2 * np.pi / av_period) ** 2
    down_reduction = np.exp(-xi * time * 2 * np.pi / av_period) ** 2 * surf_reduction
    return (up_reduction + down_reduction) / 2


def est_case_1d_millen_et_al_2019(sp, asig, depth, xi, g_mod_red=1.0, trim=False, start=False, period=0.5, exact=False,
                                   in_loc=1, g_scale_limit=1e3, nodal=True, cace=True):
    """
    Calculates the Cumulative absolute change in strain energy according to Millen et al. (2019)

    Parameters
    ----------
    sp: sfsimodels.SoilProfile object
    asig: eqsig.AccSignal object
        upward propagating wave at base of profile
    depth: the distance from the surface where energy should be estimated
    xi: float
        approximate viscous damping of whole soil profile
    g_mod_red: float
        reduction factor to apply to shear modulus
    trim: bool
        If true to return CASE with same length as asig.npts
    start: bool,
        if True then accounts for travel time from input location to depth
    period: float
        The average energy period of the ground motion
    exact: bool,
        If exact then assume wave is sine wave and reduce upward and downward components
    in_loc: int
        Location of input motion, base=1, surface=0
    g_scale_limit: int or float
        Limiting ratio that CASE can be scaled by when changing shear wave velocity
    nodal: bool
        If true then surface is a nodal (zero stress), if false then surface is anti-nodal (zero displacement)

    :return:
    """
    from scipy.interpolate import interp1d
    sp.gen_split(props=['shear_vel', 'unit_mass'], target=0.25)
    dis_depths = np.cumsum(sp.split["thickness"])
    split_depths = np.cumsum(sp.split["thickness"])

    dis_depths = np.insert(dis_depths, 0, 0)
    dis_shear_vel = sp.split["shear_vel"] * np.sqrt(g_mod_red)
    split_g_mod = dis_shear_vel ** 2 * sp.split["unit_mass"]
    travel_times = sp.split["thickness"] / dis_shear_vel
    dis_time_from_surface = np.cumsum(travel_times)
    dis_time_from_surface = np.insert(dis_time_from_surface, 0, 0)
    time_at_depth = np.interp(depth, dis_depths, dis_time_from_surface)
    total_time = dis_time_from_surface[-1]

    vs = interp1d(dis_depths[:-1], dis_shear_vel, kind='previous', fill_value='extrapolate')(depth)
    rho = interp1d(dis_depths[:-1], sp.split["unit_mass"], kind='previous', fill_value='extrapolate')(depth)

    g_mod = vs ** 2 * rho

    # tau is conserved across a boundary, so E=tau^2/G, so imp^2
    if in_loc == 0:
        in_depth = 0
    else:
        in_depth = sp.height

    if exact:
        up_red = np.exp(-xi * (total_time - time_at_depth) * 2 * np.pi / period) ** 2
        surf_reduction = np.exp(-xi * total_time * 2 * np.pi / period) ** 2
        if in_loc == 0:  # reset to 1 at surface
            soil_in = sp.get_soil_at_depth(0)
            up_red = up_red - surf_reduction + 1
            surf_reduction = 1
            stt = 0.0
        else:
            soil_in = sp.get_soil_at_depth(sp.height)
            stt = total_time
        down_red = np.exp(-xi * time_at_depth * 2 * np.pi / period) ** 2 * surf_reduction
        if cace:
            spectra_series = eqsig.surface.calc_cum_abs_surface_energy(asig, time_at_depth, up_red=up_red,
                                                                   down_red=down_red, trim=trim, nodal=nodal, stt=stt)
        else:
            spectra_series = eqsig.surface.calc_surface_energy(asig, time_at_depth, up_red=up_red,
                                                                       down_red=down_red, trim=trim, nodal=nodal,
                                                                       stt=stt)
        spectra_series = np.asarray(spectra_series)
    else:
        red_at_surf = calc_energy_ratio_w_time(xi, total_time, 0, av_period=period)
        if in_loc == 0:
            soil_in = sp.get_soil_at_depth(0)
            red_ratio = 1 + (1 - red_at_surf) / 2 * time_at_depth / total_time
            stt = 0.0
        else:
            soil_in = sp.get_soil_at_depth(sp.height)
            red_ratio = red_at_surf + (1 - red_at_surf) / 2 * time_at_depth / total_time
            stt = total_time
        if cace:
            spectra_series = eqsig.surface.calc_cum_abs_surface_energy(asig, time_at_depth, trim=trim, nodal=nodal, stt=stt,
                                                                   start=start)
        else:
            spectra_series = eqsig.surface.calc_surface_energy(asig, time_at_depth, trim=trim, nodal=nodal, stt=stt, start=start)
        if hasattr(red_ratio, '__len__'):
            spectra_series = red_ratio[:, np.newaxis] * np.asarray(spectra_series)
        else:
            spectra_series *= red_ratio
    rho_in = soil_in.unit_dry_weight / 9.8
    g_in = np.interp(in_depth, split_depths, split_g_mod)
    g_scale = (g_mod / g_in)
    if not nodal:
        g_scale_limit = 1.0
    g_scale = np.clip(g_scale, 1.0 / g_scale_limit, g_scale_limit)  # simple scaling from Millen et al. (2019)
    if hasattr(g_scale, '__len__'):
        estimate = spectra_series * rho_in / g_scale[:, np.newaxis]
    else:
        estimate = spectra_series * rho_in / g_scale
    return estimate


def calc_time_surf_to_depth(sp, depth):
    sp.gen_split(props=['shear_vel', 'unit_mass'])
    edge_depths = np.zeros(len(sp.split["thickness"]) + 1)
    edge_depths[1:] = np.cumsum(sp.split["thickness"])
    travel_times = sp.split["thickness"] / sp.split["shear_vel"]
    times_from_surface = np.zeros_like(edge_depths)
    times_from_surface[1:] = np.cumsum(travel_times)
    return np.interp(depth, edge_depths, times_from_surface)


class TimeShiftProfile(object):  # Under development
    _exact = None
    _kinetic_energy_series = None
    _strain_energy_series = None

    def __init__(self, sp, asig, ys, g_mod_red=1.0, period=0.5, exact=False, in_loc=1):
        self.sp = sp
        self.asig = asig
        self.g_mod_red = g_mod_red
        self.period = period
        self._exact = exact

        self.sp.gen_split(props=['shear_vel', 'unit_mass'])
        edge_depths = np.cumsum(sp.split["thickness"])
        self.edge_depths = np.insert(edge_depths, 0, 0)

        self.shear_vel = sp.split["shear_vel"] * np.sqrt(g_mod_red)
        self.g_mod = self.shear_vel ** 2 * sp.split["unit_mass"]
        self.rho = sp.split["unit_mass"]
        travel_times = sp.split["thickness"] / self.shear_vel
        dis_time_from_surface = np.cumsum(travel_times)
        self.times_from_surface = np.insert(dis_time_from_surface, 0, 0)
        self.total_time = self.times_from_surface[-1]
        self.in_loc = in_loc

        if self.in_loc == 0:  # reset to 1 at surface
            self.rho_in = np.interp(0, self.edge_depths[:-1], self.rho)
            self.g_mod_in = np.interp(0, self.edge_depths[:-1], self.g_mod)
        else:
            self.rho_in = np.interp(self.sp.height, self.edge_depths[:-1], self.rho)
            self.g_mod_in = np.interp(self.sp.height, self.edge_depths[:-1], self.g_mod)
        # depths = np.insert(self.edge_depths, 0, 0)
        self.time_ys = np.interp(ys, self.edge_depths, self.times_from_surface)
        self.time_shifts = 2 * self.time_ys
        self.g_mod_ys = np.interp(ys, self.edge_depths[1:], self.g_mod)

        self.shifts = np.array(self.time_shifts / asig.dt, dtype=int)
        self.up_wave = np.pad(asig.values, (0, np.max(self.shifts)), mode='constant', constant_values=0)  # 1d
        self.down_waves = eqsig.put_array_in_2d_array(asig.values, self.shifts)
        if self.in_loc == 0:
            self.start_time = -self.time_ys
        else:
            self.start_time = self.total_time - self.time_ys
        self.e_cache = {}
        # self.clear_cache()

    @property
    def exact(self):
        return self._exact

    # def clear_cache(self):
    #     self._kinetic_energy_series = None
    #     self._strain_energy_series = None

    def get_red_factors(self, xi):  # TODO: cache properties
        if self.exact:
            up_amp_red = np.exp(-xi * (self.total_time - self.time_ys) * 2 * np.pi / self.period)
            surf_reduction = np.exp(-xi * self.total_time * 2 * np.pi / self.period)
            if self.in_loc == 0:  # reset to 1 at surface
                up_amp_red = up_amp_red - surf_reduction + 1
                surf_reduction = 1
            down_amp_red = np.exp(-xi * self.time_ys * 2 * np.pi / self.period) * surf_reduction
            total_red = np.ones_like(down_amp_red)
        else:
            red_at_surf = calc_energy_ratio_w_time(xi, self.total_time, 0, av_period=self.period)
            if self.in_loc == 0:
                total_red = 1 + (1 - red_at_surf) / 2 * self.time_ys / self.total_time
            else:
                total_red = red_at_surf + (1 - red_at_surf) / 2 * self.time_ys / self.total_time
            up_amp_red = np.ones_like(total_red)
            down_amp_red = np.ones_like(total_red)
        return up_amp_red, down_amp_red, total_red

    def get_unit_energy_series(self, xi, energy):
        from scipy.integrate import cumtrapz
        e_str = "{0}-{1}".format(xi, energy)
        if e_str in self.e_cache:
            return self.e_cache[e_str]
        up_amp_red, down_amp_red, total_red = self.get_red_factors(xi)

        if energy == 'kinetic':
            acc_series = down_amp_red * self.down_waves + up_amp_red * self.up_wave
        else:
            acc_series = - down_amp_red[:, np.newaxis] * self.down_waves + up_amp_red[:, np.newaxis] * self.up_wave
        velocity = cumtrapz(acc_series, dx=self.asig.dt, initial=0)
        unit_kinetic_energy = 0.5 * velocity * np.abs(velocity) * total_red[:, np.newaxis]
        self.e_cache[e_str] = unit_kinetic_energy
        return unit_kinetic_energy

    def trim_to_length(self, out, trim=False, start=False):  # TODO: switch to use esig function
        # if not trim:
        #     return out
        total_shift = int(self.total_time / self.asig.dt)
        surf_to_depth_shift = self.shifts / 2
        surf_to_depth_shift = surf_to_depth_shift.astype(int)
        if start:  # Trim front
            if self.in_loc:
                sis = total_shift - surf_to_depth_shift
            else:
                sis = -surf_to_depth_shift
            if trim:
                npts = self.asig.npts
            else:
                extras = np.max([sis, 0]) - np.min([np.min(self.shifts), 0])
                npts = self.asig.npts + extras  # plus the zero padding in front and back
        else:
            if trim:
                sis = np.zeros_like(total_shift)
                npts = self.asig.npts
            else:  # no changes required
                return out
        print('npts: ', npts)

        outs = np.zeros((len(self.shifts), npts))
        print('len_outs: ', outs.shape)
        for i in range(len(self.shifts)):
            if sis[i] < 0:
                outs[i] = out[i, -sis[i]: npts - sis[i]]
            else:
                outs[i, sis[i]:] = out[i, : npts - sis[i]]  # zero padded
        return outs

    def get_cake(self, xi, trim=False, start=False, g_scale_limit=1.):
        kin_energy = self.get_unit_energy_series(xi, energy='kinetic')
        delta_energy = np.zeros_like(kin_energy)
        delta_energy[:, 0] = kin_energy[:, 0]
        delta_energy[:, 1:] = np.diff(kin_energy, axis=1)
        cake = self.rho_in * np.cumsum(abs(delta_energy))
        g_scale = np.clip(self.g_mod_ys / self.g_mod_in, 1. / g_scale_limit, g_scale_limit)
        return self.trim_to_length(cake, trim=trim, start=start) / g_scale

    def get_case(self, xi, trim=False, start=False, g_scale_limit=1.):
        strain_energy = self.get_unit_energy_series(xi, energy='strain')
        delta_energy = np.zeros_like(strain_energy)
        delta_energy[:, 0] = strain_energy[:, 0]
        delta_energy[:, 1:] = np.diff(strain_energy, axis=1)
        case = self.rho_in * np.cumsum(abs(delta_energy), axis=1)
        g_scale = np.clip(self.g_mod_ys / self.g_mod_in, 1. / g_scale_limit, g_scale_limit)
        return self.trim_to_length(case, trim=trim, start=start) / g_scale[:, np.newaxis]

    def get_kinetic_energy(self, xi, trim=False, start=False, g_scale_limit=1.):
        ke = self.rho_in * self.get_unit_energy_series(xi, energy='kinetic')
        g_scale = np.clip(self.g_mod_ys / self.g_mod_in, 1. / g_scale_limit, g_scale_limit)
        return self.trim_to_length(ke, trim=trim, start=start) / g_scale[:, np.newaxis]

    def get_strain_energy(self, xi, trim=False, start=False, g_scale_limit=1.):
        se = self.rho_in * self.get_unit_energy_series(xi, energy='strain')
        g_scale = np.clip(self.g_mod_ys / self.g_mod_in, 1. / g_scale_limit, g_scale_limit)
        return self.trim_to_length(se, trim=trim, start=start) / g_scale[:, np.newaxis]
