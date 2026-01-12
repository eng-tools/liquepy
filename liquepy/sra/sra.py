import numpy as np
import sfsimodels as sm

PA_TO_KPA = 0.001


def compute_pysra_tf(pysra_profile, pysra_freqs=None):
    return calc_pysra_tf(pysra_profile, pysra_freqs)


def calc_pysra_tf(
    pysra_profile, pysra_freqs=None, wave_field="outcrop", absolute=False
):
    import pysra

    if pysra_freqs is None:
        pysra_freqs = np.logspace(-0.7, 1.5, num=200)
    m = pysra.motion.Motion(freqs=pysra_freqs)
    outputs = pysra.output.OutputCollection(
        [
            pysra.output.AccelTransferFunctionOutput(
                pysra_freqs,
                pysra.output.OutputLocation(wave_field, index=-1),
                pysra.output.OutputLocation("outcrop", index=0),
                absolute=absolute,
            ),
        ]
    )
    calc = pysra.propagation.LinearElasticCalculator()
    calc(m, pysra_profile, pysra_profile.location(wave_field, index=-1))
    outputs(calc)
    out_liq_tf = outputs[0].values
    return pysra_freqs, out_liq_tf


def vardanega_2013_to_modified_hyperbolic_parameters(i_p):
    a = 0.943  # Eq 22b
    j = 3.7  # Eq 23
    gamma_ref = j * (i_p / 1000)
    return gamma_ref, a


def sm_profile_to_pysra(
    sp,
    d_inc=None,
    target_height=1.0,
    base_shear_vel=None,
    base_unit_wt=None,
    base_xi=0.01,
    vs_min=0,
    k0=0.5,
):
    """
    Converts a soil profile from sfsimodels into a soil profile for pysra

    Note: pysra uses kPa whereas sfsimodels uses Pa

    :param sp:
    :param d_inc:
    :param target_height:
    :return:
    """
    import pysra

    if d_inc is None:
        d_inc = np.ones(sp.n_layers) * target_height

    strains = np.logspace(-6, -1.0, num=30)

    layers = []
    cum_thickness = 0
    for i in range(sp.n_layers):
        if sp.layer_depth(i + 1) >= sp.height:
            break
        sl = sp.layer(i + 1)
        thickness = sp.layer_height(i + 1)

        n_slices = max(1, int(thickness / d_inc[i]))

        slice_thickness = float(thickness) / n_slices
        for j in range(n_slices):
            cum_thickness += slice_thickness
            if cum_thickness > sp.gwl:
                rho = sl.unit_sat_mass
            else:
                rho = sl.unit_dry_mass
            v_eff = sp.get_v_eff_stress_at_depth(cum_thickness)
            if hasattr(sl, "get_g_mod_at_v_eff_stress"):
                g_mod = sl.get_g_mod_at_v_eff_stress(v_eff)
            else:
                g_mod = sl.g_mod
            if hasattr(sl, "g_mod_red"):
                g_mod *= sl.g_mod_red
            vs = np.sqrt(g_mod / rho)
            if cum_thickness > sp.gwl:
                unit_wt = sl.unit_sat_weight
            else:
                unit_wt = sl.unit_dry_weight
            if not hasattr(sl, "sra_type"):
                raise ValueError(f"sra_type missing on soil ({i+1})")
            if getattr(sl, "sra_type") == "hyperbolic":
                name = "hyperbolic"
                if hasattr(sl, "strain_curvature") and hasattr(sl, "strain_ref"):
                    strain_curvature = sl.strain_curvature
                    strain_ref = sl.strain_ref
                    xi_min = sl.xi_min
                    pass
                elif hasattr(sl, "peak_strain"):
                    # Octahedral shear stress
                    p_ref = v_eff * 2.0 / 3
                    tau_f = (2 * np.sqrt(2.0) * np.sin(sl.phi_r)) / (
                        3 - np.sin(sl.phi_r)
                    ) * p_ref + 2 * np.sqrt(2.0) / 3 * sl.cohesion
                    strain_curvature = 1.0
                    strain_ref = (
                        sl.peak_strain * tau_f / (g_mod * sl.peak_strain - tau_f)
                    )
                    xi_min = sl.xi_min
                    # inputs += ['strain_curvature', 'xi_min', 'sra_type', 'strain_ref']
                elif sl.type == "pm4sand":
                    from liquepy.num.models import calc_peak_angle_for_pm4sand

                    curvature = 0.75
                    ratio = 60
                    msig = v_eff * (1 + 1 * k0) / 2
                    phi_peak = calc_peak_angle_for_pm4sand(sl.relative_density, msig)
                    g_mod_min = vs_min**2 * unit_wt / 9.8
                    g_mod0 = max(sl.get_g_mod_at_m_eff_stress(msig), g_mod_min)
                    vs = np.sqrt(g_mod0 / rho)
                    tau_max = v_eff * np.sin(np.radians(phi_peak))
                    strain_curvature = curvature
                    xi_min = 0.015
                    strain_ref = (
                        tau_max * (1 + ratio**strain_curvature) / (g_mod0 * ratio)
                    )
                elif sl.type == "pm4silt":
                    alpha = 0.919
                    ratio = 20
                    # alpha = 0.88
                    # ratio = 30
                    phi_cv = 32.0  # page 67 of manual
                    phi_rat = np.sin(np.radians(phi_cv))
                    if sl.su_rat is None:
                        if sl.su is None:
                            raise ValueError("Must set either su or su_rat")
                        sur = sl.su / v_eff
                    else:
                        sur = sl.su_rat
                    tau_max = min([phi_rat, sur]) * v_eff
                    msig = v_eff * (1 + 1 * k0) / 2
                    g_mod_min = vs_min**2 * unit_wt / 9.8
                    g_mod0 = max(sl.get_g_mod_at_m_eff_stress(msig), g_mod_min)
                    vs = np.sqrt(g_mod0 / rho)

                    xi_min = 0.015
                    strain_curvature = alpha
                    strain_ref = tau_max * (1 + ratio**alpha) / (g_mod0 * ratio)
                elif (
                    hasattr(sl, "plasticity_index")
                    and getattr(sl, "plasticity_index") is not None
                ):
                    i_p = sl.plasticity_index
                    strain_ref, strain_curvature = (
                        vardanega_2013_to_modified_hyperbolic_parameters(i_p)
                    )
                    xi_min = 0.02
                    name = "vardanega (2013) I_p = %.2f" % i_p
                else:
                    raise ValueError(
                        f"Cannot build modified hyperbolic for soil ({i+1})"
                    )
                pysra_sl = pysra.site.ModifiedHyperbolicSoilType(
                    name,
                    unit_wt,
                    strain_ref=strain_ref,
                    curvature=strain_curvature,
                    damping_min=xi_min,
                    strains=strains,
                )
            elif getattr(sl, "sra_type") == "darendeli":
                if hasattr(sl, "darendeli_sigma_m_eff"):
                    ip = sl.plasticity_index
                    if ip is None:
                        ip = 0.0
                    ip *= 100  # Input is in percentage
                    pysra_sl = pysra.site.DarendeliSoilType(
                        unit_wt,
                        plas_index=ip,
                        ocr=1,
                        stress_mean=sl.darendeli_sigma_m_eff,
                        strains=strains,
                    )
                else:
                    assert isinstance(sp, sm.SoilProfile)

                    s_v_eff = sp.get_v_eff_stress_at_depth(cum_thickness)
                    k0 = 1 - np.sin(np.radians(sl.phi))
                    darendeli_sigma_m_eff = (
                        s_v_eff * (1 + 2 * k0) / 3
                    ) * PA_TO_KPA  # Needs to be in kPa
                    # print("darendeli_sigma_m_eff: ", darendeli_sigma_m_eff)
                    ip = sl.plasticity_index
                    if ip is None:
                        ip = 0.0
                    ip *= 100  # Input is in percentage
                    pysra_sl = pysra.site.DarendeliSoilType(
                        unit_wt,
                        plas_index=ip,
                        ocr=1,
                        stress_mean=darendeli_sigma_m_eff,
                        strains=strains,
                    )

            elif getattr(sl, "sra_type") == "menq":
                k0 = 0.5
                s_v_eff = sp.vertical_effective_stress(cum_thickness)
                sigma_m_eff = (s_v_eff * (1 + 2 * k0) / 3) * PA_TO_KPA
                pysra_sl = pysra.site.MenqSoilType(
                    unit_wt, uniformity_coeff=2, diam_mean=2, stress_mean=sigma_m_eff
                )

            elif getattr(sl, "sra_type") == "linear":
                if hasattr(sl, "xi"):
                    xi = sl.xi
                elif hasattr(sl, "xi_min"):
                    xi = sl.xi_min
                else:
                    raise ValueError("must set xi_min or xi for linear analysis")
                pysra_sl = pysra.site.SoilType(sl.name, unit_wt, None, xi)
            elif getattr(sl, "sra_type") == "direct":
                if sl.type == "pm4sand":
                    name = "pm4sand-direct"
                    from liquepy.num.models import calc_peak_angle_for_pm4sand

                    curvature = 0.75
                    ratio = 60
                    msig = v_eff * (1 + 1 * k0) / 2
                    phi_peak = calc_peak_angle_for_pm4sand(sl.relative_density, msig)
                    g_mod_min = vs_min**2 * unit_wt / 9.8
                    g_mod0 = max(sl.get_g_mod_at_m_eff_stress(msig), g_mod_min)
                    vs = np.sqrt(g_mod0 / rho)
                    tau_max = v_eff * np.sin(np.radians(phi_peak))
                    strain_curvature = curvature
                    xi_min = 0.02
                    strain_ref = (
                        tau_max * (1 + ratio**strain_curvature) / (g_mod0 * ratio)
                    )
                    pysra_sl = pysra.site.ModifiedHyperbolicSoilType(
                        name,
                        unit_wt,
                        strain_ref=strain_ref,
                        curvature=strain_curvature,
                        damping_min=xi_min,
                        strains=strains,
                    )
                    pysra_d_damping = pysra.site.NonlinearProperty(
                        "", strains, 1.3 * pysra_sl.damping.values, "damping"
                    )
                    pysra_strains = pysra_sl.mod_reduc.strains
                    pysra_g_o_gmax = pysra_sl.mod_reduc
                    pysra_sl = pysra.site.SoilType(
                        name, unit_wt, pysra_g_o_gmax, pysra_d_damping
                    )
                elif sl.type == "pm4silt":
                    name = "pm4silt-direct"
                    alpha = 0.919
                    ratio = 20
                    # alpha = 0.88
                    # ratio = 30
                    phi_cv = 32.0  # page 67 of manual
                    phi_rat = np.sin(np.radians(phi_cv))
                    if sl.su_rat is None:
                        if sl.su is None:
                            raise ValueError("Must set either su or su_rat")
                        sur = sl.su / v_eff
                    else:
                        sur = sl.su_rat
                    tau_max = min([phi_rat, sur]) * v_eff
                    msig = v_eff * (1 + 1 * k0) / 2
                    g_mod_min = vs_min**2 * unit_wt / 9.8
                    g_mod0 = max(sl.get_g_mod_at_m_eff_stress(msig), g_mod_min)
                    vs = np.sqrt(g_mod0 / rho)

                    xi_min = 0.02
                    strain_curvature = alpha
                    strain_ref = tau_max * (1 + ratio**alpha) / (g_mod0 * ratio)
                    pysra_sl = pysra.site.ModifiedHyperbolicSoilType(
                        name,
                        unit_wt,
                        strain_ref=strain_ref,
                        curvature=strain_curvature,
                        damping_min=xi_min,
                        strains=strains,
                    )
                    pysra_d_damping = pysra.site.NonlinearProperty(
                        "", strains, 1.3 * pysra_sl.damping.values, "damping"
                    )
                    pysra_strains = pysra_sl.mod_reduc.strains
                    pysra_g_o_gmax = pysra_sl.mod_reduc
                    pysra_sl = pysra.site.SoilType(
                        name, unit_wt, pysra_g_o_gmax, pysra_d_damping
                    )
            else:
                raise ValueError(
                    f"sra_type = {sl.sra_type} not recognised for soil: {i+1}"
                )
            lay = pysra.site.Layer(pysra_sl, slice_thickness, vs)
            layers.append(lay)
    # add one more since it is applied at top of layer, but make it elastic
    if base_unit_wt is None:
        base_unit_wt = layers[-1].unit_wt
    if base_shear_vel is None:
        base_shear_vel = layers[-1].shear_vel
    base_soil = pysra.site.SoilType("base", base_unit_wt, None, base_xi)
    lay = pysra.site.Layer(base_soil, 0.1, base_shear_vel)
    layers.append(lay)

    profile = pysra.site.Profile(layers, wt_depth=sp.gwl)
    return profile


def compute_pysra_strain_compatible_profile(
    soil_profile, in_sig, d_inc=None, cut_time=None, target_height=1.0
):
    import pysra

    m = pysra.motion.TimeSeriesMotion(
        filename=in_sig.label,
        description=None,
        time_step=in_sig.dt,
        accels=in_sig.values / 9.8,
    )
    if d_inc is None:
        d_inc = 1.0 * np.ones(soil_profile.n_layers)
    profile = sm_profile_to_pysra(
        soil_profile, target_height=target_height, d_inc=d_inc
    )

    layers = []
    calc = pysra.propagation.EquivalentLinearCalculator()
    calc(m, profile, profile.location("outcrop", depth=soil_profile.height))
    if cut_time is not None:
        i_cut = np.where(m.times > cut_time)[0][0]

        outs = []
        for i, depth in enumerate(profile.depth):
            outs.append(
                pysra.output.StrainTSOutput(
                    pysra.output.OutputLocation("within", depth=depth), in_percent=False
                )
            )

        outputs = pysra.output.OutputCollection(*outs)
        outputs(calc)
        for i, depth in enumerate(profile.depth):
            max_strain = np.max(np.abs(outputs[i].values[:i_cut]))
            # set new strain comp
            org_layer = profile.location("outcrop", depth=depth).layer
            org_layer.strain = calc.strain_ratio * max_strain

    for depth in profile.depth:
        org_layer = profile.location("outcrop", depth=depth).layer
        shear_vel0 = org_layer.initial_shear_vel
        shear_vel = org_layer.shear_vel
        print("vs_ratio: ", shear_vel / shear_vel0)
        unit_wt = org_layer.unit_wt
        damping = org_layer.damping
        slice_thickness = org_layer.thickness
        pysra_sl = pysra.site.SoilType("soil", unit_wt, None, damping)
        lay = pysra.site.Layer(pysra_sl, slice_thickness, shear_vel)
        layers.append(lay)

    strain_comp_profile = pysra.site.Profile(layers, wt_depth=soil_profile.gwl)

    return strain_comp_profile


def update_pysra_profile(pysra_profile, depths, xis=None, shear_vels=None):
    import pysra

    layers = []
    for depth in pysra_profile.depth:
        try:
            indy = depths.index(depth)
            if xis is not None:
                damping = xis[indy]
            else:
                damping = pysra_profile.location("outcrop", depth=depth).layer.damping
            if shear_vels is not None:
                shear_vel0 = shear_vels[indy]
            else:
                shear_vel0 = pysra_profile.location(
                    "outcrop", depth=depth
                ).layer.initial_shear_vel
        except ValueError:
            shear_vel0 = pysra_profile.location(
                "outcrop", depth=depth
            ).layer.initial_shear_vel
            damping = pysra_profile.location("outcrop", depth=depth).layer.damping

        unit_wt = pysra_profile.location("outcrop", depth=depth).layer.unit_wt
        slice_thickness = pysra_profile.location("outcrop", depth=depth).layer.thickness

        pysra_sl = pysra.site.SoilType("soil", unit_wt, None, damping)
        lay = pysra.site.Layer(pysra_sl, slice_thickness, shear_vel0)
        layers.append(lay)

    new_profile = pysra.site.Profile(layers, wt_depth=pysra_profile.wt_depth)

    return new_profile


class PysraAnalysis(object):

    def __init__(
        self,
        soil_profile,
        asig,
        odepths,
        wave_field="outcrop",
        atype="eqlin",
        outs=None,
        trim=False,
        strain_ratio=None,
    ):

        import pysra

        pysra_profile = sm_profile_to_pysra(
            soil_profile, d_inc=[0.5] * soil_profile.n_layers
        )
        # Should be input in g
        pysra_m = pysra.motion.TimeSeriesMotion(
            asig.label, None, time_step=asig.dt, accels=asig.values / 9.8
        )
        kw = {}
        if strain_ratio is not None:
            kw["strain_ratio"] = strain_ratio
        if atype == "eqlin":
            calc = pysra.propagation.EquivalentLinearCalculator(**kw)
        elif atype == "fd":
            calc = pysra.propagation.FrequencyDependentEqlCalculator(
                use_smooth_spectrum=False, **kw
            )
        elif atype == "fdk":  # k=Kausel
            calc = pysra.propagation.FrequencyDependentEqlCalculator(
                use_smooth_spectrum=True, **kw
            )
        elif atype == "linear":
            calc = pysra.propagation.LinearElasticCalculator()
        else:
            raise ValueError(f'atype must: "eqlin", "fd", "fdk", "linear". Not {atype}')

        if outs is None:
            od = {"ACCX": [], "STRS": [], "TAU": []}
        else:
            od = {}
            for item in outs:
                od[item] = []
        out_holder = []
        for i, depth in enumerate(odepths):
            if "ACCX" in od:
                od["ACCX"].append(len(out_holder))
                out_holder.append(
                    pysra.output.AccelerationTSOutput(
                        pysra.output.OutputLocation("within", depth=depth)
                    )
                )
            if "ACCXup" in od:
                od["ACCXup"].append(len(out_holder))
                out_holder.append(
                    pysra.output.AccelerationTSOutput(
                        pysra.output.OutputLocation("incoming_only", depth=depth)
                    )
                )
            if "STRS" in od:
                od["STRS"].append(len(out_holder))
                out_holder.append(
                    pysra.output.StrainTSOutput(
                        pysra.output.OutputLocation("within", depth=depth),
                        in_percent=False,
                    )
                )
            if "TAU" in od:
                od["TAU"].append(len(out_holder))
                out_holder.append(
                    pysra.output.StressTSOutput(
                        pysra.output.OutputLocation("within", depth=depth),
                        normalized=False,
                    )
                )
        outputs = pysra.output.OutputCollection(out_holder)
        calc(
            pysra_m,
            pysra_profile,
            pysra_profile.location(wave_field, depth=soil_profile.height),
        )
        outputs(calc)
        if trim:
            n = asig.npts
        else:
            n = None
        out_series = {}
        for mtype in od:
            out_series[mtype] = []
            for i in range(len(od[mtype])):
                out_series[mtype].append(outputs[od[mtype][i]].values[:n])
            out_series[mtype] = np.array(out_series[mtype])
            if mtype in ["ACCX", "ACCXup"]:
                out_series[mtype] *= 9.8
        out_series["TIME"] = (
            np.arange(0, len(out_series[list(out_series)[0]][0])) * asig.dt
        )
        self.out_series = out_series
        self.pysra_profile = pysra_profile
        self.calc = calc


class PysraDeconvolutionAnalysis(object):

    def __init__(
        self,
        soil_profile,
        asig,
        odepths,
        wave_field="outcrop",
        atype="eqlin",
        outs=None,
        trim=False,
        strain_ratio=None,
    ):

        import pysra

        pysra_profile = lq.sra.sm_profile_to_pysra(
            soil_profile, d_inc=[0.5] * soil_profile.n_layers
        )
        # Should be input in g
        pysra_m = pysra.motion.TimeSeriesMotion(
            asig.label, None, time_step=asig.dt, accels=asig.values / 9.8
        )
        kw = {}
        if strain_ratio is not None:
            kw["strain_ratio"] = strain_ratio
        if atype == "eqlin":
            calc = pysra.propagation.EquivalentLinearCalculator(**kw)
        elif atype == "fd":
            calc = pysra.propagation.FrequencyDependentEqlCalculator(
                use_smooth_spectrum=False, **kw
            )
        elif atype == "fdk":
            calc = pysra.propagation.FrequencyDependentEqlCalculator(
                use_smooth_spectrum=True, **kw
            )
        elif atype == "linear":
            calc = pysra.propagation.LinearElasticCalculator()
        else:
            raise ValueError(f'atype must: "eqlin", "fd", "linear". Not {atype}')

        if outs is None:
            od = {"ACCX": [], "STRS": [], "TAU": []}
        else:
            od = {}
            for item in outs:
                od[item] = []
        out_holder = []
        for i, depth in enumerate(odepths):
            if "ACCX" in od:
                od["ACCX"].append(len(out_holder))
                out_holder.append(
                    pysra.output.AccelerationTSOutput(
                        pysra.output.OutputLocation("within", depth=depth)
                    )
                )
            if "ACCXup" in od:
                od["ACCXup"].append(len(out_holder))
                out_holder.append(
                    pysra.output.AccelerationTSOutput(
                        pysra.output.OutputLocation("incoming_only", depth=depth)
                    )
                )
            if "STRS" in od:
                od["STRS"].append(len(out_holder))
                out_holder.append(
                    pysra.output.StrainTSOutput(
                        pysra.output.OutputLocation("within", depth=depth),
                        in_percent=False,
                    )
                )
            if "TAU" in od:
                od["TAU"].append(len(out_holder))
                out_holder.append(
                    pysra.output.StressTSOutput(
                        pysra.output.OutputLocation("within", depth=depth),
                        normalized=False,
                    )
                )
        outputs = pysra.output.OutputCollection(out_holder)
        calc(pysra_m, pysra_profile, pysra_profile.location(wave_field, depth=0))
        outputs(calc)
        if trim:
            n = asig.npts
        else:
            n = None
        out_series = {}
        for mtype in od:
            out_series[mtype] = []
            for i in range(len(od[mtype])):
                out_series[mtype].append(outputs[od[mtype][i]].values[:n])
            out_series[mtype] = np.array(out_series[mtype])
            if mtype in ["ACCX", "ACCXup"]:
                out_series[mtype] *= 9.8
        out_series["TIME"] = (
            np.arange(0, len(out_series[list(out_series)[0]][0])) * asig.dt
        )
        self.out_series = out_series
        self.pysra_profile = pysra_profile
        self.calc = calc


def run_pysra(
    soil_profile,
    asig,
    odepths,
    wave_field="outcrop",
    atype="eqlin",
    outs=None,
    strain_ratio=None,
):
    """

    :param soil_profile:
    :param asig:
    :param odepths:
    :param wave_field: str
        either - 'outcrop', 'within', 'incoming_only'
    :param atype:
    :param outs:
    :return:
    """
    pa = PysraAnalysis(
        soil_profile,
        asig,
        odepths,
        wave_field=wave_field,
        atype=atype,
        outs=outs,
        strain_ratio=strain_ratio,
    )
    return pa.out_series
