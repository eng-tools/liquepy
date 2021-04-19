from collections import OrderedDict
import numpy as np
import itertools

FS_VMIN = 0.5
FS_VMAX = 1.75
FS_VMAX_W_NONLIQ = 2.5
FS_LOW_to_0p75 = (0.75, 0, 0)
FS_0p75_to_1p0 = (0.95, 0, 0)
FS_1p0_to_1p25 = (0.9, 0.4, 0.15)
FS_1p25_to_1p5 = (1, 0.65, 0.25)
FS_1p5_to_HIGH = (0.1, 0.6, 0.1)
FS_NON_LIQ = (0.4, 0.4, 0.4)
FS_COLORS = OrderedDict([
    ('FS_LOW_to_0p75', FS_LOW_to_0p75),  # dark red
    ('FS_0p75_to_1p0', FS_0p75_to_1p0),  # red
    ('FS_1p0_to_1p25', FS_1p0_to_1p25),  # dark orange
    ('FS_1p25_to_1p5', FS_1p25_to_1p5),  # orange
    ('FS_1p5_to_HIGH', FS_1p5_to_HIGH),  # green
])
FS_COLORS_W_NONLIQ = OrderedDict([
    ('FS_LOW_to_0p75', FS_LOW_to_0p75),  # dark red
    ('FS_0p75_to_1p0', FS_0p75_to_1p0),  # red
    ('FS_1p0_to_1p25', FS_1p0_to_1p25),  # dark orange
    ('FS_1p25_to_1p5', FS_1p25_to_1p5),  # orange
    ('FS_1p5_to_HIGH', FS_1p5_to_HIGH),  # green
    ('FS_1p5_to_HIGH2', FS_1p5_to_HIGH),  # green
    ('FS_1p5_to_HIGH2', FS_1p5_to_HIGH),  # green
    ('FS_NON_LIQ', FS_NON_LIQ),  # green
])


IC_VMIN = 0.0
IC_VMAX = 4.0
IC_LIMITS = [0, 1.3, 1.8, 2.1, 2.6, 4]
IC_GRAVEL = (0.98, 0.88, 0.5)
IC_GRAVEL_dark = (0.88, 0.76, 0.2)
IC_GRAVEL_light = (1.0, 1.0, 0.8)
IC_CS = (0.48, 0.98, 0.8)
IC_CS_light = (0.6, 1.0, 0.8)
IC_SwFC = (0.57, 0.97, 1.0)
IC_SwFC_dark = (0.38, 0.66, 0.8)
IC_SwFC_light = (0.8, 1.0, 1.0)
IC_NP_silt = (0.97, 0.66, 0.73)
IC_NP_silt_light = (0.97, 0.72, 0.72)
IC_P_silt = (0.65, 0.65, 0.65)
IC_P_silt_dark = (0.5, 0.5, 0.5)
IC_P_silt_light = (0.8, 0.8, 0.8)
IC_COLORS_w_CONTRAST = OrderedDict([
    ('Gravel', IC_GRAVEL_dark),  # yellow
    ('Clean sand', IC_CS_light),  # light green
    ('Sand with fines', IC_SwFC_dark),  # light blue
    ('Non-plastic silt', IC_NP_silt_light),  # light red
    ('Plastic silt', IC_P_silt_dark),  # light grey
])
IC_COLORS = OrderedDict([
    ('Gravel', IC_GRAVEL),  # yellow
    ('Clean sand', IC_CS),  # light green
    ('Sand with fines', IC_SwFC),  # light blue
    ('Non-plastic silt', IC_NP_silt),  # light red
    ('Plastic silt', IC_P_silt),  # light grey
])

# Colors
IC_ALT_GRAVEL = [154 / 255, 176 / 255, 187 / 255]
IC_ALT_CS = [255 / 255, 246 / 255, 187 / 255]
IC_ALT_SwFC = [218 / 255, 196 / 255, 77 / 255]
IC_ALT_NP_silt = [223 / 255, 220 / 255, 151 / 255]
IC_ALT_P_silt = [196 / 255, 148 / 255, 73 / 255]
IC_ALT_COLORS = OrderedDict([
    ('Gravel', IC_ALT_GRAVEL),  # yellow
    ('Clean sand', IC_ALT_CS),  # light green
    ('Sand with fines', IC_ALT_SwFC),  # light blue
    ('Non-plastic silt', IC_ALT_NP_silt),  # light red
    ('Plastic silt', IC_ALT_P_silt),  # light grey
])

FS_CMAP = None
FS_CMAP_W_NONLIQ = None
IC_CMAP = None
IC_CMAP_w_CONTRAST = None


def build_cmaps():
    from matplotlib.colors import LinearSegmentedColormap  # Importing this is slow
    global FS_CMAP
    global FS_CMAP_W_NONLIQ
    global IC_CMAP
    global IC_CMAP_w_CONTRAST
    FS_CMAP = LinearSegmentedColormap.from_list('fs', [FS_COLORS[cname] for cname in FS_COLORS], N=5)
    FS_CMAP_W_NONLIQ = LinearSegmentedColormap.from_list('fs_w_nonliq', [FS_COLORS_W_NONLIQ[cname] for cname in FS_COLORS_W_NONLIQ], N=5)
    _incs = np.round(np.diff(IC_LIMITS) * 10).astype(int)
    _clist = list(itertools.chain(*[[IC_COLORS[cname]] * _incs[x] for x, cname in enumerate(IC_COLORS)]))
    IC_CMAP = LinearSegmentedColormap.from_list('mine', _clist, N=len(_clist))
    _clist = list(itertools.chain(*[[IC_COLORS_w_CONTRAST[cname]] * _incs[x] for x, cname in enumerate(IC_COLORS_w_CONTRAST)]))
    IC_CMAP_w_CONTRAST = LinearSegmentedColormap.from_list('mine', _clist, N=len(_clist))

def get_fs_cmap():
    from matplotlib.colors import LinearSegmentedColormap  # Importing this is slow
    return LinearSegmentedColormap.from_list('fs', [FS_COLORS[cname] for cname in FS_COLORS], N=5)


def get_fs_w_nonliq_cmap():
    from matplotlib.colors import LinearSegmentedColormap  # Importing this is slow
    return LinearSegmentedColormap.from_list('fs_w_nonliq',
                                             [FS_COLORS_W_NONLIQ[cname] for cname in FS_COLORS_W_NONLIQ], N=5)

def get_ic_cmap():
    from matplotlib.colors import LinearSegmentedColormap  # Importing this is slow
    _incs = np.round(np.diff(IC_LIMITS) * 10).astype(int)
    _clist = list(itertools.chain(*[[IC_COLORS[cname]] * _incs[x] for x, cname in enumerate(IC_COLORS)]))
    return LinearSegmentedColormap.from_list('mine', _clist, N=len(_clist))


def get_ic_w_contrast_cmap():
    from matplotlib.colors import LinearSegmentedColormap  # Importing this is slow
    _incs = np.round(np.diff(IC_LIMITS) * 10).astype(int)
    _clist = list(itertools.chain(*[[IC_COLORS_w_CONTRAST[cname]] * _incs[x]
                                    for x, cname in enumerate(IC_COLORS_w_CONTRAST)]))
    return LinearSegmentedColormap.from_list('mine', _clist, N=len(_clist))


def add_ic_colours(subplot, col_dict='alt', ic_limits=None, add_legend=False):
    from matplotlib import patches as mpatches
    if ic_limits is None:
        ic_limits = IC_LIMITS

    if col_dict == 'alt':
        dd = IC_ALT_COLORS
    elif col_dict == 'main':
        dd = IC_COLORS
    elif col_dict == 'con':
        dd = IC_COLORS_w_CONTRAST
    else:
        dd = col_dict
    colours = [dd[x] for x in dd]
    colours.append((1, 1, 1))
    patches = []
    for col in dd:
        patches.append(mpatches.Patch(color=dd[col], label=col))
    if add_legend:
        subplot.legend(patches + [(patches)], list(dd), loc=0, prop={'size': 8})

    for i in range(len(ic_limits) - 1):
        subplot.axvspan(ic_limits[i], ic_limits[i + 1], alpha=1.0, color=colours[i])
    return patches

def make_ic_plot(subplot):
    add_ic_colours(subplot)
    subplot.set_xlim([0, 3])


def get_ic_legend_patches(col_dict):
    from matplotlib import patches as mpatches
    if col_dict == 'alt':
        dd = IC_ALT_COLORS
    elif col_dict == 'main':
        dd = IC_COLORS
    elif col_dict == 'con':
        dd = IC_COLORS_w_CONTRAST
    else:
        dd = col_dict
    legend_elements = []
    patches = []
    for col in dd:
        patches.append(mpatches.Patch(color=dd[col], label=col))
    return patches


def make_factor_of_safety_plot(subplot, w_ic_lim=False):
    # add the Fs = 1 line
    subplot.axvspan(0., 0.75, alpha=0.5, color=FS_LOW_to_0p75)
    subplot.axvspan(0.75, 1.0, alpha=0.5, color=FS_0p75_to_1p0)
    subplot.axvspan(1.0, 1.25, alpha=0.5, color=FS_1p0_to_1p25)
    subplot.axvspan(1.25, 1.5, alpha=0.5, color=FS_1p25_to_1p5)
    subplot.axvspan(1.5, 2.0, alpha=0.5, color=FS_1p5_to_HIGH)
    subplot.axvspan(2.0, 2.3, alpha=0.5, color=FS_NON_LIQ)
    subplot.set_xlim([0, 2.3])


def make_cpt_plots(sps, cpt, c="gray", x_origin=True, y_origin=True):

    sps[0].plot(cpt.q_c, cpt.depth, lw=1, c=c)
    sps[1].plot(cpt.f_s, cpt.depth, lw=1, c=c)
    sps[2].plot(cpt.u_2, cpt.depth, lw=1, c=c)
    sps[2].axhline(cpt.gwl, c=c, ls="--", lw=0.7)

    # Prepare y-axis
    if y_origin:
        ylim = sps[0].get_ylim()
        sps[0].set_ylim([0, ylim[1]])

    # Prepare x-axis
    if x_origin:
        xlim = sps[0].get_xlim()
        sps[0].set_xlim([0, xlim[1]])
        xlim = sps[1].get_xlim()
        sps[1].set_xlim([0, xlim[1]])
        xlim = sps[2].get_xlim()
        sps[2].set_xlim([0, xlim[1]])

    sps[0].set_ylabel("Depth [m]")
    sps[0].set_xlabel("q_c [kPa]")
    sps[1].set_xlabel("f_s [kPa]")
    sps[2].set_xlabel("u_2 [kPa]")


def make_bi2014_outputs_plot(sps, bi2014, crr_cap=0.6):

    sps[0].plot(bi2014.pore_pressure, bi2014.depth, lw=1, c="b", label="Pore pressure")
    sps[0].plot(bi2014.sigma_v, bi2014.depth, lw=1, c="r", label="$\sigma_{v}$")
    sps[0].plot(bi2014.sigma_veff, bi2014.depth, lw=1, c="gray", label="$sigma_{v,eff}$")
    sps[0].legend(prop={"size": 8}, loc="lower left")
    sps[0].set_xlabel("Stress [kPa]")

    sps[1].plot(bi2014.i_c, bi2014.depth, "o", lw=1, c="k", alpha=0.5, ms=2)
    make_ic_plot(sps[1])
    # sps[1].legend()
    sps[1].set_xlabel("$I_c$")

    sps[2].plot(np.clip(bi2014.crr_m7p5, None, crr_cap), bi2014.depth, "o", lw=1, c="k", alpha=0.5, ms=3)
    sps[2].set_xlabel("$CRR_{m7.5}$")
    sps[2].set_xlim([0, 1.])

    sps[3].plot(bi2014.factor_of_safety, bi2014.depth, "o", lw=1, c="k", alpha=0.5, ms=3)
    make_factor_of_safety_plot(sps[3])
    sps[3].axhline(bi2014.gwl, c="k", ls="--", lw=0.7)
    sps[3].set_xlabel("Factor of Safety")

    # Prepare y-axis
    ylim = sps[0].get_ylim()
    sps[0].set_ylim([0, ylim[1]])
    sps[0].invert_yaxis()

    # Prepare x-axis
    xlim = sps[0].get_xlim()
    sps[0].set_xlim([0, xlim[1]])
    xlim = sps[1].get_xlim()
    sps[1].set_xlim([0, xlim[1]])
    xlim = sps[2].get_xlim()
    sps[2].set_xlim([0, xlim[1]])
    xlim = sps[3].get_xlim()
    sps[3].set_xlim([0, xlim[1]])

    sps[0].set_ylabel("Depth [m]")
