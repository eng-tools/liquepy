from liquepy.num import mrd_curves

# THIS FILE HAS DEPRECATED, now in liquepy.num.mrd_curves


def elastic(gamma):
    return mrd_curves.elastic(gamma)


def mohr_coloumb(gamma, elastic_mod, cu):
    return mrd_curves.mohr_coloumb(gamma, elastic_mod, cu)


def flac_default_curve(gamma, l_1, l_2):
    return mrd_curves.flac_default_curve(gamma, l_1, l_2)


def ishi_mod(gamma):
    return mrd_curves.ishi_mod(gamma)


def seed_and_sun_mod(gamma):
    return mrd_curves.seed_and_sun_mod(gamma)


def vardanega_2013_mod(gamma, i_p):
    return mrd_curves.vardanega_2013_mod(gamma, i_p)
