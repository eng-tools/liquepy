
# Note that Been and Jefferies (2015) provides critical state lines for many soils but assumes ksi=1
def get_soil_params_dict():
    return {

        'toyoura_sand':
            {
                'e_min': 0.616,  # From Cubrinovski and Ishihara (2002)
                'e_max': 0.988,  # From Cubrinovski and Ishihara (2002)
                'specific_gravity': 2.68,  # Fioravante:2015ea
                'e_0': 0.934,  # fitted by Li and Wang (1998) data from Verdugo and Ishihara (1996)
                'lambda_c': 0.019,  # fitted by Li and Wang (1998) data from Verdugo and Ishihara (1996)
                'ksi': 0.7,  # fitted by Li and Wang (1998) data from Verdugo and Ishihara (1996)
                'd_50': 0.53,  # Fioravante:2015ea
                'm_c': 1.24,  # Been and Jefferies (2015) Table 2.1 (Golder Project files)
            },
        'ticino_sand':
            {
                'e_min': 0.574,  # Fioravante:2015ea
                'e_max': 0.923,  # Fioravante:2015ea
                'specific_gravity': 2.65,  # Fioravante:2015ea
                'e_0': 0.923,  # Fioravante:2015ea
                'lambda_c': 0.046,  # Fioravante:2015ea
                'ksi': 0.5,  # Fioravante:2015ea
                'd_50': 0.18,  # averaged from Fioravante:2015ea
                'm_c': 1.14,  # TRISEE tests
            },
        'nevada_sand':
            {
                'e_min': 0.516,  # Arulmoli et al. (1992)
                'e_max': 0.894,  # Arulmoli et al. (1992)
                'specific_gravity': 2.68,  # Arulmoli et al. (1992) calculated from numbers in Leshchinsky:2018ca
                'e_0': 0.843,  # Ling:2006jl
                'lambda_c': 0.0287,  # Ling:2006jl
                'ksi': 0.7,  # Ling:2006jl
                'm_c': 1.09,  # gives phi_cv = 33.0
             },
        'fuji_river_sand_loose':
            {
                'e_min': 0,
                'e_max': 0,
                'specific_gravity': 2.65,
                'e_0': 0.81,  # Ling:2006jl
                'lambda_c': 0.033,  # Ling:2006jl
                'ksi': 0.7,  # Ling:2006jl
        },
        'fuji_river_sand_dense':
            {
                'e_min': 0.5,  # default PM4Sand input
                'e_max': 0.8,  # default PM4Sand input
                'specific_gravity': 2.65,
                'e_0': 0.684,  # Ling:2006jl
                'lambda_c': 0.033,  # Ling:2006jl
                'ksi': 0.7,  # Ling:2006jl
            },
        'longstone_sand':
            {
                'a_mod': 1000,  # Tsomokos:2010cd  # TODO: confirm
                'e_min': 0.614,  # Tsomokos:2010cd
                'e_max': 0.995,  # Tsomokos:2010cd
                'specific_gravity': 2.64, # Tsomokos:2010cd
                'm_c': 1.09,  # Tsomokos:2010cd failure line at 33deg
                'd_50': 0.15, # mm, Gazetas DOI: 10.1007/s11012-014-9997-7
                'c_u': 1.4,  # d60/d10 Gazetas DOI: 10.1007/s11012-014-9997-7
                'e_0': 0.93,  # Back-calculate phi-peak from Gazetas (2014) using Manzari model
                'lambda_c': 0.1,  # Back-calculate phi-peak from Gazetas (2014) using Manzari model
                'ksi': 0.7,
            },
        'generic_sand':
            {
                'm_c': 1.25,  # Been and Jefferies (2015) range from 1.2 - 1.35
                'e_0': 0.934,  # From Toyoura Sand
                'lambda_c': 0.019,  # From Toyoura Sand
                'ksi': 0.7,  # From Toyoura Sand
                'd_50': 0.53,  # From Toyoura Sand
            },
        'generic_silt':
            {'m_c': 1.45,  # Been and Jefferies (2015) range from 1.3 - 1.6
             },
    }


def calc_g0_mod_and_a_for_nevada_sand_from_e_curr_arulmoli_et_al_1992(e_curr):
    # note that F(e) = 1 / (0.3 + 0.7 * e **2) is from Hardin (1978) - see Menq (2003) Eq. 4.3
    fe = 1 / (0.3 + 0.7 * e_curr)
    cf = 0.77  # Arulmoli et al. (1992)
    return 625 * cf * fe, 0.5


def calc_g0_mod_and_a_for_ottawa_sand_from_e_curr_hardin_and_richart_1963(e_curr):
    # from Menq (2003) table 4.1
    fe = (2.17 - e_curr) ** 2 / (1 + e_curr)
    a = 0.5
    a_g = 7000
    return a_g / (100 ** a) * fe, a


def calc_g0_mod_and_a_for_angular_crushed_quartz_from_e_curr_hardin_and_richart_1963(e_curr):
    # from Menq (2003) table 4.1
    fe = (2.97 - e_curr) ** 2 / (1 + e_curr)
    a = 0.5
    a_g = 3300
    return a_g / (100 ** a) * fe, a


def calc_g0_mod_and_a_for_clean_sand_from_e_curr_iwasaki_et_al_1978(e_curr):
    # from Menq (2003) table 4.1
    fe = (2.17 - e_curr) ** 2 / (1 + e_curr)
    a = 0.38
    a_g = 9000
    return a_g / (100 ** a) * fe, a


def calc_g0_mod_and_a_for_toyoura_sand_from_e_curr_kokusho_1980(e_curr):
    # from Menq (2003) table 4.1
    fe = (2.17 - e_curr) ** 2 / (1 + e_curr)
    a = 0.5
    a_g = 8400
    return a_g / (100 ** a) * fe, a


def calc_g0_mod_and_a_for_clean_sand_from_e_curr_yu_and_richart_1984(e_curr):
    # from Menq (2003) table 4.1
    fe = (2.17 - e_curr) ** 2 / (1 + e_curr)
    a = 0.5
    a_g = 7000
    return a_g / (100 ** a) * fe, a


def calc_g0_mod_and_a_for_crushed_rock_from_e_curr_kokusho_and_esashi_1981(e_curr):
    # from Menq (2003) table 4.1
    fe = (2.17 - e_curr) ** 2 / (1 + e_curr)
    a = 0.55
    a_g = 13000
    return a_g / (100 ** a) * fe, a


def calc_g0_mod_and_a_for_round_gravel_from_e_curr_kokusho_and_esashi_1981(e_curr):
    # from Menq (2003) table 4.1
    fe = (2.17 - e_curr) ** 2 / (1 + e_curr)
    a = 0.6
    a_g = 8400
    return a_g / (100 ** a) * fe, a


def calc_g0_mod_and_a_for_gravel_from_e_curr_tanaka_et_al_1987(e_curr):
    # from Menq (2003) table 4.1
    fe = (2.17 - e_curr) ** 2 / (1 + e_curr)
    a = 0.6
    a_g = 3080
    return a_g / (100 ** a) * fe, a


def calc_g0_mod_and_a_for_gravel_from_e_curr_goto_et_al(e_curr):
    # from Menq (2003) table 4.1
    fe = (2.17 - e_curr) ** 2 / (1 + e_curr)
    a = 0.85
    a_g = 1200
    return a_g / (100 ** a) * fe, a

