from liquepy import functions


def test_determine_t_liq_index():
    ru = [0., 0.1, 0.4, 0.5, 0.8]
    assert functions.determine_t_liq_index(ru, ru_limit=0.5) == 3
    ru = [0., 0.1, 0.4, 0.59, 0.8]
    assert functions.determine_t_liq_index(ru, ru_limit=0.5) == 3

