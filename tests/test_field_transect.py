import liquepy as lq
from tests.conftest import TEST_DATA_DIR
import sfsimodels as sm


def test_load_cpt_from_loc():
    cpt = lq.field.load_mpa_cpt_file(TEST_DATA_DIR + "standard_1.csv")
    loc = lq.field.transect.Loc(cpt=cpt)
    cpta = loc.cpt
    assert len(cpta.q_c) > 100


def test_loc_to_dict():
    cpt = lq.field.load_mpa_cpt_file(TEST_DATA_DIR + "standard_1.csv")
    sp = sm.SoilProfile()
    loc = lq.field.transect.Loc(cpt=cpt)
    loc.soil_profile = sp
    pdict = loc.to_dict()
    assert 'id' in pdict
    assert 'soil_profile_id' in pdict


if __name__ == '__main__':
    test_loc_to_dict()
