from tests.conftest import TEST_DATA_DIR
from liquepy.field import nzgd_cpt_converter as nzgd
from liquepy.field import load_mpa_cpt_file


def test_nzgd_convert_file():
    ffp = TEST_DATA_DIR + 'CPT_111111_Raw01.xlsx'
    result = nzgd.convert_file(ffp, TEST_DATA_DIR)
    assert result == 'convert_raw01_w_underscores'
    cpt = load_mpa_cpt_file(TEST_DATA_DIR + 'CPT_111111.csv')
    assert len(cpt.q_c) > 3
    # print(result)


if __name__ == '__main__':
    test_nzgd_convert_file()
