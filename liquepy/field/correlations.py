

def calc_vs_profile_mcgann_2015(cpt):
    """
    Computes the shear wave velocity profile according to :cite:`McGann:2015fd`

    Parameters
    ----------
    cpt: liquepy.field.CPT object

    Returns
    -------
    array_like
        Shear wave velocity profile corresponding to the depths in the CPT.
    """
    return 18.4 * cpt.q_c ** 0.144 * cpt.f_s ** 0.0832 * cpt.depth ** 0.278

