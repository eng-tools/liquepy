=======
History
=======

0.0.1 (2018-02-26)
------------------

* Changes to `BoulangerIdriss2014`:

    * Added `big_q` as a property
    * Modified unit weight calculation (minimum changed from 15kN/m3 to 14.715kN/m3 (lowest value in original study by
    Roberston (2010)), Added maximum unit weight 19.62kN/m3 (maximum value in study),
    * Changed atmospheric pressure from 100kPa to 101kPa (also added as an optional input)
    * Input `magnitude` -> `m_w`
    * Added support for calculation of unit weight using specific weight
    * Optional input `s_g` to override specific weight of 2.65
    * Optional input `s_g_water` to override specific weight of water


