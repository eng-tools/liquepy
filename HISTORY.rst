=======
History
=======

0.4.12 (2018-03-14)
-------------------

* Added calculation of dissipated energy and cumulative absolute change in shear stress of element tests


0.4.11 (2018-03-14)
-------------------

* Added plotting functions for CPT
* Cleaned up names of input motion saving functions, and order of args

0.4.8 - 0.4.10 (2018-03-08)
---------------------------

* Updated docstrings, readme file
* Fixed number of columns to load on CPT to be 0-3

0.4.7 (2018-02-28)
------------------

* `run_bi2014` fixed bug where water unit weight was 10 times too big

0.4.5 (2018-02-27)
------------------

* `BoulangerIdriss2014` unit weight calculation now uses the specific weight of water a gravity=9.8
* `BoulangerIdriss2014` gwl now is the ground water level during the earthquake, while cpt_gwl is the gwl at cpt measure
* Added `BoulangerIdriss2014CPT` which performs `BoulangerIdriss2014` but takes CPT as an input


0.4.4 (2018-02-27)
------------------

* Changes to `BoulangerIdriss2014`:

    * Added `big_q` as a property
    * Modified unit weight calculation (minimum changed from 15kN/m3 to 14.715kN/m3 (lowest value in original study by Roberston (2010)), Added maximum unit weight 19.62kN/m3 (maximum value in study),
    * Changed atmospheric pressure from 100kPa to 101kPa (also added as an optional input)
    * Input `magnitude` -> `m_w`
    * Added support for calculation of unit weight using specific weight
    * Optional input `s_g` to override specific weight of 2.65
    * Optional input `s_g_water` to override specific weight of water
