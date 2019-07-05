=======
History
=======

0.5.7 (2019-xx-xx)
-------------------

* Build without history in setup.py file


0.5.6 (2019-07-05)
-------------------

* Added geographiclib to install reqs.

0.5.5 (2019-07-05)
-------------------

* Added support for setting weight of the pre-drilled depth for B&I2014
* Added support for calculating Liquefaction Potential Index (LPI) and Lateral Displacement Index (LDI)
* Fixed bug where calculation of shear strain from Zhang et al. (2004) used wrong value for interpolating between
different Dr lines, such that upper or lower limit were taken not interpolated value (resulted in slight change of
values for Bray foundation settlement method)
* Added new default cpt loader `load_mpa_cpt_file`, where all values are in MPa and metres and delimeter is ','
* Fixed bug where `big_Q` in B&I2014 triggering method was calculated with the `m` from Eq. 2.15b,
now calculated using the `n` from Robertson and Wride (1997)


0.5.4 (2019-05-28)
-------------------

* Added support for computing liquefaction triggering using BI2014 for a sfsimodels.SoilProfile object
* Fixed issue in sra module where depth increment was larger than layer, then failed
* Fixed bug where 'void_ratio' method for unit weight in BI2014 triggering added dry weight twice
* Added support for loading cpt files that do not have the u2 column
* Updated factor of safety colors for plotting, added color map and colors accessible as static values
* Added npts to bi2014 object
* Added colors and color map for I_c soil values
* Minor speed improvements to the B&I2014 triggering method

0.5.3 (2019-04-08)
-------------------

* Set base layer of eqlin site response to be elastic
* Refactored crr_m7p5 function from bi2014 method
* Switched sra commands to use latests sfsimodels package

0.5.1 (2019-03-29)
-------------------

* Added more correlations
* Changed all calculation functions to start with the prefix 'calc'
* Can set cut_time for obtaining strain compatible site response profile

0.5.0 (2019-03-14)
-------------------

* Changed order of inputs in ElementTest object!
* Changed ElementTest attributes (gamma -> strain, tau -> stress)

0.4.12 (2019-03-14)
-------------------

* Added calculation of dissipated energy and cumulative absolute change in shear stress of element tests


0.4.11 (2019-03-14)
-------------------

* Added plotting functions for CPT
* Cleaned up names of input motion saving functions, and order of args

0.4.8 - 0.4.10 (2019-03-08)
---------------------------

* Updated docstrings, readme file
* Fixed number of columns to load on CPT to be 0-3

0.4.7 (2019-02-28)
------------------

* `run_bi2014` fixed bug where water unit weight was 10 times too big

0.4.5 (2019-02-27)
------------------

* `BoulangerIdriss2014` unit weight calculation now uses the specific weight of water a gravity=9.8
* `BoulangerIdriss2014` gwl now is the ground water level during the earthquake, while cpt_gwl is the gwl at cpt measure
* Added `BoulangerIdriss2014CPT` which performs `BoulangerIdriss2014` but takes CPT as an input


0.4.4 (2019-02-27)
------------------

* Changes to `BoulangerIdriss2014`:

    * Added `big_q` as a property
    * Modified unit weight calculation (minimum changed from 15kN/m3 to 14.715kN/m3 (lowest value in original study by Roberston (2010)), Added maximum unit weight 19.62kN/m3 (maximum value in study),
    * Changed atmospheric pressure from 100kPa to 101kPa (also added as an optional input)
    * Input `magnitude` -> `m_w`
    * Added support for calculation of unit weight using specific weight
    * Optional input `s_g` to override specific weight of 2.65
    * Optional input `s_g_water` to override specific weight of water
