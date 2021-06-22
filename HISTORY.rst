=======
History
=======

Pre-release
-----------

0.6.23 (2021-06-23)
-------------------
* Fixed issue with conversion of NZGD file using `convert_raw02_xlsx_or_raw_03_xls` function
* Added support for the PM4Silt constitutive model
* Added missing `e_min` and `e_max` parameters to optional parameters for PM4Sand model

0.6.22 (2021-04-19)
-------------------
* Soils now define water mass density `wmd` and `liq_sg` instead of `liq_mass_density` to allow for not water liquids.
* Can set min shear wave velocity in sra using `vs_min`.
* Default k0 conditions now calculated using poisson ratio rather than friction angle.
* Added the ELSA (Equivalent Linear Stockwell Analysis) method for performing site response analysis using time dependent transfer functions.

0.6.18 (2020-12-08)
-------------------
* Can get set liquefaction criteria for element test using double or single  amplitude strain
* Can deal with loading a cpt file when the gwl level is not set
* Added the ManzariDafaliasModel soil model
* `big_f` added as a parameter in Boulanger and Idriss CPT analysis
* fixed function to allow crr_n15 to be computed from q_c1ncs according to BI2014

0.6.15 (2020-07-22)
-------------------
* Added back-calc of n1_60 from crr_m7p5 using bi2014 `calc_n1_60cs_from_crr_m7p5`
* Added `sra.PysraAnalysis` and `sra.run_pysra` as ways to run an equivalent linear
  analysis using the PySRA package
* Can set point of liquefaction in element test using single aplitude and double amplitude strain criteria

0.6.14 (2020-04-29)
-------------------
* Added theory one and two layer transfer functions
* Minor change to LSN calculation to not use the average depth between points
* Added `run_pysra` to run a site response analysis and obtain ACCX, TAU, and STRS
* Updated logic for handling FLAC fis file inputs

0.6.13 (2020-02-24)
-------------------
* Fixed issue with creating ESP when cpt depths are irregular

0.6.12 (2020-02-21)
-------------------
* Fixed issue with fitting 3- and 5-layer ESPs now uses median to find value of least deviations

0.6.11 (2020-02-10)
-------------------
* Added MRD curves for clay
* Fixed issue with fitting 5-layer equivalent soil profile (ESP)

0.6.10 (2020-01-31)
-------------------
* Added .p_atm to Boulanger and Idriss (2014) method
* Improved plotting colors for i_c
* Changed liquefaction factor of safety for high plasticity soil to be 2.25
* Added method for detecting cyclic peaks using peak energy
* removed examples folder from pip install of package
* Added more FLAC based functions
* Significant increase in package loading speeding
* Added Equivalent Soil Profile method from Millen et al. (2020)

0.6.9 (2019-12-6)
-----------------
* Added nzgd CPT file converter
* Corrected error in `est_shear_vel_hegazy_and_mayne_2006`
* Change modules (`sra`, `fig`, `spatial`) to be optional due to heavy dependencies. Only imported at base level if all
  dependencies are satisfied.

0.6.7 (2019-11-12)
------------------
* Updated docs
* Updated names for PM4Sand FLAC
* Refactored numerical models to be in generic num.models.py file and the inherited in application specific files
  (flac, o3), where app specific details are added.
* Added permeability estimate `est_permeability_robertson_and_cabal_2012`
* Added more tests to increase coverage
* Can change empirical fitting parameter 'c0' in calculation of CRR_m7.5 from BI2016
* Added more logic for FLACSoil model to output to fis file
* Fixed issue where could not calculate transfer function in site response analysis
* `determine_t_liq_index` now uses ru greater than or equal to limit, rather than just greater than

0.6.3 (2019-10-31)
------------------
* Added missing `import eqsig` at top of lq.element.assess.py file - used in
  `calc_stored_energy_abs_incs_fd_peaks_and_indices` function
* Finalised NSES example

0.6.2 (2019-10-30)
------------------

* Fixed issue with Soil stack not working for overriding soil properties
* Added h_po to PM4Sand obj
* Added estimation of cumulative absolute (change in) strain energy (CASE) in 1D soil profile using the nodal surface
  energy spectrum (NSES) from Millen et al. (2019) as a function `lq.trigger.nses.est_case_1d_millen_et_al_2019()`
* Added example of estimation of CASE compared to linear analysis using pysra package
* Added calculation of cumulative absolute (change in) strain energy (CASE) of element test from Millen et al. (2019)
  as function `calc_case_peaks_and_indices_fd`

0.6.1 (2019-10-7)
-----------------

* Incremented numerical model inputs to make use of latest sfsimodels version

0.6.0 (2019-08-23)
------------------

* If CPT file contains key word 'Pre-Drill:' then CPT loads with pre-drilled depth removed
* Can set depth limit to LDI calculation
* Can set relative density limits to zhang_2002 relative density calculation
* Fixed issue with I_c color map for plotting
* Switched
* [Not backward compatible] Switched Zhang et al shear and volumetric strain calculations to proper naming convention
  and changed to return strain as decimal, removed old functions

0.5.7 (2019-07-05)
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
* Modified unit weight calculation (minimum changed from 15kN/m3 to 14.715kN/m3 (lowest value in original study
  by Roberston (2010)), Added maximum unit weight 19.62kN/m3 (maximum value in study),
* Changed atmospheric pressure from 100kPa to 101kPa (also added as an optional input)
* Input `magnitude` -> `m_w`
* Added support for calculation of unit weight using specific weight
* Optional input `s_g` to override specific weight of 2.65
* Optional input `s_g_water` to override specific weight of water
