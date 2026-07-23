The work is originated with goal of intelligently generating initial condition of 1-D channel modeling by Mascaret or Rezo, which is written into a file, in case of cold-start run.

The work is executed primarily by running the script "mascaret_initial_condition.py", which is dependent on api including source "sources/mascaret/API/f90/get_hydraulic_variable.f90". The script loads a json file which defines essential configurations. Keys in the json file are specified below. An example file "mascaret_initial_condition_configuration.json" is provided.

Prescribed initial condition can be provided at input model cross-section profiles within a file, which is generally written by some application program. The abscissa of cross-section profiles with prescribed initial condition must be identical to that of model input. An example file "initial_condition.txt" is provided.

If there is no prescribed initial condition or the prescribed initial condition is incomplete, the script fills initial level for input model cross-section profiles, dependending on the kernel used. The algorithm of filling initial level is specified below. Initial flowrate of input model cross-section profiles are filled with value of 0.0.

If kernel Rezo is used and type of prescribed level is used as boundary condition, the script executes the function "fill_initial_level" aiming to avoid super-critical flow and dry-zone. In case of a single reach model and prescribed level at both start and end boundaries, levels are linearly interpolated between prescribed levels at both boundaries, and then small and negative depths are fixed by input value. For ohter cases, firstly, steady kenel Sarap is executed to get level in reaches with prescribed-level boundary, defining limited flowrate by prescribed level at boundary and input limited Froude number, defining flow direction by bottom slope, and setting boundary as downstream end when zero bottom slope; then, for level in other reaches: if level of only one side has been filled, equality of conveyance is set at all cross-section profiles in a reach, and if level of both sides have been filled, use linear interpolation between two sides.
If kernel Mascaret is used or type of prescribed level is not used as boundary condition, level is fill as minimum of left and right bank level.

A pre-run could be performed, by utilizing kernel of Mascaret or Rezo, to get acceptablly quasi-steady results with steady prescribed boundary condition and initial condition. Total run period and tolerance can be defined for the pre-run. Initial condition of computational nodes in addition to input cross-section profiles is generated if pre-run is performed.

The work is on the basis of TELEMAC-MASCARET v8p5.




------------------------------------------------
------------------------------------------------
Specification of keys in configuration json file
------------------------------------------------
------------------------------------------------

-------------------------------------
-- Key "system_config/config_file" --
-------------------------------------
Full path of configuration file of TELEMAC-MASCARET system.
String type.


-------------------------------------
-- Key "system_config/config_name" --
-------------------------------------
Name of configuration of TELEMAC-MASCARET system.
String type.


-----------------------------------------------------
-- Key "models_configs/{model_name}/files/ic_cold" --
-----------------------------------------------------
Full path of prescribed initial condition file.
Use format 'RUBENS LIDOP' stated in 'mascaret.dico'.
The contents are cross-section abscissa, water level, and flowrate. They are 2D arrays, and must be 5 columns. Each element has exactly 15 characters(supporting number, letter, symbol, space), and contents other than number are parsed as 'NaN'.
The cross-section abscissa must be ordered exactly the same as in the cross-section file.
String type.


--------------------------------------------------------------
-- Key "models_configs/{model_name}/lmtd_dpth_fll_intl_lvl" --
--------------------------------------------------------------
When filling initial level for single reach and prescribed level at both boundaries, linear interpolation between level at both boundaries could cause small depth. This value is used for revising small depth.
Float type.


--------------------------------------------------------
-- Key "models_configs/{model_name}/Froude_lmtd_bdry" --
--------------------------------------------------------
Used for initial condition which could cause supercritical flow at boundary. Less Froude number leads to less limited flowrate for steady-run when filling initial level.
Float type.


---------------------------------------------------------
-- Key "models_configs/{model_name}/tlrc_fll_intl_lvl" --
---------------------------------------------------------
Tolerance used for computational iteration filling initial water level. "rltv" is relative one, and "abslt" is absolute one.
It is suggested that consider long-period averaged flow condition when setting absolute tolerance.
Less tolerance tends to output initial level with which model run executes normally, but more computational time and less chance for convergence.
Float type.


-----------------------------------------------------------
-- Key "models_configs/{model_name}/incrmt_fll_intl_lvl" --
-----------------------------------------------------------
Level increment used for computational iteration filling initial water level.
Less increment has greater chance to ensure convergence but more computational time.
Float type.


------------------------------------------
-- Key "run_configs/run_ic/run_to_stdy" --
------------------------------------------
When cold-start, whether pre-run to get initial conditions.
In pre-run, steady boundary conditions at start time are set.
When providing prescribed initial condition(key "models_configs/{model_name}/files/ic_cold") and performing pre-run, prescribed initial condition is used as internal boundary condition. However, if pre-run, it is strongly suggested not to provide prescribed initial condition, since most likely leads to computational error and failing of writing initial condition.
Bool type.


--------------------------------------
-- Key "run_configs/run_ic/run_prd" --
--------------------------------------
Maximum pre-run period in seconds.
Period which is too small may lead to failing of run-to-steady.
Float type.


---------------------------------------
-- Key "run_configs/run_ic/run_tlrc" --
---------------------------------------
Tolerance used for pre-run to get steady. "rltv_dpth" is relative one for water depth, "abslt_dpth" is absolute one for water depth, "rltv_flrt" is relative one for flowrate, and "abslt_flrt" is absolute one for flowrate.
It is suggested that consider long-period averaged flow condition when setting absolute tolerance.
Less torlerance leads to more accurate steady initial conditions, but more computational time and less chance of steady-state.
Float type.


------------------
-- Key "iprint" --
------------------
Flag of writing listing information into a file.
Bool type.


----------------------------------------------------------------------------
Other keys configuring model input files are consistent with api definition.
----------------------------------------------------------------------------