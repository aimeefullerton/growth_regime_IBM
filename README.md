
----------------------------------------------------------------------------------
Growth Regime IBM, VERSION 1.0.1

An individual based model to evaluate the contribution of seasonally warm habitat to Oncorhynchus mykiss growth, and used to generate an example presented in:


Armstrong, J.B., A.H. Fullerton, C.E. Jordan, J.L. Ebersole, J.R. Bellmore, I. Arismendi, B. Penaluna, and G.H. Reeves. The significance of warm habitat to the growth regime of coldwater fishes.

Adapted from:
1) Fullerton, A.H., B.J. Burke, J.J. Lawler, C.E. Torgersen, J.L. Ebersole, and S.G. Leibowitz. 2017. Simulated juvenile salmon growth and phenology respond to altered thermal regimes and stream network shape. Ecosphere 8(12):e02052. [Original model.]

2) Hawkins, B.L., A.H. Fulleton. B.L. Sanderson, and E.A. Steel. 2020. Individual-based simulations suggest mixed impacts of warmer temperatures and a non-native predator on Chinook salmon. Ecosphere 11(8):e03218. [Updated movement rules.]
----------------------------------------------------------------------------------

STEP I: Set up.

Follow these steps to download software, model input files, additional code, and libraries required to replicate our study.


1) Download R and RStudio.

2) Get model input files and additional code at https://github.com/aimeefullerton/growth_regime_IBM.

3) Either clone the repository or download as growth_regime_IBM.zip; when unzipped locally, this directory will serve as your R project directory.


4) Confirm that files are stored with the following structure within the growth_regime_IBM directory.

   code
	- growth_regime_IBM_v1.0.1.R - this is the main model script.
	- growth_regime_functions_v1.0.1.R - this script contains all the model functions and is sourced from the model script.
	- growth_regime_manuscript_figures.R - this script has code to create manuscript figures.
	- pre-calculate_growth.R - this script was used to create 'wt.growth.array.RData', which is also available in 'data.in'.
 
   data.in
	- wt.growth.array.RData [large file] Ð precalculated growth lookup array
	- thermal.regime.730ts.csv Ð thermal regime (temperature over time in outlet reach)
	- network-swh.ssn [folder containing network-specific files needed by the model]
 
5) Create a new R project in RStudio.

6) Install libraries in the setup section of growth_regime_IBM_v1.0.1.R.

----------------------------------------------------------------------------------


STEP II: Run simulations for four scenarios.

To run each simulation, you will need to update settings in the 'Scenarios & Startup' section

1) Run 'Baseline' scenario
      food.scenarios = "VariFood"
      mgmt.scenarios = "Base"

2) Run 'Divest in seasonally warm habitats' scenario
      food.scenarios = "VariFood"
      mgmt.scenarios = "DivestSWH"

3) Run 'Enhance perennially cold habitats' scenario
      food.scenarios = "VariFood"
      mgmt.scenarios = "EnhancePCH"

4) Run 'Constant Food' scenario
      ood.scenarios = "ConstFood"
      mgmt.scenarios = "Base"

After running all four scenarios, you should have new files in your data.out and plots folders.

   data.out
	- run.info.[scenario].txt - basic information about parameters used in this run
	- fa.[iter].steelhead.[scenario].RData - array of fish results for each time step
	- WT.[iter].steelhead.[scenario].RData - array of water temperature for each time step
	- production_[iter].csv Ð summary of fish production by habitat type

   plots
	- [iter].steelhead.[scenario].png - a quick diagnostic summary

   plots.ani - maps of each time step for one iteration



STEP III: Sensitivity analysis

To run simulations for sensitivity analysis, run the model script with parameters altered as needed (see "Scenarios.xlsx").

----------------------------------------------------------------------------------

STEP IV: Create figures and component panels for manuscript.

See code/growth_regime_manuscript_figures.R.




