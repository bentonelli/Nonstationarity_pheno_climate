Code and data for Nonstationarity in Climate x Arrival Phenology 


Please direct all comments/questions to Ben Tonelli (btonelli [at] ucla.edu)

Note that fit model object is available under a separate DOI: doi.org/10.5281/zenodo.13256390 due to its large size (~4GB).

**Associated publications:**

**Repository structure:**

General overview, code:

The code here is used to estimate arrival dates for North American migratory birds using eBird data, and then compare those arrival dates to environmental covariates derived from Daymet. This code contains three sections, the first ("1_Prepare_data") is used to filter the global eBird database to North America during the study period and according to certain checklists standards. This data is then split by year, and then used  to estimate arrival dates for species in specific locations during each year. This section also includes code to download environmental covariates via Daymet using the associated R package "daymetr", which is then combined into a single dataset. Note that to run this section, users will need to request and download the eBird dataset (see link below), and the 2010 Global Multi-resolution Terrain Elevation Data at 30 arc-second resolution (see link in data section). Therefore, this section of code is not self-contained.

Section 2 ("2_Main_analysis") includes code to format the data to be analyzed via a single Bayesian hierarchical model (also included in this section). The code to run the main analysis is included in this section as well. This section of code is self-contained, as the data used in the analysis is provided in the data folder (see below).

Section 3 ("3_Output_analysis") contains code to get posterior estimates from the model, and to additionally to reproduce figures presented in the associated publication.

* `code/` - Code to preprocess data from eBird and Daymet (section 1), run models (section 2), analyze output (section 3).
 * `0_master_script.R` - Code to run all main script files sequentially via source command. Please read comments carefully before running, as some scripts have very long runtimes, and ideally should be run in parallel.
  * `1_Prepare_data` - The code in this section is used to preprocess data from eBird and Daymet. Please note that data needs to be dowloaded before this section is run (see overview section)
  	* `1a_awk_filter_ebird.rtf` - Filter eBird data using awk and split to year-specific files (i.e. full eBird dataset available for download via https://science.ebird.org/en/use-ebird-data/download-ebird-data-products). This code needs to be run using a command-line interface. 
  	* `1b_process_GAMs.R` - Estimate arrival phenology via GAM models for each species, cell, and year.
  	* `1c_get_ebird_rngs.R` - Get species’ ranges via eBird status and trends.
	* `1d_get_spec_co_corr.R` - Get correlation between local climate and oscillation indices.
	* `1e_daymet_by_cell_par.R` - Download climate data by hexgrid cell.
	* `1f_process_sp_yr_files.R` - Combine all individual arrival (GAM) estimates.
	* `1g_combine_enviro_GAMs.R` - Combine arrival estimates with climate data.
  * `2_Main_analysis` - Code to run main analysis.
  	* `2a_comb_model_data.R` - Process data to format for running Bayesian model in Stan.
	* `2b_run_mdl.R` - Run Bayesian model in Stan
	* `pc_mdl.stan` - Bayesian model in Stan
  * `3_Output_analysis` - Code to analyze model output, create figures, etc.
	* `Figures` - All code to recreate figures in manuscript, by figure number.
	* `3a_mdl_analysis.R` - Extract parameter estimates from model
	* `Variance_explained` - Code to calculate percent variance explained by species, and by cell

General overview, data:

The data in this section includes data used to run the main analysis ("mdl_data") as well as empty containers for unprocessed data that needs to be downloaded from external sources.

* `data/` - All data used in analysis
  * `GMTED_2010` - The 2010 Global Multi-resolution Terrain Elevation Data at 30 arc-second resolution should be downloaded (via https://www.usgs.gov/centers/eros/science/usgs-eros-archive-digital-elevation-global-multi-resolution-terrain-elevation) and placed in this folder. It should be named "NA_30AS_Elev_Med.tif"
  * `climate_oscillations` - Rasters for correlations between climate oscillations and environmental variables live here. They are provided from the authors
  * `enviro_var_daily` - This folder provides the folder structure used when processing Daymet data.
  * `mdl_data` - This folder includes all data used in the main analysis.
  * `output` – Includes model diagnostic information, figures, variance explained calculations, and model parameter estimates.
  * `sp_yr_files` – This folder will fill with GAM estimates for each species and year when section 1 of the code is run.
  * `year_files` – Folder to add year-specific eBird data files into. Format as "YEAR_all.txt"
  * `spec_breed_cells` – Folder with species specific ranges calculated via script 1c.
