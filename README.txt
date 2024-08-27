Code and data for Nonstationarity in Climate x Arrival Phenology 

Note that fit model object is available under a seperate DOI: doi.org/10.5281/zenodo.13256390.

**Associated publications:**

**Repository structure:**

General outline and overview:
The code here is used to estimate arrival dates for North American migratory birds using eBird data, and than compare those arrival dates to environmental covariates derived from Daymet. This code contains three sections, the first ("1_Prepare_data") is used to filter the global eBird database to North America and according to certain checklists standards, split the data by year, and than use that data to estimate arrival dates for species in specific locations. This section also includes code to download environmental covariates via Daymet using the associated R package "daymetr", which is then combined into a single dataset. Note that to run this section, users will need to request and download the eBird dataset (see link below). Therefore, this section of code is not self-contained.

Section 2 ("2_Main_analysis") includes code to format the data to be analyzed via a single Bayesian hierarchical model (also included in this section). The code to run the main analysis is included in this section as well. This section of code is self-contained, as the data used in the analysis is included in the data folder.

Section 3 ("3_Output_analysis") contains code to get posterior estimates from the model, and to reproduce figures presented in the associated publication.

* `code/` - Code to preprocess data from eBird and Daymet (section 1), run models (section 2), analyze output (section 3).
 * `0_master_script.R` - Code to run all main script files sequentially via source command. Please read comments carefully before running, as some scripts have very long runtimes.
  * `1_Prepare_data` - The code in this section is used to preprocess data from eBird and Daymet. Please note this data needs to be dowloaded before this section is run.
  	* `1a_awk_filter_ebird.rtf` - Filter eBird data using awk (i.e. full eBird dataset available for download via https://science.ebird.org/en/use-ebird-data/download-ebird-data-products)
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
