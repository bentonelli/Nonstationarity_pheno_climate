Code and data for Nonstationarity in Climate x Arrival Phenology 

Note that fit model object is available under a seperate DOI: doi.org/10.5281/zenodo.13256390.

**Associated publications:**

**Repository structure:**

* `code/` - Code to process data, run models, analyze output
 * `0_master_script.R` - Code to run all main script files sequentially via source command.
  * `1_Prepare_data` - Code to preprocess data from eBird and Daymet. Please note this data needs to be dowloaded before this section is run.
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
