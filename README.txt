Code and data for Nonstationarity in Climate x Arrival Phenology 

**Associated publications:**

**Repository structure:**

* `code/` - Code to process data, run models, analyze output
  * `1_Prepare_data` - Code to preprocess data
  	* `1a_awk_filter_ebird.R` - Filter eBird data
  	* `1b_process_GAMs.R` - Estimate arrival phenology
  	* `1c_get_ebird_rngs.R` - Get species’ ranges
	* `1d_get_spec_co_corr.R` - Get correlation between local climate and oscillation indices
	* `1e_daymet_by_cell_par.R` - Download climate data, by hexgrid cell
	* `1f_process_sp_yr_files.R` - Combine all GAM estimates
	* `1g_combine_enviro_GAMs.R` - Combine arrival estimates with climate data
  * `2_Main_analysis` - Code to run main analysis
  	* `2a_comb_model_data.R` - Process data for Bayesian model
	* `2b_run_mdl.R` - Run Bayesian model
	* `pc_mdl.stan` - Bayesian model in Stan
  * `3_Output_analysis` - Code to analyze model output
	* `Figures` - All code to recreate figures in manuscript, by figure number.
	* `mdl_analysis.R` - Extract parameter estimates from model
	* `Variance_explained` - Code to calculate percent variance explained by species, and by cell
