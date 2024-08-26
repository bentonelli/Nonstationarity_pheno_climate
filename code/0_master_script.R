#Code to run scripts via source. Please see README before running.

# The below script is currently formatted to run an example species (Yellow-bellied Sapsucker)
# for the years 2012-2021. Due to computational requirements, it is suggested that this script be run
# in parallel split to different species and years.
source("code/1_Prepare_data/1b_process_GAMs.R")

#Script to get eBird S&T ranges.
#Access via the "ebirdst" package is needed for this to run
source("code/1_Prepare_data/1c_get_ebird_rngs.R")

#Script to get climate oscillation correlations for each speceis. 
#Access via the "ebirdst" package is needed for this to run.
source("code/1_Prepare_data/1d_get_spec_co_corr.R")

#Script to download Daymet data in parallel.
source("code/1_Prepare_data/1e_daymet_by_cell_par.R")

#Script to combine GAM estimates stored as individual files.
source("code/1_Prepare_data/1f_process_sp_yr_files.R")

#Script to combine GAM estimates now amalgamated with Daymet data.
source("code/1_Prepare_data/1g_combine_enviro_GAMs.R")

#Script to format data for model
source("code/2_Main_analysis/2a_comb_model_data.R")

#Script to run model
source("code/2_Main_analysis/2b_run_mdl.R")

#Script to run main model analysis
source("code/3_Output_analysis/3a_mdl_analysis.R")
