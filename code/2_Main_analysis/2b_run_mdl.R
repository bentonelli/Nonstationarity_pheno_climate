library(rstan)

data_in1 <- readRDS("data/mdl_data/ms_data1.rds")

fit1 <- stan(
  file = "code/2_Main_analysis/pc_mdl.stan",  # Stan program
  data = data_in1,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 9000,            # total number of iterations per chain
  cores = 4,             # number of cores (could use one per chain)
  thin = 8,
  control = list(adapt_delta = .95,
                 max_treedepth = 10,
                 stepsize=.06),
  pars = c("beta1_raw","beta2_raw","beta3_raw","phi1_raw","phi2_raw"),
  include = FALSE,
  init_r = .5,
  refresh = 50,
  save_warmup = FALSE
)

saveRDS(fit1,"data/output/PhenoClimate_mdl.rds")

