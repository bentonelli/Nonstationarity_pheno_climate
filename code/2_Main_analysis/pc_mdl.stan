// Pheno. Climate Model

data {
  int<lower=0> N;
  
  vector[N] yy;
  vector[N] yy_unc;
  
  vector[N] temp;
  vector[N] swe;
  vector[N] precip;
  
  vector[N] snow_p;
  
  vector[N] enso;
  vector[N] nao;
  
  vector[N] cell_lat;

  int<lower=1> Ncell;
  int<lower = 1, upper = Ncell> cc[N];
  
  int<lower=1> Nsp;
  int<lower = 1, upper = Nsp> ii[N];
  
  int<lower=1> N_sp_cl;
  int<lower = 1, upper = N_sp_cl> sp_cl[N];
  
  int<lower = 1, upper = Nsp> sp_cl_in[N_sp_cl];
  
  int<lower = 1, upper = Ncell> cl_in[N_sp_cl];
  
  vector[N_sp_cl] temp_niche;
  vector[N_sp_cl] swe_niche;
  vector[N_sp_cl] prcp_niche;
  
  vector[Nsp] sp_enso_pca1;
  vector[Nsp] sp_nao_pca1;
}

parameters {
  vector[N] yy_true;
  
  real gamma;
  real eta1;
  real eta2;
  real eta3;
  
  real nu1;
  real nu2;
  real nu3;
  
  real kappa1;
  real kappa2;
  
  vector[N_sp_cl] alpha;

  vector[Nsp] mu_alpha;
  
  real mu_zeta;
  vector[Nsp] zeta;

  real<lower = 0> sigma_zeta;
  
  real<lower=0> sigma_mu_alpha;
  
  real<lower=0> sigma_alpha;
  
  real<lower=0> sigma_beta1;
  real<lower=0> sigma_beta2;
  real<lower=0> sigma_beta3;

  real<lower=0>  sigma_phi1;
  real<lower=0>  sigma_phi2;
  
  real<lower=0> ln_mu;
  real<lower=0> ln_sigma;
  vector<lower=0>[Nsp] sigma;
  
  // For second part of model
  real omega1;
  real omega2;
  
  vector[Nsp] mu_phi1;
  vector[Nsp] mu_phi2;
  
  real<lower=0> sigma_mu_phi1;
  real<lower=0> sigma_mu_phi2;
  
  // For NC
  vector[N_sp_cl] beta1_raw;
  vector[N_sp_cl] beta2_raw;
  vector[N_sp_cl] beta3_raw;
  
  vector[N_sp_cl] phi1_raw;
  vector[N_sp_cl] phi2_raw;
  
  vector<lower=0>[2] sigma_b1s;
  cholesky_factor_corr[2] L_Rho_b1s;             // cholesky factor of corr matrix
  matrix[2, Nsp] z_b1s;                          // z-scores
  
  vector<lower=0>[2] sigma_b2s;
  cholesky_factor_corr[2] L_Rho_b2s;             // cholesky factor of corr matrix
  matrix[2, Nsp] z_b2s;                          // z-scores
  
  vector<lower=0>[2] sigma_b3s;
  cholesky_factor_corr[2] L_Rho_b3s;             // cholesky factor of corr matrix
  matrix[2, Nsp] z_b3s;                          // z-scores
  
}

transformed parameters {
  matrix[Nsp, 2] b1s;                 // betas & alpha: alpha, beta1, beta2, beta3
  matrix[2, 2] Rho_b1s;               // corr matrix
  
  matrix[Nsp, 2] b2s;                 // betas & alpha: alpha, beta1, beta2, beta3
  matrix[2, 2] Rho_b2s;               // corr matrix
  
  matrix[Nsp, 2] b3s;                 // betas & alpha: alpha, beta1, beta2, beta3
  matrix[2, 2] Rho_b3s;               // corr matrix
  
  vector[N_sp_cl] beta1;
  vector[N_sp_cl] beta2;
  vector[N_sp_cl] beta3;
  
  vector[Nsp] theta1;
  vector[Nsp] theta2;
  vector[Nsp] theta3;
  
  vector[Nsp] mu_beta1;
  vector[Nsp] mu_beta2;
  vector[Nsp] mu_beta3;
  
  vector[N_sp_cl] phi1;
  vector[N_sp_cl] phi2;
  
  b1s = (diag_pre_multiply(sigma_b1s, L_Rho_b1s) * z_b1s)';
  Rho_b1s = multiply_lower_tri_self_transpose(L_Rho_b1s);
  
  b2s = (diag_pre_multiply(sigma_b2s, L_Rho_b2s) * z_b2s)';
  Rho_b2s = multiply_lower_tri_self_transpose(L_Rho_b2s);
  
  b3s = (diag_pre_multiply(sigma_b3s, L_Rho_b3s) * z_b3s)';
  Rho_b3s = multiply_lower_tri_self_transpose(L_Rho_b3s);
  
  mu_beta1 = eta1 + b1s[,1];
  theta1 = nu1 + b1s[,2];
  
  mu_beta2 = eta2 + b2s[,1];
  theta2 = nu2 + b2s[,2];
  
  mu_beta3 = eta3 + b3s[,1];
  theta3 = nu3 + b3s[,2];
  
  beta1 = (mu_beta1[sp_cl_in] + theta1[sp_cl_in] .* temp_niche) + sigma_beta1 * beta1_raw;
  
  beta2 = (mu_beta2[sp_cl_in] + theta2[sp_cl_in] .* swe_niche) + sigma_beta2 * beta2_raw;
  
  beta3 = (mu_beta3[sp_cl_in] + theta3[sp_cl_in] .* prcp_niche) + sigma_beta3 * beta3_raw;
  
  phi1 = mu_phi1[sp_cl_in] + sigma_phi1 * phi1_raw;
  phi2 = mu_phi2[sp_cl_in] + sigma_phi2 * phi2_raw;
  
}

model {
  sigma_b1s ~ gamma(2,1);        // std dev beta1s, alphas
  sigma_b2s ~ gamma(2,1);        // std dev beta2s, alphas
  sigma_b3s ~ gamma(2,1);       // std dev beta3s, alphas

  gamma ~ normal(0,5);
  
  omega1 ~ std_normal();
  omega2 ~ std_normal();
  
  kappa1 ~ std_normal();
  kappa2 ~ std_normal();
  
  eta1 ~ normal(0,1);
  eta2 ~ normal(0,1);
  eta3 ~ normal(0,1);
  
  nu1 ~ normal(0,1);
  nu2 ~ normal(0,1);
  nu3 ~ normal(0,1);
  
  ln_mu ~ normal(1.35,.05);
  ln_sigma ~ normal(.3,.025);
  
  mu_zeta ~ std_normal();
  sigma_zeta ~ gamma(2,1);
  
  zeta ~ normal(mu_zeta, sigma_zeta);
  
  mu_alpha ~ normal(gamma,sigma_mu_alpha);
  
  sigma_mu_alpha ~ normal(15,5);
  
  sigma_alpha ~ normal(10,3);

  sigma_beta1 ~ gamma(3,.2);
  sigma_beta2 ~ gamma(3,.2);
  sigma_beta3 ~ gamma(3,.2);
  
  sigma_phi1 ~ gamma(3,.2);
  sigma_phi2 ~ gamma(3,.2);
  
  beta1_raw ~ std_normal();
  beta2_raw ~ std_normal();
  beta3_raw ~ std_normal();
  
  phi1_raw ~ std_normal();
  phi2_raw ~ std_normal();
  
  sigma_mu_phi1 ~ gamma(3,.3);
  sigma_mu_phi2 ~ gamma(3,.3);
  
  mu_phi1 ~ normal(omega1 + kappa1 * sp_enso_pca1,sigma_mu_phi1);
  mu_phi2 ~ normal(omega2 + kappa2 * sp_nao_pca1,sigma_mu_phi2);
  
  alpha ~ normal(mu_alpha[sp_cl_in] + zeta[sp_cl_in] .* cell_lat[cl_in],sigma_alpha);
  
  sigma ~ lognormal(ln_mu,ln_sigma);

  yy ~ normal(yy_true,yy_unc);
  yy_true ~ normal(alpha[sp_cl] + beta1[sp_cl] .* temp + snow_p .* beta2[sp_cl] .* swe + beta3[sp_cl] .* precip + phi1[sp_cl] .* enso + phi2[sp_cl] .* nao, sigma[ii]);
  
  // MVN
  to_vector(z_b1s) ~ std_normal();
  L_Rho_b1s ~ lkj_corr_cholesky(1);
  
  to_vector(z_b2s) ~ std_normal();
  L_Rho_b2s ~ lkj_corr_cholesky(1);
  
  to_vector(z_b3s) ~ std_normal();
  L_Rho_b3s ~ lkj_corr_cholesky(1);
}
