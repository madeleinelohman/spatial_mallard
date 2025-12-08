CAR_surv="
functions {
  real sparse_iar_lpdf(vector phi, real tau, array[,] int W_sparse,
                     vector D_sparse, vector lambda, int n, int W_n) {
  row_vector[n] phit_D; // phi' * D
    row_vector[n] phit_W; // phi' * W
  vector[n] ldet_terms;
  
  phit_D = (phi .* D_sparse)';
    phit_W = rep_row_vector(0, n);
    for (i in 1 : W_n) {
      phit_W[W_sparse[i, 1]] = phit_W[W_sparse[i, 1]] + phi[W_sparse[i, 2]];
      phit_W[W_sparse[i, 2]] = phit_W[W_sparse[i, 2]] + phi[W_sparse[i, 1]];
    }
    
    return 0.5 * ((n - 1) * log(tau) - tau * (phit_D * phi - (phit_W * phi)));
  }
}
data{
  int Ns; //This is the number of sites
  int Nt; //This is the number of years
  
  int<lower=0> marr_ad[Nt, Nt+1, Ns]; //Adult m-array
  int<lower=0> marr_j[Nt, Nt+1, Ns]; //Juvenile m-array
  
  matrix[Ns, Ns] W; //Neighborhood matrix
  int W_n; // number of adjacent region pairs
  
  real<lower=0, upper=1> c; // crippling rate
  
  real alpha_prior[Nt];
  real beta_prior[Nt];
  
  
  matrix[Nt, Ns] crop; //crop covariate
  matrix[Nt, Ns] fallow; //fallow covariate
  matrix[Nt, Ns] spei; //pond mean
  
}
transformed data {
  array[W_n, 2] int W_sparse; // adjacency pairs
  vector[Ns] D_sparse; // diagonal of D (number of neigbors for each site)
  vector[Ns] lambda; // eigenvalues of invsqrtD * W * invsqrtD
  
  {
    // generate sparse representation for W
    int counter;
    counter = 1;
    // loop over upper triangular part of W to identify neighbor pairs
    for (i in 1 : (Ns - 1)) {
      for (j in (i + 1) : Ns) {
        if (W[i, j] == 1) {
          W_sparse[counter, 1] = i;
          W_sparse[counter, 2] = j;
          counter = counter + 1;
        }
      }
    }
  }
  for (i in 1 : Ns) {
    D_sparse[i] = sum(W[i]);
  }
  {
    vector[Ns] invsqrtD;
    for (i in 1 : Ns) {
      invsqrtD[i] = 1 / sqrt(D_sparse[i]);
    }
    lambda = eigenvalues_sym(quad_form(W, diag_matrix(invsqrtD)));
  }
}

parameters{
  real<lower=0, upper=1> r[Nt]; // reporting rate
  
  //real<lower=0, upper=1> rho_nm; //This is the correlation coefficient for harvest mortality
  //real<lower=0, upper=1> rho_hm; //This is the correlation coefficient for harvest mortality
   
  real<lower=0> sigma_nm_j; //This is the SD on juvenile survival
  real<lower=0> sigma_k_j; //This is the SD on juvenile harvest mortality
  real<lower=0> site_sigma_j_nm;
  real<lower=0> site_sigma_j_k;
  vector[Ns] site_j_nm;
  vector[Ns] site_j_k;
  
  real<lower=0> sigma_nm_ad; //This is the SD on adult survival
  real<lower=0> sigma_k_ad; //This is the SD on adult harvest mortality
  real<lower=0> site_sigma_ad_nm;
  real<lower=0> site_sigma_ad_k;
  vector[Ns] site_ad_nm;
  vector[Ns] site_ad_k;
  
  
  vector[Nt] xi; //Adult harvest mortality untransformed
  vector[Nt] eps; //Adult natural mortality untransformed
  
  vector[Nt] xi_j; //Juvenile harvest mortality untransformed
  vector[Nt] eps_j; //Juvenile natural mortality untransformed
  
  real alpha_hm;
  real alpha_nm;
  real alpha_hm_j;
  real alpha_nm_j;
  
  
  real beta_nm[3]; //Beta coefficients for adult harvest mortality
  real beta_nm_j[3]; //Beta coefficients for adult harvest mortality
}

transformed parameters{
  //SET UP
  
  matrix[Nt, Ns] nm; //Adult natural mortality hazard rate
  matrix[Nt, Ns] hm; //Adult harvest mortality hazard rate
  matrix[Nt, Ns] nm_j; //Juvenile natural mortality hazard rate
  matrix[Nt, Ns] hm_j; //Juvenile harvest mortality hazard rate
  
  matrix<lower=0, upper=1>[Nt-1, Ns] phi; //Survival probability
  matrix<lower=0, upper=1>[Nt, Ns] kappa; //Harvest mortality probability
  matrix<lower=0, upper=1>[Nt, Ns] eta; //Natural mortality probability
  matrix<lower=0, upper=1>[Nt, Ns] f; //Direct Recovery rate
  
  matrix<lower=0, upper=1>[Nt-1, Ns] phi_j; //Juvenile Survival probability
  matrix<lower=0, upper=1>[Nt, Ns] kappa_j; //Juvenile Harvest mortality probability
  matrix<lower=0, upper=1>[Nt, Ns] eta_j; //Juvenile Natural mortality probability
  matrix<lower=0, upper=1>[Nt, Ns] f_j; //Juvenile Direct Recovery rate
  
  real<lower=0, upper=1> pr[Nt, Nt+1, Ns]; //m-array probabilities
  real<lower=0, upper=1> pr_j[Nt, Nt+1, Ns]; //Juvenile m-array probabilities
  
  
  //HAZARD RATES
  for(t in 1:Nt){
    for(s in 1:Ns){
      hm[t,s] = exp(alpha_hm + xi[t] + site_ad_k[s]);
      nm[t,s] = exp(alpha_nm + eps[t] + site_ad_nm[s] + beta_nm[1]*crop[t,s] + beta_nm[2]*spei[t,s] + beta_nm[3]*fallow[t,s]);
      hm_j[t,s] = exp(alpha_hm_j + xi_j[t] + site_j_k[s]);
      nm_j[t,s] = exp(alpha_nm_j + eps_j[t] + site_j_nm[s] + beta_nm_j[1]*crop[t,s] + beta_nm_j[2]*spei[t,s] + beta_nm_j[3]*fallow[t,s]);
      
      //PROBABILITIES
      kappa[t,s] = 1 - exp(-(hm[t,s]));
      eta[t,s] = (1-kappa[t,s]) * (1-exp(-(nm[t,s])));
      kappa_j[t,s] = 1 - exp(-(hm_j[t,s]));
      eta_j[t,s] = (1-kappa_j[t,s]) * (1-exp(-(nm_j[t,s])));
      
      f[t, s] = kappa[t, s] * r[t] * (1-c);
      f_j[t, s] = kappa_j[t, s] * r[t] * (1-c);
    }

  }
  
  for(s in 1:Ns){
    for(t in 1:(Nt-1)){
      phi[t,s] = exp(-(nm[t,s] + hm[t,s]));
      phi_j[t,s] = exp(-(nm_j[t,s] + hm_j[t,s]));
    }
  }
  
  
  //M-ARRAY PROBABILITIES
  // Adults
  for (t in 1:(Nt-1)){
    for (j in (t+1):Nt){
      for (s in 1:Ns){
        pr[t,j,s] = prod(phi[t:(j-1), s]) * f[j, s];   
      }
    }
  }
  
  for (t in 2:Nt){
    for (j in 1:(t-1)){
      for (s in 1:Ns){
        pr[t,j,s] = 0;
      }
    }
  }
  
  for (t in 1:Nt){
    for (s in 1:Ns){
      pr[t,t,s] = f[t, s];
      pr[t,(Nt+1),s] = 1 - sum(pr[t, 1:Nt, s]);
    }
  }
  
  
  // Juveniles
  for (s in 1:Ns){
    for (t in 1:(Nt-1)){
      pr_j[t,t+1,s] = phi_j[t, s] * f[t+1, s];
      if((t+2) <= Nt)
        for (j in (t+2):Nt){
          pr_j[t,j,s] = phi_j[t, s] * prod(phi[(t+1):(j-1), s]) * f[j, s]; 
        }
    }
  }
  
  for (t in 2:Nt){
    for (j in 1:(t-1)){
      for (s in 1:Ns){
        pr_j[t,j,s] = 0;
      }
    }
  }
  
  for (t in 1:Nt){
    for (s in 1:Ns){
      pr_j[t,t,s] = f_j[t, s];
      pr_j[t,(Nt+1),s] = 1 - sum(pr_j[t, 1:Nt, s]);
    }
  }
}

model{
  //PRIORS 
  alpha_hm ~ normal(-0.9,1);
  alpha_hm_j ~ normal(-1.3,1);
  alpha_nm ~ normal(-0.8,1);
  alpha_nm_j ~ normal(-1.2,1);
  
  site_sigma_ad_nm ~ uniform(0,3);
  site_sigma_j_nm ~ uniform(0,3);
  site_sigma_ad_k ~ uniform(0,3);
  site_sigma_j_k ~ uniform(0,3);

  sigma_nm_ad ~ uniform(0,3);
  sigma_nm_j ~ uniform(0,3);
  sigma_k_ad ~ uniform(0,3);
  sigma_k_j ~ uniform(0,3);
  
  // Reporting rate
  for(t in 1:Nt){
    r[t] ~ beta(alpha_prior[t], beta_prior[t]);
  }
  
  //rho_nm ~ uniform(0.9,1);
  //rho_hm ~ uniform(0.9,1);
  
  site_ad_nm ~ sparse_iar(1/(site_sigma_ad_nm^2),W_sparse, D_sparse, lambda, Ns, W_n);
  site_j_nm ~ sparse_iar(1/(site_sigma_j_nm^2), W_sparse, D_sparse, lambda, Ns, W_n);
  site_ad_k ~ sparse_iar(1/(site_sigma_ad_k^2), W_sparse, D_sparse, lambda, Ns, W_n);
  site_j_k ~ sparse_iar(1/(site_sigma_j_k^2), W_sparse, D_sparse, lambda, Ns, W_n);
    
  xi[1] ~ normal(0, sigma_k_ad);
  eps[1] ~ normal(0, sigma_nm_ad);

  xi_j[1] ~ normal(0, sigma_k_j);
  eps_j[1] ~ normal(0, sigma_nm_j);
  
  for(t in 2:Nt){
    xi[t] ~ normal(0, sigma_k_ad);
    eps[t] ~ normal(0, sigma_nm_ad);
    xi_j[t] ~ normal(0, sigma_k_j);
    eps_j[t] ~ normal(0, sigma_nm_j);
  }
  
  for(i in 1:3){
    
    beta_nm[i] ~ normal(0, 3);
  
    beta_nm_j[i] ~ normal(0, 3);
  }
  
  //LIKELIHOOD
  for (t in 1:Nt){
    for (s in 1:Ns){
      marr_ad[t,,s] ~ multinomial(to_vector(pr[t,,s]));
      marr_j[t,,s] ~ multinomial(to_vector(pr_j[t,,s]));
    }
  }
}

"


