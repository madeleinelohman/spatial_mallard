mall_noSpace_mal="
data{
  int Nt; //This is the number of years
  
  int<lower=0> marr_ad[Nt, Nt+1]; //Adult m-array
  int<lower=0> marr_j[Nt, Nt+1]; //Juvenile m-array
  
  real<lower=0, upper=1> c; // crippling rate
  
  real alpha_prior[Nt];
  real beta_prior[Nt];
  
  
  matrix[Nt] crop; //crop covariate
  matrix[Nt] fallow; //fallow covariate
  matrix[Nt] spei; //pond mean
  
}
parameters{
  real<lower=0, upper=1> r[Nt]; // reporting rate
  
  real<lower=0> sigma_nm_j; //This is the SD on juvenile survival
  real<lower=0> sigma_k_j; //This is the SD on juvenile harvest mortality
  
  real<lower=0> sigma_nm_ad; //This is the SD on adult survival
  real<lower=0> sigma_k_ad; //This is the SD on adult harvest mortality
  
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
  
  matrix[Nt] nm; //Adult natural mortality hazard rate
  matrix[Nt] hm; //Adult harvest mortality hazard rate
  matrix[Nt] nm_j; //Juvenile natural mortality hazard rate
  matrix[Nt] hm_j; //Juvenile harvest mortality hazard rate
  
  matrix<lower=0, upper=1>[Nt-1] phi; //Survival probability
  matrix<lower=0, upper=1>[Nt] kappa; //Harvest mortality probability
  matrix<lower=0, upper=1>[Nt] eta; //Natural mortality probability
  matrix<lower=0, upper=1>[Nt] f; //Direct Recovery rate
  
  matrix<lower=0, upper=1>[Nt-1] phi_j; //Juvenile Survival probability
  matrix<lower=0, upper=1>[Nt] kappa_j; //Juvenile Harvest mortality probability
  matrix<lower=0, upper=1>[Nt] eta_j; //Juvenile Natural mortality probability
  matrix<lower=0, upper=1>[Nt] f_j; //Juvenile Direct Recovery rate
  
  real<lower=0, upper=1> pr[Nt, Nt+1]; //m-array probabilities
  real<lower=0, upper=1> pr_j[Nt, Nt+1]; //Juvenile m-array probabilities
  
  
  //HAZARD RATES
  for(t in 1:Nt){
    hm[t] = exp(alpha_hm + xi[t]);
    nm[t] = exp(alpha_nm + eps[t] + beta_nm[1]*crop[t] + beta_nm[2]*spei[t] + beta_nm[3]*fallow[t]);
    hm_j[t] = exp(alpha_hm_j + xi_j[t]);
    nm_j[t] = exp(alpha_nm_j + eps_j[t] + beta_nm_j[1]*crop[t] + beta_nm_j[2]*spei[t] + beta_nm_j[3]*fallow[t]);
    
    //PROBABILITIES
    kappa[t] = 1 - exp(-(hm[t]));
    eta[t] = (1-kappa[t]) * (1-exp(-(nm[t])));
    kappa_j[t] = 1 - exp(-(hm_j[t]));
    eta_j[t] = (1-kappa_j[t]) * (1-exp(-(nm_j[t])));
    
    f[t] = kappa[t] * r[t] * (1-c);
    f_j[t] = kappa_j[t] * r[t] * (1-c);
  }

  for(t in 1:(Nt-1)){
    phi[t] = exp(-(nm[t] + hm[t]));
    phi_j[t] = exp(-(nm_j[t] + hm_j[t]));
  }
  
  
  //M-ARRAY PROBABILITIES
  // Adults
  for (t in 1:(Nt-1)){
    for (j in (t+1):Nt){
      pr[t,j] = prod(phi[t:(j-1)]) * f[j];   
    }
  }
  
  for (t in 2:Nt){
    for (j in 1:(t-1)){
      pr[t,j] = 0;
    }
  }
  
  for (t in 1:Nt){
    pr[t,t] = f[t];
    pr[t,(Nt+1)] = 1 - sum(pr[t, 1:Nt]);
  }
  
  
  // Juveniles
  for (t in 1:(Nt-1)){
    pr_j[t,t+1] = phi_j[t] * f[t+1];
    if((t+2) <= Nt)
      for (j in (t+2):Nt){
        pr_j[t,j] = phi_j[t] * prod(phi[(t+1):(j-1)]) * f[j]; 
      }
  }
  
  for (t in 2:Nt){
    for (j in 1:(t-1)){
      pr_j[t,j] = 0;
    }
  }
  
  for (t in 1:Nt){
    pr_j[t,t] = f_j[t];
    pr_j[t,(Nt+1)] = 1 - sum(pr_j[t, 1:Nt]);
  }
}

model{
  //PRIORS 
  alpha_hm ~ normal(-1.6,1);
  alpha_hm_j ~ normal(-1.6,1);
  alpha_nm ~ normal(-1.5,1);
  alpha_nm_j ~ normal(-1.5,1);
  
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
    marr_ad[t,] ~ multinomial(to_vector(pr[t,]));
    marr_j[t,] ~ multinomial(to_vector(pr_j[t,]));
  }
}

"


