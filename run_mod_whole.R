
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set up
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(beepr)
library(prioritizr)
library(rstan)
library(sf)
library(spdep)
library(tidyverse)
library(corrplot)

setwd("~/Desktop/haz")
load("data/MALL_data2026.RData")
rm(list=setdiff(ls(), c('marr.am', 'marr.jm', 'n.strata', 'n.years', 'grd', 'start', 'end')))

source("haz_iar_noSpace_male.R")

marr.a = marr.am
marr.j = marr.jm

shape = grd
n.strata = nrow(shape)

#~~~~~~~~~~~~~ 
# Reporting rates
#~~~~~~~~~~~~~ 
report <- read.csv('data/boomer_rr.csv')
#report <- report[1:n.years,]
### Moment match to get alpha and beta
y <- report$rhoB
sig <- report$rhoB.se

v <- report$rhoB.se^2

# Variance
alpha.prior <- y * (((y * (1-y)) / v) - 1)
beta.prior <- (1 - y) * (((y * (1-y)) / v) - 1)

alpha.prior <- c(rep(alpha.prior[1],2),alpha.prior,rep(alpha.prior[46],2))
beta.prior <- c(rep(beta.prior[1],2),beta.prior,rep(beta.prior[46],2))

#~~~~~~~~~~~~
# Covariates
#~~~~~~~~~~~~
load("data/spei48_whole.RData")
load("data/sum_ag_whole.RData")


#~~~~~~~~~~~~
# M-arrays
#~~~~~~~~~~~~
marr.a <- apply(marr.a,c(1,2),sum)
marr.j <- apply(marr.j,c(1,2),sum)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run model
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~
# Data
#~~~~~~~~~~~~
data=list("Nt"=n.years, 
          "marr_ad"=marr.a,
          "marr_j"=marr.j,
          "alpha_prior"=alpha.prior,
          "beta_prior"=beta.prior,
          "crop"=(crop),
          "fallow"=(fallow),
          "spei" =(spei2),
          "c"=0.2,
          "zeros" =rep(0, n.years))


#~~~~~~~~~~~~
# Initial values
#~~~~~~~~~~~~
pr = array(t(matrix(c(rep(0.015, n.years), 1-(0.015*n.years)), n.years+1, n.years)), dim=c(n.years, n.years+1))
pr.j = array(t(matrix(c(rep(0.015, n.years), 1-(0.015*n.years)), n.years+1, n.years)), dim=c(n.years, n.years+1))

# inits <- function(){list(pr_ad = pr, pr_j = pr.j)}
# initf2 <- function(chain_id=1){
#   list(pr_ad = pr, pr_j = pr.j,
#        sigma_nm_ad=0.3,
#        sigma_nm_j=0.3,
#        sigma_k_ad=0.3,
#        sigma_k_j=0.3,
#        site_sigma_j_nm=1.2,
#        site_sigma_ad_nm=1.2,
#        site_sigma_j_k=1.2,
#        site_sigma_ad_k=1.2)
# }
initf2 <- function(chain_id=1){
  list(pr_ad = pr, pr_j = pr.j,
       sigma_nm_ad=0.3,
       sigma_nm_j=0.3,
       sigma_k_ad=0.3,
       sigma_k_j=0.3,
       site_sigma_j_nm=1.2,
       site_sigma_ad_nm=1.2,
       site_sigma_j_k=1.2,
       site_sigma_ad_k=1.2,
       beta_nm=rep(0,3),
       beta_nm_j=rep(0,3))
}

params <- c("sigma_nm_ad", "sigma_k_ad", "sigma_nm_j", "sigma_k_j",
            "phi", "kappa", "eta", "phi_j", "kappa_j", "eta_j",
            'beta_nm', 'beta_nm_j', 
            'alpha_nm', 'alpha_nm_j',
            'alpha_hm', 'alpha_hm_j',
            'xi', 'xi_j',
            'eps', 'eps_j')


# params <- c("sigma_nm_ad", "sigma_k_ad", "sigma_nm_j", "sigma_k_j",
#             "phi", "kappa", "eta", "phi_j", "kappa_j", "eta_j")

#~~~~~~~~~~~~
# Run model
#~~~~~~~~~~~~

start.time = Sys.time()
m <- stan(model_code=mall_noSpace_mal, data=data, 
           chains=3, iter=2000, init=initf2,
           cores=3, pars=params)
end.time = Sys.time()
end.time - start.time

save.image("haz_iar_noSpace_male.RData")


library(beepr)
beep(9)
