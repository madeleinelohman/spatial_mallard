
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
load("data/mall_data_new_0224_07072025.RData")
rm(list=setdiff(ls(), c('marr.am', 'marr.jm', 'n.strata', 'n.years', 'grd', 'start', 'end')))

source("haz_iar_MalePriors_noBeta.R")

marr.a = marr.am
marr.j = marr.jm

shape = grd
n.strata = nrow(shape)

#~~~~~~~~~~~~~ 
# Reporting rates
#~~~~~~~~~~~~~ 
report <- read.csv('data/mall_rr_new.csv')
report <- report[1:n.years,]

### Moment match to get alpha and beta
y <- report$p
sig <- report$se

v <- report$se^2

# Variance
alpha.prior <- y * (((y * (1-y)) / v) - 1)
beta.prior <- (1 - y) * (((y * (1-y)) / v) - 1)




#~~~~~~~~~~~~
# Spatial matrices
#~~~~~~~~~~~~
W <- nb2mat(poly2nb(shape), style = "B")
D <- matrix(0, nrow(shape), nrow(shape))
diag(D) <- rowSums(W)

#~~~~~~~~~~~~
# Covariates
#~~~~~~~~~~~~

#load("data/ag.RData")
load("data/spei48.RData")
load("data/sum_ag.RData")
# load("data/spei24.RData")

# corrplot(cor(spei2, fallow))
# summary(c(cor(spei2, fallow)))
# 
# corrplot(cor(crop, fallow))
# summary(c(cor(crop, fallow)))
# 
# corrplot(cor(spei2, crop))
# summary(c(cor(spei2, crop)))
# crop[which(is.na(crop),arr.ind=T)] <- 0
# fallow[which(is.na(fallow),arr.ind=T)] <- 0


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run model
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~
# Data
#~~~~~~~~~~~~
data=list("Ns"=n.strata, 
          "Nt"=n.years, 
          "D"=D,
          "W"=W,
          "marr_ad"=marr.a,
          "marr_j"=marr.j,
          "alpha_prior"=alpha.prior,
          "beta_prior"=beta.prior,
          "W_n" = sum(W) / 2,
          "crop"=t(crop),
          "fallow"=t(fallow),
          "spei" =t(spei2),
          "c"=0.2,
          "zeros" =rep(0, n.years))


#~~~~~~~~~~~~
# Initial values
#~~~~~~~~~~~~
pr = array(t(matrix(c(rep(0.015, n.years), 1-(0.015*n.years)), n.years+1, n.years)), dim=c(n.years, n.years+1, n.strata))
pr.j = array(t(matrix(c(rep(0.015, n.years), 1-(0.015*n.years)), n.years+1, n.years)), dim=c(n.years, n.years+1, n.strata))

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
            'site_sigma_ad_nm','site_sigma_j_nm',
            'site_sigma_ad_k','site_sigma_j_k',
            'alpha_nm', 'alpha_nm_j',
            'alpha_hm', 'alpha_hm_j',
            'xi', 'xi_j',
            'eps', 'eps_j',
            'site_ad_k', 'site_j_k',
            'site_ad_nm','site_j_nm')


# params <- c("sigma_nm_ad", "sigma_k_ad", "sigma_nm_j", "sigma_k_j",
#             "phi", "kappa", "eta", "phi_j", "kappa_j", "eta_j")

#~~~~~~~~~~~~
# Run model
#~~~~~~~~~~~~

start.time = Sys.time()
m <- stan(model_code=CAR_surv, data=data, 
           chains=3, iter=2000, init=initf2,
           cores=3, pars=params)
end.time = Sys.time()
end.time - start.time

save.image("IAR_haz_male_noBetaHM_priors.RData")


library(beepr)
beep(9)
