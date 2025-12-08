library(beepr)
library(prioritizr)
library(rstan)
library(gstat)
library(sf)
library(spdep)
library(tidyverse)

post <- rstan::extract(m)
m@model_pars


#~~~~~~~~~~~~~~~~~~~
# Convergence
#~~~~~~~~~~~~~~~~~~~
conv <- summary(m)$summary
summary(conv[,10])
hist(conv[,10])
bad <- which(conv[,10] > 1.1, arr.ind=T)
conv[bad,]
summary(conv[bad,10])
hist(conv[bad,10])
length(bad) / length(conv)

#~~~~~~~~~~~~~~~~~~~
# Sigmas
#~~~~~~~~~~~~~~~~~~~
plot(post$sigma_k_ad, type='l')
plot(post$sigma_k_j, type='l')

plot(post$sigma_nm_ad, type='l')
plot(post$sigma_nm_j, type='l')


#~~~~~~~~~~~~~~~~~~~
# Survival (rough, not w concern to time or site)
#~~~~~~~~~~~~~~~~~~~
hist(post$phi)
hist(post$phi_j)
plot(apply(post$phi, 2, mean))
plot(apply(post$phi_j, 2, mean))


#~~~~~~~~~~~~~~~~~~~
# Mortality intercepts
#~~~~~~~~~~~~~~~~~~~
plot(post$alpha_hm, type='l')
plot(post$alpha_hm_j, type='l')
plot(post$alpha_nm, type='l')
plot(post$alpha_nm_j, type='l')

plot(plogis(post$alpha_hm), type='l')
plot(plogis(post$alpha_hm_j), type='l')
plot(post$alpha_nm, type='l')
plot(post$alpha_nm_j, type='l')

#~~~~~~~~~~~~~~~~~~~
# Harvest mortality (rough, not w concern to time or site)
#~~~~~~~~~~~~~~~~~~~
hist(post$kappa)
hist(post$kappa_j)

plot(apply(post$kappa,c(2),mean)[1:48])
arrows(y0=apply(post$kappa,c(2),quantile, probs=0.025)[1:48], 
       y1=apply(post$kappa,c(2),quantile, probs=0.975)[1:48],
       x0=1:48,
       length=0)
plot(apply(post$kappa_j,c(2),mean)[1:48])
arrows(y0=apply(post$kappa_j,c(2),quantile, probs=0.025)[1:48], 
       y1=apply(post$kappa_j,c(2),quantile, probs=0.975)[1:48],
       x0=1:48,
       length=0)

# for(i in 1:n.years){
#   plot(post$eps_ad[,i], type='l')
#   Sys.sleep(0.2)
# }
# 
# for(i in 1:n.years){
#   plot(post$eps_j[,i], type='l')
#   Sys.sleep(0.2)
# }

plot(apply(post$kappa, 2, mean)[1:49] + apply(post$phi, 2, mean))
plot(apply(post$kappa_j, 2, mean)[1:49] + apply(post$phi_j, 2, mean))

for(i in 1:n.strata){
  plot(apply(post$kappa_j[,,i], 2, mean)[1:43] + apply(post$phi_j[,,i], 2, mean))
  Sys.sleep(0.2)
}


#~~~~~~~~~~~~~~~~~~~
# Natural mortality (rough, not w concern to time or site)
#~~~~~~~~~~~~~~~~~~~
hist(post$eta)
hist(post$eta_j)

plot(apply(post$eta,c(2),mean))
plot(apply(post$eta_j,c(2),mean))

plot(apply(marr.j,c(1),sum))
plot(apply(marr.j,c(2),sum))

#~~~~~~~~~~~~~~~~~~~
# Beta values
#~~~~~~~~~~~~~~~~~~~
# hist(post$beta_hm[,1], main="HM Adult Land")
# abline(v=0, col="red", lwd=2)
# 
# hist(post$beta_hm[,2], main="HM Adult SPEI")
# abline(v=0, col="red", lwd=2)
# 
# hist(post$beta_hm[,3], main="HM Adult fallow")
# abline(v=0, col="red", lwd=2)


# hist(post$beta_p_ad[,1], main="Survival Adult Land")
# abline(v=0, col="red", lwd=2)
# plot(post$beta_p_ad[,1], main="Survival Adult Land", type='l')
# 
# hist(post$beta_p_ad[,2], main="Survival Adult SPEI")
# abline(v=0, col="red", lwd=2)
# 
# plot(post$beta_p_ad[,2], main="Survival Adult SPEI",type='l')
# 
# hist(post$beta_p_ad[,3], main="Survival Adult Fallow")
# abline(v=0, col="red", lwd=2)
# plot(post$beta_p_ad[,3], main="Survival Adult Fallow",type='l')


hist(post$beta_k_j[,1], main="HM Juv Land")
abline(v=0, col="red", lwd=2)

hist(post$beta_k_j[,2], main="HM Juv SPEI")
abline(v=0, col="red", lwd=2)

hist(post$beta_k_j[,3], main="HM Juv Fallow")
abline(v=0, col="red", lwd=2)
plot(post$beta_k_j[,3], main="HM Juv Fallow", type='l')

hist(post$beta_p_j[,1], main="Survival Juv Land")
abline(v=0, col="red", lwd=2)
plot(post$beta_p_j[,1], main="Survival Juv Land", type='l')

hist(post$beta_p_j[,2], main="Survival Juv SPEI")
abline(v=0, col="red", lwd=2)
plot(post$beta_p_j[,2], main="Survival Juv SPEI", type='l')

hist(post$beta_p_j[,3], main="Survival Juv Fallow")
abline(v=0, col="red", lwd=2)
plot(post$beta_p_j[,3], main="Survival Juv Fallow", type='l')



betas <- data.frame(age="Adult", rate='harvest mortality', cov='fallow', value=post$beta_hm[,3])
betas <- rbind(betas, data.frame(age="Juvenile", rate='harvest mortality', cov='fallow', value=post$beta_hm_j[,3]))
betas <- rbind(betas, data.frame(age="Juvenile", rate='harvest mortality', cov='land', value=post$beta_hm[,1]))
betas <- rbind(betas, data.frame(age="Adult", rate='harvest mortality', cov='land', value=post$beta_hm[,1]))
betas <- rbind(betas, data.frame(age="Juvenile", rate='harvest mortality', cov='spei', value=post$beta_hm_j[,2]))
betas <- rbind(betas, data.frame(age="Adult", rate='harvest mortality', cov='spei', value=post$beta_hm[,2]))

ggplot(betas, aes(x=age, y=value, fill=cov)) +
  geom_violin() +
  geom_hline(yintercept=0, color='red') +
  title('Harvest mortality')

betas %>%
  group_by(age, cov) %>%
  summarise(mu=mean(value),
            q2.5=quantile(value, probs=0.025),
            q97.5=quantile(value, probs=0.975))




betas <- data.frame(age="Adult", rate='natural mortality', cov='fallow', value=post$beta_nm[,3])
betas <- rbind(betas, data.frame(age="Juvenile", rate='natural mortality', cov='fallow', value=post$beta_nm_j[,3]))
betas <- rbind(betas, data.frame(age="Juvenile", rate='natural mortality', cov='land', value=post$beta_nm_j[,1]))
betas <- rbind(betas, data.frame(age="Adult", rate='natural mortality', cov='land', value=post$beta_nm[,1]))
betas <- rbind(betas, data.frame(age="Juvenile", rate='natural mortality', cov='spei', value=post$beta_nm_j[,2]))
betas <- rbind(betas, data.frame(age="Adult", rate='natural mortality', cov='spei', value=post$beta_nm[,2]))

ggplot(betas, aes(x=age, y=value, fill=cov)) +
  geom_violin() +
  geom_hline(yintercept=0, color='red') +
  title('Natural mortality')

betas %>%
  group_by(age, cov) %>%
  summarise(mu=mean(value),
            q2.5=quantile(value, probs=0.025),
            q97.5=quantile(value, probs=0.975))
