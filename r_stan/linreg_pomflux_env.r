####################################################################
## CODE FOR FITTING BAYESIAN HIERARCHICAL FLUX PROFILE #############
## - directory references relative to /stan_bayesian_modelling/r_stan/
####################################################################
library(rstan)
library(rstan)   #install package
options(mc.cores = parallel::detectCores()) #tell Stan to use multiple cores

############################################################
## READ DATA ###############################################
############################################################
dat <- read.csv('../data/GO_flux_env_o2_t.csv')

mod <- stan_model('../stan/linreg_pomflux_env.stan')







