rm(list=ls())

#d <- read.csv('d:/dropbox/working/stan_bayesian_modeling/data/indata.dir/pom_flux/GO_flux_env_o2_t.csv')
d <- read.csv('d:/dropbox/teaching/stan_examples/data/GO_flux_env_o2_t.csv')
d <- d[complete.cases(d$flux_poc),]

#-TRY SUBSETTING <30 days-#
d <- d[d$time < 30,]

par(mfrow=c(2,3))
plot(d$flux_poc,-d$depth,ylim=c(-6000,0))
plot(d$flux_pic,-d$depth,ylim=c(-6000,0))
plot(d$flux_psi,-d$depth,ylim=c(-6000,0))

plot(d$flux_pic/d$flux_poc,-d$depth,ylim=c(-6000,0),xlim=c(0,10))



