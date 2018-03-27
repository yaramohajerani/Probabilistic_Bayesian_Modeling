
rm(list=ls())

#d <- read.csv('d:/dropbox/working/stan_bayesian_modeling/data/indata.dir/pom_flux/GO_flux_env_o2_t.csv')
d <- read.csv('d:/dropbox/teaching/stan_examples/data/GO_flux_env_o2_t.csv')
d <- d[complete.cases(d$flux_poc),]

IDs <- unique(d$id)
ids <- c()
bs  <- c()
logJ0s <- c()
ssts <- c()
t_avgs <- c()
npps <- c()
o2_avgs <- c()


pdf('d:/dropbox/teaching/stan_examples/greg/plots/ind_regressions.pdf')
	for(i in 1:length(IDs)){
		#if(i==197 | i==202){i = i+1}
		dtmp <- d[d$id==IDs[i],]
		if(nrow(dtmp)>1){if(sd(dtmp$depth)>0){
			plot(log(dtmp$depth/100), log(dtmp$flux_poc))
			fit <- lm(log(dtmp$flux_poc) ~ log(dtmp$depth/100))
			bs     <- c(bs,summary(fit)$coefficients[2,1])
			logJ0s <- c(logJ0s,summary(fit)$coefficients[1,1])
			ids <- c(ids,IDs[i])
			ssts <- c(ssts,mean(dtmp$sst,na.rm=TRUE))
		}}
	}
dev.off()



