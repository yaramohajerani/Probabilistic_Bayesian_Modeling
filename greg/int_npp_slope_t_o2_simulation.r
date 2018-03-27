library(rstan)   #install package
library(MASS)
options(mc.cores = parallel::detectCores()) #tell Stan to use multiple cores
# setwd('c:/users/greg/dropbox/teaching/stan_examples/')  #setwd(<wd must contain .stan file>)
# setwd('d:/dropbox/teaching/stan_examples/')  #setwd(<wd must contain .stan file>)

rm(list=ls())

#######################################################
## Setup your data ####################################
#######################################################
p  <- 20 #number of locations
n  <- sample(20:100,p,replace=TRUE) #data points per location
N  <- sum(n) #total number of data points
ni <- c(0,cumsum(n)) #vector index for dataset stacked as single column
y  <- numeric(N) #allocate space for y
x  <- numeric(N) #allocate space for x
T  <- runif(p,0,3)
npp <- runif(p,1,5)
o2 <- runif(p,0,1)

beta0slp_npp <- 0.05 #slope of profile slope vs. temperature
beta0int     <- 1 #noise around the line
beta0sd      <- 0.02
beta0        <- rnorm(n=p, mean=beta0int + beta0slp_npp*npp, sd=beta0sd)

beta1slp_T   <- 0.1
beta1slp_o   <- 0.2
beta1_int    <- -0.5
beta1sd      <- 0.1
beta1        <- rnorm(n=p, mean=beta1_int + beta1slp_T*T + beta1slp_o*o2, sd=beta1sd) #slope as a function of temperature plus noise

sigma_obs    <- 10

cat('Generating simulated data \n')
for(i in 1:p){ #loop through locations
	xtmp <- runif(n[i],1,500) #synthesize independent variable
	y[(ni[i]+1):ni[i+1]] <- beta0[i] + beta1[i]*xtmp + rnorm(n=n[i],mean=0,sd=sigma_obs) #generate station-specific data as an elongating vector
    x[(ni[i]+1):ni[i+1]] <- xtmp #store the independent variable in the big vector
}

#--MAKE A PLOT--################
# par(mfrow=c(2,5),mar=c(2,2,2,2)) #2 rows of panels, 3 columns
# for(i in 1:p){ #loop through locations  
	# plot(y[(ni[i]+1):ni[i+1]],x[(ni[i]+1):ni[i+1]],ylim=rev(range(x[(ni[i]+1):ni[i+1]])),xlim=range(y)) #plot like a depth profile 
	# mtext(T[i])
# }
dat <- list(N=N, ni=ni, y=y, p=p, T=T, npp=npp, o2=o2) #package data for rstan

#######################################################
## Fit Stan model #####################################
#######################################################
cat('Compiling Stan code \n')
##--COMPILE--#########
#mod <- stan_model('int_npp_slope_t_o2.stan')  #pre-compile
#save('mod',file='int_npp_slope_t_o2.rdata') #save compiled file	
load('int_npp_slope_t_o2.rdata')

#mod2 <- stan_model('int_t_slope_npp.stan')
#save('mod2',file='int_t_slope_npp.rdata') #save compiled file	
load('int_t_slope_npp.rdata')

##--FIT--##############
cat('Fitting model \n')
fit    <- sampling(mod,  data=dat, iter=4000, chains=4, warmup=2000, init=0, control=list(adapt_delta=0.95))  #fit model
post   <- extract(fit)
beta0_post <- colMeans(post$beta0)
beta1_post <- colMeans(post$beta1)
sigma_post <- mean(post$sigma)
llpost  <- sum(colMeans(post$log_lik))

fit2   <- sampling(mod2, data=dat, iter=4000, chains=4, warmup=2000, init=0, control=list(adapt_delta=0.95))
post2  <- extract(fit2)
beta0_post2 <- colMeans(post2$beta0)
beta1_post2 <- colMeans(post2$beta1)
sigma_post2 <- mean(post2$sigma)
llpost2 <- sum(colMeans(post2$log_lik)) 

llmp=llmp2 <- 0
for(i in 1:p){
	llmp  <- llmp  + sum(log(dnorm(x=y[(ni[i]+1):ni[i+1]], mean=beta0_post[i] + beta1_post[i]*x[(ni[i]+1):ni[i+1]], sd=sigma_post)))
	llmp2 <- llmp2 + sum(log(dnorm(x=y[(ni[i]+1):ni[i+1]], mean=beta0_post2[i] + beta1_post2[i]*x[(ni[i]+1):ni[i+1]], sd=sigma_post2))) 
}

dic1 <- 2*(llmp)  - 4*(llmp - llpost)
dic2 <- 2*(llmp2) - 4*(llmp2 - llpost2)

#################################################################
## PRINT THE OUTPUT #############################################
#################################################################
print(
	data.frame(model=c('model #1','model #2'),
			   DIC  =paste(round(c(dic1,dic2),3)))
)

# parmean <- c()
# for(i in 1:length(post)){
	# if(length(dim(post[[i]])) > 1) {
		# parmean <- c(parmean,colMeans(post[[i]]))}else{
			# parmean <- c(parmean,mean(post[[i]]))}
# }


# beta0_1 <- colMeans(post1$beta0)
# beta1_1 <- colMeans(post1$beta1)
# beta0mean_1 <- mean(post1$beta0mean)
# beta0sd_1 <- mean(post1$beta0_sd)
# betaT_1 <- mean(post1$betaT)
# beta1sd_1 <- mean(post1$beta1_sd)

# beta0_2 <- colMeans(post2$beta0)
# beta1_2 <- colMeans(post2$beta1)
# beta1mean_2 <- mean(post2$beta1mean)
# beta0sd_2 <- mean(post2$beta0_sd)
# betaT_2 <- mean(post2$betaT)
# beta1sd_1 <- mean(post2$beta1_sd)

# ypred1=ypred2 <- numeric(0)
# par(mfrow=c(2,5),mar=c(2,2,2,2)) #2 rows of panels, 3 columns
# for(i in 1:p){ #loop through locations  
	# plot(y[(ni[i]+1):ni[i+1]],x[(ni[i]+1):ni[i+1]],ylim=rev(range(x[(ni[i]+1):ni[i+1]])),xlim=range(y)) #plot like a depth profile 
	# mtext(T[i])
	
	# lines(beta0_1[i] + beta1_1[i]*x[(ni[i]+1):ni[i+1]], x[(ni[i]+1):ni[i+1]],ylim=rev(range(x[(ni[i]+1):ni[i+1]])))
	# lines(beta0_2[i] + beta1_2[i]*x[(ni[i]+1):ni[i+1]], x[(ni[i]+1):ni[i+1]],ylim=rev(range(x[(ni[i]+1):ni[i+1]])),col='red')
	# ypred1 <- c(ypred1,beta0_1[i] + beta1_1[i]*x[(ni[i]+1):ni[i+1]])
	# ypred2 <- c(ypred2,beta0_2[i] + beta1_2[i]*x[(ni[i]+1):ni[i+1]])
# }

# plot(y,ypred1)
# plot(y,ypred2)


# fit_lm  <- optimizing(mod, data=dat, hessian=TRUE, 
	# init=list(beta0=colMeans(post$beta0),
			  # beta1=colMeans(post$beta1),
			  # beta0slp_npp=mean(post$beta0slp_npp),
			  # beta0int=mean(post$beta0int),
			  # beta0_sd=mean(post$beta0_sd),
			  # beta1=colMeans(post$beta1),
			  # beta1slp_T=mean(post$beta1slp_T),
			  # beta1slp_o=mean(post$beta1slp_o),
			  # beta1int=mean(post$beta1int),
			  # beta1_sd=mean(post$beta1_sd),
			  # sigma=mean(post$sigma)))
			  
# parhat <- fit_lm$par

# I       <- -fit_lm$hessian #Fisher information matrix
# S       <- solve(I) #variance covariance matrix (variances on diagonal)
# sds     <- sqrt(diag(S)) #standard deviations as the square root of the variances

# fit_lm2  <- optimizing(mod2, data=dat, hessian=TRUE, draws=0, init='random') #optimize log posterior of second model return Hessian
# I2       <- -fit_lm2$hessian #Fisher information matrix
# S2       <- solve(I2) #variance covariance matrix (variance on diagonal)
# sds2     <- sqrt(diag(S2)) #standard deviations as square root of variances

# ev1    <- fit_lm$value + log(abs(det(fit_lm$hessian))) #second-order approximation to evidence
# BIC1   <- fit_lm$value + length(fit_lm$par)*log(N) #BIC (asymptotic approximation)

# ev2    <- fit_lm2$value + log(abs(det(fit_lm2$hessian))) #second-order approximation to evidence
# BIC2   <- fit_lm2$value + length(fit_lm2$par)*log(N) #BIC (asymptotic approximation)

# data.frame(BIC=c(BIC1,BIC2),evidence_approx=c(ev1,ev2),row.names=c('mod1','mod2')) #display the results in a table

######################################################
# Analyze Stan output ################################
######################################################
# post <- extract(fit)   #extract samples
# post_lm <- fit_lm$par
# beta0_lm <- post_lm[1:10]
# beta0sds <- sds[1:10]
# beta1_lm <- post_lm[13:23]
# beta1sds <- sds[13:23]

# par(mfrow=c(4,5),mar=c(1,1,1,1),oma=c(1,1,1,1))
# for(i in 1:p){
	# hist(post$beta0[,i],main='',xlab=i)
	# abline(v=beta0t[i],lwd=2); abline(v=quantile(post$beta0[,i],c(0.025,0.975)),lty=2)

	# abline(v=beta0_lm[i],lwd=2,col='red'); #abline(v=quantile(post$beta0[,i],c(0.025,0.975)),lty=2)
	# abline(v=beta0_lm[i]+2*c(beta0sds[i],-beta0sds[i]),lwd=2,col='red',lty=2); #abline(v=quantile(post$beta0[,i],c(0.025,0.975)),lty=2)
# }

# for(i in 1:p){
	# hist(post$beta1[,i],main=''); 
	# abline(v=-beta1t[i],lwd=2); abline(v=quantile(post$beta1[,i],c(0.025,0.975)),lty=2)

	# abline(v=beta1_lm[i],lwd=2,col='red'); #abline(v=quantile(post$beta0[,i],c(0.025,0.975)),lty=2)
	# abline(v=beta1_lm[i]+2*c(beta1sds[i],-beta1sds[i]), lwd=2,col='red',lty=2)
# }	



#######################################################
## 






