library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

p  <- 10 #number of locations
n  <- sample(20:100,p,replace=TRUE) #data points per location
N  <- sum(n) #total number of data points
ni <- c(0,cumsum(n)) #vector index for dataset stacked as single column
y  <- numeric(N) #allocate space for y
x  <- numeric(N) #allocate space for x
T  <- 1:p

betaT    <- 0.05 #slope of profile slope vs. temperature
beta1sdt <- 0.01 #noise around the line
beta1t   <- betaT*T + rnorm(p,sd=beta1sdt) #slope as a function of temperature plus noise
beta0sdt <- 10 #standard deviation of intercepts for synthetic data
beta0mean<- 10 #mean of intercepts for synthetic data
beta0t   <- rnorm(p,beta0mean,beta0sdt) #slope for synthetic data
sigmat   <- 10 #standard deviation for the measurement error

for(i in 1:p){ #loop through locations
	xtmp <- runif(n[i],1,500) #synthesize independent variable
	y[(ni[i]+1):ni[i+1]] <- beta0t[i] - beta1t[i]*xtmp + rnorm(n[i],sd=sigmat) #generate station-specific data as an elongating vector
    x[(ni[i]+1):ni[i+1]] <- xtmp #store the independent variable in the big vector
}
#############################
## Make a plot ##############
#############################
par(mfrow=c(2,5),mar=c(2,2,2,2)) #2 rows of panels, 3 columns
for(i in 1:p){ #loop through locations  
	plot(y[(ni[i]+1):ni[i+1]],x[(ni[i]+1):ni[i+1]],ylim=rev(range(x[(ni[i]+1):ni[i+1]])),xlim=range(y)) #plot like a depth profile 
	mtext(T[i])
}

dat <- list(N=N,ni=ni,y=y,p=p,T=T,NI=0,NS=1,n=n,x=x,M0=cbind(rep(1,p)),M1=cbind(rep(1,p),T)) #package data for rstan

mod <- stan_model('d:/dropbox/working/stan_bayesian_modeling/stan/linreg_pomflux_env_log_lik.stan')
fit <- sampling(mod,data=dat,iter=2000,warmup=1000)


