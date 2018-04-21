##################################################################################
## 
rm(list=ls())

#d <- read.csv('d:/dropbox/working/stan_bayesian_modeling/data/indata.dir/pom_flux/GO_flux_env_o2_t.csv')
d <- read.csv('d:/dropbox/teaching/stan_examples/data/GO_flux_env_o2_t.csv')
d <- d[complete.cases(d$flux_poc),]

#-TRY SUBSETTING <30 days-#
d <- d[d$time < 30,]

d$id2 <- paste(d$id,d$month,sep='')
d$id3 <- paste(d$id,d$month,d$year,sep='')

#################################################################################
## ID: LOCATION #################################################################
#################################################################################
IDs <- unique(d$id)
X <- data.frame(b=numeric(),logJ0=numeric(),n=numeric(),id=numeric(),
				sst=numeric(),t_avg=numeric(),npp=numeric(),o2_avg=numeric(),chl=numeric(),kd=numeric(),par=numeric(),
				lat=numeric())

pdf('d:/dropbox/working/stan_bayesian_modeling/greg/plots/ind_regressions_pic.pdf')
par(mfrow=c(5,5),mar=c(1,1,1,1),cex.axis=0.7,oma=c(3,3,3,3))
	for(i in 1:length(IDs)){
		dtmp <- d[d$id==IDs[i],]
		if(nrow(dtmp)>2){if(sd(dtmp$depth)>10){
			plot(log(dtmp$depth/100), log(dtmp$flux_poc),ylim=c(-3,7),xlim=c(-1,4))
			fit <- lm(log(dtmp$flux_poc) ~ log(dtmp$depth/100))
			abline(fit)
			X <- rbind(X,data.frame(b    =summary(fit)$coefficients[2,1],
								   logJ0 =summary(fit)$coefficients[1,1],
								   id    =IDs[i],
								   sst   =mean(dtmp$sst,na.rm=TRUE),
								   t_avg =mean(dtmp$t_avg,na.rm=TRUE),
								   npp   =mean(dtmp$npp,na.rm=TRUE),
								   o2_avg=mean(dtmp$o2_avg,na.rm=TRUE),
								   chl   =mean(dtmp$chl,na.rm=TRUE),
								   kd    =mean(dtmp$kd,na.rm=TRUE),
								   par   =mean(dtmp$par,na.rm=TRUE),
								   n     =nrow(dtmp),
								   lat   =dtmp$lat[1]))
			}
		}
	}
dev.off()

par(mfrow=c(2,2),mar=c(1,1,2,1))
hist(X$b,breaks=100,main=''); mtext('bs')
	abline(v=mean(X$b),lwd=2); legend('topright',legend=paste('bavg = ',round(mean(bs),3)),bty='n')
hist(X$b,breaks=100,xlim=c(-5,5),main=''); mtext('bs')
	abline(v=mean(X$b),lwd=2); 
hist(X$logJ0,breaks=100,main=''); mtext('log J0')
	abline(v=mean(X$logJ0),lwd=2); legend('topleft',legend=paste('logJ0avg = ',round(mean(X$logJ0),3)),bty='n')
hist(X$logJ0,breaks=100,xlim=c(-10,10),main=''); mtext('log J0')
	abline(v=mean(X$logJ0),lwd=2); 
mtext('Location ID',outer=TRUE)
#################################################################################
## ID2s: LOCATION/MONTH #########################################################
#################################################################################
ID2s <- unique(d$id2)
X2 <- data.frame(b=numeric(),logJ0=numeric(),n=numeric(),id=numeric(),
				sst=numeric(),t_avg=numeric(),npp=numeric(),o2_avg=numeric(),chl=numeric(),kd=numeric(),par=numeric(),
				month=numeric(),lat=numeric())

pdf('d:/dropbox/working/stan_bayesian_modeling/greg/plots/ind_regressions_id2_pic.pdf')
par(mfrow=c(5,5),mar=c(1,1,1,1),cex.axis=0.7,oma=c(3,3,3,3))
	for(i in 1:length(ID2s)){

	dtmp <- d[d$id2==ID2s[i],]
		if(nrow(dtmp)>2){if(diff(range(dtmp$depth))>100){
			plot(log(dtmp$depth/100), log(dtmp$flux_poc),ylim=c(-3,7),xlim=c(-1,4))
			fit <- lm(log(dtmp$flux_poc) ~ log(dtmp$depth/100))
			abline(fit)
			# X2 <- rbind(X2,data.frame(b    =summary(fit)$coefficients[2,1],
								   # logJ0 =summary(fit)$coefficients[1,1],
								   # id    =IDs[i],
								   # sst   =mean(dtmp$sst,na.rm=TRUE),
								   # t_avg =mean(dtmp$t_avg,na.rm=TRUE),
								   # npp   =mean(dtmp$npp,na.rm=TRUE),
								   # o2_avg=mean(dtmp$o2_avg,na.rm=TRUE),
								   # chl   =mean(dtmp$chl,na.rm=TRUE),
								   # kd    =mean(dtmp$kd,na.rm=TRUE),
								   # par   =mean(dtmp$par,na.rm=TRUE),
								   # n     =nrow(dtmp),
								   # month =dtmp$month[1],
								   # lat   =dtmp$lat[1]))
			}
		}
	}
dev.off()

par(mfrow=c(2,2),mar=c(1,1,2,1))
hist(X2$b,breaks=100,main=''); mtext('bs')
	abline(v=mean(X2$b),lwd=2); legend('topright',legend=paste('bavg = ',round(mean(X2$b),3)),bty='n')
hist(X2$b,breaks=100,xlim=c(-5,5),main=''); mtext('bs')
	abline(v=mean(X2$b),lwd=2); 
hist(X2$logJ0,breaks=100,main=''); mtext('log J0')
	abline(v=mean(X2$logJ0),lwd=2); legend('topleft',legend=paste('logJ0avg = ',round(mean(X2$logJ0),3)),bty='n')
hist(X2$logJ0,breaks=100,xlim=c(-10,10),main=''); mtext('log J0')
	abline(v=mean(X2$logJ0),lwd=2); 
mtext('Location ID x Month',outer=TRUE)
#################################################################################
## ID2s: LOCATION/MONTH/YEAR ####################################################
#################################################################################
ID3s <- unique(d$id3)
X3 <- data.frame(b=numeric(),logJ0=numeric(),n=numeric(),id=numeric(),
				sst=numeric(),t_avg=numeric(),npp=numeric(),o2_avg=numeric(),chl=numeric(),kd=numeric(),par=numeric(),
				month=numeric(),year=numeric(),lat=numeric())

pdf('d:/dropbox/working/stan_bayesian_modeling/greg/plots/ind_regressions_id3_pic.pdf')
par(mfrow=c(5,5),mar=c(1,1,1,1),cex.axis=0.7,oma=c(3,3,3,3))
	for(i in 1:length(ID3s)){
		dtmp <- d[d$id3==ID3s[i],]
		dtmp <- dtmp[complete.cases(dtmp$flux_pic),]
		if(nrow(dtmp)> 2){if(diff(range(dtmp$depth))>500){
			#plot(log(dtmp$depth/100), log(dtmp$flux_pic + 0.0001),ylim=c(-3,7),xlim=c(-1,4))
			plot(dtmp$depth, dtmp$flux_pic,xlim=c(0,6000))

			#fit <- lm(log(dtmp$flux_pic + 0.0001) ~ log(dtmp$depth/100))
			fit <- lm(dtmp$flux_pic ~ dtmp$depth)
			abline(fit)
			# X3 <- rbind(X3,data.frame(b    =summary(fit)$coefficients[2,1],
								   # logJ0 =summary(fit)$coefficients[1,1],
								   # id    =IDs[i],
								   # sst   =mean(dtmp$sst,na.rm=TRUE),
								   # t_avg =mean(dtmp$t_avg,na.rm=TRUE),
								   # npp   =mean(dtmp$npp,na.rm=TRUE),
								   # o2_avg=mean(dtmp$o2_avg,na.rm=TRUE),
								   # chl   =mean(dtmp$chl,na.rm=TRUE),
								   # kd    =mean(dtmp$kd,na.rm=TRUE),
								   # par   =mean(dtmp$par,na.rm=TRUE),
								   # n     =nrow(dtmp),
								   # month =dtmp$month[1],
								   # year  =dtmp$year[1],
								   # lat   =dtmp$lat[1]))
			}
		}
	}
dev.off()

par(mfrow=c(2,2),mar=c(1,1,2,1))
hist(X3$b,breaks=100,main=''); mtext('bs')
	abline(v=mean(X3$b),lwd=2); legend('topright',legend=paste('bavg = ',round(mean(X3$b),3)),bty='n')
hist(X3$b,breaks=100,xlim=c(-5,5),main=''); mtext('bs')
	abline(v=mean(X3$b),lwd=2); 
hist(X3$logJ0,breaks=100,main=''); mtext('log J0')
	abline(v=mean(X3$logJ0),lwd=2); legend('topleft',legend=paste('logJ0avg = ',round(mean(X3$logJ0),3)),bty='n')
hist(X3$logJ0,breaks=100,xlim=c(-10,10),main=''); mtext('log J0')
	abline(v=mean(X3$logJ0),lwd=2); 
mtext('Location ID x Month',outer=TRUE)

################################################################################################
### SCATTERPLOT MATRICES #######################################################################
################################################################################################
X <- X[X$b<3 & X$b>-3,]
X$logkd <- log(X$kd)
X$logchl <- log(X$chl)
X <- X[,!names(X)=='n' & !names(X)=='id']
splom(X,upper.panel=NULL,
	lower.panel=function(x, y, ...) {
          panel.hexbinplot(x, y, ...)
          fit <- lm(y ~ x)
          panel.abline(fit,lty=3,lwd=1.5,col='red')
      }	,
	diag.panel = function(x, ...){
		yrng <- current.panel.limits()$ylim
		d <- density(x, na.rm=TRUE)
		d$y <- with(d, yrng[1] + 0.95 * diff(yrng) * y / max(y) )
		panel.lines(d)
		diag.panel.splom(x, ...)
	},pscale=0, varname.cex=0.7,trans=log,colramp=BTC)

	
	
X2 <- X2[X2$b<3 & X2$b>-3,]
X2$logkd <- log(X2$kd)
X2$logchl <- log(X2$chl)
X2 <- X2[,!names(X2)=='n' & !names(X2)=='id']
splom(X2,upper.panel=NULL,
	lower.panel=function(x, y, ...) {
          panel.hexbinplot(x, y, ...)
          fit <- lm(y ~ x)
          panel.abline(fit,lty=3,lwd=1.5,col='red')
      }	,
	diag.panel = function(x, ...){
		yrng <- current.panel.limits()$ylim
		d <- density(x, na.rm=TRUE)
		d$y <- with(d, yrng[1] + 0.95 * diff(yrng) * y / max(y) )
		panel.lines(d)
		diag.panel.splom(x, ...)
	},pscale=0, varname.cex=0.7,trans=log,colramp=BTC)


X3 <- X3[X3$b<3 & X3$b>-3,]
X3$logkd <- log(X3$kd)
X3$logchl <- log(X3$chl)
X3 <- X3[,!names(X3)=='n' & !names(X3)=='id']
splom(X3,upper.panel=NULL,
	lower.panel=function(x, y, ...) {
          panel.hexbinplot(x, y, ...)
          fit <- lm(y ~ x)
          panel.abline(fit,lty=3,lwd=1.5,col='red')
      }	,
	diag.panel = function(x, ...){
		yrng <- current.panel.limits()$ylim
		d <- density(x, na.rm=TRUE)
		d$y <- with(d, yrng[1] + 0.95 * diff(yrng) * y / max(y) )
		panel.lines(d)
		diag.panel.splom(x, ...)
	},pscale=0, varname.cex=0.7,trans=log,colramp=BTC)

#################################################################################
## REGRESSIONS ##################################################################
#################################################################################
library(car)

XX <- X3[complete.cases(X3),]

summary(lm(b ~ sst, data=XX))
summary(lm(b ~ t_avg, data=XX))
summary(lm(b ~ npp, data=XX))
summary(lm(b ~ o2_avg, data=XX))

fit <- lm(b ~ as.factor(month), data=XX)

coeffs <- as.numeric(summary(fit)$coefficients[,1])




fit <- lm(b ~ sst + t_avg + npp + o2_avg + par + kd + chl + logchl + logkd + as.factor(month), data=XX)
summary(step(fit,direction='both',k=log(nrow(XX))))

summary(lm(b ~ t_avg + kd, data=XX))


fit <- lm(logJ0 ~ sst + t_avg + npp + o2_avg + par + kd + chl + logchl + logkd, data=XX)
summary(step(fit,direction='both',k=log(nrow(XX))*9))


summary(lm(logJ0s ~ ssts))
summary(lm(logJ0s ~ t_avgs))
summary(lm(logJ0s ~ npps))
summary(lm(logJ0s ~ log(npps)))

summary(lm(logJ0s ~ o2_avgs))
summary(lm(logJ0s ~ chls))
summary(lm(logJ0s ~ pars))
summary(lm(logJ0s ~ kds))

fit <- lm(logJ0s ~ ssts + t_avgs + npps + o2_avgs + pars + kds + chls, data=XJ0)
summary(step(fit,direction='both',k=log(4000))


summary(lm(bs ~ ssts,weights=sqrt(ns)))
summary(lm(bs ~ t_avgs,weights=sqrt(ns)))
summary(lm(bs ~ npps,weights=sqrt(ns)))
summary(lm(bs ~ o2_avgs,weights=sqrt(ns)))
fit <- lm(bs ~ ssts + t_avgs + npps + o2_avgs + pars + kds + chls, data=X,weights=sqrt(ns))
	
