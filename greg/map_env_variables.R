
library(colorRamps)
library(rworldmap)

d <- read.csv('d:/dropbox/teaching/stan_examples/data/GO_flux_env_o2_t.csv')


d$colsst <- matlab.like(10)[as.numeric(cut(d$sst,breaks = 10))]
d$colnpp <- matlab.like(10)[as.numeric(cut(log(d$npp),breaks = 10))]
d$colchl <- matlab.like(10)[as.numeric(cut(log(d$chl),breaks = 10))]
d$colo2_avg <- matlab.like(10)[as.numeric(cut(d$o2_avg,breaks = 10))]
d$colkd     <- matlab.like(10)[as.numeric(cut(log(d$kd),breaks=10))]
d$colt_avg  <- matlab.like(10)[as.numeric(cut(d$t_avg,breaks=10))]


X <- data.frame(sst=d$sst,log_npp=log(d$npp),log_chl=log(d$chl),log_kd=log(d$kd),o2_avg=d$o2_avg,t_avg=d$t_avg)
pairs(X)



pdf('d:/dropbox/teaching/stan_examples/greg/plots/obs_maps.pdf',height=4,width=7)
map()
	points(d$lon,d$lat,col=d$colsst,pch=19,cex=0.6)
	image.plot(matrix(seq(min(d$sst,na.rm=TRUE),max(d$sst,na.rm=TRUE))),legend.only=TRUE,legend.mar=1,col=matlab.like(10))	
	title('');mtext('SST [degree C]')
map()
	points(d$lon,d$lat,col=d$colnpp,pch=19,cex=0.6)
	image.plot(matrix(seq(min(log(d$npp),na.rm=TRUE),max(log(d$npp),na.rm=TRUE))),legend.only=TRUE,legend.mar=1,col=matlab.like(10))	
	title('');mtext('log (NPP [mgC/m2/day]')
map()
	points(d$lon,d$lat,col=d$colchl,pch=19,cex=0.6)
	image.plot(matrix(seq(min(log(d$chl),na.rm=TRUE),max(log(d$chl),na.rm=TRUE))),legend.only=TRUE,legend.mar=1,col=matlab.like(10))	
	title('');mtext('log (CHL [mgCHL/m2]')
map()
	points(d$lon,d$lat,col=d$colkd,pch=19,cex=0.6)
	image.plot(matrix(seq(min(log(d$kd),na.rm=TRUE),max(log(d$kd),na.rm=TRUE))),legend.only=TRUE,legend.mar=1,col=matlab.like(10))	
	title('');mtext('log(attenuation coefficient [/m])')	
map()
	points(d$lon,d$lat,col=d$colo2_avg,pch=19,cex=0.6)
	image.plot(matrix(seq(min(d$o2_avg,na.rm=TRUE),max(d$o2_avg,na.rm=TRUE))),legend.only=TRUE,legend.mar=1,col=matlab.like(10))	
	title('');mtext('O2 (depth average) [mmol/m3]')
map()
	points(d$lon,d$lat,col=d$colt_avg,pch=19,cex=0.6)
	image.plot(matrix(seq(min(d$t_avg,na.rm=TRUE),max(d$t_avg,na.rm=TRUE))),legend.only=TRUE,legend.mar=1,col=matlab.like(10))	
	title('');mtext('temperature (depth average) [degree C]')
dev.off()


