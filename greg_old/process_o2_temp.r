
library(ncdf4)
library(colorRamps)
library(rworldmap)

######################################################################################
##-PROCESS OXYGEN-####################################################################
######################################################################################
setwd('d:/dropbox/data/woa/monthly/oxygen/one_deg/')

o2nc <- nc_open('woa13_all_o01_01.nc')
z    <- ncvar_get(o2nc,'depth_bnds')[1,]
nc_close(o2nc)

files <- list.files()
O2 <- array(NA,dim=c(360,180,57,12))
for(i in 1:12){
	o2nc       <- nc_open(files[i])
	for(j in 1:57){
		O2[,,j,i]   <- ncvar_get(o2nc,'o_an')[,,j]
	}
	nc_close(o2nc)
}

##-DOWNSCALE FROM ONE DEGREE-#########################
# O2 <- array(NA,dim=c(90,45,57,12))
# for(i in 1:12){
	# o2nc       <- nc_open(files[i])
	# for(j in 1:57){
		# O2[,,j,i]   <- resize_bilinear(z=ncvar_get(o2nc,'o_an')[,,j],xin=360,yin=180,xout=90,yout=45)
	# }
	# nc_close(o2nc)
# }

####################################################################################
##-PROCESS TEMPERATURE-#############################################################
####################################################################################

setwd('d:/dropbox/data/woa/monthly/temperature/one_deg_nc/')

Tnc <- nc_open('woa13_decav_t01_01v2.nc')
z    <- ncvar_get(o2nc,'depth_bnds')[1,]
nc_close(Tnc)

files <- list.files()
T <- array(NA,dim=c(360,180,57,12))
for(i in 1:12){
	Tnc       <- nc_open(files[i])
	for(j in 1:57){
		T[,,j,i]   <- ncvar_get(Tnc,'t_an')[,,j]
	}
	nc_close(Tnc)
}

##-DOWNSCALE FROM ONE DEGREE-#########################
# T <- array(NA,dim=c(90,45,57,12))
# for(i in 1:12){
	# Tnc       <- nc_open(files[i])
	# for(j in 1:57){
		# T[,,j,i]   <- resize_bilinear(z=ncvar_get(Tnc,'t_an')[,,j],xin=360,yin=180,xout=90,yout=45)
	# }
	# nc_close(Tnc)
# }

###################################################################################
##-PLOTS-##########################################################################
###################################################################################
# lats <- seq(-90,90,length.out=45)
# lons <- seq(-180,180,length.out=90)

#-OXYGEN-########
# pdf('d:/dropbox/teaching/stan_examples/greg/plots/o2_depths_map.pdf')
# par(mfrow=c(4,3),mar=c(2,2,2,2),cex.axis=0.8,oma=c(2,2,4,2))
# for(j in 1:12){
	# for(i in seq(1,length(z),5)){
		# image.plot(lons,lats,O2[,,i,j],zlim=c(0,12),xlab='',ylab='',col=matlab.like(20))
		# box(lwd=1); title('')
		# mtext(paste('z = ',z[i],'m'),cex=0.8)
	# }
	# mtext(month.name[j],outer=TRUE)
# }
# dev.off()
 
#-TEMPERATURE-#########
# pdf('d:/dropbox/teaching/stan_examples/greg/plots/T_depths_map.pdf')
# par(mfrow=c(4,3),mar=c(2,2,2,2),cex.axis=0.8,oma=c(2,2,4,2))
# for(j in 1:12){
	# for(i in seq(1,length(z),5)){
		# image.plot(lons,lats,T[,,i,j],zlim=c(-2,30),xlab='',ylab='',col=matlab.like(20))
		# box(lwd=1); title('')
		# mtext(paste('z = ',z[i],'m'),cex=0.8)
	# }
	# mtext(month.name[j],outer=TRUE)
# }
# dev.off()
 
 
 
 
