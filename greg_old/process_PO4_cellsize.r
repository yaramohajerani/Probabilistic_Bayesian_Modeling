library(ncdf4)
source('d:/dropbox/code/functions/resize_bilinear().R')
########################################################################
## Surface phosphate ###################################################
########################################################################
setwd('d:/dropbox/DATA/WOA/monthly/phosphate/')
files <- grep('.nc',list.files(),value=TRUE)

P <- array(NA,dim=c(360,180,12))

for(i in 1:12){
	nc     <- nc_open(files[1])
	p      <- ncvar_get(nc,'p_an')
	P[,,i] <- apply(p[,,1:3],c(1,2),mean)
	nc_close(nc)
}

########################################################################
## Cell size ###########################################################
########################################################################
setwd('d:/dropbox/DATA/kostadinov_et_al_2016/')
files <- list.files()

PICO=NANO=MICRO <- array(NA,dim=c(180,360,12))

for(i in c(1:10,10,12)){ #file #11 is corrupt, have to download the new one, for now replace month with the last month
	nc <- nc_open(files[i])
	pico_bio <- ncvar_get(nc,'C_biomass_picoplankton')
	nano_bio <- ncvar_get(nc,'C_biomass_nanoplankton')
	micro_bio<- ncvar_get(nc,'C_biomass_microplankton')
	total_bio<- pico_bio + nano_bio + micro_bio
	
	PICO[,,i]  <- resize_bilinear(xin=2160,yin=4320,xout=180,yout=360,z=pico_bio/total_bio)
	NANO[,,i]  <- resize_bilinear(xin=2160,yin=4320,xout=180,yout=360,z=nano_bio/total_bio)
	MICRO[,,i] <- resize_bilinear(xin=2160,yin=4320,xout=180,yout=360,z=micro_bio/total_bio)
}

