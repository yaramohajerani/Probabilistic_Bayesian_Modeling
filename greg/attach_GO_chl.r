######################################################################################
##-ATTACH GO_chl to GO_flux ##########################################################
######################################################################################
library(lubridate)

#-DFLUX-####################################################
dflux           <- read.delim('c:/users/greg/dropbox/teaching/stan_examples/data/GO_flux.tab',header=TRUE,skip=85) #read in the data; data start at line 86
#dflux           <- read.delim('d:/dropbox/teaching/stan_examples/data/GO_flux.tab',header=TRUE,skip=85,stringsAsFactor=FALSE) #read in the data; data start at line 86
colnames(dflux) <- c('id_ref','id','type','lat','lon','trap','depth_bathy','depth','start','end','time','area','flux_pom','flux_pom_sd','flux_c','flux_c_sd','flux_poc','flux_poc_sd','flux_pic','flux_pic_sd','flux_pon','flux_pon_sd','flux_pop','flux_pop_sd','flux_psi','flux_psi_sd','flux_psio2','flux_psioh4','flux_pai','flux_pai_sd','flux_chl','flux_pheop','flux_caco3','flux_caco3_sd','flux_fe','flux_fe_sd','flux_ba','flux_det','flux_ti','ref','url','url2') #write over the crazy variable names in the original data
dflux$year      <- substr(as.character(dflux$start),1,4); 
dflux$month     <- substr(as.character(dflux$start),6,7); 
dflux$day       <- substr(as.character(dflux$start),9,10) #extract date data from format given
dflux$dyear     <- decimal_date(ymd(paste(dflux$year,dflux$month,dflux$day,sep='-')))
dflux$id_mon    <- paste(dflux$id,dflux$month,sep='-')
dflux$id_mon_yr <- paste(dflux$id,dflux$month,dflux$year,sep='-')
#-REMOVE traps deployed for more than a month-#
dflux <- dflux[dflux$time < 32,]

#-DCHL-#####################################################
#dchl           <- read.delim('c:/users/greg/dropbox/teaching/stan_examples/data/GO_chl.tab',header=TRUE,skip=18)
dchl           <- read.delim('d:/dropbox/teaching/stan_examples/data/GO_chl.tab',header=TRUE,skip=18)
colnames(dchl) <- c('id','lat','lon','date','chl','npp','sst','kd','zeu','par')
dchl$year      <- substr(as.character(dchl$date),1,4); 
dchl$month     <- substr(as.character(dchl$date),6,7); 
dchl$day       <- substr(as.character(dchl$date),9,10) #extract date data from format given
dchl$dyear     <- decimal_date(ymd(paste(dchl$year,dchl$month,dchl$day,sep='-')))
dchl$id_mon    <- paste(dchl$id,dchl$month,sep='-') 
dchl$id_mon_yr <- paste(dchl$id,dchl$month,dchl$year,sep='-') 

#-Z-########################################################
o2nc <- nc_open('d:/dropbox/data/woa/monthly/oxygen/one_deg/woa13_all_o01_01.nc') #just to extract the depth variable, processing of temperature and oxygen is done via a script
z    <- ncvar_get(o2nc,'depth_bnds')[1,]
nc_close(o2nc)
dz <- c(diff(z),50)

#-PROCESS VARIABLES-#########################################
#source('d:/dropbox/working/stan_bayesian_modeling/greg/process_o2_temp.r') #run script to process T and O2 matrices; gives 180x360x57 arrays
#source('d:/dropbox/working/stan_bayesian_modeling/greg/process_PO4_cellsize.r')

lats <- seq(-90,90,length.out=180)
lons <- seq(-180,180,length.out=360)

###################################################################
## ID #############################################################
###################################################################
IDS <- unique(dflux$id)
for(i in 1:length(IDS)){
	
	d_tmp <- dflux[dflux$id==IDS[i],]
	if(nrow(d_tmp)>1){
		lat   <- d_tmp$lat[1]
		lon   <- d_tmp$lon[1]
		depth <- max(d_tmp$depth)
		
		lati  <- which((lats-lat)^2 == min((lats-lat)^2))[1]
		loni  <- which((lons-lon)^2 == min((lons-lon)^2))[1]
		
		dchl_tmp <- dchl[dchl$id==d_tmp$id[1],]
		t  <- rowMeans(T[loni,lati,,],na.rm=TRUE)
		o2 <- rowMeans(O2[loni,lati,,],na.rm=TRUE) 
		
		dflux$npp_id[dflux$id==IDS[i]]        <- mean(dchl_tmp$npp,na.rm=TRUE)
		dflux$sst_id[dflux$id==IDS[i]]        <- mean(dchl_tmp$sst,na.rm=TRUE)
		dflux$t_avg1000_id[dflux$id==IDS[i]]  <- sum(t[z<1000]*(dz[z<1000]/sum(dz[z<1000],na.rm=TRUE)),na.rm=TRUE)
		dflux$o2_avg1000_id[dflux$id==IDS[i]] <- sum(o2[z<1000]*(dz[z<1000]/sum(dz[z<1000],na.rm=TRUE)),na.rm=TRUE)
		dflux$t_avgZ_id[dflux$id==IDS[i]]     <- sum(t*(dz/sum(dz,na.rm=TRUE)),na.rm=TRUE)
		dflux$o2_avgZ_id[dflux$id==IDS[i]]    <- sum(o2*(dz/sum(dz,na.rm=TRUE)),na.rm=TRUE)
		dflux$PO4_id[dflux$id==IDS[i]]        <- mean(P[loni,lati,],na.rm=TRUE)
		dflux$nano_id[dflux$id==IDS[i]]       <- mean(NANO[lati,loni,],na.rm=TRUE)
		dflux$micro_id[dflux$id==IDS[i]]      <- mean(MICRO[lati,loni,],na.rm=TRUE)
		dflux$pico_id[dflux$id==IDS[i]]       <- mean(PICO[lati,loni,],na.rm=TRUE)
	}
}	

###################################################################
## IDxMONTH #######################################################
###################################################################
ID2S <- unique(dflux$id_mon)
for(i in 1:length(ID2S)){
	
	d_tmp <- dflux[dflux$id_mon==ID2S[i],]
	if(nrow(d_tmp)>1){
		lat   <- d_tmp$lat[1]
		lon   <- d_tmp$lon[1]
		depth <- max(d_tmp$depth)
		month <- as.numeric(d_tmp$month)[1]
		monthc<- d_tmp$month[1]
		
		lati  <- which((lats-lat)^2 == min((lats-lat)^2))[1]
		loni  <- which((lons-lon)^2 == min((lons-lon)^2))[1]
		
		dchl_tmp <- dchl[dchl$month==monthc & dchl$id==d_tmp$id[1],]
		t  <- T[loni,lati,,month]
		o2 <- O2[loni,lati,,month] 
		
		dflux$npp_mon[dflux$id_mon==ID2S[i]]        <- mean(dchl_tmp$npp,na.rm=TRUE)
		dflux$sst_mon[dflux$id_mon==ID2S[i]]        <- mean(dchl_tmp$sst,na.rm=TRUE)
		dflux$t_avg1000_mon[dflux$id_mon==ID2S[i]]  <- sum(t[z<1000]*(dz[z<1000]/sum(dz[z<1000],na.rm=TRUE)),na.rm=TRUE)
		dflux$o2_avg1000_mon[dflux$id_mon==ID2S[i]] <- sum(o2[z<1000]*(dz[z<1000]/sum(dz[z<1000],na.rm=TRUE)),na.rm=TRUE)
		dflux$t_avgZ_mon[dflux$id_mon==ID2S[i]]     <- sum(t*(dz/sum(dz,na.rm=TRUE)),na.rm=TRUE)
		dflux$o2_avgZ_mon[dflux$id_mon==ID2S[i]]    <- sum(o2*(dz/sum(dz,na.rm=TRUE)),na.rm=TRUE)
		dflux$PO4_mon[dflux$id_mon==ID2S[i]]        <- P[loni,lati,month]
		dflux$nano_mon[dflux$id_mon==ID2S[i]]       <- NANO[lati,loni,month]
		dflux$micro_mon[dflux$id_mon==ID2S[i]]      <- MICRO[lati,loni,month]
		dflux$pico_mon[dflux$id_mon==ID2S[i]]       <- PICO[lati,loni,month]
	}
}	

###################################################################
## IDxMONTHxYEAR ##################################################
###################################################################
ID3S <- unique(dflux$id_mon_yr)
for(i in 1:length(ID3S)){
	d_tmp <- dflux[dflux$id_mon_yr==ID3S[i],]
	if(nrow(d_tmp)>1){
		lat   <- d_tmp$lat[1]
		lon   <- d_tmp$lon[1]
		depth <- max(d_tmp$depth)
		month <- as.numeric(d_tmp$month)[1]
		monthc<- d_tmp$month[1]
		year  <- as.numeric(d_tmp$year)[1]
		yearc <- d_tmp$year[1]
		
		lati  <- which((lats-lat)^2 == min((lats-lat)^2))[1]
		loni  <- which((lons-lon)^2 == min((lons-lon)^2))[1]
		
		dchl_tmp <- dchl[dchl$month==monthc & dchl$id==d_tmp$id[1] & dchl$year==year,]

		dflux$npp_mon_yr[dflux$id_mon_yr==ID3S[i]]        <- mean(dchl_tmp$npp,na.rm=TRUE)
		dflux$sst_mon_yr[dflux$id_mon_yr==ID3S[i]]        <- mean(dchl_tmp$sst,na.rm=TRUE)	
	}
}	

#####################################################################
## WRITE CSV ########################################################
#####################################################################
write.csv(dflux,file='d:/dropbox/working/stan_bayesian_modeling/greg/GO_flux_ENV.csv')

