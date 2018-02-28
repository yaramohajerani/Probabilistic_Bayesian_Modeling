#! /usr/bin/env python
u"""
linreg_pomflux.py
by Yara Mohajerani

Update History
	02/2018	Written
"""
import os
import pystan
import pickle
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#-- directory setup
#- current directory
ddir = os.path.dirname(os.path.realpath(__file__))
#- stan code directory
stan_dir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'stan'))
#- data input
indata = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '../..', 'stan_data.dir/indata.dir'))
#- data output
outdata = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '../..', 'stan_data.dir/outdata.dir'))

#######################################################
## Setup your data ####################################
#######################################################
d = pd.read_table(os.path.join(indata,'pom_flux','GO_flux.tab'),sep='\t',skiprows=85,header=0)


#-- overwrite the column names
d.columns = ['id_ref','id','type','lat','lon','trap','depth_bathy','depth','start','end','time','area','flux_pom',\
'flux_pom_sd','flux_c','flux_c_sd','flux_poc','flux_poc_sd','flux_pic','flux_pic_sd','flux_pon','flux_pon_sd',\
'flux_pop','flux_pop_sd','flux_psi','flux_psi_sd','flux_psio2','flux_psioh4','flux_pai','flux_pai_sd','flux_chl',\
'flux_pheop','flux_caco3','flux_caco3_sd','flux_fe','flux_fe_sd','flux_ba','flux_det','flux_ti','ref','url','url2'] 

print d['start'][0]
print d['start'][0]
print d['start'][0]

#-- remove rows where the depth or desired fluxes are missing

print(d['flux_poc'].notnull())
d = d[d['flux_poc'].notnull()]
d = d[d['flux_pon'].notnull()]
d = d[d['flux_pop'].notnull()]
d = d[d['depth'].notnull()]
#-- number of filtered observations
n = len(d)
print d['start'][0]
#-- extract date from the sting format (1991-07-18T00:00:00)
print np.array([d['start'][i][:4] for i in range(5)],dtype=np.int)
#years = np.array([d['start'][i:i+1][:4] for i in range(n)],dtype=np.int)
#months = np.array([(d['start'][i:i+1][5:7]) for i in range(n)],dtype=np.int)
#days = np.array([(d['start'][i:i+1][8:10]) for i in range(n)],dtype=np.int)

#print years
print n
"""
log_fpoc <- log(d$flux_poc) #log of particulate organic carbon flux
log_fpon <- log(d$flux_pon); log_fpon[!is.finite(log_fpon)] <- NA #log of particulate organic nitrogen flux
log_fpop <- log(d$flux_pop); log_fpop[!is.finite(log_fpop)] <- NA #log of particulate organic phosphorus flux
log_z    <- log(d$depth) #log of depth
"""

