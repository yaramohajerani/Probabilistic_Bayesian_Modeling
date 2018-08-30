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
d0 = pd.read_table(os.path.join(indata,'pom_flux','GO_flux.tab'),sep='\t',skiprows=85,header=0)


#-- overwrite the column names
d0.columns = ['id_ref','id','type','lat','lon','trap','depth_bathy','depth','start','end','time','area','flux_pom',\
'flux_pom_sd','flux_c','flux_c_sd','flux_poc','flux_poc_sd','flux_pic','flux_pic_sd','flux_pon','flux_pon_sd',\
'flux_pop','flux_pop_sd','flux_psi','flux_psi_sd','flux_psio2','flux_psioh4','flux_pai','flux_pai_sd','flux_chl',\
'flux_pheop','flux_caco3','flux_caco3_sd','flux_fe','flux_fe_sd','flux_ba','flux_det','flux_ti','ref','url','url2']

#-- remove rows where the depth or desired fluxes are missing
ind1 = np.squeeze(np.nonzero(~np.isnan(d0['flux_poc'])))
ind2 = np.squeeze(np.nonzero(~np.isnan(d0['flux_pon'])))
ind3 = np.squeeze(np.nonzero(~np.isnan(d0['flux_pop'])))
#-- get intersection
ind = list(set(ind1) & set(ind2) & set(ind3))

d = {}
for c in d0.columns:
	d[c] = d0[c][ind].values

print('Reduced rows to common indices.')

n = len(d[c])
print('number of data points = %i'%n)
#-- extract date from the sting format (1991-07-18T00:00:00)
years = np.array([d['start'][i][:4] for i in range(n)],dtype=np.int)
months = np.array([(d['start'][i][5:7]) for i in range(n)],dtype=np.int)
days = np.array([(d['start'][i][8:10]) for i in range(n)],dtype=np.int)

logList = {}
logList['flux_poc'] = np.log(d['flux_poc']) #log of particulate organic carbon flux
logList['flux_pon'] = np.log(d['flux_pon']); #log of particulate organic nitrogen flux
logList['flux_pop'] = np.log(d['flux_pop']); #log of particulate organic phosphorus flux
logList['depth'] = np.log(d['depth']) #log of depth

####################################
## Plots of the data ###############
####################################
##--POC,PON,POP profiles: arithmetic and log scale--##
f,ax = plt.subplots(2,3, figsize=(12,8))
#-- original data
for i,v,s in zip([0,1,2],['flux_poc','flux_pon','flux_pop'],[1,12,106]):
	ax[0,i].plot(d[v],d['depth'],'k.')
	#ax[0,i].set_xlim([0,500./s])
	ax[0,i].set_ylim(ax[0,i].get_ylim()[::-1])
	ax[0,i].set_title(v)
	ax[0,i].set_ylabel('depth')
#-- log data
for i,v in zip([0,1,2],['flux_poc','flux_pon','flux_pop']):
	ax[1,i].plot(d[v],d['depth'],'k.')
	#ax[1,i].set_xlim([-7,6])
	ax[1,i].set_ylim(ax[1,i].get_ylim()[::-1])
	ax[1,i].set_title('log(%s)'%v)
	ax[1,i].set_ylabel('depth')
	ax[1,i].set_xlabel(v)
plt.tight_layout()
plt.savefig(os.path.join(outdata,'all_station_data.pdf'),format='pdf')
plt.close(f)
print('Plotted data.')

#########################################
## Package the data for Stan ############
#########################################
dat = dict(N=len(logList['flux_poc']),x=logList['depth'], y=logList['flux_poc'])#,zpred=np.arange(501),Npred=501)

#######################################################
## Fit Stan model #####################################
#######################################################
#-- First check if the compiled file exists. If not, compile model.
compiled_file = os.path.join(ddir,'linreg_pomflux.pkl')
if os.path.isfile(compiled_file):
	mod = pickle.load(open(compiled_file, 'rb'))
else:
	mod = pystan.StanModel(os.path.join(stan_dir,'linreg_pomflux.stan')) #pre-compile

	# save it to the file 'model.pkl' for later use
	with open(compiled_file, 'wb') as f:
	    pickle.dump(mod, f)

print('Compiled Model.')

fit = mod.sampling(data=dat, iter=2000, chains=4, warmup=1000) #fit model
print('Fit model.')

#######################################################
## Analyze Stan output ################################
#######################################################
post = fit.extract(permuted=True)   #extract samples
print(post.keys())        	    #contains lists of samples for posterior of beta0, beta1, sigma, lp

##--Model Parameters--##
post_b   =  post['beta1'] #extract posterior samples for 'b'
post_f0  = np.exp(post['beta0']) #extract posterior samples for intercept
post_sigma =  post['sigma'] #extract posterior samples for sigma

f2,axarr = plt.subplots(3,1, figsize=(8,8))
#-- beta1
yhist, xhist, _ = axarr[0].hist(post['beta1'])
axarr[0].set_title('beta1')
yline = np.arange(np.max(yhist))
axarr[0].plot(np.ones(len(yline))*np.mean(post['beta1']),yline,'k-',linewidth=2)
axarr[0].plot(np.ones(len(yline))*np.percentile(post['beta1'],2.5),yline,'k--')
axarr[0].plot(np.ones(len(yline))*np.percentile(post['beta1'],97.5),yline,'k--')
#-- beta0 mean
yhist, xhist, _  = axarr[1].hist(post_f0)
axarr[1].set_title('exp(beta0)')
yline = np.arange(np.max(yhist))
axarr[1].plot(np.ones(len(yline))*np.mean(post_f0),yline,'k-',linewidth=2)
axarr[1].plot(np.ones(len(yline))*np.percentile(post_f0,2.5),yline,'k--')
axarr[1].plot(np.ones(len(yline))*np.percentile(post_f0,97.5),yline,'k--')
#-- sigma
yhist, xhist, _  = axarr[2].hist(post_sigma)
axarr[2].set_title('post_sigma')
yline = np.arange(np.max(yhist))
axarr[2].plot(np.ones(len(yline))*np.mean(post_sigma),yline,'k-',linewidth=2)
axarr[2].plot(np.ones(len(yline))*np.percentile(post_sigma,2.5),yline,'k--')
axarr[2].plot(np.ones(len(yline))*np.percentile(post_sigma,97.5),yline,'k--')
plt.tight_layout()
plt.savefig(os.path.join(outdata,'all_stations_histogram.pdf'),format='pdf')
plt.close(f2)

##--Profiles--##
zs = np.arange(1,1001) #input depth variables for plot
fpred = np.zeros((len(zs),len(post_f0)))
for i in range(len(zs)):
	fpred[i,:] = post_f0[:]*(zs[i]/100.)**post_b[:]


f2,axarr = plt.subplots(2,1, figsize=(6,8))
ind_list = [np.arange(1000),np.arange(100,1000)]
for i,j in zip([0,1],ind_list):
	axarr[i].plot(np.mean(fpred[j],axis=1),zs[j],'k-')
	axarr[i].plot(np.percentile(fpred[j],2.5,axis=1),zs[j],color='orange',linestyle='--')
	axarr[i].plot(np.percentile(fpred[j],97.5,axis=1),zs[j],color='orange',linestyle='--')
	axarr[i].set_ylim(axarr[i].get_ylim()[::-1])
	axarr[i].set_ylabel('Depth')
	axarr[i].set_xlabel('fpred')

plt.tight_layout()
plt.savefig(os.path.join(outdata,'all_station_profiles.pdf'),format='pdf')
plt.close(f2)
