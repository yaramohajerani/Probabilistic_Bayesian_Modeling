#! /usr/bin/env python
u"""
linreg_stations_ragged.py
by Yara Mohajerani

Update History
    02/2018 Written
"""
import numpy as np
import matplotlib.pyplot as plt
import pystan
import pickle
import os

#-- directory Setup
#-- stan directory
stan_dir = '/Users/yaramohajerani/Google Drive File Stream/Team Drives/Stan_class/code/stan'
#-- current directory
ddir = os.path.dirname(os.path.realpath(__file__))

#######################################################
## Setup your data ####################################
#######################################################
p  = 6 #number of locations
n  = np.random.choice(np.arange(100,1001), size=p, replace=True) #data points per location
N  = np.sum(n) #total number of data points
ni = np.concatenate([np.array([0]),np.cumsum(n)])
y  = np.zeros(N) #allocate space for y
x  = np.zeros(N) #allocate space for x
beta1t  = 0.1 #slope for synthetic data
beta0sdt = 50 #standard deviation of intercepts for synthetic data
beta0t = np.random.normal(loc=0,scale=beta0sdt,size=p) #slope for synthetic data
sigmat = 2 #standard deviation for the measurement error

for i in range(p): # loop through locations
    xtmp = np.arange(n[i]) #synthesize independent variable
    y[ni[i]:ni[i+1]] = beta0t[i] - beta1t*xtmp + np.random.normal(loc=0,scale=sigmat,size=n[i]) #generate station-specific data as an elongating vector
    x[ni[i]:ni[i+1]] = xtmp

#############################
## Make a plot ##############
#############################
#f, axarr = plt.subplots(2, 3)
f1, axarr = plt.subplots(p,figsize=(8,8))
for i in range(p):
    axarr[i].plot(y[ni[i]+1:ni[i+1]],x[ni[i]+1:ni[i+1]])
    axarr[i].set_ylim(axarr[i].get_ylim()[::-1])
plt.tight_layout()
plt.savefig(os.path.join(ddir,'station_data_points.pdf'),format='pdf')
plt.close(f1)

dat = dict(N=N, ni=ni, n=n, y=y, p=p, x=x) #package data for pystan
#######################################################
## Fit Stan model #####################################
#######################################################
#-- First check if the compiled file exists. If not, compile model.
compiled_file = os.path.join(ddir,'linreg_stations_ragged_compiled.pkl')
if os.path.isfile(compiled_file):
	mod = pickle.load(open(compiled_file, 'rb'))
else:
	mod = pystan.StanModel(os.path.join(stan_dir,'linreg_stations_ragged.stan')) #pre-compile

	# save it to the file 'model.pkl' for later use
	with open(compiled_file, 'wb') as f:
	    pickle.dump(mod, f)

fit = mod.sampling(data=dat, iter=2000, chains=4, warmup=1000) #fit model

#######################################################
## Analyze Stan output ################################
#######################################################
post = fit.extract(permuted=True)   #extract samples
print(post.keys())        	    #contains lists of samples for posterior of beta0, beta1, sigma, lp

print post['beta0'].shape
f2, axarr = plt.subplots(p, 1, figsize=(8,8))
f2.suptitle('beta0')
for i in range(p):
    yhist, xhist, _ = axarr[i].hist(post['beta0'][:,i])
    #-- make y axis for plotting vertical lines
    yline = np.arange(np.max(yhist))
    axarr[i].plot(np.ones(len(yline))*beta0t[i],yline,'k-',linewidth=2)
    axarr[i].plot(np.ones(len(yline))*np.percentile(post['beta0'][:,i],2.5),yline,'k--')
    axarr[i].plot(np.ones(len(yline))*np.percentile(post['beta0'][:,i],97.5),yline,'k--')
    axarr[i].set_title = 'beta0 (%i)'%i
plt.tight_layout()
plt.savefig(os.path.join(ddir,'beta0_histogram.pdf'),format='pdf')
plt.close(f2)

f3,axarr = plt.subplots(3,1, figsize=(6,8))
#-- beta1
yhist, xhist, _ = axarr[0].hist(post['beta1'])
axarr[0].set_title('beta1')
yline = np.arange(np.max(yhist))
axarr[0].plot(-np.ones(len(yline))*beta1t,yline,'k-',linewidth=2)
axarr[0].plot(np.ones(len(yline))*np.percentile(post['beta1'],2.5),yline,'k--')
axarr[0].plot(np.ones(len(yline))*np.percentile(post['beta1'],97.5),yline,'k--')
#-- beta0 mean
yhist, xhist, _  = axarr[1].hist(post['beta0mean'])
axarr[1].set_title('beta0mean')
yline = np.arange(np.max(yhist))
axarr[1].plot(np.zeros(len(yline)),yline,'k-',linewidth=2)
axarr[1].plot(np.ones(len(yline))*np.percentile(post['beta0mean'],2.5),yline,'k--')
axarr[1].plot(np.ones(len(yline))*np.percentile(post['beta0mean'],97.5),yline,'k--')
#-- beta0 standard deviation
axarr[2].hist(post['beta0_sd'])
axarr[2].set_title('beta0_sd')
plt.tight_layout()
plt.savefig(os.path.join(ddir,'beta1_beta0mean_histogram.pdf'),format='pdf')
plt.close(f3)
