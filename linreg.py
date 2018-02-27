#!/usr/bin/env python
u"""
linreg.py
by Yara Moherani

Update history:
	02/2018	Update Directory strucutre, automatically check if precompiled model exists.
	01/2018	Written
"""
import os
import pystan
import pickle
import numpy as np
import matplotlib.pyplot as plt

#-- directory setup
ddir = os.path.dirname(os.path.realpath(__file__))
stan_dir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'stan'))

#######################################################
## Setup your data ####################################
#######################################################
N = 100 #number of synthetic data points
beta0t = 2 #true value of beta0 for synthetic data
beta1t = 1 #true value of beta1 for synthetic data
sigmat = 1 #true value for measurement error of synthetic data

x = np.random.normal(0, 1, N) #generate normally distributed independent variables
y = beta0t + beta1t*x + np.random.normal(0,sigmat,100) #generate the linear regression according to parameters
fig1 = plt.figure(1)
plt.scatter(x,y) #plot synthetic data

dat = dict(N=N, x=x, y=y) #package data into format for Rstan

#######################################################
## Fit Stan model #####################################
#######################################################
#-- First check if the compiled file exists. If not, compile model.
compiled_file = os.path.join(ddir,'linreg_compiled.pkl')
if os.path.isfile(compiled_file):
	mod = pickle.load(open(compiled_file, 'rb'))
else:
	mod = pystan.StanModel(os.path.join(stan_dir,'linreg.stan')) #pre-compile

	# save it to the file 'model.pkl' for later use
	with open(os.path.join(ddir,'linreg_compiled.pkl'), 'wb') as f:
	    pickle.dump(mod, f)

fit = mod.sampling(data=dat, iter=2000, chains=4, warmup=1000) #fit model

print(fit)

#######################################################
## Analyze Stan output ################################
#######################################################
post = fit.extract(permuted=True)   #extract samples
print(post.keys())        	    #contains lists of samples for posterior of beta0, beta1, sigma, lp

#--Characterize inference on the parameters--#
f2, axarr = plt.subplots(3, 1, figsize=(8,8))
axarr[0].hist(post['beta0'])
axarr[0].set_title('beta0')
axarr[1].hist(post['beta1'])
axarr[1].set_title('beta1')
axarr[2].hist(post['sigma'])
axarr[2].set_title('sigma')

#--Characterize the fit to the data--#
j = np.argsort(x) #index of x in order or increasing x
b0p_mean = np.mean(post['beta0']) #mean of posterior samples for b0
b1p_mean = np.mean(post['beta1']) #mean of posterior samples for b1
sp_mean  = np.mean(post['sigma']) #mean of posterior samples for the standard deviation of the errors

#Option #1: the distribution of y
fig3 = plt.figure(3)
plt.scatter(x,y)
plt.plot(x[j],b0p_mean + b1p_mean*x[j], color = 'orange') #line of best fit
plt.plot(x[j],b0p_mean + b1p_mean*x[j] + sp_mean, color='r') #line of best fit plus measurement error
plt.plot(x[j],b0p_mean + b1p_mean*x[j] - sp_mean, color='r') #line of best fit minus measurement error
plt.title('the distribution of y')

#Option #2: the 'true value' of y given x; analytical
hypf = ((x[j] - np.mean(x))**2)/np.sum((x[j]-np.mean(x))**2) #factor that enters analytical formula for 'prediction interval'

fig4 = plt.figure(4)
plt.scatter(x,y) #plot data
plt.plot(x[j],b0p_mean + b1p_mean*x[j], color = 'orange') #line of best fit
plt.plot(x[j],b0p_mean + b1p_mean*x[j] + np.sqrt((sp_mean**2)*(1/N + hypf)),color='r') #analytical upper confidence interval
plt.plot(x[j],b0p_mean + b1p_mean*x[j] - np.sqrt((sp_mean**2)*(1/N + hypf)),color='r') #analytical lower confidence interval
plt.title("the 'true value' of y given x; analytical")

#Option #3: true value via MCMC
fig5 = plt.figure(5)
plt.scatter(x,y) #plot data
plt.scatter(x[j],np.mean(post['y_pred'],axis=0)[j],color='orange') #plot mean of posterior for y(t)
plt.plot(x[j],np.mean(post['y_pred'],axis=0)[j]+np.std(post['y_pred'],axis=0)[j],color='r')
plt.plot(x[j],np.mean(post['y_pred'],axis=0)[j]-np.std(post['y_pred'],axis=0)[j],color='r') 
plt.title('True value via MCMC')

plt.show()
