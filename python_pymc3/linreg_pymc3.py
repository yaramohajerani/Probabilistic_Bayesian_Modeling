#!/usr/bin/env python
u"""
linreg_pymc3.py
by Yara Moherani

same as linreg.py and linreg.stan but implemented in pymc3

Update history:
	04/2018	Written
"""
import pymc3 as pm
import numpy as np
import matplotlib.pyplot as plt
from theano import shared

#######################################################
## Setup your data ####################################
#######################################################
N = 100 #number of synthetic data points
beta0t = 2 #true value of beta0 for syntetic data
beta1t = 1 #true value of beta1 for synthetic data
sigmat = 1 #true value for measurement error of synthetic data

x = np.random.normal(0, 1, N) #generate normally distributed independent variables
y = beta0t + beta1t*x + np.random.normal(0,sigmat,N) #generate the linear regression according to parameters
fig1 = plt.figure(1)
plt.scatter(x,y) #plot synthetic data

#######################################################
#-- set up model and perform inteference
#######################################################
with pm.Model() as model:
	# Define priors
	sigma = pm.HalfCauchy('sigma', beta=10, testval=1.)
	beta0 = pm.Normal('beta0', 10, sd=20)
	beta1 = pm.Normal('beta1', 0, sd=20)

	# Define likelihood
	pm.Normal('y', mu=beta1*x + beta0,sd=sigma, observed=y)

	# Inference
	post = pm.sample(3000, cores=2) # draw 3000 posterior samples
	#-- also generate prediction values (200 datasets containing N points)
	predictions = pm.sample_ppc(post, model=model)
	#-- extract
	y_pred = predictions['y']


#######################################################
#-- Get the posterior and analyze output
#######################################################
print(pm.summary(post))

#fig2, ax = plt.subplots(3, 2, )
pm.traceplot(post[100:])
plt.tight_layout();

#fig3 = plt.figure(3)
pm.plot_posterior(post)
plt.tight_layout();


#-- extract parameters
#--Characterize the fit to the data--#
j = np.argsort(x) #index of x in order or increasing x
b0p_mean = np.mean(post['beta0']) #mean of posterior samples for b0
b1p_mean = np.mean(post['beta1']) #mean of posterior samples for b1
sp_mean  = np.mean(post['sigma']) #mean of posterior samples for the standard deviation of the errors

#Option #1: the distribution of y
fig3 = plt.figure(4)
plt.scatter(x,y)
plt.plot(x[j],b0p_mean + b1p_mean*x[j], color = 'orange') #line of best fit
plt.plot(x[j],b0p_mean + b1p_mean*x[j] + sp_mean, color='r') #line of best fit plus measurement error
plt.plot(x[j],b0p_mean + b1p_mean*x[j] - sp_mean, color='r') #line of best fit minus measurement error
plt.title('the distribution of y')
plt.tight_layout()

#Option #2: the 'true value' of y given x; analytical
hypf = ((x[j] - np.mean(x))**2)/np.sum((x[j]-np.mean(x))**2) #factor that enters analytical formula for 'prediction interval'

fig4 = plt.figure(5)
plt.scatter(x,y) #plot data
plt.plot(x[j],b0p_mean + b1p_mean*x[j], color = 'orange') #line of best fit
plt.plot(x[j],b0p_mean + b1p_mean*x[j] + np.sqrt((sp_mean**2)*(1/N + hypf)),color='r') #analytical upper confidence interval
plt.plot(x[j],b0p_mean + b1p_mean*x[j] - np.sqrt((sp_mean**2)*(1/N + hypf)),color='r') #analytical lower confidence interval
plt.title("the 'true value' of y given x; analytical")
plt.tight_layout()

#Option #3: true value via MCMC
fig5 = plt.figure(6)
plt.scatter(x,y) #plot data
plt.scatter(x[j],np.mean(y_pred,axis=0)[j],color='orange') #plot mean of posterior for y(t)
plt.plot(x[j],np.mean(y_pred,axis=0)[j]+np.std(y_pred,axis=0)[j],color='r')
plt.plot(x[j],np.mean(y_pred,axis=0)[j]-np.std(y_pred,axis=0)[j],color='r')
plt.title('True value via MCMC')
plt.tight_layout()

plt.show()
