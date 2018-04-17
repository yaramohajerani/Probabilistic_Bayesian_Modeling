#!/usr/bin/env python
u"""
linreg_edward.py
by Yara Moherani

same as linreg.py and linreg.stan but implemented in edward

Update history:
	04/2018	Written
"""
import edward as ed
from edward.models import Normal,Empirical,Uniform
import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt

#######################################################
## Setup your data ####################################
#######################################################
N = 100 #number of synthetic data points
T = 500 #-- number of posterior samples
D = 1 #-- number of dimensions ("features")
beta0t = 2 #true value of beta0 for syntetic data
beta1t = 1 #true value of beta1 for synthetic data
sigmat = 1 #true value for measurement error of synthetic data

x = np.random.normal(0, 1, N) #generate normally distributed independent variables
y = beta0t + beta1t*x + np.random.normal(0,sigmat,N) #generate the linear regression according to parameters
fig1 = plt.figure(1)
plt.scatter(x,y) #plot synthetic data

#######################################################
#-- Make placeholders for inference and set up model
#######################################################
x_holder = tf.placeholder(tf.float32, [N,D])
#-- priors on slope and intercept
beta0 = Normal(loc=tf.ones(1)*10, scale=tf.ones(1)*100)
beta1 = Normal(loc=tf.zeros(D), scale=tf.ones(D)*100)
sigma = Uniform(low=tf.zeros(N),high=tf.ones(N)*10.)
y_model = Normal(loc= ed.dot(x_holder,beta1) + beta0, scale=sigma)

#######################################################
#-- Inference
#######################################################
qbeta1 = Empirical(params=tf.get_variable("qbeta0/params", [T,D]))
qbeta0 = Empirical(params=tf.get_variable("qbeta1/params", [T,1]))
qsigma = Empirical(params=tf.get_variable("qsigma/params", [T,N]))
#-- Stochastic gradient Hamiltonian Monte Carlo
inference = ed.SGHMC({beta0: qbeta0, beta1: qbeta1, sigma: qsigma},
    data={x_holder: x.reshape((N,D)), y_model: y})
inference.run(step_size=1e-3)

#######################################################
#-- Get the posterior and analyze output
#######################################################
y_post = ed.copy(y_model, {beta1: qbeta1, beta0: qbeta0, sigma: qsigma})
print('Mean Square Error: %.4f'%\
    ed.evaluate('mean_squared_error', data={x_holder: x.reshape((N,D)), y_model: y}))

#-- convert posterior tensors to numpy arrays
n_post = 10
post = {}
post['beta1'] = qbeta1.sample(n_post).eval()
post['beta0'] = qbeta0.sample(n_post).eval()
post['sigma'] = qsigma.sample(n_post).eval()
post['y_pred'] = np.zeros((n_post,N))
for i in range(n_post):
    post['y_pred'][i,:] =  post['beta1'][i]*x + post['beta0'][i]

# post = {}
# post['beta1'] = beta1.eval()
# post['beta0'] = beta0.eval()
# post['sigma'] = sigma.eval()
# post['y_pred'] = np.zeros((T,N))
# for i in range(T):
#     post['y_pred'][i,:] =  post['beta1'][i]*x + post['beta0'][i]

print beta1.eval()
print qbeta1.eval()

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
