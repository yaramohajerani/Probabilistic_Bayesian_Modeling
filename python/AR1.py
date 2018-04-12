#!/usr/bin/env python
u"""
AR1.py
by Yara Mohajerani
"""
import pickle
import pystan
import os
import matplotlib.pyplot as plt

#-- directory setup
ddir = os.path.dirname(os.path.realpath(__file__))
pkl_dir = os.path.join(ddir,'compiled_models.dir')
stan_dir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'stan'))

T = 100
y0 = 1
phi = 0.95
sigma = 1
dat = dict(T=T,y0=y0,phi=phi,sigma=sigma)


#######################################################
## Fit Stan model #####################################
#######################################################
#-- First check if the compiled file exists. If not, compile model.
compiled_file = os.path.join(pkl_dir,'ar1_compiled.pkl')
if os.path.isfile(compiled_file):
        mod = pickle.load(open(compiled_file, 'rb'))
else:
        mod = pystan.StanModel(os.path.join(stan_dir,'ar1.stan')) #pre-compile

        # save it to the file 'model.pkl' for later use
        with open(os.path.join(pkl_dir,'ar1_compiled.pkl'), 'wb') as f:
            pickle.dump(mod, f)

fit = mod.sampling(data=dat, iter=2000, chains=4, warmup=0,algorithm='Fixed_param')

post = fit.extract(permuted=True)   #extract samples
print post['y_hat'].shape
for i in range(6):
	plt.plot(post['y_hat'][i])
plt.show()
