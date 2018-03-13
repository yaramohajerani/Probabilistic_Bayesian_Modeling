#!/usr/bin/env python
u"""
linreg_pomflux_ragged_env.py
by Yara Mohajerani (03/2018)
"""
import os
import sys
import pystan
import pickle
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import getopt

#-- directory setup
#- current directory
ddir = os.path.dirname(os.path.realpath(__file__))
#- stan code directory
stan_dir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'stan'))
#- data input
indata = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '../..', 'stan_data.dir/indata.dir'))
#- data output
outdata = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '../..', 'stan_data.dir/outdata.dir'))


def fit_var(var,min_stn,max_stn,niter,nchains,nwarm,PLOT):
    print('looking at stations %i to %i.'%(min_stn,max_stn))
    print('Fit configurations: iterations=%i , chain=%i, warmup=%i'%(niter,nchains,nwarm))

    #######################################################
    ## Setup your data ####################################
    #######################################################
    d0 = pd.read_csv(os.path.join(indata,'pom_flux','GO_flux_env.csv'))

    #-- remove rows where the depth or desired fluxes are missing
    ind1 = np.squeeze(np.nonzero(~np.isnan(d0['flux_poc'])))
    ind2 = np.squeeze(np.nonzero(~np.isnan(d0[var])))
    #-- get intersection
    ind = list(set(ind1) & set(ind2))

    #-- remove nans and convert to dictionarties with python arrays (easier to work with)
    d = {}
    for c in ['id','depth','flux_poc',var]:
    	d[c] = d0[c][ind].values

    #-- only look at the first 10 stations to beging with
    indstn = []
    for i in range(min_stn,max_stn+1):
        if np.count_nonzero(d['id']==i) > 1:
            indstn += list(np.squeeze(np.nonzero(d['id']==i)))


    #-- Set up data for stan. we want each station to be separate
    #-- number of unique stations
    station_nums = np.unique(d['id'][indstn])
    p = len(station_nums)
    #-- number of data points per station
    n = np.zeros(p)
    for i in range(p):
        #-- only count stations that have more than 1 value
        if np.count_nonzero(d['id'][indstn]==station_nums[i]) > 1:
            n[i] = np.count_nonzero(d['id'][indstn]==station_nums[i])
    N = np.int(np.sum(n)) #-- total number of points

    #vector index for dataset stacked as single column
    ni = np.concatenate([np.array([0]),np.cumsum(n)])
    ni = np.array(ni,dtype=np.int)

    #-- average all the environmental values for each station
    var_avg = np.zeros(p)
    for i in range(p):
        var_avg[i] = np.mean(d[var][indstn][ni[i]:ni[i+1]])

    #-- package data for Stan
    dat = dict(N=N, ni=ni, x = d['depth'][indstn], y=d['flux_poc'][indstn], p=p, v=var_avg)

    #######################################################
    ## Fit Stan model #####################################
    #######################################################
    #-- First check if the compiled file exists. If not, compile model.
    compiled_file = os.path.join(ddir,'linreg_pomflux_ragged_env.pkl')
    if os.path.isfile(compiled_file):
    	mod = pickle.load(open(compiled_file, 'rb'))
    else:
    	mod = pystan.StanModel(os.path.join(stan_dir,'linreg_pomflux_ragged_env.stan')) #pre-compile

    	# save it to the file 'model.pkl' for later use
    	with open(compiled_file, 'wb') as f:
    	    pickle.dump(mod, f)

    print('Compiled Model.')

    fit = mod.sampling(data=dat, iter=niter, chains=nchains, warmup=nwarm) #fit model
    print('Fit model.')

    #######################################################
    ## Analyze Stan output ################################
    #######################################################
    post = fit.extract(permuted=True)   #extract samples
    print(post.keys())        	    #contains lists of samples for posterior of beta0, beta1, sigma, lp

    if PLOT in ['y','Y']:
        f1, axarr = plt.subplots(p, 1, figsize=(8,8))
        f1.suptitle('beta0')
        for i in range(p):
            yhist, xhist, _ = axarr[i].hist(post['beta0'][:,i])
            #-- make y axis for plotting vertical lines
            yline = np.arange(np.max(yhist))
            axarr[i].plot(np.ones(len(yline))*post['beta0mean'][i],yline,'k-',linewidth=2)
            axarr[i].plot(np.ones(len(yline))*np.percentile(post['beta0'][:,i],2.5),yline,'k--')
            axarr[i].plot(np.ones(len(yline))*np.percentile(post['beta0'][:,i],97.5),yline,'k--')
            axarr[i].set_title = 'beta0 (%i)'%i
        plt.tight_layout()
        plt.savefig(os.path.join(outdata,'%s_beta0_histogram_stn%i-%i_%iiter_%ichains_%iwarmup.pdf'\
            %(var,min_stn,max_stn,niter,nchains,nwarm)),format='pdf')
        plt.close(f1)

        f2, axarr = plt.subplots(p, 1, figsize=(8,8))
        f2.suptitle('beta1')
        for i in range(p):
            yhist, xhist, _ = axarr[i].hist(post['beta1'][:,i])
            #-- make y axis for plotting vertical lines
            yline = np.arange(np.max(yhist))
            axarr[i].plot(np.ones(len(yline))*np.mean(post['beta1'][:,i]),yline,'k-',linewidth=2)
            axarr[i].plot(np.ones(len(yline))*np.percentile(post['beta1'][:,i],2.5),yline,'k--')
            axarr[i].plot(np.ones(len(yline))*np.percentile(post['beta1'][:,i],97.5),yline,'k--')
            axarr[i].set_title = 'beta1 (%i)'%i
        plt.tight_layout()
        plt.savefig(os.path.join(outdata,'%s_beta1_histogram_stn%i-%i_%iiter_%ichains_%iwarmup.pdf'\
            %(var,min_stn,max_stn,niter,nchains,nwarm)),format='pdf')
        plt.close(f2)

    #-- write results to file
    f = open(os.path.join(outdata,'%s_stn%i-%i_%iiter_%ichains_%iwarmup.txt'\
        %(var,min_stn,max_stn,niter,nchains,nwarm)),'w')
    for i in range(p):
        f.write("Stn %i: slope = %.4f (%.4f-%.4f) ; intercept =  %.4f (%.4f-%.4f)\n"%(station_nums[i],\
            np.mean(post['beta1'][:,i]),np.percentile(post['beta1'][:,i],2.5),np.percentile(post['beta1'][:,i],97.5),\
            np.mean(post['beta0'][:,i]),np.percentile(post['beta0'][:,i],2.5),np.percentile(post['beta0'][:,i],97.5)))

#-- help function for usage
def usage_info():
    print("--help or -h to display help information.")
    print("--variable=X or -V:X to choose environmental variable for fit. Default npp.")
    print("--min_stn=X or -i:X to set intial station number X. Default 1.")
    print("--max_stn=X or -f:X to set intial station number X. Default 840.")
    print("--iter=X or -I:X to set number of interations X. Default 2000.")
    print("--chains=X or -C:X to set number of chains X. Default 4.")
    print("--warmup=X or -W:X to set number of warmup interations X. Default 1000.")
    print("--PLOT=Y or -P:Y to plot histograms. Otherwise set to 'N'. Default 'N'.")

#-- main function to get parameters and pass them along to fitting function
def main():
    #-- Read the system arguments listed after the program
    long_options = ['help','variable=','min_stn=','max_stn=','iter=','chains=','warmup=','PLOT=']
    optlist,arglist = getopt.getopt(sys.argv[1:],'hV:i:f:I:C:W:P:',long_options)

    #-- set defaults
    var = 'npp'
    min_stn = 1
    max_stn = 840
    niter = 2000
    nchains = 4
    nwarm = 1000
    PLOT = 'N'
    #-- set parameters
    for opt, arg in optlist:
        if opt in ("-h","--help"):
            usage_info()
            sys.exit()
        elif opt in ("-V","--variable"):
            var = arg
        elif opt in ("-i","--min_stn"):
            min_stn = np.int(arg)
        elif opt in ("-f","--max_stn"):
            max_stn = np.int(arg)
        elif opt in ("-I","--iter"):
            niter = np.int(arg)
        elif opt in ("-C","--chains"):
            nchains = np.int(arg)
        elif opt in ("-W","--warmup"):
            nwarm = np.int(arg)
        elif opt in ("-P","--PLOT"):
            PLOT = arg

    #-- pass parameters to fitting function
    fit_var(var,min_stn,max_stn,niter,nchains,nwarm,PLOT)

#-- run main program
if __name__ == '__main__':
        main()
