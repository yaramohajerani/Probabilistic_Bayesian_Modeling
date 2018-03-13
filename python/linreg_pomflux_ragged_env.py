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


def fit_var(var,env,min_stn,max_stn,niter,nchains,nwarm,PLOT,NP):
    print('looking at stations %i to %i.'%(min_stn,max_stn))
    print('Plotting: %s'%PLOT)
    if env == 0:
        print('Intercept is based on %s. Slope is sampled from normal distribution.'%var)
        env_config = 'intercept'
    elif env == 1:
        print('Slope is based on %s. Intercept is sampled from normal distribution.'%var)
        env_config = 'slope'
    print('Fit configurations: iterations=%i , chain=%i, warmup=%i'%(niter,nchains,nwarm))
    print('Number of parallel processes = %i'%NP)

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
    X = np.log(d['depth'][indstn]/100.)
    Y = np.log(d['flux_poc'][indstn])
    dat = dict(N=N, ni=ni, x = X, y= Y, p=p, v=var_avg)

    #######################################################
    ## Fit Stan model #####################################
    #######################################################
    #-- First check if the compiled file exists. If not, compile model.
    model_name = 'linreg_pomflux_ragged_env-%s'%env_config
    compiled_file = os.path.join(ddir,'%s.pkl'%model_name)
    if os.path.isfile(compiled_file):
    	mod = pickle.load(open(compiled_file, 'rb'))
    else:
    	mod = pystan.StanModel(os.path.join(stan_dir,'%s.stan'%model_name)) #pre-compile

    	# save it to the file 'model.pkl' for later use
    	with open(compiled_file, 'wb') as f:
    	    pickle.dump(mod, f)

    print('Compiled Model.')

    fit = mod.sampling(data=dat, iter=niter, chains=nchains, warmup=nwarm,n_jobs=NP) #fit model
    print('Fit model.')
    print(fit)

    #######################################################
    ## Analyze Stan output ################################
    #######################################################
    post = fit.extract(permuted=True)   #extract samples

    if PLOT in ['y','Y']:
        f1, axarr = plt.subplots(p, 1, figsize=(8,10))
        #f1.suptitle('beta0')
        for i in range(p):
            yhist, xhist, _ = axarr[i].hist(post['beta0'][:,i],bins=np.int(np.sqrt(len(post['beta0'][:,i])))\
                ,color='grey')
            #-- make y axis for plotting vertical lines
            yline = np.arange(np.max(yhist))
            if env==1:
                axarr[i].plot(np.ones(len(yline))*post['beta0mean'][i],yline,'r-',linewidth=2)
            axarr[i].plot(np.ones(len(yline))*np.mean(post['beta0'][:,i]),yline,'k-',linewidth=2)
            axarr[i].plot(np.ones(len(yline))*np.percentile(post['beta0'][:,i],2.5),yline,'k--')
            axarr[i].plot(np.ones(len(yline))*np.percentile(post['beta0'][:,i],97.5),yline,'k--')
            axarr[i].set_title('beta0, stn #%i'%station_nums[i])
        plt.tight_layout()
        plt.savefig(os.path.join(outdata,'%s_beta0_histogram_stn%i-%i_%iiter_%ichains_%iwarmup.pdf'\
            %(var,min_stn,max_stn,niter,nchains,nwarm)),format='pdf')
        plt.close(f1)

        f2, axarr = plt.subplots(p, 1, figsize=(8,10))
        #f2.suptitle('beta1')
        for i in range(p):
            yhist, xhist, _ = axarr[i].hist(post['beta1'][:,i],bins=np.int(np.sqrt(len(post['beta1'][:,i]))),\
                color='grey')
            #-- make y axis for plotting vertical lines
            yline = np.arange(np.max(yhist))
            if env==0:
                axarr[i].plot(np.ones(len(yline))*post['beta1mean'][i],yline,'r-',linewidth=2)
            axarr[i].plot(np.ones(len(yline))*np.mean(post['beta1'][:,i]),yline,'k-',linewidth=2)
            axarr[i].plot(np.ones(len(yline))*np.percentile(post['beta1'][:,i],2.5),yline,'k--')
            axarr[i].plot(np.ones(len(yline))*np.percentile(post['beta1'][:,i],97.5),yline,'k--')
            axarr[i].set_title('beta1, stn #%i'%station_nums[i])
        plt.tight_layout()
        plt.savefig(os.path.join(outdata,'%s_beta1_histogram_stn%i-%i_%iiter_%ichains_%iwarmup.pdf'\
            %(var,min_stn,max_stn,niter,nchains,nwarm)),format='pdf')
        plt.close(f2)

    #-- the environmental dependence histograms are always plotted (since the # of graphs doesn't
    #-- increase with the number of stations.)
    f3, axarr = plt.subplots(2, 1, figsize=(8,8))
    f3.suptitle("Dependence of %s on %s"%(env_config,var))
    for i in range(2):
        yhist, xhist, _ = axarr[i].hist(post['betaV%i'%i],bins=np.int(np.sqrt(len(post['betaV%i'%i]))),\
            color='grey')
        #-- make y axis for plotting vertical lines
        yline = np.arange(np.max(yhist))
        axarr[i].plot(np.ones(len(yline))*np.mean(post['betaV%i'%i]),yline,'k-',linewidth=2)
        axarr[i].plot(np.ones(len(yline))*np.percentile(post['betaV%i'%i],2.5),yline,'k--')
        axarr[i].plot(np.ones(len(yline))*np.percentile(post['betaV%i'%i],97.5),yline,'k--')
        axarr[i].set_title("betaV%i"%i)

    #plt.tight_layout()
    plt.savefig(os.path.join(outdata,'%s_environmental-dependence_histogram_stn%i-%i_%iiter_%ichains_%iwarmup.pdf'\
        %(var,min_stn,max_stn,niter,nchains,nwarm)),format='pdf')
    plt.close(f3)

    #-- write results to file
    f = open(os.path.join(outdata,'%s_stn%i-%i_%iiter_%ichains_%iwarmup.txt'\
        %(var,min_stn,max_stn,niter,nchains,nwarm)),'w')
    for i in range(p):
        f.write("Stn %i: slope = %.4f (%.4f : %.4f) std=%.4f; intercept =  %.4f (%.4f : %.4f) std=%.4f\n"%(station_nums[i],\
            np.mean(post['beta1'][:,i]),np.percentile(post['beta1'][:,i],2.5),np.percentile(post['beta1'][:,i],97.5),\
            np.std(post['beta1'][:,i]),np.mean(post['beta0'][:,i]),np.percentile(post['beta0'][:,i],2.5),\
            np.percentile(post['beta0'][:,i],97.5),np.std(post['beta0'][:,i])))

    f.write("\nSlope between %s and %ss = %.4f (%.4f : %.4f) std=%.4f\n"%(var,env_config,np.mean(post['betaV1'])\
        ,np.percentile(post['betaV1'],2.5),np.percentile(post['betaV1'],97.5),np.std(post['betaV1'])))
    f.write("Intercept between %s and %ss = %.4f (%.4f : %.4f) std=%.4f\n"%(var,env_config,np.mean(post['betaV0'])\
        ,np.percentile(post['betaV0'],2.5),np.percentile(post['betaV0'],97.5),np.std(post['betaV0'])))

    f.close()

#-- help function for usage
def usage_info():
    print("--help or -h to display help information.")
    print("--variable=X or -V:X to choose environmental variable for fit. Default npp.")
    print("--env=X or -E:X to choose if slope or intercept is based on environment.\n\t '0' for intercept and '1' for slope. Default = 0")
    print("--min_stn=X or -i:X to set intial station number X. Default 1.")
    print("--max_stn=X or -f:X to set intial station number X. Default 840.")
    print("--iter=X or -I:X to set number of interations X. Default 2000.")
    print("--chains=X or -C:X to set number of chains X. Default 4.")
    print("--warmup=X or -W:X to set number of warmup interations X. Default 1000.")
    print("--PLOT=Y or -P:Y to plot histograms. Otherwise set to 'N'. Default 'N'.")
    print("--NP=X or -n:X to set number of parallel jobs. Defaulr is the the number of chasins.")

#-- main function to get parameters and pass them along to fitting function
def main():
    #-- Read the system arguments listed after the program
    long_options = ['help','variable=','min_stn=','max_stn=','iter=','chains=','warmup=','PLOT=',\
        'env=','NP=']
    optlist,arglist = getopt.getopt(sys.argv[1:],'hV:i:f:I:C:W:P:E:n:',long_options)

    #-- set defaults
    var = 'npp'
    env = 0
    min_stn = 1
    max_stn = 840
    niter = 2000
    nchains = 4
    nwarm = 1000
    PLOT = 'N'
    NP = np.copy(nchains)
    #-- set parameters
    for opt, arg in optlist:
        if opt in ("-h","--help"):
            usage_info()
            sys.exit()
        elif opt in ("-V","--variable"):
            var = arg
        elif opt in ("-E","--env"):
            env = np.int(arg)
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
        elif opt in ("-n","--NP"):
            NP = np.int(arg)

    #-- pass parameters to fitting function
    fit_var(var,env,min_stn,max_stn,niter,nchains,nwarm,PLOT,NP)

#-- run main program
if __name__ == '__main__':
        main()
