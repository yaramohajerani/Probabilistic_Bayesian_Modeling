#!/usr/bin/env python
u"""
linreg_pomflux_ragged_env.py
by Yara Mohajerani (04/2018)

Update History
    04/2018   Debug and make change directory stucture so data is outside of the
                github repository folder.
    03/26/18  Add option to take log of input envionemntal variables
                ** NOTE ** the code assumes for now the log is the only
                option so if you add anything after the name it won't work
    03/22/18  Add option for a general number of env variables
                using linreg_pomflux_env.stan
    03/2018 Written
"""
from __future__ import print_function

import os
import sys
import pystan
import pickle
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import getopt
import copy

#-- directory setup
#- current directory
ddir = os.path.dirname(os.path.realpath(__file__))
#- stan code directory
stan_dir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'stan'))
#- data input
indata = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..','..','poc_data/indata.dir'))
#- data output
outdata = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..','..','poc_data/outdata.dir'))


def fit_var(parameters):
    slope_vars = parameters['slope_variables'].split(',')
    intercept_vars = parameters['intercept_variables'].split(',')
    min_stn = np.int(parameters['min_stn'])
    max_stn = np.int(parameters['max_stn'])
    niter = np.int(parameters['iterations'])
    nchains = np.int(parameters['chains'])
    nwarm = np.int(parameters['warmup'])
    PLOT = parameters['PLOT'].upper()
    NP = np.int(parameters['NP'])
    #-- if number of processes was never set, set it equal to # of chains
    if NP == 0:
        NP = np.copy(nchains)

    print('looking at stations %i to %i.'%(min_stn,max_stn))
    print('Plotting: %s'%PLOT)
    intercept_str = ''
    print('Intercept Dependence:')
    for v in intercept_vars:
        print(v)
        intercept_str += '%s,'%v
    slope_str = ''
    print('Slope Dependence:')
    for v in slope_vars:
        print(v)
        slope_str += '%s,'%v
    title_str = 'intercept_%s-dependent_slope_%s-dependent'%(intercept_str[:-1],slope_str[:-1])

    print('Fit configurations: iterations=%i , chain=%i, warmup=%i'%(niter,nchains,nwarm))
    print('Number of parallel processes = %i'%NP)

    #-- extract names and operations from the given variable names (e.g) lognpp --> log of npp
    intercept_varNames = copy.copy(intercept_vars)
    intercept_ops = copy.copy(intercept_vars)
    slope_varNames = copy.copy(slope_vars)
    slope_ops = copy.copy(slope_vars)
    for i,v in enumerate(intercept_vars):
        if v.startswith('log'):
            #-- for now we're assuming the only operation is log
            intercept_varNames[i] = v[3:]
            intercept_ops[i] = 'log'
        else:
            intercept_varNames[i] = v[:]
            intercept_ops[i] = '1'
    for i,v in enumerate(slope_vars):
        if v.startswith('log'):
            #-- for now we're assuming the only operation is log
            slope_varNames[i] = v[3:]
            slope_ops[i] = 'log'
        else:
            intercept_varNames[i] = v[:]
            slope_ops[i] = '1'

    #######################################################
    ## Setup your data ####################################
    #######################################################
    d0 = pd.read_csv(os.path.join(indata,'pom_flux','GO_flux_env.csv'))

    #-- remove rows where the depth or desired fluxes are missing
    #-- strat by just the poc indices and take intersections iteratively
    ind = np.squeeze(np.nonzero(~np.isnan(d0['flux_poc'])))
    for v in intercept_varNames:
        ind_temp = np.squeeze(np.nonzero(~np.isnan(d0[v])))
        ind = list(set(ind) & set(ind_temp))
    for v in slope_varNames:
        ind_temp = np.squeeze(np.nonzero(~np.isnan(d0[v])))
        ind = list(set(ind) & set(ind_temp))

    #-- remove nans and convert to dictionarties with python arrays (easier to work with)
    d = {}
    for v in ['id','depth','flux_poc']:
    	d[v] = d0[v][ind].values
    #-- also apply operations to the environmental variables
    for i,v in enumerate(intercept_varNames):
        if v != None:
            if intercept_ops[i] == 'log':
                d[v] = np.log(d0[v][ind].values)
            else:
                d[v] = d0[v][ind].values
    for i,v in enumerate(slope_varNames):
        if v != None:
            if slope_ops[i] == 'log':
                d[v] = np.log(d0[v][ind].values)
            else:
                d[v] = d0[v][ind].values

    #-- pick a subset of stations
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
    intercept_avg = {}
    for v in intercept_varNames:
        intercept_avg[v] = np.zeros(p)
        for i in range(p):
            intercept_avg[v][i] = np.mean(d[v][indstn][ni[i]:ni[i+1]])
    slope_avg = {}
    for v in slope_varNames:
        slope_avg[v] = np.zeros(p)
        for i in range(p):
            slope_avg[v][i] = np.mean(d[v][indstn][ni[i]:ni[i+1]])

    #-- package data for Stan
    X = np.log(d['depth'][indstn]/100.)
    Y = np.log(d['flux_poc'][indstn])
    NI = len(intercept_vars)
    NS = len(slope_vars)
    #-- make design matrices for slope and intercept
    M0 = np.ones((p,NI+1))
    M1 = np.ones((p,NS+1))
    for i,v in enumerate(intercept_varNames):
        M0[:,i+1] = intercept_avg[v][:]
    for i,v in enumerate(slope_varNames):
        M1[:,i+1] = slope_avg[v][:]

    dat = dict(p=p, N=N, NI=NI, NS=NS, ni=ni, x=X, y=Y, M0=M0, M1=M1)


    #######################################################
    ## Fit Stan model #####################################
    #######################################################
    #-- First check if the compiled file exists. If not, compile model.
    model_name = 'linreg_pomflux_env'
    compiled_file = os.path.join(ddir,'compiled_models.dir','%s.pkl'%model_name)
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
        #f1.suptitle('logJ0')
        for i in range(p):
            yhist, xhist, _ = axarr[i].hist(post['logJ0'][:,i],bins=np.int(np.sqrt(len(post['logJ0'][:,i])))\
                ,color='grey')
            #-- make y axis for plotting vertical lines
            yline = np.arange(np.max(yhist))
            axarr[i].plot(np.ones(len(yline))*np.mean(post['logJ0'][:,i]),yline,'k-',linewidth=2)
            axarr[i].plot(np.ones(len(yline))*np.percentile(post['logJ0'][:,i],2.5),yline,'k--')
            axarr[i].plot(np.ones(len(yline))*np.percentile(post['logJ0'][:,i],97.5),yline,'k--')
            axarr[i].set_title('logJ0, stn #%i'%station_nums[i])
        plt.tight_layout()
        plt.savefig(os.path.join(outdata,'%s_logJ0_histogram_stn%i-%i_%iiter_%ichains_%iwarmup.pdf'\
            %(title_str,min_stn,max_stn,niter,nchains,nwarm)),format='pdf')
        plt.close(f1)

        f2, axarr = plt.subplots(p, 1, figsize=(8,10))
        #f2.suptitle('b')
        for i in range(p):
            yhist, xhist, _ = axarr[i].hist(post['b'][:,i],bins=np.int(np.sqrt(len(post['b'][:,i]))),\
                color='grey')
            #-- make y axis for plotting vertical lines
            yline = np.arange(np.max(yhist))
            axarr[i].plot(np.ones(len(yline))*np.mean(post['b'][:,i]),yline,'k-',linewidth=2)
            axarr[i].plot(np.ones(len(yline))*np.percentile(post['b'][:,i],2.5),yline,'k--')
            axarr[i].plot(np.ones(len(yline))*np.percentile(post['b'][:,i],97.5),yline,'k--')
            axarr[i].set_title('b, stn #%i'%station_nums[i])
        plt.tight_layout()
        plt.savefig(os.path.join(outdata,'%s_b_histogram_stn%i-%i_%iiter_%ichains_%iwarmup.pdf'\
            %(title_str,min_stn,max_stn,niter,nchains,nwarm)),format='pdf')
        plt.close(f2)

    #-- the environmental dependence histograms are always plotted (since the # of graphs doesn't
    #-- increase with the number of stations.)
    f3, axarr = plt.subplots(NI+1, 1, figsize=(8,6))
    f3.suptitle("Dependence of Intercept on %s"%intercept_str)
    for i in range(NI+1):
        yhist, xhist, _ = axarr[i].hist(post['beta0'][:,i],bins=np.int(np.sqrt(len(post['beta0'][:,i]))),\
            color='grey')
        #-- make y axis for plotting vertical lines
        yline = np.arange(np.max(yhist))
        axarr[i].plot(np.ones(len(yline))*np.mean(post['beta0'][:,i]),yline,'k-',linewidth=2)
        axarr[i].plot(np.ones(len(yline))*np.percentile(post['beta0'][:,i],2.5),yline,'k--')
        axarr[i].plot(np.ones(len(yline))*np.percentile(post['beta0'][:,i],97.5),yline,'k--')
        if i ==0:
            axarr[i].set_title("Constant")
        else:
            axarr[i].set_title(intercept_vars[i-1])
    #plt.tight_layout()
    plt.savefig(os.path.join(outdata,'%s_intercept_histogram_stn%i-%i_%iiter_%ichains_%iwarmup.pdf'\
        %(title_str,min_stn,max_stn,niter,nchains,nwarm)),format='pdf')
    plt.close(f3)

    f4, axarr = plt.subplots(NS+1, 1, figsize=(8,6))
    f4.suptitle("Dependence of Slope on %s"%slope_str)
    for i in range(NI+1):
        yhist, xhist, _ = axarr[i].hist(post['beta1'][:,i],bins=np.int(np.sqrt(len(post['beta1'][:,i]))),\
            color='grey')
        #-- make y axis for plotting vertical lines
        yline = np.arange(np.max(yhist))
        axarr[i].plot(np.ones(len(yline))*np.mean(post['beta1'][:,i]),yline,'k-',linewidth=2)
        axarr[i].plot(np.ones(len(yline))*np.percentile(post['beta1'][:,i],2.5),yline,'k--')
        axarr[i].plot(np.ones(len(yline))*np.percentile(post['beta1'][:,i],97.5),yline,'k--')
        if i ==0:
            axarr[i].set_title("Constant")
        else:
            axarr[i].set_title(slope_vars[i-1])
    #plt.tight_layout()
    plt.savefig(os.path.join(outdata,'%s_slope_histogram_stn%i-%i_%iiter_%ichains_%iwarmup.pdf'\
        %(title_str,min_stn,max_stn,niter,nchains,nwarm)),format='pdf')
    plt.close(f4)


    #-- write results to file
    f = open(os.path.join(outdata,'%s_stn%i-%i_%iiter_%ichains_%iwarmup.txt'\
        %(title_str,min_stn,max_stn,niter,nchains,nwarm)),'w')
    for i in range(p):
        f.write("Stn %i: slope = %.4f (%.4f : %.4f) std=%.4f; intercept =  %.4f (%.4f : %.4f) std=%.4f\n"%(station_nums[i],\
            np.mean(post['b'][:,i]),np.percentile(post['b'][:,i],2.5),np.percentile(post['b'][:,i],97.5),\
            np.std(post['b'][:,i]),np.mean(post['logJ0'][:,i]),np.percentile(post['logJ0'][:,i],2.5),\
            np.percentile(post['logJ0'][:,i],97.5),np.std(post['logJ0'][:,i])))


    f.write("\n\nIntercept Dependence:\n")
    for i in range(NI+1):
        if i==0:
            vname = 'Constant'
        else:
            vname = intercept_vars[i-1]
        f.write("%s: %.4f (%.4f : %.4f) std=%.4f\n"%(vname,np.mean(post['beta0'][:,i]),np.percentile(post['beta0'][:,i],2.5)\
            ,np.percentile(post['beta0'][:,i],97.5),np.std(post['beta0'][:,i])))
    f.write("\n\nSlope Dependence:\n")
    for i in range(NS+1):
        if i==0:
            vname = 'Constant'
        else:
            vname = slope_vars[i-1]
        f.write("%s: %.4f (%.4f : %.4f) std=%.4f\n"%(vname,np.mean(post['beta1'][:,i]),np.percentile(post['beta1'][:,i],2.5)\
            ,np.percentile(post['beta1'][:,i],97.5),np.std(post['beta1'][:,i])))

    f.close()

#-- help function for usage
def usage_info():
    print("--help or -h to display help information.")
    print("Need to input ascii (txt) parameter files with the run confugrations:")
    print("slope_variables  v1,v2,...   # envionemntal dependence if slope")
    print("intercept_variables  v1,v2,...   #envionemntal dependence if intercept")
    print("min_stn      X   # set intial station number X")
    print("max_stn      X   # set final station number X")
    print("iterations   X   # set number of iterations")
    print("chains       X   # set number of chains")
    print("warmup       X   # set number of warmup cycles")
    print("PLOT         Y/N # Y make plots of slopes and intercepts. N no plots except for environmental histograms")
    print("NP           X   # Number of parallel jobs for chains. 0: same as # of chains.")

#-- main function to get parameters and pass them along to fitting function
def main():
    if (len(sys.argv) == 1):
        sys.exit('You need to input at least one parameter file to set run configurations.')
    else:
        if sys.argv[1] in ['-h','--help']:
            usage_info()
            sys.exit()
        else:
            #-- Input Parameter Files (sys.argv[0] is the python code)
            input_files = sys.argv[1:]
            #-- for each input parameter file
            for file in input_files:
                #-- keep track of progress
                print(os.path.basename(file))
                #-- variable with parameter definitions
                parameters = {}
                #-- Opening parameter file and assigning file ID number (fid)
                fid = open(file, 'r')
                #-- for each line in the file will extract the parameter (name and value)
                for fileline in fid:
                    #-- Splitting the input line between parameter name and value
                    part = fileline.split()
                    #-- filling the parameter definition variable
                    parameters[part[0]] = part[1]
                #-- close the parameter file
                fid.close()

                #-- pass parameters to fitting function
                fit_var(parameters)

#-- run main program
if __name__ == '__main__':
        main()
