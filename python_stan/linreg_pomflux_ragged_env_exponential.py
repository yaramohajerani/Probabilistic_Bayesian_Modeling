#!/usr/bin/env python
u"""
linreg_pomflux_ragged_env_exponential.py
by Yara Mohajerani (04/2018)

NOT COMPLETED YET; DON'T USE UNTIL PUSHED TO MASTER.

Update History
    04/2018   Forked from linreg_pomflux_ragged_env.py
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
    bInf_vars = parameters['bInf_variables'].split(',')
    b1_vars = parameters['b1_variables'].split(',')
    b2_vars = parameters['b2_variables'].split(',')
    min_stn = np.int(parameters['min_stn'])
    max_stn = parameters['max_stn']
    niter = np.int(parameters['iterations'])
    nchains = np.int(parameters['chains'])
    nwarm = np.int(parameters['warmup'])
    NP = np.int(parameters['NP'])
    #-- if number of processes was never set, set it equal to # of chains
    if NP == 0:
        NP = np.copy(nchains)

    print('looking at stations %s to %s.'%(min_stn,max_stn))
    bInf_str = ''
    print('bInf Dependence:')
    for v in bInf_vars:
        print(v)
        bInf_str += '%s,'%v
    b1_str = ''
    print('b1 Dependence:')
    for v in b1_vars:
        print(v)
        b1_str += '%s,'%v
    b2_str = ''
    print('b2 Dependence:')
    for v in b2_vars:
        print(v)
        b2_str += '%s,'%v
    title_str = 'bInf_%s-dependent_b1_%s-dependent_b2_%s-dependent'%(bInf_str[:-1],b1_str[:-1],b2_str[:-1])

    print('Fit configurations: iterations=%i , chain=%i, warmup=%i'%(niter,nchains,nwarm))
    print('Number of parallel processes = %i'%NP)

    #-- extract names and operations from the given variable names (e.g) lognpp --> log of npp
    bInf_varNames = copy.copy(bInf_vars)
    bInf_ops = copy.copy(bInf_vars)
    b1_varNames = copy.copy(b1_vars)
    b1_ops = copy.copy(b1_vars)
    b2_varNames = copy.copy(b2_vars)
    b2_ops = copy.copy(b2_vars)
    for i,v in enumerate(bInf_vars):
        if v.startswith('log'):
            #-- for now we're assuming the only operation is log
            bInf_varNames[i] = v[3:]
            bInf_ops[i] = 'log'
        else:
            bInf_varNames[i] = v[:]
            bInf_ops[i] = '1'
    for i,v in enumerate(b1_vars):
        if v.startswith('log'):
            #-- for now we're assuming the only operation is log
            b1_varNames[i] = v[3:]
            b1_ops[i] = 'log'
        else:
            b1_varNames[i] = v[:]
            b1_ops[i] = '1'
    for i,v in enumerate(b2_vars):
        if v.startswith('log'):
            #-- for now we're assuming the only operation is log
            b2_varNames[i] = v[3:]
            b2_ops[i] = 'log'
        else:
            b2_varNames[i] = v[:]
            b2_ops[i] = '1'

    #######################################################
    ## Setup your data ####################################
    #######################################################
    d0 = pd.read_csv(os.path.join(indata,'pom_flux','GO_flux_ENV.csv'))

    #-- remove rows where the depth or desired fluxes are missing
    #-- strat by just the poc indices and take intersections iteratively
    ind = np.squeeze(np.nonzero(~np.isnan(d0['flux_poc'])))
    for v in bInf_varNames:
        ind_temp = np.squeeze(np.nonzero(~np.isnan(d0[v])))
        ind = list(set(ind) & set(ind_temp))
    for v in b1_varNames:
        ind_temp = np.squeeze(np.nonzero(~np.isnan(d0[v])))
        ind = list(set(ind) & set(ind_temp))
    for v in b2_varNames:
        ind_temp = np.squeeze(np.nonzero(~np.isnan(d0[v])))
        ind = list(set(ind) & set(ind_temp))

    #-- remove nans and convert to dictionarties with python arrays (easier to work with)
    d = {}
    for v in ['id_mon_yr','depth','flux_poc']:
    	d[v] = d0[v][ind].values
    #-- also apply operations to the environmental variables
    for i,v in enumerate(bInf_varNames):
        if v != None:
            if bInf_ops[i] == 'log':
                d[v] = np.log(d0[v][ind].values)
            else:
                d[v] = d0[v][ind].values
    for i,v in enumerate(b1_varNames):
        if v != None:
            if b1_ops[i] == 'log':
                d[v] = np.log(d0[v][ind].values)
            else:
                d[v] = d0[v][ind].values
    for i,v in enumerate(b2_varNames):
        if v != None:
            if b2_ops[i] == 'log':
                d[v] = np.log(d0[v][ind].values)
            else:
                d[v] = d0[v][ind].values

    #-- go through station numbers and make sure they are all numbers
    #-- remove hyphens if necessary
    for i in range(len(d['id_mon_yr'])):
        try:
            d['id_mon_yr'][i] = int(float(d['id_mon_yr'][i]))
        except:
            d['id_mon_yr'][i] = d['id_mon_yr'][i].replace('-','')
            d['id_mon_yr'][i] = int(float(d['id_mon_yr'][i]))

    if max_stn in ['none','NONE','None']:
        indstn = np.arange(len(d['id_mon_yr']))
    else:
        indstn = []
        for i in range(min_stn,np.int(max_stn)+1):
            #-- when doing a dubset of stations get stations with more than
            #-- one poitns since we just want a small subset anyways
            if np.count_nonzero(d['id_mon_yr']==i) > 1:
                indstn += list(np.squeeze(np.nonzero(d['id_mon_yr']==i)))

    #-- Set up data for stan. we want each station to be separate
    #-- number of unique stations
    station_nums = np.unique(d['id_mon_yr'][indstn])
    p = len(station_nums)
    #-- number of data points per station
    n = np.zeros(p)
    for i in range(p):
        n[i] = np.count_nonzero(d['id_mon_yr'][indstn]==station_nums[i])
    N = np.int(np.sum(n)) #-- total number of points
    #vector index for dataset stacked as single column
    ni = np.concatenate([np.array([0]),np.cumsum(n)])
    ni = np.array(ni,dtype=np.int)

    #-- average all the environmental values for each station (the should all
    #-- be the same when using the id_mon_yr ID)
    bInf_avg = {}
    for v in bInf_varNames:
        bInf_avg[v] = np.zeros(p)
        for i in range(p):
            bInf_avg[v][i] = np.mean(d[v][indstn][ni[i]:ni[i+1]])
    b1_avg = {}
    for v in b1_varNames:
        b1_avg[v] = np.zeros(p)
        for i in range(p):
            b1_avg[v][i] = np.mean(d[v][indstn][ni[i]:ni[i+1]])
    b2_avg = {}
    for v in b2_varNames:
        b2_avg[v] = np.zeros(p)
        for i in range(p):
            b2_avg[v][i] = np.mean(d[v][indstn][ni[i]:ni[i+1]])

    #-- package data for Stan
    X = d['depth'][indstn]-100.
    Y = np.log(d['flux_poc'][indstn])
    Ninf = len(bInf_vars)
    N1 = len(b1_vars)
    N2 = len(b2_vars)
    #-- make design matrices for slope and intercept
    MInf = np.ones((p,Ninf+1))
    M1 = np.ones((p,N1+1))
    M2 = np.ones((p,N2+1))
    for i,v in enumerate(bInf_varNames):
        MInf[:,i+1] = bInf_avg[v][:]
    for i,v in enumerate(b1_varNames):
        M1[:,i+1] = b1_avg[v][:]
    for i,v in enumerate(b2_varNames):
        M2[:,i+1] = b2_avg[v][:]

    dat = dict(p=p, N=N, Ninf=Ninf, N1=N1, N2=N2, ni=ni, x=X, y=Y,\
        MInf=MInf, M1=M1, M2=M2)


    #######################################################
    ## Fit Stan model #####################################
    #######################################################
    #-- First check if the compiled file exists. If not, compile model.
    model_name = 'linreg_pomflux_env_exponential'
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


    #-- plot the environmental dependence histograms
    f1, axarr = plt.subplots(Ninf+1, 1, figsize=(10,10))
    f1.suptitle("Dependence of bInf on %s"%bInf_str)
    for i in range(Ninf+1):
        yhist, xhist, _ = axarr[i].hist(post['betaInf'][:,i],bins=np.int(np.sqrt(len(post['betaInf'][:,i]))),\
            color='grey')
        #-- make y axis for plotting vertical lines
        yline = np.arange(np.max(yhist))
        axarr[i].plot(np.ones(len(yline))*np.mean(post['betaInf'][:,i]),yline,'k-',linewidth=2)
        axarr[i].plot(np.ones(len(yline))*np.percentile(post['betaInf'][:,i],2.5),yline,'k--')
        axarr[i].plot(np.ones(len(yline))*np.percentile(post['betaInf'][:,i],97.5),yline,'k--')
        if i ==0:
            axarr[i].set_title("Constant")
        else:
            axarr[i].set_title(bInf_vars[i-1])
    plt.tight_layout()
    plt.subplots_adjust(top=0.90)
    plt.savefig(os.path.join(outdata,'%s_bInf_histogram_stn%s-%s_%iiter_%ichains_%iwarmup.pdf'\
        %(title_str,min_stn,max_stn,niter,nchains,nwarm)),format='pdf')
    plt.close(f1)

    f2, axarr = plt.subplots(N1+1, 1, figsize=(10,10))
    f2.suptitle("Dependence of b1 on %s"%b1_str)
    for i in range(N1+1):
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
            axarr[i].set_title(b1_vars[i-1])
    plt.tight_layout()
    plt.subplots_adjust(top=0.90)
    plt.savefig(os.path.join(outdata,'%s_b1_histogram_stn%s-%s_%iiter_%ichains_%iwarmup.pdf'\
        %(title_str,min_stn,max_stn,niter,nchains,nwarm)),format='pdf')
    plt.close(f2)

    f3, axarr = plt.subplots(N2+1, 1, figsize=(10,10))
    f3.suptitle("Dependence of b2 on %s"%b2_str)
    for i in range(N1+1):
        yhist, xhist, _ = axarr[i].hist(post['beta2'][:,i],bins=np.int(np.sqrt(len(post['beta2'][:,i]))),\
            color='grey')
        #-- make y axis for plotting vertical lines
        yline = np.arange(np.max(yhist))
        axarr[i].plot(np.ones(len(yline))*np.mean(post['beta2'][:,i]),yline,'k-',linewidth=2)
        axarr[i].plot(np.ones(len(yline))*np.percentile(post['beta2'][:,i],2.5),yline,'k--')
        axarr[i].plot(np.ones(len(yline))*np.percentile(post['beta2'][:,i],97.5),yline,'k--')
        if i ==0:
            axarr[i].set_title("Constant")
        else:
            axarr[i].set_title(b2_vars[i-1])
    plt.tight_layout()
    plt.subplots_adjust(top=0.90)
    plt.savefig(os.path.join(outdata,'%s_b2_histogram_stn%s-%s_%iiter_%ichains_%iwarmup.pdf'\
        %(title_str,min_stn,max_stn,niter,nchains,nwarm)),format='pdf')
    plt.close(f3)


    #-- write results to file
    f = open(os.path.join(outdata,'%s_stn%s-%s_%iiter_%ichains_%iwarmup.txt'\
        %(title_str,min_stn,max_stn,niter,nchains,nwarm)),'w')
    for i in range(p):
        f.write("Stn %i: bInf = %.4f (%.4f : %.4f) std=%.4f; b1 =  %.4f (%.4f : %.4f) std=%.4f; \
            b2 =  %.4f (%.4f : %.4f) std=%.4f\n"%(station_nums[i],np.mean(post['bInf'][:,i]),\
            np.percentile(post['bInf'][:,i],2.5),np.percentile(post['bInf'][:,i],97.5),\
            np.std(post['bInf'][:,i]),np.mean(post['b1'][:,i]),np.percentile(post['b1'][:,i],2.5),\
            np.percentile(post['b1'][:,i],97.5),np.std(post['b1'][:,i]),np.mean(post['b2'][:,i]),\
            np.percentile(post['b2'][:,i],2.5),np.percentile(post['b2'][:,i],97.5),np.std(post['b2'][:,i])))


    f.write("\n\nbInf Dependence:\n")
    for i in range(Ninf+1):
        if i==0:
            vname = 'Constant'
        else:
            vname = bInf_vars[i-1]
        f.write("%s: %.4f (%.4f : %.4f) std=%.4f\n"%(vname,np.mean(post['betaInf'][:,i]),\
            np.percentile(post['betaInf'][:,i],2.5),np.percentile(post['betaInf'][:,i],97.5)\
            ,np.std(post['betaInf'][:,i])))
    f.write("\n\nb1 Dependence:\n")
    for i in range(N1+1):
        if i==0:
            vname = 'Constant'
        else:
            vname = b1_vars[i-1]
        f.write("%s: %.4f (%.4f : %.4f) std=%.4f\n"%(vname,np.mean(post['beta1'][:,i]),np.percentile(post['beta1'][:,i],2.5)\
            ,np.percentile(post['beta1'][:,i],97.5),np.std(post['beta1'][:,i])))
    f.write("\n\nb2 Dependence:\n")
    for i in range(N2+1):
        if i==0:
            vname = 'Constant'
        else:
            vname = b2_vars[i-1]
        f.write("%s: %.4f (%.4f : %.4f) std=%.4f\n"%(vname,np.mean(post['beta2'][:,i]),np.percentile(post['beta2'][:,i],2.5)\
            ,np.percentile(post['beta2'][:,i],97.5),np.std(post['beta2'][:,i])))

    f.close()

#-- help function for usage
def usage_info():
    print("--help or -h to display help information.")
    print("Need to input ascii (txt) parameter files with the run confugrations:")
    print("bInf_variables  v1,v2,...   # envionemntal dependence of bInf")
    print("b1_variables  v1,v2,...   #envionemntal dependence of b1")
    print("b2_variables  v1,v2,...   #envionemntal dependence of b2")
    print("min_stn      X   # set intial station number X")
    print("max_stn      X   # set final station number X")
    print("iterations   X   # set number of iterations")
    print("chains       X   # set number of chains")
    print("warmup       X   # set number of warmup cycles")
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
