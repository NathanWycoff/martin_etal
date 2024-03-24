#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#  symptom_prevalence.py Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 02.09.2024

## Fully Bayesian version of the results in the Prevalence of Symptoms section.

## 1.

## And also 

import arviz as az
import numpy as np
from python.repro.COVIDInfectionSurvey import Constant, PowerTerm 
import pyjags
import matplotlib.pyplot as plt
import pandas as pd

week = np.array([0, 5, 12])
prevalence = np.array([0.5, 0.2, 0.1])

week_pred = np.linspace(0,15,num=100)

theirpred = Constant*np.exp(PowerTerm*week_pred)

inputs = {}

#model = 'linU'
#model = 'linN'
#model = 'exp'
#model = 'beta'
for model in ['linU','linN','exp','beta']:

    if model[:3]=='lin':
        spec = ''' 
        model {
          ## DRY1
          for (i in 1:length(log_prevalence)) {
            log_prevalence[i] ~ dnorm(c + p*week[i], phi)
          }
          ## DRY1
          ## DRY1
          for (i in 1:length(week_pred)) {
            pred_prevalence[i] ~ dnorm(c + p*week_pred[i], phi)
          }
          ## DRY1
          phi ~ dt(0,1,1)T(0,)
          '''
        if model[3] == 'U':
            spec += '''
            c ~ dunif(0,1)
            p ~ dunif(-1,0)
            }'''
        elif model[3]=='N':
            spec += '''
            c ~ dnorm(0,1e-3)
            p ~ dnorm(0,1e-3)
            }'''
        else: 
            raise Exception("Unknown prior for loglinear model.")
        data = {'log_prevalence':np.log(prevalence), 'week' : week, 'week_pred' : week_pred} 
    elif model=='exp':
        spec = ''' 
        model {
          ## DRY2
          for (i in 1:length(prevalence)) {
            prevalence[i] ~ dnorm(c * exp(p*week[i]), phi)T(0,)
          }
          ## DRY2
          ## DRY2
          for (i in 1:length(week_pred)) {
            pred_prevalence[i] ~ dnorm(c * exp(p*week_pred[i]), phi)T(0,)
          }
          ## DRY2

          c ~ dunif(0,1)
          p ~ dunif(-1,0)
          phi ~ dt(0,1,1)T(0,)
        }
        '''
        data = {'prevalence':prevalence, 'week' : week, 'week_pred' : week_pred}
    elif model=='beta':
        spec = ''' 
        model {
          ## DRY3
          for (i in 1:length(prevalence)) {
            mu[i] <- c*exp(p*week[i])
            alpha[i] <- mu[i] * phi
            beta[i] <- (1-mu[i]) * phi
            prevalence[i] ~ dbeta(alpha[i], beta[i])
          }
          ## DRY3

          ## DRY3
          for (i in 1:length(week_pred)) {
            mu_pred[i] <- c*exp(p*week_pred[i])
            alpha_pred[i] <- mu_pred[i] * phi
            beta_pred[i] <- (1-mu_pred[i]) * phi
            pred_prevalence[i] ~ dbeta(alpha_pred[i], beta_pred[i])
          }
          ## DRY3

          c ~ dunif(0,1)
          p ~ dunif(-1,0)
          phi ~ dt(0,1,1)T(0,)
        }
        '''
        data = {'prevalence':prevalence, 'week' : week, 'week_pred' : week_pred}
    else:
        raise Exception("model not known.")


    jm = pyjags.Model(code=spec, data=data, chains=3)

    #iters = 100000
    iters = 1000
    burnin = iters//2
    params = ['c','p','phi']
    samp = jm.sample(iters, vars = params+['pred_prevalence'])
    param_samp = dict([(v,samp[v]) for v in params]) 
    pred_samp = samp['pred_prevalence']

    plt.figure()
    az.plot_trace(param_samp)
    plt.tight_layout()
    plt.savefig("images/prev_trace.pdf")
    plt.close()

    mea = {'c':Constant,'p':PowerTerm}

    ## Marginal histograms
    fig = plt.figure()
    for i,p in enumerate(['c','p']):
        plt.subplot(1,2,i+1)
        keepsamp = param_samp[p][:,:burnin,:].flatten()

        if model=='lin':
            alpha = 0.02
            bounds = np.quantile(keepsamp,[alpha/2,1-alpha/2])
            #plt.xlim(np.quantile(keepsamp,0.01), np.quantile(keepsamp,0.99))
            keepsamp = keepsamp[np.logical_and(keepsamp>bounds[0], keepsamp<bounds[1])]

        if p=='c' and model=='lin':
            keepsamp = np.exp(keepsamp)
        plt.hist(keepsamp, bins = 100)
        ax = plt.gca()
        ll,ul = ax.get_ylim()
        plt.vlines(np.mean(keepsamp), ymin = ll, ymax = ul, color = 'orange', label = 'Mean')
        plt.vlines(mea[p], ymin = ll, ymax = ul, color = 'red', label = 'Martin et al')
        plt.title(model + " " + p + " Estimate")
        if i==0:
            #plt.legend(prop={'shrink':0.5})
            plt.legend()

    plt.savefig("images/prev_"+model+".pdf")
    plt.close()

    ## Posterior Predictive.
    fig = plt.figure()
    pred = pred_samp[:,:burnin,:].reshape((len(week_pred),-1))

    if model[:3]=='lin':
        pred = np.exp(pred)

    med = np.median(pred, axis = 1)
    lb = np.quantile(pred, 0.95, axis = 1)
    ub = np.quantile(pred, 0.05, axis = 1)

    plt.plot(week_pred, med, color = 'blue')
    plt.plot(week_pred, lb, color = 'blue', linestyle = '--')
    plt.plot(week_pred, ub, color = 'blue', linestyle = '--')
    plt.plot(week_pred, theirpred, color = 'orange')
    plt.scatter(week, prevalence)

    plt.title(model)

    plt.savefig("images/pred_"+model+".pdf")
    plt.close()

    ## Answer 3.iv
    eightweeks = np.argmin(np.abs(8-week_pred))
    inputs['8wsymp_nonhosp']['samp'] = pred[eightweeks,:]
    inputs['8wsymp_nonhosp']['est'] = np.mean(pred[eightweeks,:])

    ## Multivariate shit
    samp_flat = dict([(v,param_samp[v][:,burnin:,:].flatten()) for v in param_samp])
    sdf = pd.DataFrame(samp_flat)

    sdfc = sdf - np.array(np.mean(sdf, axis = 0))[np.newaxis,:]
    ndfc = (sdf - np.array(np.mean(sdf, axis = 0))[np.newaxis,:])/np.array(np.std(sdf, axis = 0))[np.newaxis,:]
    S = sdfc.T @ sdfc
    sd = np.sqrt(np.diag(S))
    C = np.diag(1/sd) @ S @ np.diag(1/sd)

    np.linalg.eigh(S)
    np.linalg.eigh(C)

    for pair in [1,2]:
        fig = plt.figure()
        if pair==2:
            std = np.sqrt(1/samp_flat['phi'])
            #std = np.minimum(np.sqrt(1/samp_flat['phi']),1)
            if model=='exp':
                std = np.minimum(std, 1)
            plt.hist2d(samp_flat['p'], std, bins = 30)
            plt.xlabel("Power")
            if model=='beta':
                plt.ylabel("Scale")
            else:
                plt.ylabel("Standard Deviation")
        elif pair==1:
            plt.hist2d(samp_flat['p'], samp_flat['c'], bins = 30)
            plt.xlabel("Power")
            plt.ylabel("Constant")
        else:
            raise Exception("Unknown plot.")
        Ns = 100
        sind = np.random.choice(len(samp_flat['p']),Ns,replace=False)
        plt.scatter(samp_flat['p'][sind], np.minimum(np.sqrt(1/samp_flat['phi']),1)[sind], color='orange')
        pn = 'phi' if pair==2 else 'c'
        plt.savefig("images/biv_"+pn+"_"+model+".pdf")
        plt.close()

    pnorm = (samp_flat['p'] - np.mean(samp_flat['p'])) / np.std(samp_flat['p'])
    cnorm = (samp_flat['c'] - np.mean(samp_flat['c'])) / np.std(samp_flat['c'])
    np.mean(pnorm*cnorm)
