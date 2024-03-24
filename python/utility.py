#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#  utility.py Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 03.24.2024

##4) Utility
import numpy as np
from python.lib import get_mean_post_prop
from scipy.optimize import minimize_scalar
from scipy.stats import norm, gamma
from python.group_dist import inputs as gd_inputs

nsamp = 100

inputs = {}

##i) Proportion of people with Fatigue about 48 days after hospitalization.
nward = 68
xward = 41
nicu_fatigue = 32
xicu_fatigue = 23

inputs['6w_pward_fatigue'] = get_mean_post_prop(xward_fatigue, nward)
inputs['6w_picu_fatigue'] = get_mean_post_prop(xicu_fatigue, nicu)

##ii) Average loss of utility for hospitalised COVID- 19 patients, split by general ward-care (-6.1%), and care on ITU treatment (-15.5%) was reported in the UK after a mean of 48 days [15].
# https://onlinelibrary.wiley.com/doi/epdf/10.1002/jmv.26368
# The study does not report varition. However, they do report the number of patients with a value < -0.05. 
xbar_ward = -0.061
xbar_icu = -0.155
nl05_ward = 31
pl05_ward = nl05_ward / nward
nl05_icu = 22
pl05_icu = nl05_icu / nicu

#def cost_ward(s):
#    return np.abs(norm.cdf(-0.05, loc = xbar_ward, scale = s) - pl05_ward)
#
#minimize_scalar(cost_ward, [1e-4, 1])
#cost_ward(121312312)
#cost_ward(0.1)

#mu = 0.1
#sig2 = 0.4
#a = mu*mu/sig2
#s = sig2/mu
#gamma.mean(a=a,scale=s)
gamma.var(a=a,scale=s)

def cost_ward(sig2):
    mu = -xbar_ward
    a = mu*mu/sig2
    s = sig2/mu
    return np.abs(gamma.cdf(0.05, a = a, scale = s) - pl05_ward)
def cost_icu(sig2):
    mu = -xbar_icu
    a = mu*mu/sig2
    s = sig2/mu
    return np.abs(gamma.cdf(0.05, a = a, scale = s) - pl05_icu)

some_sigs = np.logspace(-8,0,num=100)

costs = [cost_ward(sig2) for sig2 in some_sigs]
imin = np.argmin(costs)
opt_ward = minimize_scalar(cost_ward, [some_sigs[imin-1], some_sigs[imin+1]])
sig2 = opt_ward.x / np.sqrt(nward)
mu = -xbar_ward
a_ward = mu*mu/sig2
s_ward = sig2/mu

costs = [cost_icu(sig2) for sig2 in some_sigs]
imin = np.argmin(costs)
opt_icu = minimize_scalar(cost_icu, [some_sigs[imin-1], some_sigs[imin+1]])
sig2 = opt_icu.x / np.sqrt(nicu)
mu = -xbar_icu
a_icu = mu*mu/sig2
s_icu = sig2/mu

inputs['utility_change_ward'] = {}
inputs['utility_change_ward']['est'] = -gamma.mean(a=a_ward, scale=s_ward)
inputs['utility_change_ward']['samp'] = -gamma.rvs(a=a_ward, scale=s_ward, size = nsamp)
inputs['utility_change_icu'] = {}
inputs['utility_change_icu']['est'] = -gamma.mean(a=a_icu, scale=s_icu)
inputs['utility_change_icu']['samp'] = -gamma.rvs(a=a_icu, scale=s_icu, size = nsamp)

##iii) COVID-19 symptomatic utility change of -6.1%/60% = -10% for ward patients and -15.5%/72% = -22% for ITU patients symptomatic at 48 days post COVID-19

inputs['SUC_ward'] = {}
inputs['SUC_ward']['est'] = inputs['utility_change_ward']['est'] / inputs['6w_pward_fatigue']['est']
inputs['SUC_ward']['samp'] = inputs['utility_change_ward']['samp'] / inputs['6w_pward_fatigue']['samp']
inputs['SUC_icu'] = {}
inputs['SUC_icu']['est'] = inputs['utility_change_icu']['est'] / inputs['6w_picu_fatigue']['est']
inputs['SUC_icu']['samp'] = inputs['utility_change_icu']['samp'] / inputs['6w_picu_fatigue']['samp']

## Assume nonhosp is ward.
inputs['SUC_nonhosp'] = inputs['SUC_ward']

##iv) UK [21]. At one-year post-discharge, mean utility was 0.58, both for those above and below the age of 65. Taking into account the reference population utility for the UK (0.856) [22], the ARDS specific utility at 1-year was calculated at 0.58/0.856 = 0.68.

mean = 0.58
sd = 0.35
N = 795

reference_qual = 0.85

inputs['util_permin'] = {}
inputs['util_permin']['est'] = mean / reference_qual
inputs['util_permin']['sample'] = norm.rvs(loc=mean,scale=sd/np.sqrt(N),size=nsamp) / reference_qual

##v) Aggregate COVID-specific utility for the symptomatic COVID cohort was calculated at -11%, using a weighted sum of the utilities for the three treatment groups (ward, ITU and population).
## They report 11% in paper but I see 10% in the Excel sheet as well as from these calculations.

inputs['agg_util'] = {}
inputs['agg_util']['est'] = gd_inputs['prop_icu']['est'] * inputs['SUC_icu']['est'] + gd_inputs['prop_ward']['est'] * inputs['SUC_ward']['est'] + (1-gd_inputs['prop_hosp']['est']) * inputs['SUC_nonhosp']['est']
inputs['agg_util']['samp'] = gd_inputs['prop_icu']['samp'] * inputs['SUC_icu']['samp'] + gd_inputs['prop_ward']['samp'] * inputs['SUC_ward']['samp'] + (1-gd_inputs['prop_hosp']['samp']) * inputs['SUC_nonhosp']['samp']
