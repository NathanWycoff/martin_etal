#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#  prev_perm.py Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 03.24.2024

import numpy as np
from python.lib import get_mean_post_prop, nsamp
from python.group_dist import inputs as gd_inputs

inputs = {}

# Question 3

## i) 12% of the 66.6 million UK population having been infected
##     a) ONS. Coronavirus (COVID-19) Infection Survey: antibody data for the UK, January 2021. 2021.

pops = {'england' : 56536000, 'wales' : 3105000, 'ni' : 1905000, 'scotland' : 5480000}
#pops_vec = np.array([56536000, 3105000, 1905000, 5480000])
pops_vec = np.array([pops[v] for v in pops])
ukpop = np.sum(pops_vec)
pop_weights = dict([(v, pops[v] / ukpop) for v in pops])

# Looks similar.
# https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/datasets/coronaviruscovid19antibodydatafortheuk/2021
# Feb 21
#props = np.array([0.153, 0.112, 0.092, 0.101])

# https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/datasets/coronaviruscovid19infectionsinthecommunityinengland
# Jan 19
#sampsize = np.array([20878,539,207,1241])
#nposanti = np.array([2205, 50, 11,114])

# Since they use weights but do not tell us how they were calculated, it is not clear how to do an exact Bayesian version of their procedure.
# We will use the success count necessary to match their estimate in expectation
#sampsize = np.array([20878,539,207,1241])
#props = np.array([0.121, 0.098, 0.078, 0.089])
sampsize = {'england' : 20878, 'wales' : 539, 'ni' : 207, 'scotland' : 1241}
props = {'england' : 0.121, 'wales' : 0.098, 'ni' : 0.078, 'scotland' : 0.089}
implicit_successes = dict([(v,sampsize[v]*props[v]) for v in props])

inputs['prop_infected_by_country'] = {}
for v in pops:
    inputs['prop_infected_by_country'][v] = get_mean_post_prop(implicit_successes[v], sampsize[v], nsamp = nsamp)

inputs['prop_infected'] = {}
inputs['prop_infected']['est'] = np.sum([inputs['prop_infected_by_country'][v]['est'] * pop_weights[v] for v in pops])
inputs['prop_infected']['samp'] = np.sum([inputs['prop_infected_by_country'][v]['samp'] * pop_weights[v] for v in pops], axis = 0)

## ii) 33% of SARS-Cov-2 infections resulted in a positive test (assuming no reinfections).
##     a) GOV.UK Coronavirus (COVID-19) in the UK. 2020
numpos = 2.66e6
inputs['prob_test_given_pos'] = {}
inputs['prob_test_given_pos']['est'] = numpos / (inputs['prop_infected']['est'] * ukpop)
inputs['prob_test_given_pos']['samp'] = numpos / (inputs['prop_infected']['samp'] * ukpop)

## iii) Symptom prevalence at 6 weeks was set at 72% for ITU and 60% for ward patients
##     a) https://doi.org/10.1002/jmv.26368
nicu = 32
nward = 68
picu = 0.72
pward = 0.603
xicu = int(np.round(nicu * picu))
xward = int(np.round(nward * pward))

inputs['6wsymp_icu'] = get_mean_post_prop(xicu, nicu, nsamp = nsamp)
inputs['6wsymp_ward'] = get_mean_post_prop(xward, nward, nsamp = nsamp)

## iv) 16% for nonhosp
##     a) Office for National Statistics. The prevalence of long COVID symptoms and COVID-19 complications.  2020
## Answered in symptom_prevalence.py

## v) "Reflecting this uplift, here it was estimated that 50% of the
# ITU survivors with symptoms at 6 weeks post-COVID will be left permanently injured. For
# the ward-based group and the non-hospitalised group, here it was assumed 5% permanently
# injured and 0.5% permanently injured amongst survivors symptomatic at 6 weeks post-
# COVID, respectively (10% and 1% of the ITU rate, respectively"
# These numbers kind just get made out of thin air.

inputs['perminj_icu'] = {}
inputs['perminj_icu']['est'] = 0.5
inputs['perminj_icu']['samp'] = 0.5*np.ones(nsamp)
inputs['perminj_ward'] = {}
inputs['perminj_ward']['est'] = 0.05
inputs['perminj_ward']['samp'] = 0.05*np.ones(nsamp)
inputs['perminj_nonhosp'] = {}
inputs['perminj_nonhosp']['est'] = 0.005
inputs['perminj_nonhosp']['samp'] = 0.005*np.ones(nsamp)


## vi) Go from those assumption to final prevalence.

icu_survivors_est = gd_inputs['prop_icu']['est']*(1-gd_inputs['prop_critdied']['est'])
ward_survivors_est = gd_inputs['prop_ward']['est'] *(1- gd_inputs['prop_warddied']['est'])
icu_survivors_samp = gd_inputs['prop_icu']['samp']*(1-gd_inputs['prop_critdied']['samp'])
ward_survivors_samp = gd_inputs['prop_ward']['samp'] *(1- gd_inputs['prop_warddied']['samp'])

undetect_correction = 1./3.
overall_permin_est = inputs['perminj_icu']['est'] * icu_survivors_est  + inputs['perminj_ward']['est'] * ward_survivors_est  + inputs['perminj_ward']['est'] * (1-gd_inputs['prop_hosp']['est']) *undetect_correction
overall_permin_samp = inputs['perminj_icu']['samp'] * icu_survivors_samp  + inputs['perminj_ward']['samp'] * ward_survivors_samp  + inputs['perminj_ward']['samp'] * (1-gd_inputs['prop_hosp']['samp']) *undetect_correction
inputs['overall_permin'] = {}
inputs['overall_permin']['est'] = overall_permin_est
inputs['overall_permin']['samp'] = overall_permin_samp
