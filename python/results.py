#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#  results.py Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 03.26.2024

from python.utility import inputs as u_inputs
from python.repro.proportionLong import population
from python.prev_perm import inputs as  pp_inputs
from python.symptom_prevalence import inputs as  sp_inputs
from python.group_dist import inputs as  gd_inputs
import numpy as np
import matplotlib.pyplot as plt
from python.repro.sim_settings import modeltype

outputs = {}

################################################################################################################################################
### Copied from repor:
daysperyear = 365.25

# In B20-B24
### TODO: This one is wrong.
FAR = 0.6 # "Ultimate Attack Rate (proportion of UK infected)"
totInfected = population*FAR if modeltype=='Projection - no vaccine' else NumCases

discount = 0.015
discountDaily = np.power(1+discount, 1/daysperyear)-1

## Table starting at B31-J31
n_days = int(np.ceil(daysperyear * 10))+1
Day = np.arange(n_days)
theFloor = pp_inputs['overall_permin']['samp']

## TODO: Posteriori predictive
print("Now using post predictive, just post condition exp.")
Constant = sp_inputs['c']['samp']
PowerTerm = sp_inputs['p']['samp']
phi = sp_inputs['phi']['samp']

## First index gives MC iteration, second day of study.
prop_unrecovered = Constant[:,np.newaxis]*np.exp(PowerTerm[:,np.newaxis]*Day[np.newaxis,:]/7) # B column. "proportion of infections un-recovered"
prop_unrecovered_nofloor = (1-theFloor[:,np.newaxis])*prop_unrecovered # C column. "proportion of infections un-recovered (without the floor)"
cum_day_npi = np.cumsum(prop_unrecovered_nofloor, axis = 1) # E column. "Cumulative days non-permanent injury"
cum_day_inj = theFloor[:,np.newaxis] * np.arange(1,n_days+1)[np.newaxis,:] # F column. "Cumulative days of injury"

# NOTE: We are using lostQALYs_Sympt to calc qaly_perm and lostQALYs_Injured to calc qaly_symp. This is not a mistake; we are following Martin et al's naming convention.
qaly_perm = u_inputs['lostQALYs_Sympt']['samp'][:,np.newaxis]*cum_day_inj*totInfected / daysperyear # G column. "QALYs lost to permanent injury"
qaly_symp = u_inputs['lostQALYs_Injured']['samp'][:,np.newaxis] * cum_day_npi * totInfected / daysperyear # H column. "QALYs lost symptomatic COVID"

qaly_total = qaly_perm + qaly_symp # I column "Total undiscounted QALYs lost due to post-COVID illness"

qaly_total_disc = qaly_total[:,0][:,np.newaxis] + np.cumsum(np.diff(qaly_total)/np.power(1+discountDaily, Day[1:]), axis = 1)
qaly_total_disc = np.concatenate([qaly_total[:,0][:,np.newaxis], qaly_total_disc], axis = 1) # J column "Total discounted QALYs lost due to post-COVID illness"

### TODO: Hosp/nonhosp Split
###  Table L31-N31
#HospFrac = (ITU_injured + Ward_injured)/(ITU_injured + Ward_injured + NonHosp_injured)#
#NonHospFrac = (NonHosp_injured)/(ITU_injured + Ward_injured + NonHosp_injured)#
#
#hosp_undisc = (qaly_perm + qaly_symp) * HospFrac # L Column "Hospitalized"
#nonhosp_undisc = (qaly_perm + qaly_symp) * NonHospFrac # M Column "Non-hospitalized"
#both_undisc = hosp_undisc + nonhosp_undisc 
#
###  Table L31-N31
#hosp_disc = hosp_undisc[0] + np.cumsum(np.diff(hosp_undisc)/np.power(1+discountDaily, Day[1:]))
#hosp_disc = np.concatenate([[hosp_undisc[0]], hosp_disc]) # P column  "Hospitalized"
## NOTE: Why not next line?
##nonhosp_disc = nonhosp_undisc[0] + np.cumsum(np.diff(nonhosp_undisc)/np.power(1+discountDaily, Day[1:]))
#nonhosp_disc = nonhosp_undisc[0] + np.cumsum(np.diff(nonhosp_undisc)/np.power(1+discountDaily, prop_unrecovered[1:]))
#nonhosp_disc = np.concatenate([[nonhosp_undisc[0]], nonhosp_disc]) # Q column  "Non-hospitalized"
#both_disc = hosp_disc + nonhosp_disc
### TODO: Hosp/nonhosp Split

# In J13-J25
WTPpQ =  60000. #J13 "Life expectancy for hospitalised patients"
timeHorizon = 10. # J15 "Time horizon (life expectancy of survivors, years)"
# NOTE: on next line, EXCEL round 3652.5 up, while python rounds it down. Add 1 to compensate.
offset1 = round(timeHorizon*daysperyear)+1
# NOTE: But sometimes it's not rounded, and seems to truncate.
offset2 = round(timeHorizon*daysperyear) 

outputs['inj_qaly_lost'] = qaly_perm[:,offset1]  # J16 "QALYs lost to permanent injury from COVID"
outputs['symp_qaly_lost'] = qaly_symp[:,offset1]  # J17 "QALYs lost to symptomatic COVID"
outputs['total_qaly_lost'] = qaly_total[:,offset1]  # J18 "Total QALYs lost with 10 year time horizon"
outputs['lostQ_Tot_Undiscounted'] = qaly_total_disc[:,offset1] #  J21 "Total discounted QALYs lost with 10 year time horizon"
### TODO: Hosp/nonhosp Split
#hosp_qaly_lost = hosp_undisc[offset2] # J22 "QALYs lost in the hospitalised"
#nonhosp_qaly_lost = nonhosp_undisc[offset2] # J23 " QALYs lost in the non-hospitalised" 
#hosp_qaly_lost_disc = hosp_disc[offset2] # K22 "QALYs lost in the hospitalised" Discounted
#nonhosp_qaly_lost_disc = nonhosp_disc[offset2] # K23 " QALYs lost in the non-hospitalised" Discounted
### TODO: Hosp/nonhosp Split

outputs['NICE_TWP'] = outputs['lostQ_Tot_Undiscounted']*WTPpQ # J24 " HM Treasury value of QALY loss"
outputs['num_covid_injured'] = totInfected * theFloor # J25 "Estimate of people living with COVID injury"

from python.repro.theModel import inj_qaly_lost, symp_qaly_lost, total_qaly_lost, lostQ_Tot_Undiscounted, num_covid_injured, NICE_TWP

martin = {}
martin['inj_qaly_lost'] = inj_qaly_lost 
martin['symp_qaly_lost'] = symp_qaly_lost
martin['total_qaly_lost'] = total_qaly_lost
martin['lostQ_Tot_Undiscounted'] = lostQ_Tot_Undiscounted
martin['num_covid_injured'] = num_covid_injured
martin['NICE_TWP'] = NICE_TWP

fig = plt.figure(figsize=[8,6])
trans = lambda x: np.log10(x)
for vi, v in enumerate(outputs):
    plt.subplot(2,3,1+vi)
    plt.hist(trans(outputs[v]))
    ax = plt.gca()
    ul,ll = ax.get_ylim()
    plt.vlines(trans(martin[v]),ul,ll, color = 'red', linestyle = '--')
    plt.title(v)
plt.tight_layout()
plt.savefig("output.pdf")
plt.close()


