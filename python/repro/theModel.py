#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

from python.repro.proportionLong import theFloor,population
from python.repro.Utility import Ugeneral,UInjury, ITU_injured, Ward_injured, NonHosp_injured
from python.repro.sim_settings import modeltype
from python.repro.UK_COVIDCases import NumCases
from python.repro.COVIDInfectionSurvey import Constant, PowerTerm
from python.repro.lifeExpectancy import lifeExpectancyAdjusted

daysperyear = 365.25

# In B20-B24
lostQALYs_Injured = 1-Ugeneral #
lostQALYs_Sympt = 1-UInjury # ""QoL lost for the permanently injured (COVID injured)"
FAR = 0.6 # "Ultimate Attack Rate (proportion of UK infected)"
totInfected = population*FAR if modeltype=='Projection - no vaccine' else NumCases

# 
discount = 0.015
discountDaily = np.power(1+discount, 1/daysperyear)-1

## Table starting at B31-J31
n_days = int(np.ceil(daysperyear * 10))+1
Day = np.arange(n_days)

prop_unrecovered = Constant*np.exp(PowerTerm*Day/7) # B column. "proportion of infections un-recovered"
prop_unrecovered_nofloor = (1-theFloor)*prop_unrecovered # C column. "proportion of infections un-recovered (without the floor)"
cum_day_npi = np.cumsum(prop_unrecovered_nofloor) # E column. "Cumulative days non-permanent injury"
cum_day_inj = theFloor * np.arange(1,n_days+1) # F column. "Cumulative days of injury"
qaly_perm = (1-UInjury)*cum_day_inj*totInfected / daysperyear # G column. "QALYs lost to permanent injury"
qaly_symp = lostQALYs_Injured * cum_day_npi * totInfected / daysperyear # H column. "QALYs lost symptomatic COVID"
qaly_total = qaly_perm + qaly_symp # I column "Total undiscounted QALYs lost due to post-COVID illness"

qaly_total_disc = qaly_total[0] + np.cumsum(np.diff(qaly_total)/np.power(1+discountDaily, Day[1:]))
qaly_total_disc = np.concatenate([[qaly_total[0]], qaly_total_disc]) # J column "Total discounted QALYs lost due to post-COVID illness"

##  Table L31-N31
HospFrac = (ITU_injured + Ward_injured)/(ITU_injured + Ward_injured + NonHosp_injured)#
NonHospFrac = (NonHosp_injured)/(ITU_injured + Ward_injured + NonHosp_injured)#

hosp_undisc = (qaly_perm + qaly_symp) * HospFrac # L Column "Hospitalized"
nonhosp_undisc = (qaly_perm + qaly_symp) * NonHospFrac # M Column "Non-hospitalized"
both_undisc = hosp_undisc + nonhosp_undisc 

##  Table L31-N31
hosp_disc = hosp_undisc[0] + np.cumsum(np.diff(hosp_undisc)/np.power(1+discountDaily, Day[1:]))
hosp_disc = np.concatenate([[hosp_undisc[0]], hosp_disc]) # P column  "Hospitalized"
# NOTE: Why not next line?
#nonhosp_disc = nonhosp_undisc[0] + np.cumsum(np.diff(nonhosp_undisc)/np.power(1+discountDaily, Day[1:]))
nonhosp_disc = nonhosp_undisc[0] + np.cumsum(np.diff(nonhosp_undisc)/np.power(1+discountDaily, prop_unrecovered[1:]))
nonhosp_disc = np.concatenate([[nonhosp_undisc[0]], nonhosp_disc]) # Q column  "Non-hospitalized"
both_disc = hosp_disc + nonhosp_disc

# In J13-J25
WTPpQ =  60000. #J13 "Life expectancy for hospitalised patients"
timeHorizon = 10. # J15 "Time horizon (life expectancy of survivors, years)"
# NOTE: on next line, EXCEL round 3652.5 up, while python rounds it down. Add 1 to compensate.
offset1 = round(timeHorizon*daysperyear)+1
# NOTE: But sometimes it's not rounded, and seems to truncate.
offset2 = round(timeHorizon*daysperyear) 

inj_qaly_lost = qaly_perm[offset1]  # J16 "QALYs lost to permanent injury from COVID"
symp_qaly_lost = qaly_symp[offset1]  # J17 "QALYs lost to symptomatic COVID"
total_qaly_lost = qaly_total[offset1]  # J18 "Total QALYs lost with 10 year time horizon"
lostQ_Tot_Undiscounted = qaly_total_disc[offset1] #  J21 "Total discounted QALYs lost with 10 year time horizon"
hosp_qaly_lost = hosp_undisc[offset2] # J22 "QALYs lost in the hospitalised"
nonhosp_qaly_lost = nonhosp_undisc[offset2] # J23 " QALYs lost in the non-hospitalised" 
hosp_qaly_lost_disc = hosp_disc[offset2] # K22 "QALYs lost in the hospitalised" Discounted
nonhosp_qaly_lost_disc = nonhosp_disc[offset2] # K23 " QALYs lost in the non-hospitalised" Discounted

NICE_TWP = lostQ_Tot_Undiscounted*WTPpQ # J24 " HM Treasury value of QALY loss"
num_covid_injured = totInfected * theFloor # J25 "Estimate of people living with COVID injury"
