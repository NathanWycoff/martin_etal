#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from python.repro.WardSurvival import ward_deathrate
from python.repro.sim_settings import modeltype

#1 "To identify the proportion of patients hospitalised"
# UK-wide covid data.
positive_tests = 2657305 
hospitalized_uk = 287662 

not_hosp = positive_tests - hospitalized_uk #D13
prop_hosp = hospitalized_uk / positive_tests #E13

#5 "To estimate the proportion of undetected cases" (out of order because is necessary to compute 2.1).
population = 6.66e7 # UK Population	
prop_infected = 0.12 # Proportion of UK population have been infected
covid_infections = population * prop_infected
all_case_correction = positive_tests / covid_infections
puduh = 1. #Proportion of undetected cases are non-hospitalised	

#2 "To estimate the proportion of surviving hospitalised patients of all known cases"
#2.1 "To identify the proportion of hospitalised patients that are admitted to ITU, we used Docherty et al.This was a study of 20,133 patients admitted to hospital with COVID-19 in the UK using the ISARIC data set."
# From the study of Docherty et al
hospitalized_dochert = 18183
critical_care = 3001
prop_itu = critical_care / hospitalized_dochert  #D20

# 2.2 estimate the proportion of positive cases in the Sudre(2020) population that were addmitted to critical care:
prop_itu = prop_hosp * prop_itu # D24
prop_ward = prop_hosp - prop_itu # Proportion of cases admitted to ward only
prop_itu_corrected = prop_ward * all_case_correction # TODO: Why isn't this prop_itu * acc, instead of prop_ward * acc?

# 2.3
prop_itu_dies = 0.38 #"Proportion of ITU patients died"
prop_itu_survivors = (1-prop_itu_dies) * prop_itu # Proportion of known cases who become ITU survivors.
prop_itu_survivors_corrected = prop_itu_survivors  * all_case_correction # "Proportion of all cases suriviving ITU care"

# 2.4 "To estimate the surviving rate of patients admitted to ward only among all postivie cases"
prop_ward_dies = 0.232
prop_ward_survivors = (1-ward_deathrate)*prop_ward
prop_ward_survivors_corrected = prop_ward_survivors*all_case_correction

# 3 "To estimate the proportion of surviving patients not admitted to hospital of all positive cases"
prop_nonhosp_dies = 0.
prop_nonhosp_survivors = (1-prop_nonhosp_dies)*(1-prop_hosp)
#prop_nonhosp_survivors_corrected = prop_nonhosp_survivors*all_case_correction ?
prop_nonhosp_survivors_corrected = (covid_infections - hospitalized_uk) / (covid_infections * (1-prop_nonhosp_dies))

# Compute proportions from 2.4-3
denom = prop_itu_survivors_corrected + prop_ward_survivors_corrected + prop_nonhosp_survivors_corrected
pITU_survAlln = prop_itu_survivors_corrected / denom # /G29 "Proportion of all survivors suriviving ITU care"
pWard_SurvAlln = prop_ward_survivors_corrected / denom # /G35 "Proportion of all survivors surviving ward"
pNonHosp_SurvAlln = prop_nonhosp_survivors_corrected / denom # /G39 "Proportion of all survivors that are non-hospitalised and survive"

# 4.1 "From Halpin, taken as the prevalence of the most common symptom (fatigue) at 4-10 weeks post-discharge (mean of 48 days)."
p48_ITU = 0.72     # p48_ITU "Prevalence of symptoms at 48 days for ITU survivors"
p48_Ward = 0.6     # p48_Ward "Prevalence of symptoms at 48 days for ward survivors"
p48_NonHosp = 0.16 # p48_NonHosp "Prevalence of symptoms at 56 days for non-hospitalised"

# 4.2 is omitted in the Excel sheet.

# 4.3 "To estimate symptomatic-COVID at 6 weeks of all known cases"
symp_oak_itu =  p48_ITU * prop_itu_survivors # p48_ITU_allCases 
symp_oak_ward =  p48_Ward * prop_ward_survivors # p48Ward_allCases 
symp_oak_nonhosp =  p48_NonHosp * prop_nonhosp_survivors # p48NonHosp_allCases 

# 4.4 "Assumption for proportion of symptomatic-COVID permanently injured"
# ppi - proportion of permanent injury
injuryRate_ITU = 0.5 # "Proportion of symptomatic-COVID at 6 weeks that are permanently injured in ICU survivors."
injuryRate_Ward = injuryRate_ITU / 10 # "Proportion of symptomatic-COVID at 6 weeks that are permanently injured in ward survivors."
injuryRate_NonHosp = injuryRate_Ward / 10 # "Proportion of symptomatic-COVID at 6 weeks that are permanently injured in non-hospitalised survivors."

# 4.5 "To estimate permanently injured after 6 weeks."
ppi_oak_hosp = symp_oak_itu * injuryRate_ITU + symp_oak_ward * injuryRate_Ward # D59 "Estimated proportion of all known cases permanently injured after hospitalisation.	"
ppi_oak_nonhosp = symp_oak_nonhosp * injuryRate_NonHosp # D60 "Estimated proportion of all known cases permanently injured in the non-hospitalised	"

# 6 "Adjust prevalence of permanent injury for all infections, known and unknown"
theFloor = ppi_oak_hosp + ppi_oak_nonhosp  # "Adjusted prevalence of permanent injury for all infections, known and unknown"
if modeltype != 'To date':
    theFloor *= all_case_correction
