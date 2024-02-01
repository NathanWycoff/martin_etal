#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd

icu_admissions = 10910
icu_deaths = 4293

hospital_outcomes = pd.DataFrame({
    'admissions' : [113816,14029,1678],
    'deaths' : [29617,1598,561]
    }).T
hospital_outcomes.columns = ['England','Wales','NI']

hosp_admissions = np.sum(hospital_outcomes.loc['admissions',:])
hosp_deaths = np.sum(hospital_outcomes.loc['deaths',:])

ward_admissions = hosp_admissions - icu_admissions
ward_deaths = hosp_deaths - icu_deaths

ward_deathrate = ward_deaths / ward_admissions
