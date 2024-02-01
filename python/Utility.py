#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
from python.proportionLong import p48_ITU, p48_Ward, pITU_SurvAlln, pWard_SurvAlln, pNonHosp_SurvAlln

eq_df = pd.DataFrame({
            'itu_pre' : [1.6, 1.2, 1.3, 1.4, 1.4],
            'itu_post' : [2.2, 1.3, 2.2, 1.8, 1.9],
            'nonitu_pre' : [1.7, 1.3, 1.5, 1.8, 1.3],
            'nonitu_post' : [2.1, 1.5, 2.0, 1.9, 1.5]
        })
eq_df.index = ['Mobility','Self-care','Usual activites','Pain','Anxiety/depression']

covid_change_ITU = 1-0.155
covid_change_Ward = 1-0.061

j14 = 1. # "Adjustment for lost QALYS for non-hospitalised versus ward care."

covid_change_ITU_adjusted = 1-(1-covid_change_ITU)/p48_ITU
covid_change_Ward_adjusted = 1-(1-covid_change_Ward)/p48_Ward
covid_change_NonHosp_adjusted = 1-(1-covid_change_Ward_adjusted )/j14

denom = pITU_SurvAlln + pWard_SurvAlln + pNonHosp_SurvAlln
Ugeneral = (covid_change_ITU_adjusted*pITU_SurvAlln + covid_change_Ward_adjusted*pWard_SurvAlln + covid_change_NonHosp_adjusted*pNonHosp_SurvAlln) / denom

Uards = 0.58 # "The EQ-5D index at 1 year post ITU discharge for ARDS in those under 65 in Marti et al "
Uengland = 0.85 # "Population reference"
UInjury = Uards / Uengland
