#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#  group_dist.py Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 02.16.2024

## Question 2.

import numpy as np
from python.settings import nsamp
from python.lib import get_mean_post_prop

inputs = {}

# 2.i from Martin et al.
ntests = 2657305
nhosp = 287662
inputs['prop_hosp'] = get_mean_post_prop(nhosp, ntests)

# 2.ii https://www.bmj.com/content/369/bmj.m1985
ncrit = 3001
nhosp = 18183
inputs['prop_crit'] = get_mean_post_prop(ncrit, nhosp)

inputs['prop_icu'] = {}
inputs['prop_icu']['est'] = inputs['prop_hosp']['est'] * inputs['prop_crit']['est']
inputs['prop_icu']['samp'] = inputs['prop_hosp']['samp'] * inputs['prop_crit']['samp']

inputs['prop_ward'] = {}
inputs['prop_ward']['est'] = inputs['prop_hosp']['est'] - inputs['prop_icu']['est']
inputs['prop_ward']['samp'] = inputs['prop_hosp']['samp'] - inputs['prop_icu']['samp']

# 2.iii 
# Could not find their dataset, so using data from 
# https://www.icnarc.org/DataServices/Attachments/Download/aa75698e-6dde-eb11-9132-00505601089b Page 54.
ncrit = 26550
ndied = 9956
inputs['prop_critdied'] = get_mean_post_prop(ndied, ncrit)

# 2.iv from Martin et al.
nward = 118613
ndied = 27483
inputs['prop_warddied'] = get_mean_post_prop(ndied, nward)

for v in inputs:
    print("Their est:")
    print(inputs[v]['est'])
    print("Post mean")
    print(np.mean(inputs[v]['samp']))
