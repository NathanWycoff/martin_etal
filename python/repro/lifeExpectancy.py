#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np

mle = pd.read_csv('data/male_life_expectancy.csv')
fle = pd.read_csv('data/female_life_expectancy.csv')

weighted_pop_ex = (mle['ex'] * mle['Pop 2019'] + fle['ex'] * fle['Pop 2019']) / (mle['Pop 2019'] + fle['Pop 2019'])

ages = pd.read_csv('data/age_distribution.csv')
ages['Age distribution'][0]

lifeExpectancy = np.sum(ages['Age distribution'] * weighted_pop_ex.iloc[:ages.shape[0]]) # "Life expectancy of those admitted."
LE_RedFact = 0.5 #"Reduction factor for the admitted population "
lifeExpectancyAdjusted = lifeExpectancy*LE_RedFact
