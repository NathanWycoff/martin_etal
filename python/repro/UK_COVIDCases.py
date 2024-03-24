import pandas as pd
import numpy as np

df = pd.read_csv("./data/uk_covid_cases.csv")
df = df.T
df.columns = df.iloc[0,:]
df = df.iloc[1:,:]

NumCases = np.sum(df['Cases'])
numHosp = np.sum(df['Admissions'])
numNonHosp = np.sum(df['Non-hospitalised cases'])
