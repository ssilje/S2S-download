import pandas as pd
import json

from S2S.local_configuration import config

path     = config['EIDE']
filename = 'eidedata.xlsx'

df = pd.read_excel(path+filename)
print(df)
