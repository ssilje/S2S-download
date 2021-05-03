import numpy as np
import csv
import json

path =         '/home/heau/Norce/S2S-download/data/BW/BW_temperature'
bw_fname =     'BW_temperature_export.csv'
m_data_fname = 'metadata_BW_sites.json'

with open('./BW_temperature_export.csv','r') as file:
    data = csv.DictReader(file)

    for line in data:
        print(line)
        exit()
