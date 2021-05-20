from ftpretty import ftpretty
from S2S.local_configuration import config
import xarray as xr

local_path = config['IMR']

host = 'ftp.nmdc.no'
path = '/nmdc/IMR/fastestasjoner'

user     = 'anonymous'
password = ''

filenames = []

f = ftpretty(host, user, password)

for long_path in f.list(path):
    for file in f.list(long_path):
        if file[-3:]=='.nc':
            local_name = '/'+file.split('/')[-1]
            f.get(file,local_path+local_name)
            filenames.append(local_name)

for filename in filenames:
    print(xr.open_dataset(local_path+filename))
