import numpy as np
import csv
import json
import pandas as pd
import xarray as xr
import os

from S2S.local_configuration import config

def make_dir(path):
    """
    Creates directory if it does not exist.

    args:
        path: str
    """
    if not os.path.exists(path):
        os.makedirs(path)

def name_from_loc(loc):
    with open(config['SITES'], 'r') as file:
        data = json.load(file)
        for line in data:
            if line["localityNo"]==int(loc):
                return line['name']

def organize_files():

    path   = config['BW']
    f_path = path + 'BW_temperature_export.csv'
    m_path = path + 'metadata_BW_sites.json'

    make_dir(path)

    ##############################################
    #### Create localityNo - meta data parser ####
    ##############################################
    with open(m_path,'r') as file:
        raw_site = json.load(file)

    site = {}
    for instance in raw_site:
        locNo = str(instance['localityNo'])
        site[locNo]                  = instance
        site[locNo]['var_name']      = 'sst'
        site[locNo]['date']          = []
        site[locNo]['active_status'] = []
        site[locNo]['value']         = []

    ##############################
    #### Append data to sites ####
    ##############################
    with open(f_path,'r') as file:
        data = csv.DictReader(file)
        for line in data:
            
            locNo = str(line['locNo'])
            try:
                site[locNo]['date'].append(line['date'])
                site[locNo]['active_status'].append(line['active_status'])
                site[locNo]['value'].append(line['temperature'])
            except KeyError:
                pass

    #########################################
    #### Select only locations with data ####
    #########################################
    end_date   = pd.to_datetime('1900-01-01')
    start_date = pd.to_datetime('2100-01-01')

    final_site = {}
    for locNo in site:
        if site[locNo]['value']:
            final_site[locNo] = site[locNo]

            t_start_date = pd.to_datetime(final_site[locNo]['date'][0])
            t_end_date   = pd.to_datetime(final_site[locNo]['date'][-1])

            if start_date > t_start_date:
                start_date = t_start_date

            if end_date < t_end_date:
                end_date = t_end_date

    with open(path+'temp_BW.json', 'w') as outfile:
        json.dump(final_site, outfile)

    ###############################
    #### Final sites to xarray ####
    ###############################
    print('Location','\t\t','obs not nan','\t','total obs')
    for key in final_site:

        filename  = '_'.join(['barentswatch',key])+'.nc'

        sst  = np.array(final_site[key]['value'],dtype='float32')
        time = pd.to_datetime(final_site[key]['date'])
        print(name_from_loc(key),'\t\t',np.count_nonzero(~np.isnan(sst)),'\t',len(sst))
        df = pd.DataFrame(
                    {
                        'sst' :sst,
                        'time':time
                    }
                ).set_index(['time'])

        df = df.groupby(df.index).mean()
        # df = df.reindex(pd.date_range(start=start_date,end=end_date,freq='W'))
        df.index.name = 'time'
        # df['lon'] = np.stack([final_site[key]['lon']]*df.shape[0])
        # df['lat'] = np.stack([final_site[key]['lat']]*df.shape[0])
        df['location'] = np.stack([key]*df.shape[0])
        df.reset_index().set_index(['location','time']).to_xarray()\
            .assign_coords(
                {
                    'lon':('location',[final_site[key]['lon']]),
                    'lat':('location',[final_site[key]['lat']])
                }
            ).to_netcdf(path+filename)
