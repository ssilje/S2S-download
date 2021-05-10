#!/usr/bin/env python
from ecmwfapi import *
import os,sys
import pandas as pd
from datetime import datetime


server = ECMWFService("mars")

product = 'forecast' # forecast, hincast
dirbase = '/nird/projects/nird/NS9001K/sso102/S2S/DATA/grib'
dir = '%s/%s/%s/'%(dirbase,product,'/ECMWF/sfc')

forcastcycle = 'CY46R1_CY47R1_raw'



if product == 'hindcast':
    STREAM =  'enfh' 
if product == 'forecast':
    STREAM = 'enfo' 
  



basedict = {
    'class': 'od',
    'expver': '1',
    'stream': STREAM ,
    'time': '00:00:00',
    'grid': '0.38/0.38, 
}

l = range(0,1128,24)
paired = ['-'.join([str(v) for v in l[i:i + 2]]) for i in range(len(l))]
final = '/'.join(paired[0:-1])

meta = {
    'tp': {
        'param': '228228',  
        'levtype': 'sfc',
        #'step': '/'.join(['%i'%i for i in range(0,1128,24)]) 
        'step': '0/to/1104/by/24'
    },
    
     't2m': {
        'param': '167',  
        'levtype': 'sfc',
        'step': '/'.join([final]) 
    },
    
     'sst': {
        'param': '34.128',  
        'levtype': 'sfc',
        'step': '0/to/1104/by/24'
       # 'step': '/'.join([final]) 
    },
    
     'u10': {
        'param': '165',  
        'levtype': 'sfc',
        'step': '0/to/1104/by/24'
    },
    
    'v10': {
        'param': '166',  
        'levtype': 'sfc',
        'step': '0/to/1104/by/24'
    },
  
    'mslp': {
        'param': '151',  
        'levtype': 'sfc',
        'step': '0/to/1104/by/24'
    },
     
    'sal': {
        'param': '151175',  
        'levtype': 'o2d',
        'step': '/'.join([final]) 
    }
}

dates_monday = pd.date_range("20200123", periods=52, freq="7D") # forecasts start Thursday
dates_thursday = pd.date_range("20200127", periods=52, freq="7D") # forecasts start Monday
dates_fcycle = dates_monday.union(dates_thursday) 
    
   # Program start
for filename in (
    'sst',
    
):
    for prefix in (
        'pf',
        'cf',
    ):
        for dates in dates_fcycle:
            d = dates.strftime('%Y-%m-%d')
            refyear = int(d[:4])
            datadir = '%s/%s'%(dir,filename)
            if not os.path.exists(datadir)  :
                os.makedirs(datadir)
            hdate = '/'.join([d.replace('%i'%refyear,'%i'%i) for i in range(refyear-20,refyear)])
            target = '%s/%s_%s_%s_%s_%s.grb'%(datadir,filename,forcastcycle,d,prefix,product)
            if not os.path.isfile(target):
               dic = basedict.copy()
               for k,v in meta[filename].items():
                   dic[k] = v
               dic['date'] = d
               dic['type'] = prefix
               if ( product == 'hindcast' ):
                   dic['hdate'] = hdate
                   if prefix == 'pf':
                       dic['number'] =  '1/to/10'
               if ( product == 'forecast' ):
                   if prefix == 'pf':
                       dic['number'] =  '1/to/50'
              # out = '%s,%s=%s'%(dic, '"target"',target)
               out= '%s=%s'%("'target'", target) 
               print(out)
               if server is not None:
                   server.execute(dic,out)
                    
print('DONE')
