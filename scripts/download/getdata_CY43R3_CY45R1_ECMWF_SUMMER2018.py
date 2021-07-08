#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
import os,sys
import pandas as pd
from datetime import datetime

server = ECMWFDataServer()

product = 'forecast' # forecast, hincast
dirbase = '/nird/projects/nird/NS9853K/DATA/S2S/SUMMER2018/'
dir = '%s/%s/%s/'%(dirbase,product,'/ECMWF/sfc')

forecastcycle = 'CY43R3_CY45R1'

if product == 'hindcast':
    STREAM =  'enfh', 
if product == 'forecast':
    STREAM = 'enfo', 
        
basedict = {
    'class': 's2',
    'dataset': 's2s',
    'expver': 'prod',
    'model': 'glob',
    'origin': 'ecmf',
    'stream': STREAM ,
    'time': '00:00:00'
}

l = range(0,1128,24)
paired = ['-'.join([str(v) for v in l[i:i + 2]]) for i in range(len(l))]
final = '/'.join(paired[0:-1])

meta = {
    'tp': {
        'param': '228228',  
        'levtype': 'sfc',
        'step': '/'.join(['%i'%i for i in range(0,1128,24)]) 
    },
    
     't2m': {
        'param': '167',  
        'levtype': 'sfc',
        'step': '/'.join([final]) 
    },
    
     'sst': {
        'param': '34',  
        'levtype': 'sfc',
        'step': '/'.join([final]) 
    },
    
     'u10': {
        'param': '165',  
        'levtype': 'sfc',
        'step': '/'.join(['%i'%i for i in range(0,1128,24)]) 
    },
    
    'v10': {
        'param': '166',  
        'levtype': 'sfc',
        'step': '/'.join(['%i'%i for i in range(0,1128,24)]) 
    },
    
    'sm20': {
        'param': '228086',  
        'levtype': 'sfc',
        'step': '/'.join([final]) 
    }
}

#2018-05-03
dates_monday = pd.date_range("20180301", periods=20, freq="7D") # forecasts start Monday
dates_thursday = pd.date_range("20180305", periods=20, freq="7D") # forecasts start Thursday
dates_fcycle = dates_monday.union(dates_thursday)   
    
   # Program start
for filename in (
    'tp',
    't2m',
    'sm20',
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
            target = '%s/%s_%s_%s_%s.grb'%(datadir,filename,forecastcycle,d,prefix)
            if not os.path.isfile(target):
               dic = basedict.copy()
               for k,v in meta[filename].items():
                   dic[k] = v
               dic['date'] = d 
               dic['hdate'] = hdate
               #dic['number'] =  '0/to/10'
               dic['type'] = prefix
               if ( product == 'hindcast' ):
                   dic['hdate'] = hdate
                   if prefix == 'pf':
                       dic['number'] =  '1/to/10'
               if ( product == 'forecast' ):
                   if prefix == 'pf':
                       dic['number'] =  '1/to/50'
               dic['target'] = target    
               print(dic)
               if server is not None:
                   server.retrieve(dic)
