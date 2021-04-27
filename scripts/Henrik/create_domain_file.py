import json
import os

# local dependencies
from S2S.local_configuration_H import config

def make_dir(path):
    """
    Creates directory if it does not exist.

    args:
        path: str
    """
    if not os.path.exists(path):
        os.makedirs(path)

def site2domain(site):
    return {
        'name':site['name'],
        'localityNo':site['localityNo'],
        'bounds':(site['lon'],site['lat']),
        'bounds_fmt':'(lon,lat)'
        }

domain_dict = {
            'full_grid':{
                'name':'FullGrid',
                'bounds':(0,360,-90,90),
                'bounds_fmt':'(min lon, max lon, min lat, max lat)'
                    },
            'NVK':{
                'name':'NorVestKyst',
                'bounds':(4,6,58,62),
                'bounds_fmt':'(min lon, max lon, min lat, max lat)'
                    },
            'NVS':{
                'name':'NorVestSjo',
                'bounds':(3,5,58,62),
                'bounds_fmt':'(min lon, max lon, min lat, max lat)'
                    }
                }

with open(config['SITES'],'r') as file:
    sites = json.load(file)

for site in sites:
    domain_dict[site['name']] = site2domain(site)

make_dir('/'.join(config['DOMAINS'].split('/')[:-1]))

with open(config['DOMAINS'], 'w') as outfile:
    json.dump(domain_dict, outfile)
