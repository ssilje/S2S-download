import json
import os

# local dependencies
from S2S.local_configuration import config

def make_dir(path):
    """
    Creates directory if it does not exist.

    args:
        path: str
    """
    if not os.path.exists(path):
        os.makedirs(path)

def bump2map(domain):
    map_domain = domain.copy()
    map_domain['bounds'] = (
                        map_domain['bounds'][0]-3,map_domain['bounds'][1]+3,
                        map_domain['bounds'][2]-1,map_domain['bounds'][3]+1
                        )
    map_domain['name'] = map_domain['name']+'_map'
    return map_domain


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
                'bounds':(2,6,57,63),
                'bounds_fmt':'(min lon, max lon, min lat, max lat)'
                    },
            'NVS':{
                'name':'NorVestSjo',
                'bounds':(3,5,58,62),
                'bounds_fmt':'(min lon, max lon, min lat, max lat)'
                    },
            'norwegian_coast':{
                'name':'Norwegian Coast',
                'bounds':(0,28,55,75),
                'bounds_fmt':'(min lon, max lon, min lat, max lat)'
                    }
                }

for domainID in list(domain_dict):
    domain_dict[domainID+'_map'] = bump2map(domain_dict[domainID])

with open(config['SITES'],'r') as file:
    sites = json.load(file)

for site in sites:
    domain_dict[site['name']] = site2domain(site)

make_dir('/'.join(config['DOMAINS'].split('/')[:-1]))

with open(config['DOMAINS'], 'w') as outfile:
    json.dump(domain_dict, outfile)
