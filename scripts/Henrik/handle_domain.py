import json

from S2S.local_configuration import config

def get_bounds(domainID):
    with open(config['DOMAINS'], 'r') as file:
        domain_dict = json.load(file)
    return domain_dict[domainID]['bounds']

def get_domain(domainID):
    with open(config['DOMAINS'], 'r') as file:
        domain_dict = json.load(file)
    return domain_dict[domainID]

### FUNCTIONS NOT IN USE ###
def try_keys(dictionary,list_of_keys):
    for n,key in enumerate(list_of_keys):
        try:
            if n<len(list_of_keys)-1:
                dictionary = dictionary[key]
            else:
                return n,key
        except KeyError:
            return n,key

def update_storage_dict(
                        domainID,
                        ftype,
                        model,
                        var_abb,
                        t_res,
                        i_time,
                        e_time,
                        dtype,
                        store_path
                        ):

    with open('./data/EIDE/domains.json', 'r') as file:
        domain_dict = json.load(file)

    list_of_keys = [ftype,model,var_abb,t_res]

    basedict = {
                model:{
                    var_abb:{
                        t_res:{
                            'i_time':i_time,
                            'e_time':e_time,
                            'dtype':dtype,
                            'store_path':store_path
                                }
                            }
                        }
                    }

    n,key = try_keys(domain_dict[domainID]['storage'],list_of_keys)

    if n==0:
        domain_dict[domainID]['storage'][ftype] = basedict
    if n==1:
        domain_dict[domainID]['storage'][ftype][model] = basedict[model]
    if n==2:
        domain_dict[domainID]['storage'][ftype][model][var_abb] =\
            basedict[model][var_abb]
    if n==3:
        domain_dict[domainID]['storage'][ftype][model][var_abb][t_res] =\
            basedict[model][var_abb][t_res]

    with open('./data/EIDE/domains.json', 'w') as outfile:
        json.dump(domain_dict,outfile)
