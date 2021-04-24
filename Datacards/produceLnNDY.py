import os
import json
import yaml

with open("CompareDY/factors.json","r") as handle:
    factors = json.load(handle)
d = {}
for cat,val in factors.items():
    u = min(0.99,max(0,abs(1-val)))
    era = None
    if '2016' in cat:
        era = '2016'
    if '2017' in cat:
        era = '2017'
    if '2018' in cat:
        era = '2018'
    entry = {'hist':cat,'proc':'DY','val':[round(1-u,3),round(1+u,3)]}
    if era is not None:
        entry['era'] = era
        entry['hist'] = entry['hist'].replace('_'+era,'')
    d[cat.replace('HH_','dy_nonclosure_')] = [entry]
with open("config/lnNDY.yml",'w') as handle:
    yaml.dump(d,handle)

from pprint import pprint
pprint (d)
