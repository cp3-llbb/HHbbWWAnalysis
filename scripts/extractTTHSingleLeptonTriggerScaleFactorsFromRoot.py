#! /bin/env python

import math
import argparse
import ROOT
import re
import json


parser = argparse.ArgumentParser()
parser.add_argument('--SF_file', help='ROOT file containing nominal scale factors')
parser.add_argument('-s', '--suffix', help='Suffix to append at the end of the output filename', required=True)

args = parser.parse_args()


f = ROOT.TFile.Open(args.SF_file)
json_content = {'dimension': 2, 'variables': ['AbsEta', 'Pt'], 'binning': {}, 'data': [], 'error_type': 'absolute'}
pt_binning = []
eta_binning = []

for key in f.GetListOfKeys():
    name = key.GetName()
    if not 'ZMass' in name:
        continue
    if not 'MC' in name:
        continue

    obj_MC = f.Get(name)
    obj_Data = f.Get(name.replace('MC','Data'))

    eta_val = re.findall(r"\d+p\d+",name)
    eta_val = [float(eta.replace('p','.')) for eta in eta_val]
    if len(eta_val) == 1:
        if 'Gt' in name:
            eta_val = eta_val + [2.5]
        elif 'Lt' in name:
            eta_val = [0.] + eta_val
        else:
            raise RuntimeError("No bound for eta")

    if eta_val[0] not in eta_binning:
        eta_binning.append(eta_val[0]) 
    if eta_val[1] not in eta_binning:
        eta_binning.append(eta_val[1]) 

    eta_dict = {'bin':eta_val,'values':[]}
    
    x = obj_MC.GetX()
    eff_MC = obj_MC.GetY()
    eff_Data = obj_Data.GetY()
    for i in range(0,obj_MC.GetN()):
        pt_center = x[i]
        pt_bin = [x[i]-obj_MC.GetErrorXlow(i),x[i]+obj_MC.GetErrorXhigh(i)]
        if pt_bin[0] not in pt_binning:
            pt_binning.append(pt_bin[0]) 
        if pt_bin[1] not in pt_binning:
            pt_binning.append(pt_bin[1]) 
        value = min(eff_Data[i]/max(1.e-6,eff_MC[i]),1.)
        eta_dict['values'].append({'bin':pt_bin,'value':value,'error_low':0.,'error_high':0})

    json_content['data'].append(eta_dict)

json_content['binning'] = {'x':sorted(eta_binning),'y':sorted(pt_binning)}


filename = 'TTH_trigger_%s.json' % (args.suffix)
with open(filename, 'w') as j:
    json.dump(json_content, j, indent=2)
print ("Saved json as %s"%filename)

f.Close()













