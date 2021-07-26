#! /bin/env python

import math
import argparse
import json
import csv
import sys
import copy

def operating_point_to_string(operating_point):
    if operating_point == 0:
        return "loose"
    elif operating_point == 1:
        return "medium"
    elif operating_point == 2:
        return "tight"
    elif operating_point == 3:
        return "discr_reshaping"
    else:
        return "unknown"

def jet_flavor_to_string(jet_flavor):
    if jet_flavor == 0:
        return "bjets"
    elif jet_flavor == 1:
        return "cjets"
    elif jet_flavor == 2:
        return "lightjets"
    else:
        return "unknown"

valid_syst_types = {0 : ["jes", "lf", "hf", "hfstats1", "hfstats2", "lfstats1" , "lfstats2"],
                    1 : ["cferr1" , "cferr2"],
                    2 : ["jes", "lf", "hf", "hfstats1", "hfstats2", "lfstats1" , "lfstats2"]}
           

parser = argparse.ArgumentParser()
parser.add_argument('file', help='CSV file containing b-tagging scale factors')
parser.add_argument('-s', '--suffix', help='Suffix to append at the end of the output filename', required=True)

args = parser.parse_args()

# See https://twiki.cern.ch/twiki/bin/view/CMS/BTagCalibration for CSV format definition

all_content = {}

def get_token(syst_type, jet_flavor):
    return (syst_type, jet_flavor)

with open(args.file, 'r') as f:
    data = csv.reader(f, skipinitialspace=True)

    operating_point = None # Can be loose = 0, medium = 1, tight = 2, discriminator reshaping = 3
    measurement_type = None  # A string representing the type of measurement done to obtain the scale factors
    jet_flavor = None
    for row in data:
        if 'formula' in row:
            continue # Avoid potential header
        if len(row) != 11:
            continue
        operating_point = int(row[0])
        measurement_type = row[1]
        syst_type = row[2]
        jet_flavor = int(row[3])

        
        if syst_type != 'central' and not syst_type.replace('up_','').replace('down_','') in valid_syst_types[jet_flavor]:
            continue
        if measurement_type != 'iterativefit':
            continue

        token = get_token(syst_type, jet_flavor)


        eta_bin = (float(row[4]), float(row[5]))
        pt_bin = (float(row[6]), float(row[7]))
        discr_bin = (float(row[8]), float(row[9]))
        formula = row[10].strip()

        if token not in all_content.keys():
            all_content[token] = {'Eta':[eta_bin],'Pt':[pt_bin],'BTagDiscri':[discr_bin],'formula':[formula]}
        else:
            all_content[token]['Eta'].append(eta_bin)
            all_content[token]['Pt'].append(pt_bin)
            all_content[token]['BTagDiscri'].append(discr_bin)
            all_content[token]['formula'].append(formula)

def makeJson(jet_flavor,syst_type,tokens,all_content):
    json_content = {'dimension': 3, 'variables': ['Eta', 'Pt', 'BTagDiscri'], 'binning': {'x': [], 'y': [], 'z': []}, 'data': [], 'error_type': 'variated', 'formula': True, 'variable': 'z'}
    
    # Check if up and down have same binning #
    for key in ['Eta', 'Pt', 'BTagDiscri']:
        for up_bin, down_bin in zip(all_content[tokens['up']][key],all_content[tokens['down']][key]):
            if up_bin!=down_bin:
                raise RuntimeError("Up and down fluctuations of syst %s do not have same binning"%syst_type)

    # Check if central as same binning or only one large bin containing all #
    central_one_bin = len(all_content[tokens['central']]['Eta']) == 1 
    if not central_one_bin:
        for key in ['Eta', 'Pt', 'BTagDiscri']:
            for up_bin, down_bin in zip(all_content[tokens['up']][key],all_content[tokens['central']][key]):
                if up_bin!=down_bin:
                    raise RuntimeError("Up and down fluctuations of syst %s do not have same binning as central"%syst_type)

    xbins = []
    ybins = []
    zbins = []
    data = {}
    for idx in range(len(all_content[tokens['up']]['Eta'])):
        eta_bin   = all_content[tokens['up']]['Eta'][idx]
        pt_bin    = all_content[tokens['up']]['Pt'][idx]
        discr_bin = all_content[tokens['up']]['BTagDiscri'][idx]
        formula_central = all_content[tokens['central']]['formula'][idx] if not central_one_bin else  all_content[tokens['central']]['formula'][0]
        formula_up      = all_content[tokens['up']]['formula'][idx]
        formula_down    = all_content[tokens['down']]['formula'][idx]

        if eta_bin[0] not in xbins:
            xbins.append(eta_bin[0])
        if eta_bin[1] not in xbins:
            xbins.append(eta_bin[1])
        if pt_bin[0] not in ybins:
            ybins.append(pt_bin[0])
        if pt_bin[1] not in ybins:
            ybins.append(pt_bin[1])
        if discr_bin[0] not in zbins:
            zbins.append(discr_bin[0])
        if discr_bin[1] not in zbins:
            zbins.append(discr_bin[1])

        form_dict = {'value':formula_central,'error_high':formula_up,'error_low':formula_down}
        if eta_bin not in data.keys():
            data[eta_bin] = {}
            if pt_bin not in data[eta_bin].keys():
                data[eta_bin][pt_bin] = {}
                if discr_bin not in data[eta_bin][pt_bin].keys():
                    data[eta_bin][pt_bin][discr_bin] = form_dict
                else:
                    raise RuntimeError("Two formulas for same bin") 
        else:
            if pt_bin not in data[eta_bin].keys():
                data[eta_bin][pt_bin] = {}
            else:
                if discr_bin not in data[eta_bin][pt_bin].keys():
                    data[eta_bin][pt_bin][discr_bin] = form_dict
                else:
                    raise RuntimeError("Two formulas for same bin") 

    xbins = sorted(xbins)
    ybins = sorted(ybins)
    zbins = sorted(zbins)
    if all(x>=0 for x in xbins):
        json_content['variables'] = ['AbsEta', 'Pt', 'BTagDiscri']
        
    for eta_bin, pt_disc in data.items():
        eta_list = []
        for ptbin, disc in pt_disc.items():
            pt_dict = {'bins':ptbin,'values':[{'bins':disc_bin,**disc_val} for disc_bin,disc_val in disc.items()]}
            eta_list.append(pt_dict)
        json_content['data'].append({'bin':eta_bin, 'values':eta_list})

    json_content['binning']['x'] = xbins 
    json_content['binning']['y'] = ybins 
    json_content['binning']['z'] = zbins 
                    
    return json_content


for jet_flavor, syst_types in valid_syst_types.items():
    token_central = ('central',jet_flavor)
    for syst_type in syst_types:
        token_up = ('up_'+syst_type,jet_flavor)
        token_down = ('down_'+syst_type,jet_flavor)
        token_comb = {'central':token_central , 'up':token_up , 'down':token_down}
        json_content = makeJson(jet_flavor,syst_type,token_comb,all_content)

        filename = 'BTagging_%s_%s_%s_%s.json'%(operating_point_to_string(operating_point), jet_flavor_to_string(jet_flavor), syst_type , args.suffix)
        with open(filename, 'w') as j:
            json.dump(json_content, j, indent=2)
        print ("Saved SF as %s"%filename)
