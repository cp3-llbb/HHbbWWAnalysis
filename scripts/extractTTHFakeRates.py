#! /bin/env python

import math
import argparse
import ROOT
import json

parser = argparse.ArgumentParser()
parser.add_argument('--file', required=True, type=str, help='ROOT file containing fake rates')
parser.add_argument('--era', required=True, type=str, help='Era')
parser.add_argument('--wp', required=True, type=str,choices=['Loose','Tight'], help='WP for the lepton selection : Loose | Tight')
parser.add_argument('--type', required=False, default='data', type=str, choices=['data','mc'], help='Type : data (actual datadriven) [default] | mc (derived from QCD, used to MC closure -> only for SL)')
parser.add_argument('--channel', required=True, type=str, choices=['DL','SL'], help='Lepton channel : DL | SL')

args = parser.parse_args()

assert args.wp == "Loose" or args.wp == "Tight"
assert args.channel == "DL" or args.channel == "SL"

f = ROOT.TFile.Open(args.file)

if args.wp == "Tight": 
    el_nom = f.Get('FR_mva080_el_data_comb_NC')
    el_list_up   = ['FR_mva080_el_data_comb_NC_pt1','FR_mva080_el_data_comb_NC_up','FR_mva080_el_data_comb_NC_be1']
    el_list_down = ['FR_mva080_el_data_comb_NC_pt2','FR_mva080_el_data_comb_NC_down','FR_mva080_el_data_comb_NC_be2']
    mu_nom = f.Get('FR_mva085_mu_data_comb')
    mu_list_up   = ['FR_mva085_mu_data_comb_pt1','FR_mva085_mu_data_comb_up','FR_mva085_mu_data_comb_be1']
    mu_list_down = ['FR_mva085_mu_data_comb_pt2','FR_mva085_mu_data_comb_down','FR_mva085_mu_data_comb_be2']
    syst_suffixes = ['ptSyst','normSyst','barrelSyst']
    if args.type == 'mc':
        raise RuntimeError('Tight WP with mc closure not implemented')
if args.wp == "Loose": 
    if args.type == 'data':
        el_nom = f.Get('FR_mva030_el_data_comb')
        el_list_up   = ['FR_mva030_el_data_comb_pt1','FR_mva030_el_data_comb_up','FR_mva030_el_data_comb_be1']
        el_list_down = ['FR_mva030_el_data_comb_pt2','FR_mva030_el_data_comb_down','FR_mva030_el_data_comb_be2']
        mu_nom = f.Get('FR_mva050_mu_data_comb')
        mu_list_up   = ['FR_mva050_mu_data_comb_pt1','FR_mva050_mu_data_comb_up','FR_mva050_mu_data_comb_be1']
        mu_list_down = ['FR_mva050_mu_data_comb_pt2','FR_mva050_mu_data_comb_down','FR_mva050_mu_data_comb_be2']
        syst_suffixes = ['ptSyst','normSyst','barrelSyst']
    if args.type == 'mc':
        el_nom = f.Get('FR_mva030_el_data_comb_QCD_fakes')
        el_list_up   = ['FR_mva030_el_data_comb_QCD_fakes']
        el_list_down = ['FR_mva030_el_data_comb_QCD_fakes']
        mu_nom = f.Get('FR_mva050_mu_data_comb_QCD_fakes')
        mu_list_up   = ['FR_mva050_mu_data_comb_QCD_fakes']
        mu_list_down = ['FR_mva050_mu_data_comb_QCD_fakes']
        syst_suffixes = ['QCD']
            
        

mu_hist_up = [f.Get(lup) for lup in mu_list_up]
mu_hist_down = [f.Get(ldown) for ldown in mu_list_down]
el_hist_up = [f.Get(lup) for lup in el_list_up]
el_hist_down = [f.Get(ldown) for ldown in el_list_down]

for cat,h_nom,h_ups,h_downs in zip(["Electron","Muon"],[el_nom,mu_nom],[el_hist_up,mu_hist_up],[el_hist_down,mu_hist_down]):
    xAxis = h_nom.GetXaxis()
    yAxis = h_nom.GetYaxis()
    for h_up, h_down, syst_suffix in zip(h_ups,h_downs,syst_suffixes):
        xbinning = []
        ybinning = []
        data = []
        for y in range(1,h_nom.GetNbinsY()+1):
            yMin = yAxis.GetBinLowEdge(y)
            yMax = yAxis.GetBinUpEdge(y)
            if yMin not in ybinning:
                ybinning.append(yMin)
            if yMax not in ybinning:
                ybinning.append(yMax)
            ydict = {'bin':[yMin,yMax],'values':[]}
            for x in range(1,h_nom.GetNbinsX()+1):
                nom = h_nom.GetBinContent(x,y)
                if args.channel == "SL":
                    nom = min(nom,0.8) # clipping
                xMin = xAxis.GetBinLowEdge(x)
                xMax = xAxis.GetBinUpEdge(x)
                if xMin not in xbinning:
                    xbinning.append(xMin)
                if xMax not in xbinning:
                    xbinning.append(xMax)
                up = h_up.GetBinContent(x,y)
                down = h_down.GetBinContent(x,y)
                err_up = abs(up-nom)/nom
                err_down = abs(nom-down)/nom
                ydict['values'].append({'bin':[xMin,xMax],
                                        'value':nom,
                                        'error_low':err_down,
                                        'error_high':err_up})
            data.append(ydict)

        json_content = {'dimension': 2, 'variables': ['AbsEta', 'Pt'], 'binning': {'x': ybinning, 'y': xbinning}, 'data': data, 'error_type': 'relative'}
            # Inverted binning because it is clearer to have eta bins containing conept bins

        filename = f'TTHFakeRates_{args.wp}MVA_{args.channel}_{cat}_{args.era}_{syst_suffix}.json'
        with open(filename, 'w') as j:                                                                   
            json.dump(json_content, j, indent=2)
        print ("Created file %s"%filename)



            

            
