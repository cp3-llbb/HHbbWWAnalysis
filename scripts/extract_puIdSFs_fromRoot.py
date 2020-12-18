#! /bin/env python

import math
import argparse
import ROOT
import json

IGNORE_LAST_PT_BIN = False

def format_eta_bin(eta_bin):
    return 'ptabseta<%.1f' % (eta_bin[1]) if (eta_bin[0] == 0) else 'ptabseta%.1f-%.1f' % (eta_bin[0], eta_bin[1])

def getBinning(h):
    eta_binning = []
    for i in range(1, h.GetYaxis().GetNbins() + 1):
        if len(eta_binning) == 0:
            eta_binning.append(h.GetYaxis().GetBinLowEdge(i))
            eta_binning.append(h.GetYaxis().GetBinUpEdge(i))
        else:
            eta_binning.append(h.GetYaxis().GetBinUpEdge(i))

    pt_binning = []
    for i in range(1, h.GetXaxis().GetNbins() + 1):
        if len(pt_binning) == 0:
            pt_binning.append(h.GetXaxis().GetBinLowEdge(i))
            pt_binning.append(h.GetXaxis().GetBinUpEdge(i))
        else:
            pt_binning.append(h.GetXaxis().GetBinUpEdge(i))

    if IGNORE_LAST_PT_BIN and len(pt_binning) > 2:
        pt_binning.pop()

    return eta_binning,pt_binning

parser = argparse.ArgumentParser()
parser.add_argument('--eff', required=True, help='ROOT file containing efficiency and mistag of the PU jet ID')
parser.add_argument('--sf', required=True, help='ROOT file containing SF for efficiency and mistag + systematics')

args = parser.parse_args()



for fargs,fName in zip([args.eff,args.sf],['EFF','SF']):
    f = ROOT.TFile.Open(fargs)
    keys = [key.GetName() for key in f.GetListOfKeys()]
    for key in keys:
        if 'Systuncty' in key:
            continue

        print (key)
        h = f.Get(key)
        h_syst = None
        if key+'_Systuncty' in keys:
            h_syst = f.Get(key+'_Systuncty')

        eta_binning,pt_binning = getBinning(h)
        if h_syst:
            eta_binning_syst,pt_binning_syst = getBinning(h_syst)
            assert eta_binning_syst == eta_binning
            assert pt_binning_syst == pt_binning

        eta = 'Eta' if eta_binning[0] < 0 else 'AbsEta'
        json_content = {'dimension': 2, 'variables': [eta, 'Pt'], 'binning': {'x': eta_binning, 'y': pt_binning}, 'data': [], 'error_type': 'absolute'}
        json_content_data = json_content['data']

        for y in range(0, len(eta_binning) - 1):
            eta_data = {'bin': [eta_binning[y], eta_binning[y + 1]], 'values': []}
            mean_eta = (eta_binning[y] + eta_binning[y + 1]) / 2.
            for x in range(0, len(pt_binning) - 1):
                mean_pt = (pt_binning[x] + pt_binning[x + 1]) / 2.
                bin = h.FindBin(mean_pt, mean_eta)
                if h_syst:
                    error = h_syst.GetBinContent(bin)
                else:
                    error = 0.
                value = h.GetBinContent(bin)
                if value == 0:
                    print ('value 0',mean_eta,mean_pt)

                pt_data = {'bin': [pt_binning[x], pt_binning[x + 1]], 'value': value, 'error_low': error, 'error_high': error}

                eta_data['values'].append(pt_data)

            json_content_data.append(eta_data)

      
        filename = 'PUID_%s_%s.json' % (fName,key)
        print (filename)
        with open(filename, 'w') as j:
            json.dump(json_content, j, indent=2)
    
    f.Close()
