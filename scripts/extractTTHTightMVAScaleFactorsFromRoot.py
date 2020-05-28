#! /bin/env python

import math
import argparse
import ROOT

def format_eta_bin(eta_bin):
    return 'ptabseta<%.1f' % (eta_bin[1]) if (eta_bin[0] == 0) else 'ptabseta%.1f-%.1f' % (eta_bin[0], eta_bin[1])

parser = argparse.ArgumentParser()
#parser.add_argument('file', help='ROOT file containing electron scale factors')
parser.add_argument('--SF_file', help='ROOT file containing nominal scale factors')
parser.add_argument('--err_eta', help='ROOT file containing systematic errors in eta')
parser.add_argument('--err_pt', help='ROOT file containing systematic errors in pt')
parser.add_argument('-s', '--suffix', help='Suffix to append at the end of the output filename', required=True)

args = parser.parse_args()

# Moriond17: last pt bin is the same as the previous one, but with 100% error and should be ignored
IGNORE_LAST_PT_BIN = False # Not needed for ttH

f = ROOT.TFile.Open(args.SF_file)

f_eta = ROOT.TFile.Open(args.err_eta)
h_eta = f_eta.Get("histo_eff_data")
f_pt = ROOT.TFile.Open(args.err_pt)
h_pt = f_pt.Get("histo_eff_data")

for key in f.GetListOfKeys():
    wp = key.GetName()

    if not 'SF' in wp:
        continue

    if not '2D' in wp:
        continue

    h = f.Get(wp)

    # Get binning
    eta_binning = []
    for i in range(1, h.GetXaxis().GetNbins() + 1):
        if len(eta_binning) == 0:
            eta_binning.append(h.GetXaxis().GetBinLowEdge(i))
            eta_binning.append(h.GetXaxis().GetBinUpEdge(i))
        else:
            eta_binning.append(h.GetXaxis().GetBinUpEdge(i))

    pt_binning = []
    for i in range(1, h.GetYaxis().GetNbins() + 1):
        if len(pt_binning) == 0:
            pt_binning.append(h.GetYaxis().GetBinLowEdge(i))
            pt_binning.append(h.GetYaxis().GetBinUpEdge(i))
        else:
            pt_binning.append(h.GetYaxis().GetBinUpEdge(i))

    if IGNORE_LAST_PT_BIN and len(pt_binning) > 2:
        pt_binning.pop()

    eta = 'Eta' if eta_binning[0] < 0 else 'AbsEta'
    json_content = {'dimension': 2, 'variables': [eta, 'Pt'], 'binning': {'x': eta_binning, 'y': pt_binning}, 'data': [], 'error_type': 'absolute'}
    json_content_data = json_content['data']

    for i in range(0, len(eta_binning) - 1):
        eta_data = {'bin': [eta_binning[i], eta_binning[i + 1]], 'values': []}
        mean_eta = (eta_binning[i] + eta_binning[i + 1]) / 2.
        for j in range(0, len(pt_binning) - 1):
            mean_pt = (pt_binning[j] + pt_binning[j + 1]) / 2.
            bin_SF = h.FindBin(mean_eta, mean_pt)
            val_SF = h.GetBinContent(bin_SF)
            bin_eta = h_eta.FindBin(mean_eta)
            bin_pt = h_pt.FindBin(mean_pt)
            shift_eta = abs(1-h_eta.GetBinContent(bin_eta))
            shift_pt = abs(1-h_pt.GetBinContent(bin_pt))
            val_SF_up = max(1+shift_eta,1+shift_pt)*val_SF
            val_SF_down = min(1-shift_eta,1-shift_pt)*val_SF
            #error_tot = math.sqrt(error_eta**2+error_pt**2)
            pt_data = {'bin': [pt_binning[j], pt_binning[j + 1]], 'value': val_SF, 'error_low': val_SF-val_SF_down, 'error_high': val_SF_up-val_SF}

            eta_data['values'].append(pt_data)

        json_content_data.append(eta_data)

    # Save JSON file
    filename = 'TTHSF_%s_%s.json' % (wp, args.suffix)
    with open(filename, 'w') as j:
        import json
        json.dump(json_content, j, indent=2)
