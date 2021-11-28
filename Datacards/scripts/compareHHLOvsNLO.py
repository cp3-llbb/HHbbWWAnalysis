import os
import sys
import math
import ROOT
import ctypes

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

exact = '/nfs/scratch/fynu/fbury/Datacards/bbww_dl/NonResonant/datacard_fit_NonResonant_DNN11_GGF_raw_FR2' # Exact GGF NLO samples
GGFLO = '/nfs/scratch/fynu/fbury/Datacards/bbww_dl/NonResonant/datacard_fit_NonResonant_DNN11_GGFLO_raw_FR2' # GGF LO -> GGF NLO samples
GGFNLO = '/nfs/scratch/fynu/fbury/Datacards/bbww_dl/NonResonant/datacard_fit_NonResonant_DNN11_GGFNLO_raw_FR2' # GGF NLO -> GGF NLO samples

categories = [
 'HH_resolved1b_GGF',
 'HH_resolved1b_VBF',
 'HH_resolved2b_GGF',
 'HH_resolved2b_VBF',
 'HH_boosted_GGF',
 'HH_boosted_VBF',
 'HH_inclusive_DY_VVV',
 'HH_resolved_other',
 'HH_boosted_other',
# 'HH_resolved_l1_pt',
# 'HH_boosted_l1_pt',
# 'HH_resolved_met_pt',
# 'HH_boosted_met_pt',
# 'HH_resolved_b1_pt',
# 'HH_resolved_M_ll',
# 'HH_resolved_M_ll',
# 'HH_resolved_DR_ll',
# 'HH_boosted_DR_ll',
# 'HH_resolved_DR_bb',
# 'HH_resolved_DR_llbb',
# 'HH_resolved_HT',
# 'HH_resolved_M_HH',
# 'HH_boosted_M_HH',
# 'HH_boosted_B_pt',
# 'HH_boosted_B_M',
# 'HH_boosted_DR_llB',
]

eras = ['2016','2017','2018']

samples = [
  'ggHH_kl_0_kt_1_hbbhwwdl',
  'ggHH_kl_1_kt_1_hbbhwwdl',
  'ggHH_kl_2p45_kt_1_hbbhwwdl',
  'ggHH_kl_5_kt_1_hbbhwwdl',
  'ggHH_kl_0_kt_1_hbbhtt',
  'ggHH_kl_1_kt_1_hbbhtt',
  'ggHH_kl_2p45_kt_1_hbbhtt',
  'ggHH_kl_5_kt_1_hbbhtt',
]

rebin = 4

yields = {era:{sample:{name: [0.,0.] for name in ['exact','GGFLO','GGFNLO']} for sample in samples} for era in eras}

C = ROOT.TCanvas("C","C",800,600)
pdfName = f"NLO_samples_FR2.pdf"
C.Print(pdfName+'[')

for sample in samples:
    for cat in categories:
        for era in eras:
            catEra = f'{cat}_{era}.root'
            exactFile = ROOT.TFile(os.path.join(exact,catEra),'READ')
            GGFLOFile = ROOT.TFile(os.path.join(GGFLO,catEra),'READ')
            GGFNLOFile = ROOT.TFile(os.path.join(GGFNLO,catEra),'READ')

            h_exact = exactFile.Get(sample)
            h_GGFLO = GGFLOFile.Get(sample)
            h_GGFNLO = GGFNLOFile.Get(sample)

            if h_exact.GetNbinsX() > 100:
                h_exact.Rebin(rebin)
                h_GGFLO.Rebin(rebin)
                h_GGFNLO.Rebin(rebin)

            err = ctypes.c_double(0.)
            yields[era][sample]['exact'][0] += h_exact.IntegralAndError(1,h_exact.GetNbinsX(),err)
            yields[era][sample]['exact'][1] += err.value**2
            yields[era][sample]['GGFLO'][0] += h_GGFLO.IntegralAndError(1,h_GGFLO.GetNbinsX(),err)
            yields[era][sample]['GGFLO'][1] += err.value**2
            yields[era][sample]['GGFNLO'][0] += h_GGFNLO.IntegralAndError(1,h_GGFNLO.GetNbinsX(),err)
            yields[era][sample]['GGFNLO'][1] += err.value**2

            #h_exact.Scale(1/h_exact.Integral())
            #h_GGFLO.Scale(1/h_GGFLO.Integral())
            #h_GGFNLO.Scale(1/h_GGFNLO.Integral())

            h_exact.SetLineWidth(2)
            h_GGFLO.SetLineWidth(2)
            h_GGFNLO.SetLineWidth(2)
            h_exact.SetLineColor(ROOT.kGreen+1)
            h_GGFLO.SetLineColor(ROOT.kRed+1)
            h_GGFNLO.SetLineColor(ROOT.kBlue+1)
            h_exact.SetMarkerColor(ROOT.kGreen+1)
            h_GGFLO.SetMarkerColor(ROOT.kRed+1)
            h_GGFNLO.SetMarkerColor(ROOT.kBlue+1)
            h_exact.SetMarkerStyle(25)
            h_GGFLO.SetMarkerStyle(8)
            h_GGFNLO.SetMarkerStyle(4)
            

            C.Clear()
            pad1 = ROOT.TPad("pad1", "pad1", 0.0,0.3,1.0,1.0)
            pad2 = ROOT.TPad("pad2", "pad2", 0.0,0.0,1.0,0.3)
            
            pad1.SetRightMargin(0.05)
            pad1.SetLeftMargin(0.15)
            pad1.SetTopMargin(0.12)
            pad1.SetBottomMargin(0.02)
            #pad1.SetLogy()

            pad2.SetRightMargin(0.05)
            pad2.SetLeftMargin(0.15)
            pad2.SetTopMargin(0.02)
            pad2.SetBottomMargin(0.30)

            
            pad1.Draw()
            pad2.Draw()

            # PAD 1 #
            pad1.cd()

            h_exact.Draw("P0")
            h_GGFLO.Draw("P0 same")
            h_GGFNLO.Draw("P0 same")

            ymax = max([h.GetMaximum() for h in [h_exact,h_GGFLO,h_GGFNLO]]) * 1.5
            h_exact.SetTitle(f'Sample = {sample}, category = {cat}, era = {era}')
            h_exact.GetYaxis().SetTitle('Yield')
            h_exact.GetYaxis().SetRangeUser(0.,ymax)
            h_exact.GetYaxis().SetLabelFont(43)
            h_exact.GetYaxis().SetLabelSize(14)
            h_exact.GetYaxis().SetTitleOffset(1.5)
            h_exact.GetYaxis().SetTitleFont(43)
            h_exact.GetYaxis().SetTitleSize(24)
            h_exact.GetXaxis().SetLabelFont(43)
            h_exact.GetXaxis().SetLabelSize(0)
            h_exact.GetXaxis().SetTitleFont(43)
            h_exact.GetXaxis().SetTitleSize(0)


            leg_up = ROOT.TLegend(0.60,0.6,0.92,0.85)
            leg_up.AddEntry(h_exact,f'#splitline{{Exact NLO sample}}{{Yield = {h_exact.Integral():6.2f}}}')
            leg_up.AddEntry(h_GGFLO,f'#splitline{{Reweighted from LO samples}}{{Yield = {h_GGFLO.Integral():6.2f}}}')
            leg_up.AddEntry(h_GGFNLO,f'#splitline{{Reweighted from NLO samples}}{{Yield = {h_GGFNLO.Integral():6.2f}}}')
            leg_up.Draw()

            # PAD 2 #
            pad2.cd()
            pad2.SetGridy()

            ratio_GGFLO = h_GGFLO.Clone()
            ratio_GGFLO.Add(h_exact,-1)
            ratio_GGFLO.Divide(h_exact)

            ratio_GGFNLO = h_GGFNLO.Clone()
            ratio_GGFNLO.Add(h_exact,-1)
            ratio_GGFNLO.Divide(h_exact)

            ratio_GGFLO.SetLineWidth(2)
            ratio_GGFNLO.SetLineWidth(2)
            ratio_GGFNLO.SetLineColor(ROOT.kRed+1)
            ratio_GGFLO.SetLineColor(ROOT.kBlue+1)

            ratio_GGFLO.SetTitle('')
            ratio_GGFLO.GetYaxis().SetTitle('#frac{h^{reweighted}-h^{exact}}{h^{exact}}')
            ratio_GGFLO.GetYaxis().SetRangeUser(-2,2)
            ratio_GGFLO.GetYaxis().SetLabelFont(43)
            ratio_GGFLO.GetYaxis().SetLabelSize(14)
            ratio_GGFLO.GetYaxis().SetTitleOffset(1.2)
            ratio_GGFLO.GetYaxis().SetTitleFont(43)
            ratio_GGFLO.GetYaxis().SetTitleSize(24)
            ratio_GGFLO.GetXaxis().SetLabelFont(43)
            ratio_GGFLO.GetXaxis().SetLabelSize(14)
            ratio_GGFLO.GetXaxis().SetTitleOffset(3)
            ratio_GGFLO.GetXaxis().SetTitleFont(43)
            ratio_GGFLO.GetXaxis().SetTitleSize(24)
#            from IPython import embed
#            embed()

            ratio_GGFLO.Draw("C")
            ratio_GGFNLO.Draw("C same")

            leg_down = ROOT.TLegend(0.6,0.75,0.9,0.95)
            leg_down.AddEntry(h_GGFLO,'Reweighted from LO samples')
            leg_down.AddEntry(h_GGFNLO,'Reweighted from NLO samples')
            leg_down.Draw()

            C.Print(pdfName,f'Title:{sample}_{cat}_{era}')

            exactFile.Close()
            GGFLOFile.Close()
            GGFNLOFile.Close()

C.Print(pdfName+']')

for era in yields.keys():
    print (f'Era {era}')
    for sample in yields[era].keys():
        print (f'... Sample : {sample}')
        for name in yields[era][sample].keys():
            yields[era][sample][name][1] = math.sqrt(yields[era][sample][name][1])
            print (f'       -> {name:10s} inclusive yield = {yields[era][sample][name][0]:6.2f} +/- {yields[era][sample][name][1]:6.2f}')

