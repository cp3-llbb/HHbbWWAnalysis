import os
import ROOT
import glob

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

masses = [260,270,280,300,350,400,450,500,550,600,650,700,750,800,850,900]
pdfName = 'ResFix.pdf'

C = ROOT.TCanvas('C','C',2400,900)

C.Print (pdfName+'[')

path = "/nfs/scratch/fynu/fbury/Datacards/bbww_dl/TestInterpolation/ComparisonFix/datacard_fit_Resonant_1D_syst_M_{mass}_{type}_noSyst_FR2"
for mass in masses:

    masspath = path.replace('{mass}',str(mass))
    categories = [f'HH_DL_{mass}_resolved1b_GGF',
                  f'HH_DL_{mass}_resolved2b_GGF',
                  f'HH_DL_{mass}_boosted_GGF',
                  f'HH_DL_{mass}_boosted_other',
                  f'HH_DL_{mass}_resolved_other',
                  f'HH_DL_{mass}_inclusive_DY_VVV']

    samples = [f'ggHH_Graviton_M_{mass}_hbbhwwdl',
               f'ggHH_Radion_M_{mass}_hbbhwwdl']

    eras = ['2016','2017','2018']
    for cat in categories:
        for sample in samples:
            C.Clear()
            C.SetTopMargin(0.2)
            title = ROOT.TPaveText(0.15,1.,0.75,0.8)
            title.AddText(f'Sample {sample}, category {cat}')
            title.SetFillStyle(4000)
            title.SetBorderSize(0)
            C.Divide(3,1)
            Fs = []
            eras = ['2016','2017','2018']
            hs = {era:{era:None for era in eras} for era in eras} # first key is test era, second is sample era
            for test_era in ['2016','2017','2018']: # test era 
                for sample_era in ['2016','2017','2018']: # exact era
                    filepath = os.path.join(masspath,f'{cat}_{sample_era}.root').replace('{type}',f'fix_{test_era}')
                    if os.path.exists(filepath):
                        F = ROOT.TFile(filepath,'READ')
                        h = F.Get(sample)
                    else:
                        h = None
                    if h is not None:
                        #h.SetTitle(f'Sample era {x_era}, DNN input era {y_era}')
                        h.SetLineWidth(2)
                        h.Scale(1./h.Integral())
                        h.Rebin(5)
                        if sample_era == '2016':
                            h.SetLineColor(ROOT.TColor.GetColor('#648FFF'))
                        if sample_era == '2017':
                            h.SetLineColor(ROOT.TColor.GetColor('#785EF0'))
                        if sample_era == '2018':
                            h.SetLineColor(ROOT.TColor.GetColor('#DC267F'))
                    Fs.append(F)
                    hs[test_era][sample_era] = h


            leg = [None]*3
            for i,test_era in enumerate(eras):
                pad = C.cd(i+1)
                pad.SetTopMargin(0.15)
                if hs[test_era]['2016'] is not None and hs[test_era]['2017'] is not None and hs[test_era]['2018'] is not None:
                    leg[i] = ROOT.TLegend(0.15,0.6,0.45,0.82)
                    hmax = max([h.GetMaximum() for h in hs[test_era].values()])*1.1
                    hs[test_era]['2016'].SetTitle(f'Using {test_era} as DNN input')
                    hs[test_era]['2016'].SetMaximum(hmax)
                    hs[test_era]['2016'].Draw("hist")
                    hs[test_era]['2017'].Draw("hist same")
                    hs[test_era]['2018'].Draw("hist same")
                    for era in eras:
                        leg[i].AddEntry(hs['2016'][era],f'{era} sample')
                    leg[i].Draw()
 
            C.cd()
            title.Draw()
            C.Print (pdfName,f'Title:{sample}_{cat}')

            for F in Fs:
                F.Close()
                 

C.Print (pdfName+']')
