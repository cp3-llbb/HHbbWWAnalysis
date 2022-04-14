import os
import sys
import copy
import ROOT
import ctypes

from IPython import embed

sys.path.append('../BambooDatacardProducer/')

from context import TFileOpen
from interpolation import Interpolation

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)


main_path = '/nfs/scratch/fynu/fbury/Datacards/bbww_dl/InterpolationClosure/test_true_{branch}_{mass}_2018'

def makePlots(canvas,pdfName,m1,m2,m3,branches,cat,samples,xlabel):
    assert float(m1) < float(m2)
    assert float(m2) < float(m3)
#    if all([float(m)>=500 for m in [m1,m2,m3]]):
#        branch = 'HighMass'
#    elif all([float(m)<=500 for m in [m1,m2,m3]]):
#        branch = 'LowMass'
#    else:
#        raise RuntimeError(f'No mix in branch: m1={m1}, m2={m2}, m3={m3}')

    f1 = os.path.join(main_path.format(branch=branches[0],mass=m1),cat.format(mass=m1))
    f2 = os.path.join(main_path.format(branch=branches[1],mass=m2),cat.format(mass=m2))
    f3 = os.path.join(main_path.format(branch=branches[2],mass=m3),cat.format(mass=m3))

    interpolate = Interpolation(float(m1),float(m3),float(m2))
    
    hs = None
    for sample in samples:
        with TFileOpen(f1,'r') as F1, \
             TFileOpen(f2,'r') as F2, \
             TFileOpen(f3,'r') as F3:
            h_tmp = [
                 copy.deepcopy(F1.Get(sample.format(mass=m1))),
                 copy.deepcopy(F2.Get(sample.format(mass=m2))),
                 copy.deepcopy(F3.Get(sample.format(mass=m3))),
            ]
            if hs is None:
                hs = h_tmp
            else:
                for i in range(len(hs)):
                    hs[i].Add(h_tmp[i])
                
    names = [
        f'True M_{{HH}} = {m1}',
        f'True M_{{HH}} = {m2}',
        f'True M_{{HH}} = {m3}',
    ]

    hs[0].SetLineColor(ROOT.kBlue+1)
    hs[1].SetLineColor(ROOT.kOrange+1)
    hs[2].SetLineColor(ROOT.kGreen+1)

    # Interpolation #
    hs.insert(2,interpolate(hs[0],hs[2],'interp'))
    names.insert(2,f'Interpolated M_{{HH}} = {m2}')
    hs[2].SetLineColor(ROOT.kRed+1)
    hs[2].SetLineStyle(2)


    ints = []
    errs = []
    cdfs = []

    hmax = 0.
    for h in hs:
        h.SetLineWidth(2)
        h.SetTitle(f"Shapes")
        h.GetXaxis().SetTitle(xlabel)
        cdf = copy.deepcopy(h)
        cdf.Scale(1./cdf.Integral())
        cdf = cdf.GetCumulative()
        cdf.SetTitle(f"CDF")
        cdf.GetXaxis().SetTitle(xlabel)
        cdf.GetYaxis().SetRangeUser(0,1)
        cdfs.append(cdf)
        err = ctypes.c_double(0.)
        ints.append(h.IntegralAndError(1,h.GetNbinsX(),err))
        errs.append(err.value)
        if h.GetNbinsX() > 200:
            h.Rebin(5)
        else:
            h.Rebin(1)

        hmax = max(hmax,h.GetMaximum())

    hs[0].SetMaximum(hmax*1.5)

    canvas.Clear()
    canvas.Divide(2)

    canvas.cd(1)
    #opt = 'HCP'
    opt = 'hist'
    for h in hs:
        h.Draw(opt)
        if 'same' not in opt:
            opt += ' same'
    leg1 = ROOT.TLegend(0.15,0.6,0.5,0.9)
    for i in range(0,len(hs)):
        leg1.AddEntry(hs[i],f"#splitline{{{names[i]}}}{{Yield = {ints[i]:5.2f} +- {errs[i]:5.2f}}}")
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0)
    leg1.Draw()

    canvas.cd(2)
    opt = ''
    leg2 = ROOT.TLegend(0.15,0.6,0.55,0.85)
    for i in range(0,len(hs)):
        cdfs[i].Draw(opt)
        if 'same' not in opt:
            opt += ' same'
        leg2.AddEntry(cdfs[i],names[i])
    leg2.SetBorderSize(0)
    leg2.SetFillStyle(0)
    leg2.Draw()

    canvas.Print(pdfName)

        
categories = ['HH_DL_{mass}_GGF_2018.root','HH_DL_{mass}_HME_2018.root']
xlabels = ['GGF output node','HME']
signals = ['signal_ggf_spin2_{mass}_hbbhww','signal_ggf_spin2_{mass}_hbbhtt']

masses = [
    ('260','270','280'),
    ('270','280','300'),
    ('280','300','350'),
    ('300','350','400'),
    ('350','400','450'),
    ('400','450','500'),
    ('500','550','600'),
    ('550','600','650'),
    ('600','650','700'),
    ('650','700','750'),
    ('700','750','800'),
    ('750','800','850'),
    ('800','850','900'),
]

branches = [
    ('LowMass','LowMass','LowMass'),
    ('LowMass','LowMass','LowMass'),
    ('LowMass','LowMass','LowMass'),
    ('LowMass','LowMass','LowMass'),
    ('LowMass','LowMass','LowMass'),
    ('LowMass','LowMass','LowMass'),
    ('HighMass','HighMass','HighMass'),
    ('HighMass','HighMass','HighMass'),
    ('HighMass','HighMass','HighMass'),
    ('HighMass','HighMass','HighMass'),
    ('HighMass','HighMass','HighMass'),
    ('HighMass','HighMass','HighMass'),
    ('HighMass','HighMass','HighMass'),
]

C = ROOT.TCanvas("C","C",1200,600)
pdfName = 'closure_interpolation.pdf'

C.Print(pdfName+'[')

for (m1,m2,m3),(b1,b2,b3) in zip(masses,branches):
    print ()
    print (m1,b1)
    print (m2,b2)
    print (m3,b3)
    for cat,xlabel in zip(categories,xlabels):
        makePlots(C,pdfName,m1,m2,m3,(b1,b2,b3),cat,signals,xlabel)

C.Print(pdfName+']')


