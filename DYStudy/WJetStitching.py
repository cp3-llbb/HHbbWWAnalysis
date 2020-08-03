import os
import sys
import ROOT
import pandas as pd
from pprint import pprint

path = sys.argv[1]

f_inc = ROOT.TFile(os.path.join(path,'WJetsToLNu.root'))
f_1j  = ROOT.TFile(os.path.join(path,'W1JetsToLNu.root'))
f_2j  = ROOT.TFile(os.path.join(path,'W2JetsToLNu.root'))
f_3j  = ROOT.TFile(os.path.join(path,'W3JetsToLNu.root'))
f_4j  = ROOT.TFile(os.path.join(path,'W4JetsToLNu.root'))
f_70HT = ROOT.TFile(os.path.join(path,'WJetsToLNu_HT-70To100.root'))
f_100HT = ROOT.TFile(os.path.join(path,'WJetsToLNu_HT-100To200.root'))
f_200HT = ROOT.TFile(os.path.join(path,'WJetsToLNu_HT-200To400.root'))
f_400HT = ROOT.TFile(os.path.join(path,'WJetsToLNu_HT-400To600.root'))
f_600HT = ROOT.TFile(os.path.join(path,'WJetsToLNu_HT-600To800.root'))
f_800HT = ROOT.TFile(os.path.join(path,'WJetsToLNu_HT-800To1200.root'))
f_1200HT = ROOT.TFile(os.path.join(path,'WJetsToLNu_HT-1200To2500.root'))
f_2500HT = ROOT.TFile(os.path.join(path,'WJetsToLNu_HT-2500ToInf.root'))

h_inc   = f_inc.Get("LHE_HTVSNjets")
h_1j    = f_1j.Get("LHE_HTVSNjets")
h_2j    = f_2j.Get("LHE_HTVSNjets")
h_3j    = f_3j.Get("LHE_HTVSNjets")
h_4j    = f_4j.Get("LHE_HTVSNjets")
h_70HT  = f_70HT.Get("LHE_HTVSNjets")
h_100HT = f_100HT.Get("LHE_HTVSNjets")
h_200HT = f_200HT.Get("LHE_HTVSNjets")
h_400HT = f_400HT.Get("LHE_HTVSNjets")
h_600HT = f_600HT.Get("LHE_HTVSNjets")
h_800HT = f_800HT.Get("LHE_HTVSNjets")
h_1200HT= f_1200HT.Get("LHE_HTVSNjets")
h_2500HT= f_2500HT.Get("LHE_HTVSNjets")

def returnContent(h):
    return {'0j0HT' :    h.GetBinContent(1,1),
            '0j70HT' :   h.GetBinContent(2,1),
            '0j100HT' :  h.GetBinContent(3,1),
            '0j200HT' :  h.GetBinContent(4,1),
            '0j400HT' :  h.GetBinContent(5,1),
            '0j600HT' :  h.GetBinContent(6,1),
            '0j800HT' :  h.GetBinContent(7,1),
            '0j1200HT' : h.GetBinContent(8,1),
            '0j2500HT' : h.GetBinContent(9,1),

            '1j0HT' :    h.GetBinContent(1,2),
            '1j70HT' :   h.GetBinContent(2,2),
            '1j100HT' :  h.GetBinContent(3,2),
            '1j200HT' :  h.GetBinContent(4,2),
            '1j400HT' :  h.GetBinContent(5,2),
            '1j600HT' :  h.GetBinContent(6,2),
            '1j800HT' :  h.GetBinContent(7,2),
            '1j1200HT' : h.GetBinContent(8,2),
            '1j2500HT' : h.GetBinContent(9,2),

            '2j0HT' :    h.GetBinContent(1,3),
            '2j70HT' :   h.GetBinContent(2,3),
            '2j100HT' :  h.GetBinContent(3,3),
            '2j200HT' :  h.GetBinContent(4,3),
            '2j400HT' :  h.GetBinContent(5,3),
            '2j600HT' :  h.GetBinContent(6,3),
            '2j800HT' :  h.GetBinContent(7,3),
            '2j1200HT' : h.GetBinContent(8,3),
            '2j2500HT' : h.GetBinContent(9,3),

            '3j0HT' :    h.GetBinContent(1,4),
            '3j70HT' :   h.GetBinContent(2,4),
            '3j100HT' :  h.GetBinContent(3,4),
            '3j200HT' :  h.GetBinContent(4,4),
            '3j400HT' :  h.GetBinContent(5,4),
            '3j600HT' :  h.GetBinContent(6,4),
            '3j800HT' :  h.GetBinContent(7,4),
            '3j1200HT' : h.GetBinContent(8,4),
            '3j2500HT' : h.GetBinContent(9,4),

            '4j0HT' :    h.GetBinContent(1,5),
            '4j70HT' :   h.GetBinContent(2,5),
            '4j100HT' :  h.GetBinContent(3,5),
            '4j200HT' :  h.GetBinContent(4,5),
            '4j400HT' :  h.GetBinContent(5,5),
            '4j600HT' :  h.GetBinContent(6,5),
            '4j800HT' :  h.GetBinContent(7,5),
            '4j1200HT' : h.GetBinContent(8,5),
            '4j2500HT' : h.GetBinContent(9,5),

            '5j0HT' :    h.GetBinContent(1,6),
            '5j70HT' :   h.GetBinContent(2,6),
            '5j100HT' :  h.GetBinContent(3,6),
            '5j200HT' :  h.GetBinContent(4,6),
            '5j400HT' :  h.GetBinContent(5,6),
            '5j600HT' :  h.GetBinContent(6,6),
            '5j800HT' :  h.GetBinContent(7,6),
            '5j1200HT' : h.GetBinContent(8,6),
            '5j2500HT' : h.GetBinContent(9,6),
            }

params = list(returnContent(h_inc).keys())

block0j = ['0j0HT','0j70HT','0j100HT','0j200HT','0j400HT','0j600HT','0j800HT','0j1200HT','0j2500HT']
block1j = ['1j0HT','1j70HT','1j100HT','1j200HT','1j400HT','1j600HT','1j800HT','1j1200HT','1j2500HT']
block2j = ['2j0HT','2j70HT','2j100HT','2j200HT','2j400HT','2j600HT','2j800HT','2j1200HT','2j2500HT']
block3j = ['3j0HT','3j70HT','3j100HT','3j200HT','3j400HT','3j600HT','3j800HT','3j1200HT','3j2500HT']
block4j = ['4j0HT','4j70HT','4j100HT','4j200HT','4j400HT','4j600HT','4j800HT','4j1200HT','4j2500HT']
block5j = ['5j0HT','5j70HT','5j100HT','5j200HT','5j400HT','5j600HT','5j800HT','5j1200HT','5j2500HT']
block0HT = ['0j0HT','1j0HT','2j0HT','3j0HT','4j0HT','5j0HT']
block70HT = ['0j70HT','1j70HT','2j70HT','3j70HT','4j70HT','5j70HT']
block100HT = ['0j100HT','1j100HT','2j100HT','3j100HT','4j100HT','5j100HT']
block200HT = ['0j200HT','1j200HT','2j200HT','3j200HT','4j200HT','5j200HT']
block400HT = ['0j400HT','1j400HT','2j400HT','3j400HT','4j400HT','5j400HT']
block600HT = ['0j600HT','1j600HT','2j600HT','3j600HT','4j600HT','5j600HT']
block800HT = ['0j800HT','1j800HT','2j800HT','3j800HT','4j800HT','5j800HT']
block1200HT = ['0j1200HT','1j1200HT','2j1200HT','3j1200HT','4j1200HT','5j1200HT']
block2500HT = ['0j2500HT','1j2500HT','2j2500HT','3j2500HT','4j2500HT','5j2500HT']

blocks = [block0j,block1j,block2j,block3j,block4j,block5j,block0HT,block70HT,block100HT,block200HT,block400HT,block600HT,block800HT,block1200HT,block2500HT]

dict_content = {'Inclusive' : returnContent(h_inc),
                'Exclusive 1J' : returnContent(h_1j),
                'Exclusive 2J' : returnContent(h_2j),
                'Exclusive 3J' : returnContent(h_3j),
                'Exclusive 4J' : returnContent(h_4j),
                'Exclusive HT [70,100]' : returnContent(h_70HT),
                'Exclusive HT [100,200]' : returnContent(h_100HT),
                'Exclusive HT [200,400]' : returnContent(h_200HT),
                'Exclusive HT [400,600]' : returnContent(h_400HT),
                'Exclusive HT [600,800]' : returnContent(h_600HT),
                'Exclusive HT [800,1200]' : returnContent(h_800HT),
                'Exclusive HT [1200,2500]' : returnContent(h_1200HT),
                'Exclusive HT [2500,inf[' : returnContent(h_2500HT)}

dict_xsec = {'Inclusive' : 61526.7,
             'Exclusive 1J' : 9442.49,
             'Exclusive 2J' : 3252.49,
             'Exclusive 3J' : 1153.42,
             'Exclusive 4J' : 634.05,
             'Exclusive HT [70,100]' : 1504.92,
             'Exclusive HT [100,200]' : 1625.08,
             'Exclusive HT [200,400]' : 477.96,
             'Exclusive HT [400,600]' : 67.441,
             'Exclusive HT [600,800]' : 15.096,
             'Exclusive HT [800,1200]' : 6.3626,
             'Exclusive HT [1200,2500]' : 1.2658,
             'Exclusive HT [2500,inf[' : 0.009405}

df = pd.DataFrame() 
for name, content in dict_content.items():
    df = df.append(pd.DataFrame([list(content.values())+[dict_xsec[name]]],columns=list(content.keys())+['xsec'],index=[name]))

df['total'] = df[params].sum(axis=1)

print ('Event weight sum per sample per bins')
for block in blocks:
    print (df[block+['xsec','total']])

print ('N per p parameter')
sums_params = {col:df[col].sum() for col in params}
pprint (sums_params)

print ('Xsec per p parameter')
xsec_params = {p:(df[p]*df['xsec']/df['total']).sum() for p in params}
for p in xsec_params.keys():
    if ('0j' in p and not '0HT' in p) or (not '0j' in p and '0HT' in p):
        xsec_params[p] /= 2
    elif not '0j' in p and not '0HT' in p:
        xsec_params[p] /= 3
pprint (xsec_params)

exclusiveParams = ['0j','1j','2j','3j','4j','0HT','70HT','100HT','200HT','400HT','600HT','800HT','1200HT','2500HT']
xsecExclu = {exp:sum([x for p,x in xsec_params.items() if exp in p]) for exp in exclusiveParams}
for p,x in xsecExclu.items():
    df_inc = df[[col for col in df.columns if p in col]].loc['Inclusive']
    correct_xsec = dict_xsec['Inclusive'] * df[[col for col in df.columns if p in col]].loc['Inclusive'].sum() / df['total']['Inclusive']
    print ("Exclusive sample %s, Xsec %0.3f - Correct %0.3f"%(p,x,correct_xsec))
    
print ('Total Xsec per param %0.3f - Inclusive Xsec %0.3f'%(sum(xsec_params.values()),dict_xsec['Inclusive']))


print ('Lumiscale weight per sample')
df_w_sample = df['xsec']/df['total']
df_w_sample.columns = ['weight_per_sample']
print (df_w_sample)

print ('Lumiscale weights per jet multiplicity + HT bins')
w_params = {p:xsec_params[p]/sums_params[p] if sums_params[p]!=0 else 0. for p in sums_params.keys()}
df_w_params = pd.DataFrame(list(w_params.values()), columns=['weight_per_parameter'],index=list(w_params.keys()))
print (df_w_params)

print ("Stitching weights")
df_w_stitching = pd.DataFrame([[0]*df_w_sample.shape[0]]*df_w_params.shape[0],index=df_w_params.index,columns=df_w_sample.index)
for col in df_w_stitching.columns:
    for row in df_w_stitching.index:
        wp = df_w_params.loc[row]
        ws = df_w_sample.loc[col]
        df_w_stitching[col][row] = wp/ws if ws != 0 else 0
print (df_w_stitching)
print (df_w_sample.sum(),df_w_params.sum())
