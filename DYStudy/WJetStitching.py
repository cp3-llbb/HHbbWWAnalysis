import os
import sys
import ROOT
import pandas as pd
import json
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
            }

params = list(returnContent(h_inc).keys())

block0j = ['0j0HT','0j70HT','0j100HT','0j200HT','0j400HT','0j600HT','0j800HT','0j1200HT','0j2500HT']
block1j = ['1j0HT','1j70HT','1j100HT','1j200HT','1j400HT','1j600HT','1j800HT','1j1200HT','1j2500HT']
block2j = ['2j0HT','2j70HT','2j100HT','2j200HT','2j400HT','2j600HT','2j800HT','2j1200HT','2j2500HT']
block3j = ['3j0HT','3j70HT','3j100HT','3j200HT','3j400HT','3j600HT','3j800HT','3j1200HT','3j2500HT']
block4j = ['4j0HT','4j70HT','4j100HT','4j200HT','4j400HT','4j600HT','4j800HT','4j1200HT','4j2500HT']
block0HT = ['0j0HT','1j0HT','2j0HT','3j0HT','4j0HT']
block70HT = ['0j70HT','1j70HT','2j70HT','3j70HT','4j70HT']
block100HT = ['0j100HT','1j100HT','2j100HT','3j100HT','4j100HT']
block200HT = ['0j200HT','1j200HT','2j200HT','3j200HT','4j200HT']
block400HT = ['0j400HT','1j400HT','2j400HT','3j400HT','4j400HT']
block600HT = ['0j600HT','1j600HT','2j600HT','3j600HT','4j600HT']
block800HT = ['0j800HT','1j800HT','2j800HT','3j800HT','4j800HT']
block1200HT = ['0j1200HT','1j1200HT','2j1200HT','3j1200HT','4j1200HT']
block2500HT = ['0j2500HT','1j2500HT','2j2500HT','3j2500HT','4j2500HT']

blocks = [block0j,block1j,block2j,block3j,block4j,block0HT,block70HT,block100HT,block200HT,block400HT,block600HT,block800HT,block1200HT,block2500HT]

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

dict_samples = {'Inclusive' : 'WJetsToLNu',
                'Exclusive 1J' : 'W1JetsToLNu',
                'Exclusive 2J' : 'W2JetsToLNu',
                'Exclusive 3J' : 'W3JetsToLNu',
                'Exclusive 4J' : 'W4JetsToLNu',
                'Exclusive HT [70,100]' : 'WJetsToLNu_HT-70To100',
                'Exclusive HT [100,200]' : 'WJetsToLNu_HT-100To200',
                'Exclusive HT [200,400]' : 'WJetsToLNu_HT-200To400',
                'Exclusive HT [400,600]' : 'WJetsToLNu_HT-400To600',
                'Exclusive HT [600,800]' : 'WJetsToLNu_HT-600To800',
                'Exclusive HT [800,1200]' : 'WJetsToLNu_HT-800To1200',
                'Exclusive HT [1200,2500]' : 'WJetsToLNu_HT-1200To2500',
                'Exclusive HT [2500,inf[' : 'WJetsToLNu_HT-2500ToInf'}

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

samples = dict_xsec.keys()

df = pd.DataFrame() 
for name, content in dict_content.items():
    df = df.append(pd.DataFrame([list(content.values())+[dict_xsec[name]]],columns=list(content.keys())+['xsec'],index=[name]))

df['total'] = df[params].sum(axis=1)

if len(sys.argv)>2:
    df.to_excel(sys.argv[2]+'_N.xls',float_format='%0.10f')

print ('Event weight sum per sample per bins')
for block in blocks:
    print (df[block+['xsec','total']])

print ('N per p parameter')
sums_params = {col:df[col].sum(axis=0) for col in params}

print ('Xsec per parameter')
xsec_params = {p:(df[p]*df['xsec']/df['total']).sum() for p in params}
for p in xsec_params.keys():
    pjet = p[:2]
    pHT  = p[2:]
    if ('0j' == pjet and not '0HT' == pHT) or (not '0j' == pjet and '0HT' == pHT):
        xsec_params[p] /= 2
    elif not '0j' == pjet and not '0HT' == pHT:
        xsec_params[p] /= 3
df_xsec_params = pd.DataFrame(list(xsec_params.values()), index=xsec_params.keys(), columns=['Xsec per parameter'])
for p in params:
    print ("Parameter %s, xsec = %0.3f"%(p,xsec_params[p]))
if len(sys.argv)>2:
    df_xsec_params.to_excel(sys.argv[2]+'_xsec.xls',float_format='%0.10f')

print ('Xsec per exclusif parameter')
xsec_test = {}
exclusiveParams = ['0j','1j','2j','3j','4j','0HT','70HT','100HT','200HT','400HT','600HT','800HT','1200HT','2500HT']
for p in exclusiveParams:
    col = [c for c in df.columns if p == c[:2] or p == c[2:]]
    xsec_test[p] = df[col].loc['Inclusive'].sum()/df['total'].loc['Inclusive']*dict_xsec['Inclusive']
print ("Xsec sum over jet bins :",sum([v for k,v in xsec_test.items() if 'j' in k]))
print ("Xsec sum over HT  bins :",sum([v for k,v in xsec_test.items() if 'HT' in k]))
print ("Inclusive xsec",dict_xsec['Inclusive'])
    
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
df_w_stitching = pd.DataFrame([[0.]*df_w_sample.shape[0]]*df_w_params.shape[0],index=df_w_params.index,columns=df_w_sample.index)
for col in df_w_stitching.columns:
    for row in df_w_stitching.index:
        wp = df_w_params.loc[row]
        ws = df_w_sample.loc[col]
        df_w_stitching[col][row] = float(wp.values[0])/ws if ws != 0 else 0
print (df_w_stitching)
conv_jet = {'0j':0,'1j':1,'2j':2,'3j':3,'4j':4,'5j':5}
conv_HT  = {'0HT':(0,70),'70HT':(70,100),'100HT':(100,200),'200HT':(200,400),'400HT':(400,600),'600HT':(600,800),'800HT':(800,1200),'1200HT':(1200,2500),'2500HT':(2500,1000000)}
dict_weights = {}
for col in df_w_stitching.columns:
    sample = dict_samples[col]
    dict_weights[sample] = {}
    for index in df_w_stitching.index:
        jetBin = conv_jet[index[:2]]
        HTBin = conv_HT[index[2:]]
        if not jetBin in dict_weights[sample].keys():
            dict_weights[sample][jetBin] = []
        dict_weights[sample][jetBin].append({'low':HTBin[0],'up':HTBin[1] if len(HTBin)>1 else HTBin[0],'value':df_w_stitching[col].loc[index]})
if len(sys.argv)>2:
    with open(sys.argv[2]+'.json','w') as f:
        json.dump(dict_weights,f,indent=4)
    df_w_sample.to_excel(sys.argv[2]+'_weight_sample.xls',float_format='%0.10f')
    df_w_params.to_excel(sys.argv[2]+'_weight_param.xls',float_format='%0.10f')
    df_w_stitching.to_excel(sys.argv[2]+'_weight_stitch.xls',float_format='%0.10f')


print ('Fractions inclusive')
f_inc = {p:df.loc['Inclusive'][p]/df.loc['Inclusive']['total'] for p in params}
pprint (f_inc)
print ("Sum %0.5f"%(sum(f_inc.values())))
print ('Fractions exclusive')
rows = [s for s in samples if s != 'Inclusive']
f_exc = {p:(df[p].loc[rows]*df['xsec'].loc[rows]/df['total'].loc[rows]/dict_xsec['Inclusive']).sum() for p in params}
print ("Sum %0.5f"%(sum(f_exc.values())))
print ('Fractions stitching')
f_stitch = {p:xsec_params[p]/dict_xsec['Inclusive'] for p in params}
pprint (f_stitch)
print ("Sum %0.5f"%(sum(f_stitch.values())))
for p in params:
    print ('Param %10s fractions : Incl = %10.10f, Excl = %10.10f -> Mean = %10.10f VS Stitch = %10.10f'
            %(p,f_inc[p],f_exc[p],abs((f_inc[p]+f_exc[p])/2),f_stitch[p]))

df_fractions = pd.DataFrame([[a,b,c] for a,b,c in zip(list(f_inc.values()),list(f_exc.values()),list(f_stitch.values()))], columns=['Inclusive fraction','Exclusive fraction','Stitch fraction'],index=list(f_inc.keys()))
df_fractions['Mean(inc,exc)'] = (df_fractions['Inclusive fraction']+df_fractions['Exclusive fraction'])/2
if len(sys.argv)>2:
    df_fractions.to_excel(sys.argv[2]+'_fraction.xls',float_format='%0.10f')
