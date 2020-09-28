import os
import sys
import ROOT


path = sys.argv[1]

f_inc = ROOT.TFile(os.path.join(path,'DYJetsToLL_M-50.root'))
f_0j  = ROOT.TFile(os.path.join(path,'DYToLL_0J.root'))
f_1j  = ROOT.TFile(os.path.join(path,'DYToLL_1J.root'))
f_2j  = ROOT.TFile(os.path.join(path,'DYToLL_2J.root'))

h_inc = f_inc.Get('LHE_Njets')
h_0j  = f_0j.Get('LHE_Njets')
h_1j  = f_1j.Get('LHE_Njets')
h_2j  = f_2j.Get('LHE_Njets')

xsec_inc = 6077.2
xsec_0j  = 4843.6
xsec_1j  = 897.8
xsec_2j  = 335.8

N0_inc = h_inc.GetBinContent(1)
N1_inc = h_inc.GetBinContent(2)
N2_inc = h_inc.GetBinContent(3)
N3_inc = h_inc.GetBinContent(4)
N_inc = N0_inc + N1_inc + N2_inc + N3_inc

N0_0j = h_0j.GetBinContent(1)
N1_0j = h_0j.GetBinContent(2)
N2_0j = h_0j.GetBinContent(3)
N3_0j = h_0j.GetBinContent(4)
N_0j = N0_0j + N1_0j + N2_0j + N3_0j

N0_1j = h_1j.GetBinContent(1)
N1_1j = h_1j.GetBinContent(2)
N2_1j = h_1j.GetBinContent(3)
N3_1j = h_1j.GetBinContent(4)
N_1j = N0_1j + N1_1j + N2_1j + N3_1j

N0_2j = h_2j.GetBinContent(1)
N1_2j = h_2j.GetBinContent(2)
N2_2j = h_2j.GetBinContent(3)
N3_2j = h_2j.GetBinContent(4)
N_2j = N0_2j + N1_2j + N2_2j + N3_2j


print ('Inclusive sample :')
print ('...','N_inc^0     = %0.2E'%N0_inc)
print ('...','N_inc^1     = %0.2E'%N1_inc)
print ('...','N_inc^2     = %0.2E'%N2_inc)
print ('...','N_inc^3     = %0.2E'%N3_inc)
print ('...','Total N_inc = %0.2E'%N_inc)

print ('Exclusive 0J sample :')
print ('...','N_0J^0     = %0.2E'%N0_0j)
print ('...','N_0J^1     = %0.2E'%N1_0j)
print ('...','N_0J^2     = %0.2E'%N2_0j)
print ('...','N_0J^3     = %0.2E'%N3_0j)
print ('...','Total N_0J = %0.2E'%N_0j)

print ('Exclusive 1J sample :')
print ('...','N_1J^0     = %0.2E'%N0_1j)
print ('...','N_1J^1     = %0.2E'%N1_1j)
print ('...','N_1J^2     = %0.2E'%N2_1j)
print ('...','N_1J^3     = %0.2E'%N3_1j)
print ('...','Total N_1J = %0.2E'%N_1j)

print ('Exclusive 2J sample :')
print ('...','N_2J^0     = %0.2E'%N0_2j)
print ('...','N_2J^1     = %0.2E'%N1_2j)
print ('...','N_2J^2     = %0.2E'%N2_2j)
print ('...','N_2J^3     = %0.2E'%N3_2j)
print ('...','Total N_2J = %0.2E'%N_2j)

N0 = N0_inc+N0_0j+N0_1j+N0_2j
N1 = N1_inc+N1_0j+N1_1j+N1_2j
N2 = N2_inc+N2_0j+N2_1j+N2_2j
N3 = N3_inc+N3_0j+N3_1j+N3_2j

print ("Number per jet multiplicity")
print ("...","N^0 = %0.2E"%N0)
print ("...","N^1 = %0.2E"%N1)
print ("...","N^2 = %0.2E"%N2)
print ("...","N^3 = %0.2E"%N3)

xsec0 = 0.5*(xsec_0j*N0_0j/N_0j + xsec_1j*N0_1j/N_1j + xsec_2j*N0_2j/N_2j + xsec_inc*N0_inc/N_inc)
xsec1 = 0.5*(xsec_0j*N1_0j/N_0j + xsec_1j*N1_1j/N_1j + xsec_2j*N1_2j/N_2j + xsec_inc*N1_inc/N_inc)
xsec2 = 0.5*(xsec_0j*N2_0j/N_0j + xsec_1j*N2_1j/N_1j + xsec_2j*N2_2j/N_2j + xsec_inc*N2_inc/N_inc)
xsec3 = 0.5*(xsec_0j*N3_0j/N_0j + xsec_1j*N3_1j/N_1j + xsec_2j*N3_2j/N_2j + xsec_inc*N3_inc/N_inc)

print ('Cross-section per jet multiplicity')
print ("...",'Xsec^0 :',xsec0)
print ("...",'Xsec^1 :',xsec1)
print ("...",'Xsec^2 :',xsec2)
print ("...",'Xsec^3 :',xsec3)
print ('Total cross-section (sum over jet mult) :',xsec0+xsec1+xsec2+xsec3)
print ('Inclusive cross section :',xsec_inc)

w_inc = xsec_inc/N_inc
w_0j  = xsec_0j/N_0j
w_1j  = xsec_1j/N_1j
w_2j  = xsec_2j/N_2j

print ('Lumiscale weight per sample')
print ("...","Inclusive    = %0.2E"%w_inc)
print ("...","Exclusive 0J = %0.2E"%w_0j)
print ("...","Exclusive 1J = %0.2E"%w_1j)
print ("...","Exclusive 2J = %0.2E"%w_2j)

w0 = xsec0/N0
w1 = xsec1/N1
w2 = xsec2/N2
w3 = xsec3/N3

print ('Lumiscale weights per jet multiplicity')
print ("...","w0 = %0.2E"%w0)
print ("...","w1 = %0.2E"%w1)
print ("...","w2 = %0.2E"%w2)
print ("...","w3 = %0.2E"%w3)

w0_inc = w0/w_inc
w1_inc = w1/w_inc
w2_inc = w2/w_inc
w3_inc = w3/w_inc

w0_0j = w0/w_0j
w1_0j = w1/w_0j
w2_0j = w2/w_0j
w3_0j = w3/w_0j

w0_1j = w0/w_1j
w1_1j = w1/w_1j
w2_1j = w2/w_1j
w3_1j = w3/w_1j

w0_2j = w0/w_2j
w1_2j = w1/w_2j
w2_2j = w2/w_2j
w3_2j = w3/w_2j

if len(sys.argv)>2:
    import json
    dict_weights = {'DYJetsToLL_M-50':{0:w0_inc,1:w1_inc,2:w2_inc,3:w3_inc},
                    'DYToLL_0J':{0:w0_0j,1:w1_0j,2:w2_0j,3:w3_0j},
                    'DYToLL_1J':{0:w0_1j,1:w1_1j,2:w2_1j,3:w3_1j},
                    'DYToLL_2J':{0:w0_2j,1:w1_2j,2:w2_2j,3:w3_2j}}
    with open(sys.argv[2],'w') as f:
        json.dump(dict_weights,f,indent=4)


print ('Stitching weights')
print ('...','weight on inclusive sample      in bin 0 (w^0/w_inc) = %0.3f'%w0_inc)
print ('...','weight on inclusive sample      in bin 1 (w^1/w_inc) = %0.3f'%w1_inc)
print ('...','weight on inclusive sample      in bin 2 (w^2/w_inc) = %0.3f'%w2_inc)
print ('...','weight on inclusive sample      in bin 3 (w^3/w_inc) = %0.3f'%w3_inc)
print ('...','weight on exclusive sample (0J) in bin 0 (w^0/w_0)   = %0.3f'%w0_0j)
print ('...','weight on exclusive sample (0J) in bin 1 (w^1/w_0)   = %0.3f'%w1_0j)
print ('...','weight on exclusive sample (0J) in bin 2 (w^2/w_0)   = %0.3f'%w2_0j)
print ('...','weight on exclusive sample (0J) in bin 3 (w^3/w_0)   = %0.3f'%w3_0j)
print ('...','weight on exclusive sample (1J) in bin 0 (w^0/w_1)   = %0.3f'%w0_1j)
print ('...','weight on exclusive sample (1J) in bin 1 (w^1/w_1)   = %0.3f'%w1_1j)
print ('...','weight on exclusive sample (1J) in bin 2 (w^2/w_1)   = %0.3f'%w2_1j)
print ('...','weight on exclusive sample (1J) in bin 3 (w^3/w_1)   = %0.3f'%w3_1j)
print ('...','weight on exclusive sample (2J) in bin 0 (w^0/w_2)   = %0.3f'%w0_2j)
print ('...','weight on exclusive sample (2J) in bin 1 (w^1/w_2)   = %0.3f'%w1_2j)
print ('...','weight on exclusive sample (2J) in bin 2 (w^2/w_2)   = %0.3f'%w2_2j)
print ('...','weight on exclusive sample (2J) in bin 3 (w^3/w_2)   = %0.3f'%w3_2j)

f0_inc = N0_inc/N_inc
f1_inc = N1_inc/N_inc
f2_inc = N2_inc/N_inc
f3_inc = N3_inc/N_inc

f0_exc = N0_0j/N_0j * xsec_0j/xsec_inc + N0_1j/N_1j * xsec_1j/xsec_inc + N0_2j/N_2j * xsec_2j/xsec_inc
f1_exc = N1_0j/N_0j * xsec_0j/xsec_inc + N1_1j/N_1j * xsec_1j/xsec_inc + N1_2j/N_2j * xsec_2j/xsec_inc
f2_exc = N2_0j/N_0j * xsec_0j/xsec_inc + N2_1j/N_1j * xsec_1j/xsec_inc + N2_2j/N_2j * xsec_2j/xsec_inc
f3_exc = N3_0j/N_0j * xsec_0j/xsec_inc + N3_1j/N_1j * xsec_1j/xsec_inc + N3_2j/N_2j * xsec_2j/xsec_inc

f0_stitch = xsec0/xsec_inc
f1_stitch = xsec1/xsec_inc
f2_stitch = xsec2/xsec_inc
f3_stitch = xsec3/xsec_inc

print ('Fractions')
print ('... 0J : inclusive (%0.3f), exclusive (%0.3f) -> Mean (%0.3f) VS stiched (%0.3f)'%(f0_inc,f0_exc,(f0_inc+f0_exc)/2,f0_stitch))
print ('... 1J : inclusive (%0.3f), exclusive (%0.3f) -> Mean (%0.3f) VS stiched (%0.3f)'%(f1_inc,f1_exc,(f1_inc+f1_exc)/2,f1_stitch))
print ('... 2J : inclusive (%0.3f), exclusive (%0.3f) -> Mean (%0.3f) VS stiched (%0.3f)'%(f2_inc,f2_exc,(f2_inc+f2_exc)/2,f2_stitch))
print ('... 3J : inclusive (%0.3f), exclusive (%0.3f) -> Mean (%0.3f) VS stiched (%0.3f)'%(f3_inc,f3_exc,(f3_inc+f3_exc)/2,f3_stitch))
