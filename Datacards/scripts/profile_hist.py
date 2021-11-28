import os
import sys
import ROOT
from time import perf_counter
from array import array
import numpy as np

from hist_interface import PythonInterface, CppInterface

from ROOT import gErrorIgnoreLevel
ROOT.TH1.AddDirectory(False)

from IPython import embed

TIMER_EVAL = 10
N_BINS = 400

h1 = ROOT.TH1F('h1','h1',N_BINS,0,100)
h2 = ROOT.TH2F('h2','h2',N_BINS,0.,100.,N_BINS,0.,200.)

rg = ROOT.TRandom3()
for _ in range(100000):
    h1.Fill(rg.Gaus(50,10))
    h2.Fill(rg.Gaus(50,10),rg.Gaus(50,10))

def timer(func):
    def wrap_func(*args, **kwargs):
        times = []
        for i in range(TIMER_EVAL+1):
            t1 = perf_counter()
            result = func(*args, **kwargs)
            t2 = perf_counter()
            if i > 0:  # Burn in time
                times.append(t2-t1)
        times = np.array(times)
        print(f'Function {func.__qualname__!r:40s} executed in {times.mean():.6f}s +/- {times.std():.6f}s ({TIMER_EVAL} calls)')
        return result
    return wrap_func
                                                  

def assertTH1(h1,h2,atol=0.,rtol=1e-7):
    if h1.GetNbinsX() != h2.GetNbinsX():
        raise RuntimeError(f'N h1 [= {h1.GetNbinsX()}] != N h2 [= {h2.GetNbinsX()}]')
    for i in range(1,h1.GetNbinsX()+1):
        a = h1.GetBinContent(i)
        b = h2.GetBinContent(i)
        if abs(a-b) > atol + rtol * abs(b):
            raise RuntimeError(f'Position {i} : {a} not close to {b}')

def assertTH2(h1,h2,atol=0.,rtol=1e-7):
    if h1.GetNbinsX() != h2.GetNbinsX():
        raise RuntimeError(f'Nx h1 [= {h1.GetNbinsX()}] != Nx h2 [= {h2.GetNbinsX()}]')
    if h1.GetNbinsY() != h2.GetNbinsY():
        raise RuntimeError(f'Ny h1 [= {h1.GetNbinsY()}] != Ny h2 [= {h2.GetNbinsY()}]')
    for x in range(1,h1.GetNbinsX()+1):
        for y in range(1,h1.GetNbinsY()+1):
            a = h1.GetBinContent(x,y)
            b = h2.GetBinContent(x,y)
            if abs(a-b) > atol + rtol * abs(b):
                raise RuntimeError(f'Position ({x},{y}) : {a} not close to {b}')

## Extracting #
e1p,w1p,s1p = timer(PythonInterface.getContent1D)(h1)
e1c,w1c,s1c = timer(CppInterface.getContent1D)(h1)

e2p,w2p,s2p = timer(PythonInterface.getContent2D)(h2)
e2c,w2c,s2c = timer(CppInterface.getContent2D)(h2)
# Testing #
np.testing.assert_allclose(e1p,e1c)
np.testing.assert_allclose(w1p,w1c)
np.testing.assert_allclose(s1p,s1c)

np.testing.assert_allclose(e2p[0],e2c[0])
np.testing.assert_allclose(e2p[1],e2c[1])
np.testing.assert_allclose(w2p,w2c)
np.testing.assert_allclose(s2p,s2c)

## Filling #
h1p = timer(PythonInterface.fillHistogram1D)(e1p,w1p,s1p,'h1p')
h1c = timer(CppInterface.fillHistogram1D)(e1c,w1c,s1c,'h1c')

assertTH1(h1p,h1c)

h2p = timer(PythonInterface.fillHistogram2D)(e2p,w2p,s2p,'h2p')
h2c = timer(CppInterface.fillHistogram2D)(e2c,w2c,s2c,'h2c')

assertTH2(h2p,h2c)

