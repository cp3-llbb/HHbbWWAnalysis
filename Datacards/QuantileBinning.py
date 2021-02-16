import numpy as np
import scipy
import scipy.interpolate
from array import array
import ROOT

class Quantile:
    def __init__(self,h,q):
        if not isinstance(q,np.ndarray):
            q = np.array(q)
        if isinstance(h,ROOT.TH1F) or isinstance(h,ROOT.TH1D):
            e,w,s = self.getContent(h) 
        elif isinstance(h,list):
            e = None
            w = None
            s = None
            for hi in h:
                if not isinstance(hi,ROOT.TH1F) and not isinstance(hi,ROOT.TH1D):
                    raise RuntimeError('Not a ROOT histogram')
                ei,wi,si = self.getContent(hi)
                if e is None and w is None:
                    w = wi
                    e = ei
                    s = si**2
                else:
                    assert np.array_equal(e,ei)
                    w += wi
                    s += si**2
            s = np.sqrt(s)
        else:
            raise RuntimeError('Not a ROOT histogram')
        x = (e[:-1]+e[1:])/2
        if w[w>0].shape[0] >= q.shape[0]:
            nx = self.wquant(x[w>0],w[w>0],q)
            idx = np.digitize(nx, e) - 1
            self.ne = np.r_[e[0], e[idx], e[-1]]
            self.ne = np.unique(self.ne)
        elif w[w>0].shape[0] == 0:   
            self.ne = np.array([e[0],e[-1]])
        else:
            idx = np.digitize(x[w>0], e) - 1
            self.ne = np.r_[e[0], e[idx], e[-1]]
            self.ne = np.unique(self.ne)

    def __call__(self,h):
        e,w,s = self.getContent(h) 
        if np.isnan(w).any():
            print ('Warning : nan found in hist %s'%h.GetName())
            return None
        nw,ns = self.rebin(e,w,s,self.ne)
        return self.fillHistogram(self.ne,nw,ns,h.GetName()+'q')
        
    @staticmethod
    def getContent(h):
        e = [h.GetXaxis().GetBinUpEdge(0)]
        w = []
        s = []
        for i in range(1,h.GetNbinsX()+1):
            e.append(h.GetXaxis().GetBinUpEdge(i))
            w.append(h.GetBinContent(i))
            s.append(h.GetBinError(i))
        return np.array(e),np.array(w),np.array(s)

    @staticmethod
    def wquant(x, w, q):
        """
        x: bin centers
        w: bin heights (bin content)
        q: quantiles
        """
        assert x.shape == w.shape
        assert x.ndim == 1
        assert q.ndim == 1
        assert np.all((0 <= q) & (q <= 1))
        i = np.argsort(x)
        x = x[i]
        w = w[i]
        c = np.cumsum(w)
        inter = scipy.interpolate.interp1d(c, x, kind="nearest",fill_value="extrapolate")
        return inter(q[1:-1] * c[-1])

    @staticmethod
    def fillHistogram(e,w,s,name=""):
        h = ROOT.TH1F(name,name,e.shape[0]-1,array('d',e))
        for i in range(w.shape[0]):
            h.SetBinContent(i+1,w[i])
            h.SetBinError(i+1,s[i])
        return h

    @staticmethod
    def rebin(e,w,s,ne):
        x = [(e[i]+e[i])/2 for i in range(e.shape[0]-1)] 
        idx = np.digitize(x,ne) - 1
        nw = np.zeros(ne.shape[0]-1)
        ns = np.zeros(ne.shape[0]-1)
        for i in range(nw.shape[0]):
            nw[i] += w[idx == i].sum()
            ns[i] += (s[idx == i]**2).sum()
        ns = np.sqrt(ns)
        assert w.sum()==0. or abs(w.sum()-nw.sum())/w.sum() < 1e-9
        assert s.sum()==0. or abs((s**2).sum()-(ns**2).sum())/(s**2).sum() < 1e-9
        return nw,ns



