import numpy as np
import scipy
import scipy.interpolate
from array import array
import ROOT



class Rebin:
    """
        Base rebin method
        includes common methods to rebin, extract and fill histograms
    """
    @staticmethod
    def getContent(h):
        """
            From TH1 extract : 
                e : edges including lower and upper 
                w : content (GetBinContent)
                s : errors  (GetBinError)
            return [e,w,s]
        """
        e = [h.GetXaxis().GetBinUpEdge(0)]
        w = []
        s = []
        for i in range(1,h.GetNbinsX()+1):
            e.append(h.GetXaxis().GetBinUpEdge(i))
            w.append(h.GetBinContent(i))
            s.append(h.GetBinError(i))
        return np.array(e).round(9),np.array(w),np.array(s)

    @staticmethod
    def rebin(e,w,s,ne):
        """
            Given content of histogram 
                e  : edges of initial histogram
                w  : content of initial histogram 
                s  : errors of initial histogram 
                ne : new bin edges
            returns : [nw,ns]
                nw : new bin content
                ns : new bin errors
        """
        if e[0] != ne[0]:
            raise RuntimeError("X axis first edge not matching : {} != {}".format(e[0],ne[0]))
        if e[-1] != ne[-1]:
            raise RuntimeError("X axis last edge not matching : {} != {}".format(e[-1],ne[-1]))
        if np.any(ne[:-1] > ne[1:]):
            raise RuntimeError("X axis new edges not increasing : ["+",".join([str(ne[i]) for i in range(ne[0].shape[0])])+"]")
        if not np.isin(ne,e).all():
            wrong_edges = ne[~np.isin(ne,e)]
            raise RuntimeError("New X axis edges not contained in initial ones : ["+",".join([str(wrong_edges[i]) for i in range(wrong_edges.shape[0])])+"]")

        x = (e[1:]+e[:-1])/2
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

    @staticmethod
    def fillHistogram(e,w,s,name=""):
        """
            Inputs : 
            e : bin edges
            w : bin content
            s : bin errors
            name : name to be used in the TH1 instantiation (default = '')
            return : TH1
        """
        h = ROOT.TH1F(name,name,e.shape[0]-1,array('d',e))
        for i in range(w.shape[0]):
            h.SetBinContent(i+1,w[i])
            h.SetBinError(i+1,s[i])
        return h

    def __call__(self,h):
        """
            input : initial histogram TH1
            return : rebinned TH1
        """
        e,w,s = self.getContent(h) 
        if np.isnan(w).any():
            print ('Warning : nan found in hist %s'%h.GetName())
            return None
        if not hasattr(self,'ne'):
            raise RuntimeError('New bin edges have not been computed, is the rebin_method() not implemented')
        nw,ns = self.rebin(e,w,s,self.ne)
        return self.fillHistogram(self.ne,nw,ns,h.GetName()+'q')
        



class Quantile(Rebin):
    """
        Applies quantile binning
        -> provide quantile boundaries and bin histogram accordingly
        list of quantile values need to be optimized
    """
    def __init__(self,h,q):
        """
            h : either TH1 or list of TH1
            q : quantiles list [0 , ... 1]
        """
        if not isinstance(q,np.ndarray):
            q = np.array(q)
        if q[0] != 0. or q[-1] != 1.:
            raise RuntimeError("Invalid quantiles boundaries ["+",".join([str(q[i]) for i in range(q.shape[0])])+"]")
        if np.any(q[:-1] > q[1:]):
            raise RuntimeError("Quantile edges not increasing ["+",".join([str(q[i]) for i in range(q.shape[0])])+"]")
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
            nx = self.rebin_method(x[w>0],w[w>0],q)
            idx = np.digitize(nx, e) - 1
            self.ne = np.r_[e[0], e[idx], e[-1]]
            self.ne = np.unique(self.ne)
        elif w[w>0].shape[0] == 0:   
            self.ne = np.array([e[0],e[-1]])
        else:
            idx = np.digitize(x[w>0], e) - 1
            self.ne = np.r_[e[0], e[idx], e[-1]]
            self.ne = np.unique(self.ne)

    @staticmethod
    def rebin_method(x, w, q):
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


class Threshold(Rebin):
    """
        Applied threshold binning 
        -> content is scanned from the right of the histogram 
           until enough stats are accumulated so that 
                 bincontent - unc > threshold
           then move to next threshold
        list of threshold values need to be optimized
    """
    def __init__(self, h, thresh, use_var=False, extra=None, rsut=None):
        """
            thresh : thresholds to remain above
            h : either TH1 or list of TH1
            use_var : if True will take the variance from the bin errors, else zeroes 
            extra : list of hists to remain above 0
        """
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
        if extra is not None:
            extra = [self.getContent(extra_i)[1] for extra_i in extra]
            for idx in range(len(extra)):
                if extra[idx].shape[0] != w.shape[0]:
                    raise RuntimeError("Extra contribution index {} does not have the same number of bins".format(idx))
            extra = np.c_[tuple(extra)]
        if not isinstance(thresh,np.ndarray):
            thresh = np.array(thresh)
        if w.shape[0] < thresh.shape[0]:
            raise RuntimeError("Fewer bins than thresholds")
        if not use_var:
            s = np.zeros_like(w)
        idx = self.rebin_method(thresh=thresh,val=w,var=np.power(s,2),extra=extra,rsut=rsut)

        # fuse left-most two bins if they are rising in content
        if len(idx) > 0 and w[0 : idx[0]].sum() < w[idx[0] : (idx[1] if len(idx) > 1 else None)].sum():
            idx = idx[1:]     
        self.ne = np.r_[e[0], e[idx] , e[-1]]

    @staticmethod
    def rebin_method(thresh,val,var,extra,rsut=np.inf):
        """
            thresh : array of threshold values from right to left
            val    : array of bin height
            var    : variance of each bin
            extra  : additional contributions to be kept above threshold
            rsut   : relative stat. unc. threshold
        """
        assert thresh.shape[0] < val.shape[0]     
        assert np.all(thresh >= 0)     
        assert rsut >= 0.0 

        sum_val = 0
        sum_var = 0
        sum_extra = np.zeros_like(extra[0])

        val = val[::-1]
        var = var[::-1]
        extra = extra[::-1]

        la = len(val)
        lt = len(thresh)
        idx = np.zeros(lt, dtype=np.intp)
        tidx = -1
        for i in range(la):
            if tidx < -lt:
                break
            sum_val += val[i]
            sum_var += var[i]
            sum_extra += extra[i]
            unc = np.sqrt(sum_var)

            if (sum_val - unc) >= thresh[tidx] and np.all(sum_extra) and np.abs((unc / sum_val)) <= rsut:
                idx[tidx] = la - 1 - i
                sum_val = 0.0
                sum_var = 0.0
                sum_extra[:] = 0.0
                tidx -= 1
        return idx[1 + tidx :]     

class Boundary(Rebin):
    """
        Applied boundary binning 
        Rebin with the given boundaries for the bin edges
    """
    def __init__(self, boundaries):
        """
            boundaries : list of bin edges
        """
        self.ne = np.array(boundaries)

class Rebin2D(Rebin):
    @staticmethod
    def getContent(h):
        """
            From TH2 extract : 
                e : edges including lower and upper 
                w : content (GetBinContent)
                s : errors  (GetBinError)
            return [e,w,s]
        """
        xAxis = h.GetXaxis()
        yAxis = h.GetYaxis()
        Nx = xAxis.GetNbins() 
        Ny = yAxis.GetNbins()
        w = np.zeros((Nx,Ny))
        s = np.zeros((Nx,Ny))
        for ix in range(0,Nx):
            for iy in range(0,Ny):
                w[ix,iy] = h.GetBinContent(ix+1,iy+1)
                s[ix,iy] = h.GetBinError(ix+1,iy+1)
        e = [np.array([xAxis.GetBinLowEdge(i) for i in range(1,Nx+2)]).round(9),
             np.array([yAxis.GetBinLowEdge(i) for i in range(1,Ny+2)]).round(9)]
        return e,w,s

    @staticmethod
    def rebin(e,w,s,ne):
        """
            Given content of histogram 
                e  : edges of initial histogram
                w  : content of initial histogram 
                s  : errors of initial histogram 
                ne : new bin edges
            returns : [nw,ns]
                nw : new bin content
                ns : new bin errors
        """
        if e[0][0] != ne[0][0]:
            raise RuntimeError("X axis first edge not matching : {} != {}".format(e[0][0],ne[0][0]))
        if e[1][0] != ne[1][0]:
            raise RuntimeError("Y axis first edge not matching : {} != {}".format(e[1][0],ne[1][0]))
        if e[0][-1] != ne[0][-1]:
            raise RuntimeError("X axis last edge not matching : {} != {}".format(e[0][-1],ne[0][-1]))
        if e[1][-1] != ne[1][-1]:
            raise RuntimeError("Y axis last edge not matching : {} != {}".format(e[1][-1],ne[1][-1]))
        if np.any(ne[0][:-1] > ne[0][1:]):
            raise RuntimeError("X axis new edges not increasing : ["+",".join([str(ne[0][i]) for i in range(ne[0].shape[0])])+"]")
        if np.any(ne[1][:-1] > ne[1][1:]):
            raise RuntimeError("Y axis new edges not increasing : ["+",".join([str(ne[1][i]) for i in range(ne[1].shape[0])])+"]")
        if not np.isin(ne[0],e[0]).all():
            wrong_edges = ne[0][~np.isin(ne[0],e[0])]
            raise RuntimeError("New X axis edges not contained in initial ones : ["+",".join([str(wrong_edges[i]) for i in range(wrong_edges.shape[0])])+"]")
        if not np.isin(ne[1],e[1]).all():
            wrong_edges = ne[1][~np.isin(ne[1],e[1])]
            raise RuntimeError("New Y axis edges not contained in initial ones : ["+",".join([str(wrong_edges[i]) for i in range(wrong_edges.shape[0])])+"]")

        x = (e[0][1:]+e[0][:-1])/2
        y = (e[1][1:]+e[1][:-1])/2
        idx = np.digitize(x,ne[0]) - 1
        idy = np.digitize(y,ne[1]) - 1
        nw = np.zeros((ne[0].shape[0]-1,ne[1].shape[0]-1))
        ns = np.zeros((ne[0].shape[0]-1,ne[1].shape[0]-1))
        for ix in range(0,nw.shape[0]):
            for iy in range(0,nw.shape[1]):
                nw[ix,iy] += w[np.ix_(idx==ix,idy==iy)].sum()
                ns[ix,iy] += (s[np.ix_(idx==ix,idy==iy)]**2).sum()
        ns = np.sqrt(ns)
        assert w.sum()==0. or abs(w.sum()-nw.sum())/w.sum() < 1e-9
        assert s.sum()==0. or abs((s**2).sum()-(ns**2).sum())/(s**2).sum() < 1e-9
        return nw,ns

    @staticmethod
    def fillHistogram(e,w,s,name=""):
        """
            Inputs : 
            e : bin edges
            w : bin content
            s : bin errors
            name : name to be used in the TH2 instantiation (default = '')
            return : TH2
        """
        assert len(e) == 2
        h = ROOT.TH2F(name,name,e[0].shape[0]-1,array('d',e[0]),e[1].shape[0]-1,array('d',e[1]))
        for ix in range(w.shape[0]):
            for iy in range(0,w.shape[1]):
                h.SetBinContent(ix+1,iy+1,w[ix,iy])
                h.SetBinError(ix+1,iy+1,s[ix,iy])
        return h

class Boundary2D(Rebin2D):
    """
        Applied boundary binning 
        Rebin with the given boundaries for the bin edges
    """
    def __init__(self,bx,by):
        """
            bx : list of edges for the x axis
            by : list of edges for the y axis
        """
        assert isinstance(bx,list)
        assert isinstance(by,list)
        self.ne = [np.array(bx),np.array(by)]

class Quantile2D(Rebin2D):
    """
        Applies quantile binning
        -> provide quantile boundaries and bin histogram accordingly
        list of quantile values need to be optimized
    """
    def __init__(self,h,qx,qy):
        """
            h : either TH1 or list of TH1
            qx : quantiles list [0 , ... 1]
            qy : quantiles list [0 , ... 1]
        """
        if isinstance(h,ROOT.TH2F) or isinstance(h,ROOT.TH2D):
            h2D = h 
        elif isinstance(h,list):
            h2D = None
            for hi in h:
                if not isinstance(hi,ROOT.TH2F) and not isinstance(hi,ROOT.TH2D):
                    raise RuntimeError('Not a ROOT histogram')
                if h2D is None:
                    h2D = hi
                else:
                    h2D.Add(hi)
        else:
            raise RuntimeError('Not a ROOT histogram')

        hx = h2D.ProjectionX()
        hy = h2D.ProjectionY()

        qObjx = Quantile(hx,qx)
        qObjy = Quantile(hy,qy)

        self.ne = [qObjx.ne,qObjy.ne]



class Threshold2D(Rebin2D):
    """
        Applied threshold binning 
        -> content is scanned from the right of the histogram 
           until enough stats are accumulated so that 
                 bincontent - unc > threshold
           then move to next threshold
        list of threshold values need to be optimized
    """
    def __init__(self, h, threshx, threshy, use_var=False, extra=None, rsut=None):
        """
            thresh : thresholds to remain above
            h : either TH1 or list of TH1
            use_var : if True will take the variance from the bin errors, else zeroes 
            extra : list of hists to remain above 0
        """
        if isinstance(h,ROOT.TH2F) or isinstance(h,ROOT.TH2D):
            h2D = h 
        elif isinstance(h,list):
            h2D = None
            for hi in h:
                if not isinstance(hi,ROOT.TH2F) and not isinstance(hi,ROOT.TH2D):
                    raise RuntimeError('Not a ROOT histogram')
                if h2D is None:
                    h2D = hi
                else:
                    h2D.Add(hi)
        else:
            raise RuntimeError('Not a ROOT histogram')

        if not isinstance(threshx,np.ndarray):
            threshx = np.array(threshx)
        if not isinstance(threshy,np.ndarray):
            threshy = np.array(threshy)

        hx = h2D.ProjectionX()
        hy = h2D.ProjectionY()

        qObjx = Threshold(hx,threshx,use_var,extra,rsut)
        qObjy = Threshold(hy,threshyuse_var,extra,rsut)

        self.ne = [qObjx.ne,qObjy.ne]


class Linearize2D(Rebin2D):
    def __init__(self,major='y'):
        if major != 'x' and major != 'y':
            raise RuntimeError(f'Major {major} not understood')
        self.major = major

    def __call__(self,h):
        """
            input : initial histogram TH2
            return : linearized TH1
        """
        e,w,s = self.getContent(h) 
        if self.major == 'x':
            flatten_arg = "C"
            self.eminor = e[1]
            self.emajor = e[0]
        else:
            flatten_arg = "F"
            self.eminor = e[0]
            self.emajor = e[1]
        flatten_arg = "C" if self.major == 'x' else "F"
        nw = w.flatten(flatten_arg)
        ns = s.flatten(flatten_arg)
        ne = np.unique(np.array([self.eminor + i * self.eminor[-1] for i in range(len(self.emajor)-1)]).flatten('C'))
        self.majorsep = [i*self.eminor[-1] for i in range(1,len(self.emajor))]
        nh = ROOT.TH1F("","",ne.shape[0]-1,array('d',ne))
        for i in range(nw.shape[0]):
            nh.SetBinContent(i+1,nw[i])
            nh.SetBinError(i+1,ns[i])
            nh.GetXaxis().SetBinLabel(i+1,'#{}'.format(i%(self.eminor.shape[0]-1)+1))
        nh.GetXaxis().SetNdivisions(self.emajor.shape[0]+self.eminor.shape[0]*100)
        return nh

    def getPlotData(self):
        labels = []
        lines  = []
        ylabel = 'y'
        for i in range(len(self.majorsep)):
            label = f"{self.emajor[i]} < {ylabel} < {self.emajor[i+1]}"
            poslabel = [float(i+0.2)/len(self.majorsep),0.02]
            labels.append({'text':label,'position':poslabel})
            lines.append(float(self.majorsep[i]))

        return {'labels':labels,'vertical-lines':lines}


