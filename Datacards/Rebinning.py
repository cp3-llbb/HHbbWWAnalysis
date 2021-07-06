import sys
import numpy as np
import scipy
import logging
import scipy.interpolate
from array import array
import ROOT

def compareTwoAxes(x1,x2):
    if x1[0] != x2[0]:
        raise RuntimeError("Axis first edge not matching : {} != {}".format(x1[0],x2[0]))
    if x1[-1] != x2[-1]:
        raise RuntimeError("Axis last edge not matching : {} != {}".format(x1[-1],x2[-1]))
    if np.any(x1[:-1] > x1[1:]):
        raise RuntimeError("Axis old edges not increasing : ["+",".join([str(x1[i]) for i in range(x1.shape[0])])+"]")
    if np.any(x2[:-1] > x2[1:]):
        raise RuntimeError("Axis new edges not increasing : ["+",".join([str(x2[i]) for i in range(x2.shape[0])])+"]")
    if not np.isin(x2,x1).all():
        wrong_edges = x2[~np.isin(x2,x1)]
        raise RuntimeError("New X axis edges not contained in initial ones : ["+",".join([str(wrong_edges[i]) for i in range(wrong_edges.shape[0])])+"]")

def processHistogramToROOT(h,htype):
    if isinstance(h,htype):
        hist = h 
    elif isinstance(h,list):
        hist = None
        for hi in h:
            if not isinstance(hi,htype):
                raise RuntimeError('Not a ROOT histogram')
            if hist is None:
                hist = hi
            else:
                hist.Add(hi)
    else:
        raise RuntimeError('Not a ROOT histogram nor a list of ROOT histograms : '+str(type(h)))
    return hist

def processHistogramToNumpy(h,htype,getter):
        if isinstance(h,htype):
            e,w,s = getter(h) 
        elif isinstance(h,list):
            e = None
            w = None
            s = None
            for hi in h:
                if not isinstance(hi,htype):
                    raise RuntimeError('Not a ROOT histogram')
                ei,wi,si = getter(hi)
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
            raise RuntimeError('Not a ROOT histogram nor a list of ROOT histograms : '+str(type(h)))
        return e,w,s



class Rebin:
    """
        Base rebin method
        includes common methods to rebin, extract and fill histograms
    """
    @staticmethod
    def getContent1D(h):
        """
            From TH1 extract : 
                e : edges including lower and upper 
                w : content (GetBinContent)
                s : errors  (GetBinError)
            return [e,w,s]
        """
        assert isinstance(h,ROOT.TH1)
        e = [h.GetXaxis().GetBinUpEdge(0)]
        w = []
        s = []
        for i in range(1,h.GetNbinsX()+1):
            e.append(h.GetXaxis().GetBinUpEdge(i))
            w.append(h.GetBinContent(i))
            s.append(h.GetBinError(i))
        return np.array(e).round(9),np.array(w),np.array(s)

    @staticmethod
    def rebin1D(e,w,s,ne):
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
        compareTwoAxes(e,ne)

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
    def fillHistogram1D(e,w,s,name=""):
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
        e,w,s = self.getContent1D(h) 
        if np.isnan(w).any():
            logging.warning('Warning : nan found in hist %s'%h.GetName())
            return None
        if not hasattr(self,'ne'):
            raise RuntimeError('New bin edges have not been computed, is the rebin_method() not implemented')
        nw,ns = self.rebin1D(e,w,s,self.ne)
        return self.fillHistogram1D(self.ne,nw,ns,h.GetName()+'rebin')
        



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
        e,w,s = processHistogramToNumpy(h,ROOT.TH1,self.getContent1D)
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
    def __init__(self, h, thresh, extra=None, nbins=None, rsut=None):
        """
            thresh : thresholds to remain above
            h : either TH1 or list of TH1
            extra : list of hists to remain above 0
        """
        use_var = True 
        e,w,s = processHistogramToNumpy(h,ROOT.TH1,self.getContent1D)
        if extra is not None:
            extra = [processHistogramToNumpy(extra_i,ROOT.TH1,self.getContent1D)[1] for extra_i in extra]
            for i in range(len(extra)):
                if extra[i].shape[0] != w.shape[0]:
                    raise RuntimeError("Extra contribution index {} does not have the same number of bins".format(i))
            extra = np.c_[tuple(extra)]
        if not isinstance(thresh,np.ndarray):
            thresh = np.array(thresh)
        if w.shape[0] < thresh.shape[0]:
            raise RuntimeError("Fewer bins than thresholds")
        if not use_var:
            s = np.zeros_like(w)
#        if factor is not None:
#            thresh *= factor * w.sum() / thresh.max() 
#        idx = self.rebin_method(thresh=thresh,val=w,var=np.power(s,2),extra=extra,rsut=rsut)
#        # fuse left-most two bins if they are rising in content
#        if len(idx) > 0 and w[0 : idx[0]].sum() < w[idx[0] : (idx[1] if len(idx) > 1 else None)].sum():
#            idx = idx[1:]     
#        self.ne = np.unique(np.r_[e[0], e[idx] , e[-1]])
        if nbins is None:
            factors = [1.]
        else:
            nbins = min(nbins,thresh.shape[0])
            factors = list(np.linspace(1,0,101).round(4))
        for factor in factors:
            thresh_test = thresh * factor * w.sum() / thresh.max() 
            idx = self.rebin_method(thresh=thresh_test,val=w,var=np.power(s,2),extra=extra,rsut=rsut)
            if len(idx) > 0 and w[0 : idx[0]].sum() < w[idx[0] : (idx[1] if len(idx) > 1 else None)].sum(): # merge two first bins in case rising in content
                idx = idx[1:]     
            self.ne = np.unique(np.r_[e[0], e[idx] , e[-1]])
#            nw,ns = self.rebin1D(e,w,s,self.ne)
#            import IPython
#            IPython.embed()
#            inc_bins = ((nw[1:] > nw[:-1])*1).sum()
#            print (self.ne.shape[0],inc_bins,0.9*self.ne.shape[0],inc_bins > 0.9*self.ne.shape[0])
            if nbins is not None and self.ne.shape[0] > nbins :#and inc_bins > 0.9*self.ne.shape[0]:
                break

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

            #if (sum_val - unc) >= thresh[tidx] and np.all(sum_extra>3.) and np.abs((unc / sum_val)) <= rsut:
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
        h not used (kept for compatibility)
    """
    def __init__(self, h=None, boundaries=None):
        """
            boundaries : list of bin edges
        """
        if boundaries is None:
            raise RuntimeError("Boundary method requires to set the boundaries")
        self.ne = np.array(boundaries)

class Rebin2D(Rebin):
    @staticmethod
    def getContent2D(h):
        """
            From TH2 extract : 
                e : edges including lower and upper 
                w : content (GetBinContent)
                s : errors  (GetBinError)
            return [e,w,s]
        """
        assert isinstance(h,ROOT.TH2)
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

    def __call__(self,h):
        """
            input : initial histogram TH2
            return : rebinned TH2
        """
        e,w,s = self.getContent2D(h) 
        if np.isnan(w).any():
            logging.warning('Warning : nan found in hist %s'%h.GetName())
            return None
        if not hasattr(self,'ne'):
            raise RuntimeError('New bin edges have not been computed, is the rebin_method() not implemented')
        nw,ns = self.rebin2D(e,w,s,self.ne)
        return self.fillHistogram2D(self.ne,nw,ns,h.GetName()+'rebin')
        


    @staticmethod
    def rebin2D(e,w,s,ne):
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
        compareTwoAxes(e[0],ne[0])
        compareTwoAxes(e[1],ne[1])

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
    def fillHistogram2D(e,w,s,name=""):
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
        h2D = processHistogramToROOT(h,ROOT.TH2)
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
    def __init__(self, h, threshx, threshy, extra=None, rsut=None):
        """
            thresh : thresholds to remain above
            h : either TH1 or list of TH1
            use_var : if True will take the variance from the bin errors, else zeroes 
            extra : list of hists to remain above 0
        """
        h2D = processHistogramToROOT(h,ROOT.TH2)
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

        self.nemajor = None
        self.neminor = None

    def __call__(self,h):
        """
            input : initial histogram TH2
            return : linearized TH1
        """
        e,w,s = self.getContent2D(h) 
        if np.isnan(w).any():
            logging.warning('Warning : nan found in hist %s'%h.GetName())
            return None
        nw,ns,eminor,emajor = self.split(e,w,s)
        if self.neminor is not None:
            if len(nw) != len(self.neminor):
                raise RuntimeError(f'Number of major bins content {len(nw)} != major bins rebinning {len(self.neminor)}')
            for i in range(len(nw)):
                nw[i],ns[i] = self.rebin1D(eminor,nw[i],ns[i],self.neminor[i])
            eminor = self.neminor
        else:
            eminor = [eminor for _ in range(len(nw))] 

        ne,nw,ns = self.linearize(nw,ns,eminor,emajor)
        nh = self.fillHistogram1D(ne,nw,ns,h.GetName()+"lin")
        self.savePlotData(eminor,emajor)
        return nh

    def linearize(self,w,s,eminor,emajor):
        nw = np.concatenate(w)
        ns = np.concatenate(s)
        ne = np.unique(np.concatenate([eminor[i]+i*eminor[i-1][-1] for i in range(len(eminor))]))
        return ne,nw,ns

    def split(self,e,w,s):
        if self.major == 'x':
            emajor = e[0]
            eminor = e[1]
        else:
            emajor = e[1]
            eminor = e[0]
            w = w.T            
            s = s.T            

        if self.nemajor is not None:
            compareTwoAxes(emajor,self.nemajor)
            xmajor = (emajor[1:]+emajor[:-1])/2
            idx = np.digitize(xmajor,self.nemajor) - 1
            emajor = self.nemajor
            nw = [np.sum(w[np.where(idx==i)],axis=0) for i in np.unique(idx)] 
            ns = [np.sum(s[np.where(idx==i)],axis=0) for i in np.unique(idx)] 
        else:
            nw = [arr.ravel() for arr in np.vsplit(w,w.shape[0])]
            ns = [arr.ravel() for arr in np.vsplit(s,s.shape[0])]

        return nw,ns,eminor,emajor

    def savePlotData(self,eminor,emajor):
        labels = []
        ylabel = 'HME'
        xpos = [0.] + [float(eminor[i][-1] + i * eminor[i-1][-1]) for i in range(len(eminor))]
        lines = xpos[1:-1]
        for i in range(len(xpos)-1):
            label = f"{emajor[i]} < {ylabel} < {emajor[i+1]}"
            poslabel = [(xpos[i] + 0.1 * (xpos[i+1]-xpos[i]))/xpos[-1],0.75]
            labels.append({'text':label,'position':poslabel,'size':26-2*len(eminor)})
        self.plotData = {'labels':labels,'lines':lines}

    def getPlotData(self):
        return self.plotData

class LinearizeSplitRebin2D(Linearize2D):
    def __init__(self,h,major='y',major_class=None,major_params=None,minor_class=None,minor_params=None):
        h2D = processHistogramToROOT(h,ROOT.TH2)
        super(LinearizeSplitRebin2D,self).__init__(major)
        if major == 'y':
            hmajor = h2D.ProjectionY()
            if major_class == 'Threshold':
                major_params[1] = [h.ProjectionY() for h in major_params[1]]
        if major == 'x':
            hmajor = h2D.ProjectionX()
            if major_class == 'Threshold':
                major_params[1] = [h.ProjectionX() for h in major_params[1]]
        
        majorObj = getattr(sys.modules[__name__], major_class)(hmajor,*major_params)

        self.nemajor = majorObj.ne
        
        nw,ns,eminor,emajor = self.split(*self.getContent2D(h2D))
        if minor_class == 'Threshold':
            extra_content = [self.split(*self.getContent2D(h)) for h in minor_params[1]]
            extra_content_nw = [extra[0] for extra in extra_content]
            extra_content_ns = [extra[1] for extra in extra_content]
            extra_content_eminor = [extra[2] for extra in extra_content]
                # Component per extra contribution

        self.neminor = []
        for i,(w_tmp,s_tmp) in enumerate(zip(nw,ns)):
            h_tmp = self.fillHistogram1D(eminor,w_tmp,s_tmp,'tmp{}'.format(i))
            if minor_class == 'Threshold':
                minor_params[1] = [self.fillHistogram1D(extra_content_eminor[j],extra_content_nw[j][0],extra_content_ns[j][0]) for j in range(len(extra_content))]
            minorObj = getattr(sys.modules[__name__], minor_class)(h_tmp,*minor_params)
            self.neminor.append(minorObj.ne)

    

