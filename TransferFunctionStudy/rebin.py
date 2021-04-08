import numpy as np
from array import array
import ROOT

def getArray(h):
    xAxis = h.GetXaxis()
    yAxis = h.GetYaxis()
    Nx = xAxis.GetNbins()
    Ny = yAxis.GetNbins()
    content = np.zeros((Nx,Ny))
    error   = np.zeros((Nx,Ny))
    for ix in range(0,Nx):
        for iy in range(0,Ny):
            content[ix,iy] = h.GetBinContent(ix+1,iy+1)
            error[ix,iy] = h.GetBinError(ix+1,iy+1)
    edges = [np.array([xAxis.GetBinLowEdge(i) for i in range(1,Nx+2)]),
             np.array([yAxis.GetBinLowEdge(i) for i in range(1,Ny+2)])]
    return content,error,edges


def fillHist(content,error,edges,name='h'):
    assert len(edges) == 2
    h = ROOT.TH2F(name,name,
                  edges[0].shape[0]-1,array('d',edges[0]),
                  edges[1].shape[0]-1,array('d',edges[1]))
    for ix in range(edges[0].shape[0]-1):
        for iy in range(edges[1].shape[0]-1):
            h.SetBinContent(ix+1,iy+1,content[ix,iy])
            h.SetBinError(ix+1,iy+1,error[ix,iy])

    return h



def rebin2D(h,xnew,ynew,name='hist',capMin=None,capMax=None):
    array,error,edges = getArray(h)
    x = edges[0]
    y = edges[1]

    xnew = np.array(xnew)
    ynew = np.array(ynew)
    edgesnew = [xnew,ynew]

    if xnew[0] != x[0]:
        raise RuntimeError("New x axis does not start at same point as old one")
    if ynew[0] != y[0]:
        raise RuntimeError("New y axis does not start at same point as old one")
    if xnew[-1] != x[-1]:
        raise RuntimeError("New x axis does not end at same point as old one")
    if ynew[-1] != y[-1]:
        raise RuntimeError("New y axis does not end at same point as old one")
    if not np.all(xnew[:-1] <= xnew[1:]):
        raise RuntimeError("New x axis values not increasing")
    if not np.all(ynew[:-1] <= ynew[1:]):
        raise RuntimeError("New y axis values not increasing")
    if not np.isin(xnew,x).all():
        raise RuntimeError('New x bin edges do not match the ones in the histogram')
    if not np.isin(ynew,y).all():
        raise RuntimeError('New y bin edges do not match the ones in the histogram')

    xc = (x[1:]+x[:-1])/2
    yc = (y[1:]+y[:-1])/2

    idx = np.digitize(xc,xnew) - 1
    idy = np.digitize(yc,ynew) - 1

    arraynew = np.zeros((xnew.shape[0]-1,ynew.shape[0]-1))
    errornew = np.zeros((xnew.shape[0]-1,ynew.shape[0]-1))
    for ix in range(arraynew.shape[0]):
        for iy in range(arraynew.shape[1]):
            arraynew[ix,iy] = array[np.ix_(idx==ix,idy==iy)].sum()
            errornew[ix,iy] = (error[np.ix_(idx==ix,idy==iy)]**2).sum()
    errornew = np.sqrt(errornew)

    if capMin is not None:
        arraynew[arraynew < capMin] = capMin
        errornew[errornew < capMin] = capMin
    if capMax is not None: 
        arraynew[arraynew > capMax] = capMax
        errornew[errornew > capMax] = capMax
    #import IPython
    #IPython.embed()
    print("Original bin content sum = {:0.6e}, rebinned content sum = {:0.6e} -> relative difference = {:0.6e}".format(array.sum(),arraynew.sum(),abs(arraynew.sum()-array.sum())/array.sum()))
    print("Original bin error squared sum = {:0.6e}, rebinned error squared sum = {:0.6e} -> relative difference = {:0.6e}".format((error**2).sum(),(errornew**2).sum(),abs((errornew**2).sum()-(error**2).sum())/(error**2).sum()))
    if abs(arraynew.sum()-array.sum())/array.sum() > 1e-5:
        raise RuntimeError("Original bin content sum = {:0.6e}, rebinned content sum = {:0.6e} -> relative difference = {:0.6e} > 1e-5".format(array.sum(),arraynew.sum(),abs(arraynew.sum()-array.sum())/array.sum()))
    if abs((errornew**2).sum()-(error**2).sum())/(error**2).sum() > 1e-5:
        raise RuntimeError("Original bin error squared sum = {:0.6e}, rebinned error squared sum = {:0.6e} -> relative difference = {:0.6e} > 1e-5".format((error**2).sum(),(errornew**2).sum(),abs((errornew**2).sum()-(error**2).sum())/(error**2).sum()))

    del array,error,edges
    hnew = fillHist(arraynew,errornew,edgesnew,name)
    del arraynew,errornew,edgesnew

    return hnew

def rebin1D(h,xn,name='hist'):
    hn = h.Rebin(len(xn)-1,name,np.array(xn))
    return hn

#h = ROOT.TH2F("t","t",10,0,10,10,0,10)
#xAxis = h.GetXaxis()
#yAxis = h.GetYaxis()
#Nx = xAxis.GetNbins()
#Ny = yAxis.GetNbins()
#iN = 0
#for ix in range(1,Nx+1):
#    for iy in range(1,Ny+1):
#        h.SetBinContent(ix,iy,iN)
#        iN += 1
#h.SetContour(100)
#
#
#xn = np.r_[0:3:4j, 4:8:3j, 10]
#yn = np.r_[0:6:3j,10]
#
#hnew = rebin2D(h,xn,yn)
#
#C1 = ROOT.TCanvas()
#C1.cd()
#h.Draw("colz")
#C2 = ROOT.TCanvas()
#C2.cd()
#hnew.Draw("colz")


