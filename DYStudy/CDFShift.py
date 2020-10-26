import os
import sys
import math
from copy import copy,deepcopy
import ROOT
from pprint import pprint
from array import array

#####
# Loosely based on 
# https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit/blob/102x/interface/th1fmorph.h
#####

class CDFShift:
    def __init__(self,hist1,hist2,hist_int,norm=1.,additional_hist={},verbose=False):
        self.h1 = copy(hist1)
        self.h2 = copy(hist2)
        self.h3 = copy(hist_int)
        self.padHistNegativeBins(self.h1)
        self.padHistNegativeBins(self.h2)
        self.padHistNegativeBins(self.h3)
        self.norm = norm
        self.verbose = verbose
        print ('Output histogram will be scaled by %0.5f'%self.norm)

        # Normalize #
        if self.h1.Integral() != 0:
            self.h1.Scale(1./self.h1.Integral())
        if self.h2.Integral() != 0:
            self.h2.Scale(1./self.h2.Integral())
        if self.h3.Integral() != 0:
            self.h3.Scale(1./self.h3.Integral())

        # Get CDF #
        self.cdf1 = self.produceCDF(self.h1,'cdf1') 
        self.cdf2 = self.produceCDF(self.h2,'cdf2') 
        self.cdf3 = self.produceCDF(self.h3,'cdf3') 

        self.additional_hist = deepcopy(additional_hist)
        for name in self.additional_hist:
            self.padHistNegativeBins(self.additional_hist[name])
            if self.additional_hist[name].Integral() != 0:
                self.additional_hist[name].Scale(1./self.additional_hist[name].Integral())
        self.additional_cdf = {'cdf_'+name:self.produceCDF(hist,'cdf_'+name) for name,hist in self.additional_hist.items()}

        self.rootsave = {'h1'  : self.h1,
                         'h2'  : self.h2,
                         'h3'  : self.h3,
                         'cdf1': self.cdf1,
                         'cdf2': self.cdf2,
                         'cdf3': self.cdf3}
        self.rootsave.update(self.additional_cdf)
        self.rootsave.update(self.additional_hist)

        # Get shift for nominal, up and down #
        edges1,cont1,err1 = self.getArray(self.cdf1)
        edges2,cont2,err2 = self.getArray(self.cdf2)
        edges3,cont3,err3 = self.getArray(self.cdf3)
        cont1_up          = [cont1[i]+err1[i] for i in range(len(cont1))]
        cont2_up          = [cont2[i]+err2[i] for i in range(len(cont2))]
        cont1_down        = [max(cont1[i]-err1[i],0.) for i in range(len(cont1))]
        cont2_down        = [max(cont2[i]-err2[i],0.) for i in range(len(cont2))]
        if not cont1_up == sorted(cont1_up):
            print (cont1_up)
            cont1_up.sort()
            print ('[Warning] Had to sort cont1_up')
        if not cont1_down == sorted(cont1_down):
            print (cont1_down)
            cont1_down.sort()
            print ('[Warning] Had to sort cont1_down')
        if not cont2_up == sorted(cont2_up):
            print (cont2_up)
            cont2_up.sort()
            print ('[Warning] Had to sort cont2_up')
        if not cont2_down == sorted(cont2_down):
            print (cont2_down)
            cont2_down.sort()
            print ('[Warning] Had to sort cont2_down')

        if cont1_down[-1]<2:
            cont1_down[-1] = 1
            print ('[Warning] cdf 1 down variation does not reach 1, artificially fix that')
        if cont2_down[-1]<1:
            cont2_down[-1] = 1
            print ('[Warning] cdf 2 down variation does not reach 1, artificially fix that')

        if self.verbose:
            print ("Producing the nominal shift")
        shift, ypos             = self.getShiftFromBinning(edges1,cont1,edges2,cont2)
        if self.verbose:
            print ("Producing the up shift")
        shift_up, ypos_up       = self.getShiftFromBinning(edges1,cont1_up,edges2,cont2_up)
        if self.verbose:
            print ("Producing the down shift")
        shift_down, ypos_down   = self.getShiftFromBinning(edges1,cont1_down,edges2,cont2_down)

        self.g,self.ginv              = self.plotShift(shift,ypos)
        self.g_up,self.ginv_up        = self.plotShift(shift_up,ypos_up)
        self.g_down,self.ginv_down    = self.plotShift(shift_down,ypos_down)
        
        self.rootsave['g']          = self.g
        self.rootsave['ginv']       = self.ginv
        self.rootsave['g_up']       = self.g_up
        self.rootsave['ginv_up']    = self.ginv_up
        self.rootsave['g_down']     = self.g_down
        self.rootsave['ginv_down']  = self.ginv_down

        nedges3,ncont3,nerr3 = self.applyShift(edges3,cont3,err3,shift,ypos,shift_up,ypos_up,shift_down,ypos_down)

        nh3raw,ncdf3raw = self.makeHistAndCDF(nedges3,ncont3,nerr3,'nh3raw')
        
        self.rootsave['nh3raw']     = nh3raw
        self.rootsave['ncdf3raw']   = ncdf3raw

        nh3edges, nh3cont, nh3err = self.getArray(nh3raw)
        nh3edges, nh3cont, nh3err = self.rebin(nh3edges, nh3cont, nh3err, edges3)

        self.nh3 = ROOT.TH1F("nh3","nh3",len(nh3edges)-1,array('d',nh3edges))
        for i in range(1,len(nh3cont)+1):
            self.nh3.SetBinContent(i,nh3cont[i-1])
            self.nh3.SetBinError(i,nh3err[i-1])
        self.ncdf3 = self.produceCDF(self.nh3,'ncdf3')
        self.nh3.Scale(self.norm)

        self.rootsave['nh3']        = self.nh3
        self.rootsave['ncdf3']      = self.ncdf3

    @staticmethod
    def padHistNegativeBins(h):
        for i in range(1,h.GetNbinsX()+1):
            if h.GetBinContent(i)<0:
                h.SetBinContent(i,0)
                h.SetBinError(i,0)

    @staticmethod
    def produceCDF(h,name):
        cdfcont = [0]
        cdferr = []
        edges = []
        for i in range(1,h.GetNbinsX()+1):
            edges.append((h.GetXaxis().GetBinLowEdge(i),h.GetXaxis().GetBinUpEdge(i)))
            if h.GetBinContent(i) < 0:
                cdfcont.append(0. + cdfcont[-1])
            else:
                cdfcont.append(h.GetBinContent(i) + cdfcont[-1])
            cdferr.append(math.sqrt(sum([h.GetBinError(j)**2 for j in range(1,i+1)])))
        del cdfcont[0]
        edges = [e[0] for e in edges] + [edges[-1][1]]
        cdf = ROOT.TH1F(name,name,len(edges)-1,array('d',edges))
        for i in range(1,cdf.GetNbinsX()+1):
            cdf.SetBinContent(i,cdfcont[i-1])
            cdf.SetBinError(i,cdferr[i-1])
        cdf.SetTitle(h.GetTitle())
        cdf.GetXaxis().SetTitle(h.GetXaxis().GetTitle())
        cdf.GetYaxis().SetTitle(h.GetYaxis().GetTitle())
        return cdf

    @staticmethod
    def getArray(h):
        edges = [h.GetXaxis().GetBinLowEdge(i) for i in range(1,h.GetNbinsX()+2)]
        content = [h.GetBinContent(i) for i in range(1,h.GetNbinsX()+1)]
        error = [h.GetBinError(i) for i in range(1,h.GetNbinsX()+1)]
        if any([c<0 for c in content]):
            print ('Negative values in array, will crop to 0')
            for i in range(len(content)):
                if content[i] < 0:
                    print ('... idx %d -> %0.10f in [%.10f,%.10f]'%(i,content[i],edges[i],edges[i+1]))
        content = [c if c>=0 else 0 for c in content] # Crop negative values
        return edges,content,error

    @staticmethod
    def getFirstAndLastBins(cdf):
        # Find last bin that still increases #
        iil = len(cdf)-1
        while(cdf[iil-1] >= cdf[iil]):
            iil -= 1
        # Find beginning of curves (non-zero) #
        iif = 0
        while (cdf[iif+1] <= cdf[iif]):
            iif += 1
        return iif,iil


    def getShiftFromBinning(self,e1,c1,e2,c2):
        self.assertMonoticallyIncreasing(c1,"cdf1")
        self.assertMonoticallyIncreasing(c2,"cdf2")

        i1f,i1l = self.getFirstAndLastBins(c1)
        i2f,i2l = self.getFirstAndLastBins(c2)

        if self.verbose:
            print ("Last bins i1l %d i2l %d"%(i1l,i2l))
            print ("First bins i1f %d i2f %d"%(i1f,i2f))

        i1 = i1f
        i2 = i2f
        shift = []
        ypos = []

        select12 = -1
        while (i1 < i1l or i2 < i2l):
            if self.verbose:
                print ('-'*20)
                print ("Index i1 %d/%d [c1 = %0.10f] i2 %d/%d [c2 = %0.10f]"%(i1,i1l,c1[i1],i2,i2l,c2[i2]))
            # Check which hist to be select y
            if (c1[i1] <= c2[i2] or i2 == i2l) and i1 < i1l :
                while c1[i1+1] <= c1[i1] and i1 < i1l: # Avoid flat region
                    i1 += 1
                select12 = 1
            elif i2<i2l:
                while c2[i2+1] <= c2[i2] and i2 < i2l: # Avoid flat region
                    i2 += 1
                select12 = 2
            else:
                raise RuntimeError('Something fishy')
            if self.verbose:
                print ("Looking at h%d"%select12)
            # Take point in h1, find corresponding in h2 #
            if select12 == 1:
                x1 = e1[i1+1]
                y = c1[i1]
                if self.verbose:
                    print ('x1 = %0.10f -> y = %.10f'%(x1,y))
                if i2 == 0: # CDF does not start at 0
                    x2 = e2[0]
                    if self.verbose:
                        print ('First bin non zero, taking its x')
                elif i2 == i2l: 
                    #break # TODO : test
                    x2 = e2[i2l+1]
                    if self.verbose:
                        print ('Last bin, taking last edge')
                else:
                    x20 = e2[i2]
                    x21 = e2[i2+1]
                    y20 = c2[i2-1]
                    y21 = c2[i2]
                    if self.verbose:
                        print ('-> Found x2 in range [%.10f,%.10f] : y2 in range [%.10f,%.10f]'%(x20,x21,y20,y21))
                    assert y >= y20 and y <= y21
                    if y21>y20:
                        x2 = x20 + (x21-x20) * (y-y20)/(y21-y20)
                    else:
                        x2 = x21
                    if self.verbose:
                        print ('-> Interpolated x2 = %.10f'%x2)

                i1 += 1

            # Tak e point in h2, find corresponding in h1 #
            elif select12 == 2:
                x2 = e2[i2+1]
                y = c2[i2]
                if self.verbose:
                    print ('x2 = %0.10f -> y = %0.10f'%(x2,y))
                if i1 == 0: # CDF does not start at 0
                    x1 = e1[0]
                    if self.verbose:
                        print ('First bin non zero, taking its x')
                elif i1 == i1l:
                    #break # TODO : test
                    x1 = e1[i1l+1]
                    if self.verbose:
                        print ('Last bin, taking last edge')
                else:
                    x10 = e1[i1]
                    x11 = e1[i1+1]
                    y10 = c1[i1-1]
                    y11 = c1[i1]
                    if self.verbose:
                        print ('-> Found x1 in range [%.10f,%.10f] : y1 in range [%.10f,%.10f]'%(x10,x11,y10,y11))
                    assert y >= y10 and y <= y11
                    if y11>y10:
                        x1 = x10 + (x11-x10) * (y-y10)/(y11-y10)
                    else:
                        x1 = x11
                    if self.verbose:
                        print ('-> Interpolated x1 = %.10f'%x1)
                i2 += 1
            else:
                raise RuntimeError("Something went wrong")
            ypos.append(y)
            shift.append(x2-x1)
            if self.verbose:
                print ("x1 = %.10f, x2 = %.10f -> shift = %.10f at y = %.10f"%(x1,x2,x2-x1,y))
        if self.verbose:
            print ('Shift values :')
            for i in range(len(ypos)):
                print ('  y = %.10f -> shift = %.10f'%(ypos[i],shift[i]))
        if ypos[0] != 0.:
            print ('Warning : shift y position did not start at 0, will add that to shift')
            ypos.insert(0,0.)
            shift.insert(0,shift[0])
        if ypos[-1] < 1.:
            print ('Warning : shift y position does not stop after 1, will add that to shift')
            ypos.insert(len(ypos),1.)
            shift.insert(len(shift),shift[-1])
            
        return shift,ypos

    @staticmethod
    def assertMonoticallyIncreasing(cdf,name):
        for i in range(len(cdf)-1):
            if cdf[i]>cdf[i+1]+1e-6 and cdf[i]<1: # Sometimes errors at edge
                print (cdf)
                raise RuntimeError("CDF %s not monotically increasing : cdf[%d]=%.10f -> cdf[%d]=%.10f"%(name,i,cdf[i],i+1,cdf[i+1]))

    def applyShift(self,e3,c3,err3,sn,yn,su,yu,sd,yd):
        self.assertMonoticallyIncreasing(c3,"cdf3")
        if self.verbose:
            print ('='*30)
            print ('Applying shift')

        inf,inl = self.getFirstAndLastBins(yn)
        iuf,iul = self.getFirstAndLastBins(yu)
        idf,idl = self.getFirstAndLastBins(yd)

        ne3 = []
        nc3 = []
        nerr3 = []
        nge3 = []
        ngc3 = []
        iin = inf
        iiu = iuf
        iid = idf
        shift_var = [0]*len(c3)
        kept_idx = []
        for i in range(len(c3)):
            y = c3[i]
            if self.verbose:
                print ("Looking at y = %.10f"%y)
            while y > yn[iin+1] and iin < len(yn)-2:
                iin += 1
            y1 = yn[iin+1] # upper y of sn bin
            y0 = yn[iin]   # lower y of sn bin
            x0 = sn[iin]  # lower sn
            x1 = sn[iin+1]# upper sn

            while y > yu[iiu+1]:
                iiu += 1
            while y> yd[iid+1]:
                iid += 1

            up_var = su[iiu] + (y-yu[iiu])/(yu[iiu+1]-yu[iiu]) * (su[iiu+1]-su[iiu])
            down_var = sd[iid] + (y-yd[iid])/(yd[iid+1]-yd[iid]) * (sd[iid+1]-sd[iid])
            shift_var[i] = max(shift_var[i-1],math.sqrt(0.5*((sn[iin]-up_var)**2+(sn[iin]-down_var)**2)))
            #shift_var[i] = math.sqrt(0.5*((sn[iin]-up_var)**2+(sn[iin]-down_var)**2))

            assert y>=yu[iiu] and y<=yu[iiu+1]
            assert y>=yd[iid] and y<=yd[iid+1]
            
            if y > 1 and y>y1:
                y0 = y
                y1 = y
                x0 = sn[iin+1]
                x1 = sn[iin+1]
                if self.verbose:
                    print('Warning : y > 1  and reached the end of shift array, will use the last value %.10f -> %.10f'%(yn[iin+1],sn[iin+1]))

            if self.verbose:
                print ("-> sn at y0 = %.10f : %.10f\n-> sn at y1 = %.10f : %.10f"%(y0,x0,y1,x1))
            assert y >= y0 and y <= y1
            if y1-y0 == 0:
                x = 0
            else:
                x = x0 + (x1-x0) * (y-y0)/(y1-y0) # sn to be applied
            if self.verbose:
                print ("-> Shift is %0.10f"%x)
            
            if len(ne3)>0 and e3[i]+x < ne3[-1]:
                print ('Warning : idx %d %.10f [y = %.10f] < %.10f]'%(i,ne3[-1],y,e3[i]+x))
            elif e3[i]+x > e3[-1]:
                print ('Warning : idx %d %.10f after last bin %.10f'%(i,e3[i]+x,e3[-1]))
            else:
                kept_idx.append(i)
                nc3.append(y)
                ne3.append(e3[i]+x)
                #nerr3.append(math.sqrt(err3[i]**2+shift_var[i]**2))
                nerr3.append(err3[i])
            ngc3.append(y)
            nge3.append(e3[i]+x)
        ne3.append(e3[-1])
        nge3.append(e3[-1])
        self.assertMonoticallyIncreasing(nc3,"ncdf3")
        if self.verbose:
            print ('New cdf3 content :')
            for i in range(len(nc3)):
                print (' bin [%.10f,%.10f] -> %.10f'%(ne3[i],ne3[i+1],nc3[i]))

        self.gncdf3 = ROOT.TGraph(len(ngc3))
        for i in range(len(ngc3)):
            self.gncdf3.SetPoint(i,nge3[i+1],ngc3[i])
        self.rootsave['gncdf3'] = self.gncdf3
    
        self.cdf_bin_err = ROOT.TGraph(len(kept_idx))
        for i,ik in enumerate(kept_idx):
            self.cdf_bin_err.SetPoint(i,ne3[i+1],err3[ik])
        self.rootsave['cdf_bin_err'] = self.cdf_bin_err

        self.cdf_shift_err = ROOT.TGraph(len(kept_idx))
        for i,ik in enumerate(kept_idx):
            self.cdf_shift_err.SetPoint(i,ne3[i+1],shift_var[ik])
        self.rootsave['cdf_shift_err'] = self.cdf_shift_err

        self.cdf_tot_err = ROOT.TGraph(len(kept_idx))
        for i in range(len(kept_idx)):
            self.cdf_tot_err.SetPoint(i,ne3[i+1],nerr3[i])
        self.rootsave['cdf_tot_err'] = self.cdf_tot_err

        return ne3,nc3,nerr3

    @staticmethod
    def makeHistAndCDF(edges,cdfcont,cdferr,name):
        hcont = []
        herr = []
        for i in range(len(cdfcont)-2,-1,-1):
            hcont.append(cdfcont[i+1]-cdfcont[i])
            if hcont != 0:
                herr.append(math.sqrt(cdferr[i+1]**2-cdferr[i]**2))
            else:
                herr.append(0)
        hcont.append(cdfcont[0])
        herr.append(cdferr[0])
        hcont.reverse()
        herr.reverse()
    
        h = ROOT.TH1F(name,name,len(edges)-1,array('d',edges))
        cdf = ROOT.TH1F(name+'cdf',name+'cdf',len(edges)-1,array('d',edges))
        for i in range(1,len(hcont)+1):
            h.SetBinContent(i,hcont[i-1])
            h.SetBinError(i,herr[i-1])
        for i in range(1,len(cdfcont)+1):
            cdf.SetBinContent(i,cdfcont[i-1])
            cdf.SetBinError(i,cdferr[i-1])
        return h,cdf
            
    @staticmethod
    def plotShift(shift,ypos):
        gshift = ROOT.TGraph(len(shift))
        gshift_inv = ROOT.TGraph(len(shift))
        for i in range(len(shift)):
            gshift.SetPoint(i,shift[i],ypos[i])
            gshift_inv.SetPoint(i,ypos[i],shift[i])
        gshift.SetTitle(';shift;CDF y')
        gshift_inv.SetTitle(';CDF y;shift')
        return gshift,gshift_inv

    def saveToRoot(self,rootname):
        if not rootname.endswith('.root'):
            rootname += '.root'
        f = ROOT.TFile(rootname,"RECREATE")
        for name,obj in self.rootsave.items():
            obj.SetName(name)
            obj.Write()
        print ("CDFShift : saved root in %s"%rootname)
        f.Close()

    def rebin(self,edges1,cont1,err1,edges2):
        assert len(cont1) == len(edges1)-1
        assert len(err1) == len(edges1)-1
        cont2 = [0.]*(len(edges2)-1)
        err2 = [0.]*(len(edges2)-1)
        i1 = 0
        i2 = 0
        if self.verbose:
            print ('='*30)
            print ('Rebin')
        while (i1 < len(edges1)-1 and i2 < len(edges2)-1):
            if self.verbose:
                print ('-'*20)
                print ('Looking at bin i1 %d/%d edges = [%.3f,%.3f] -> %.5f (cont), %.5f (err squared), compared to bin i2 %d/%d edges = [%.3f,%.3f] -> %.5f (cont), %.5f (err squared)'                            %(i1,len(edges1)-1,edges1[i1],edges1[i1+1],cont1[i1],err1[i1]**2,i2,len(edges2)-1,edges2[i2],edges2[i2+1],cont2[i2],err2[i2]**2))
            if edges1[i1+1] <= edges2[i2+1]: 
                if self.verbose:
                    print ('edges1 <= edges2 : adding %.5f (cont) and %.5f (err squared) to i2 %d edges = [%.3f,%.3f]'%(cont1[i1],err1[i1]**2,i2,edges2[i2],edges2[i2+1]))
                cont2[i2] += cont1[i1]
                err2[i2] += err1[i1]**2
            else:
                pr = cont1[i1]
                sve = edges1[i1]
                si2 = i2
                cont_frac = []
                nums = []
                while edges1[i1+1] >= edges2[i2+1] and i2 < len(edges2)-2:
                    pl = (edges2[i2+1]-sve)/(edges1[i1+1]-sve) * pr
                    nums.append(edges2[i2+1]-sve)
                    sve = edges2[i2+1]
                    pr -= pl
                    assert pl >= 0
                    assert pr >= 0
                    cont_frac.append(pl)
                    i2 += 1
                cont_frac.append(pr)
                nums.append((edges1[i1+1]-edges1[i1])-sum(nums))
                for cf,num in zip(cont_frac,nums):
                    cont2[si2] += cf
                    add_err = num**2/sum([n**2 for n in nums]) * err1[i1]**2
                    err2[si2] += add_err
                    if self.verbose:
                        print ('edges1 > edges2 : adding %.5f (cont) and %.5f (err squared) to i2 %d edges = [%.3f,%.3f]'%(cf,add_err,si2,edges2[si2],edges2[si2+1]))
                    si2 += 1
            i1 += 1
        assert abs(sum(cont1)-sum(cont2))<1e-9
        assert abs(sum([e**2 for e in err1])-sum(err2))<1e-9
        err2 = [math.sqrt(e) for e in err2]
        return edges2,cont2,err2


    def returnHist(self):
        print ("\nWeight 2b integral : %0.5f\n"%self.weight_2b.Integral())
        return self.nh3
        
