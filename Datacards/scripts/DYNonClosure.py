import os
import sys
import json
import numpy as np
import ROOT

class NonClosureSystematic:
    def __init__(self,factors):
        self.factors        = factors
        self.coefficients   = None
        self.covariance     = None

    def SetCategory(self,cat):
        if cat not in self.factors.keys():
            raise RuntimeError("Category {} not in the json file".format(cat))
        self.cat = cat
        self.coefficients = self.factors[cat]['coefficients']
        self.covariance = self.factors[cat]['covariance']
        self._fitComputations()

    def _fitComputations(self):
        self.fit = lambda x : self.coefficients[0] + self.coefficients[1] * x

        cov = np.array(self.covariance)
        assert cov.shape == (2,2)
        eigenValues , eigenVectors = np.linalg.eigh(cov)

        self.shape1_up   = lambda x : (self.coefficients[0] + np.sqrt(eigenValues[0]) * eigenVectors[0][0]) \
                                    + (self.coefficients[1] + np.sqrt(eigenValues[0]) * eigenVectors[1][0]) * x
        self.shape1_down = lambda x : (self.coefficients[0] - np.sqrt(eigenValues[0]) * eigenVectors[0][0]) \
                                    + (self.coefficients[1] - np.sqrt(eigenValues[0]) * eigenVectors[1][0]) * x
        self.shape2_up   = lambda x : (self.coefficients[0] + np.sqrt(eigenValues[1]) * eigenVectors[0][1]) \
                                    + (self.coefficients[1] + np.sqrt(eigenValues[1]) * eigenVectors[1][1]) * x
        self.shape2_down = lambda x : (self.coefficients[0] - np.sqrt(eigenValues[1]) * eigenVectors[0][1]) \
                                    + (self.coefficients[1] - np.sqrt(eigenValues[1]) * eigenVectors[1][1]) * x

    def _makeTGraphs(self):
        if not hasattr(self,'fit'):
            raise RuntimeError("Something is wrong, did you initialize the category (with SetCategory()) ?")
        x = np.linspace(0,1,101)
        return {'nominal'     : ROOT.TGraph(x.shape[0],x,self.fit(x)),
                'shape1_up'   : ROOT.TGraph(x.shape[0],x,self.shape1_up(x)),
                'shape1_down' : ROOT.TGraph(x.shape[0],x,self.shape1_down(x)),
                'shape2_up'   : ROOT.TGraph(x.shape[0],x,self.shape2_up(x)),
                'shape2_down' : ROOT.TGraph(x.shape[0],x,self.shape2_down(x))}

    def SaveFitToRoot(self,path='shapes.root'):
        if not path.endswith('.root'):
            path += '.root'

        F = ROOT.TFile(path,'RECREATE')
        for name , g in self._makeTGraphs().items():
            g.Write(name) 
        print ("Saved shapes to {}".format(path))

    def PlotFitToPdf(self,path='shapes.pdf'):
        if not path.endswith('.pdf'):
            path += '.pdf'
        
        C = ROOT.TCanvas("","",800,600)
        legend = ROOT.TLegend(0.65,0.65,0.85,0.85)
        graphDict = self._makeTGraphs()
        opt = ""
        colors = [1,2,2,4,4]
        styles = [1,2,3,2,3]
        minVal = min([g.GetHistogram().GetMinimum() for g in graphDict.values()])
        maxVal = max([g.GetHistogram().GetMaximum() for g in graphDict.values()])
        for i, (name, g) in enumerate(graphDict.items()):
            legend.AddEntry(g,name)
            g.SetLineWidth(2)
            g.SetLineStyle(styles[i])
            g.SetLineColor(colors[i])
            g.SetTitle("{};DNN output;Correction".format(self.cat))
            g.GetHistogram().GetYaxis().SetRangeUser(minVal,maxVal)
            g.GetHistogram().GetXaxis().SetRangeUser(0.,1.)
            g.Draw(opt)
            if not "same" in opt:
                opt += 'same'
        legend.Draw()
        C.Print(path)

        print ("Saved shapes to {}".format(path))


    def GetNonClosureShape(self,h):
        if not isinstance(h,ROOT.TH1):
            raise RuntimeError('Only works with ROOT TH1')
        if not hasattr(self,'fit'):
            raise RuntimeError("Something is wrong, did you initialize the category (with SetCategory()) ?")
    
        h_shape1_up   = h.Clone("shape1_up")
        h_shape1_down = h.Clone("shape1_down")
        h_shape2_up   = h.Clone("shape2_up")
        h_shape2_down = h.Clone("shape2_down")

        for i in range(1,h.GetNbinsX()+1):
            x = h.GetXaxis().GetBinCenter(i)
            y = h.GetBinContent(i)
            h_shape1_up.SetBinContent(i, y * self.shape1_up(x))
            h_shape1_down.SetBinContent(i, y * self.shape1_down(x))
            h_shape2_up.SetBinContent(i, y * self.shape2_up(x))
            h_shape2_down.SetBinContent(i, y * self.shape2_down(x))

        return ((h_shape1_up,h_shape1_down),(h_shape2_up,h_shape2_down))

    def __call__(self,h):
        assert isinstance(h,ROOT.TH1)
        if not hasattr(self,'fit'):
            raise RuntimeError("Something is wrong, did you initialize the category (with SetCategory()) ?")
        for i in range(1,h.GetNbinsX()+1):
            x = h.GetXaxis().GetBinCenter(i)
            y = h.GetBinContent(i)
            h.SetBinContent(i,y*self.fit(x))

    @classmethod 
    def load_from_json(cls,path_to_json):
        with open(path_to_json,'r') as handle:
            factors = json.load(handle)
        return cls(factors)



# HOW TO USE 
# 1) Load json file into the object
#
# objDY = NonClosureSystematic.load_from_json(path_to_json)
#
# 2) Select correct category ( == key name in the json), eg
#
# objDY.SetCategory("HH_inclusive_DY_VVV_2016")
#
# Note : you can modify the keys in the json file to fit your FW, use the modified ones here then
#
# 3) Get the non closure systematic shapes
#
# h_up, h_down = objDY.GetNonClosureShape(h)
#
# where h is a ROOT TH1 histogram and is the nominal histogram of the category
#
# 4) Scale each DY histogram (both nominal and systematics)
#
# objDY(h) where h is any histogram (either the nominal or a systematic)
#
# Note : the script direcly changes the bin content, so the input histogram will have changed 
# (python passes the object by reference), if you want to return the modified histogram but 
# keep the initial histogram constant, consider using either the ROOT Clone or the copy 
# python package
#
# Note * Important * : GetNonClosureShape() produces the up and down fluctuations where the 
# slope corrections are already taken into account. Because of that, make sure that :
#    - the nominal histogram has not been slope corrected already 
#       -> GetNonClosureShape() needs to be called before you scale the nominal and the other shapes
#    - the two (up and down) non closure fluctuations from GetNonClosureShape() are not corrected afterwards 
#           (eg because you ran in a loop over all the systematics)
#
# 5) profit
#
# BONUS) two methods have been used to visualize the linear fit and associated uncertainties 
#
#   - objDY.SaveFitToRoot(path_to_root) will save the nominal fit and two the shapes to a ROOT file
#   - objDY.PlotFitToPdf(path_to_pdf) will plot the nominal and the two shapes into a pdf
#
