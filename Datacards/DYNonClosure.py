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
        self.coefficients = self.factors[cat]['coefficients']
        self.covariance = self.factors[cat]['covariance']
        self.formula = lambda x : self.coefficients[0] + self.coefficients[1] * x

    def GetNonClosureShape(self,h):
        if not isinstance(h,ROOT.TH1):
            raise RuntimeError('Only works with ROOT TH1')
        if self.coefficients is None or self.covariance is None:
            raise RuntimeError("Something is wrong, did you initialize the category (with SetCategory()) ?")
        cov = np.array(self.covariance)
        eigenValues , eigenVectors = np.linalg.eig(cov)
    
        formula_up   = lambda x : self.coefficients[0]* \
                                    (1+eigenVectors[:,0].dot(np.sqrt(eigenValues))) \
                                + self.coefficients[1]* \
                                    (1+eigenVectors[:,1].dot(np.sqrt(eigenValues))) * x

        formula_down = lambda x : self.coefficients[0]* \
                                    (1-eigenVectors[:,0].dot(np.sqrt(eigenValues))) \
                                + self.coefficients[1]* \
                                    (1-eigenVectors[:,1].dot(np.sqrt(eigenValues))) * x

        h_up   = h.Clone("up")
        h_down = h.Clone("down")

        for i in range(1,h.GetNbinsX()+1):
            x = h.GetXaxis().GetBinCenter(i)
            y = h.GetBinContent(i)
            h_up.SetBinContent(i,y*(1+(formula_up(x))))
            h_down.SetBinContent(i,y*(1+formula_down(x)))

        return h_up,h_down 

    def __call__(self,h):
        assert isinstance(h,ROOT.TH1)
        fit = lambda x : self.coefficients[0] + self.coefficients[1] * x
        for i in range(1,h.GetNbinsX()+1):
            x = h.GetXaxis().GetBinCenter(i)
            y = h.GetBinContent(i)
            h.SetBinContent(i,y*(1+self.formula(x)))

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


