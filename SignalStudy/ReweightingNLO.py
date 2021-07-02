import os
import sys
import json
import copy
import ROOT 
import glob
import argparse

import  NonResonantModelNLO

class ReweightingNLO:
    def __init__(self,filePath,histName):
        self.filePath = filePath
        inputFile = ROOT.TFile(self.filePath,"READ")
        if histName not in [key.GetName() for key in inputFile.GetListOfKeys()]:
            raise RuntimeError("Could not find histogram {} in file {}".format(histName,filePath))
        self.hist = copy.deepcopy(inputFile.Get(histName))
        inputFile.Close()
        self.Nevtot = self.hist.Integral()
        self.model = NonResonantModelNLO.NonResonantModelNLO()
        self.model.ReadCoefficients("pm_pw_NLO_Ais_13TeV_V2.txt")
        self.couplings = None

        self.mHHEdges,self.cosThetaEdges = self._getBinEdges(self.hist)

    @staticmethod
    def _getBinEdges(hist):
        xEdges = [(hist.GetXaxis().GetBinLowEdge(x),hist.GetXaxis().GetBinUpEdge(x)) for x in range(1,hist.GetNbinsX()+1)]
        yEdges = [(hist.GetYaxis().GetBinLowEdge(y),hist.GetYaxis().GetBinUpEdge(y)) for y in range(1,hist.GetNbinsY()+1)]

        return xEdges,yEdges

    def _loopOverBins(self):
        content = {}
        sumWeight = 0.
        for (mHHLow,mHHUp) in self.mHHEdges:
            content[(mHHLow,mHHUp)] = {}
            for (cosThetaLow,cosThetaUp) in self.cosThetaEdges:
                # Get Bin center #
                mHH = (mHHUp+mHHLow)/2
                cosTheta = (cosThetaUp+cosThetaLow)/2
                # Get Event content #
                Nev = self.hist.GetBinContent(self.hist.FindBin(mHH,cosTheta))
                if self.couplings is None:
                    raise RuntimeError("A set of couplings or a benchmark must have been chosen")
                # Get differential Xsec #
                XS = self.model.getDifferentialXS2D(mHH, cosTheta, *self.couplings)
                # Scale by bin width #
                Noutputev = XS * self.model.getmHHbinwidth(mHH) * self.model.getcosthetabinwidth(cosTheta) 
                # Get total Xsec #
                XStot = self.model.getTotalXS(*self.couplings)
                # Get Weight #
                if Nev != 0.:
                    weight = Noutputev/Nev * self.Nevtot/XStot
                else: 
                    weight = 0.
                # Record #
                content[(mHHLow,mHHUp)][(cosThetaLow,cosThetaUp)] = weight
                # Get new sum of weight for sanity check #
                sumWeight += Nev * weight
        if abs(sumWeight-self.Nevtot)/self.Nevtot > 0.05:
            raise RuntimeError("File {} with set of couplings [{}] :\nSum weights passed from {:.2f} to {:.2f} ({:3.2f}% change)".format(
                                            self.filePath,
                                            ['{:.1f}'.format(x) for x in self.couplings],
                                            self.Nevtot,
                                            sumWeight,
                                            abs(sumWeight-self.Nevtot)/self.Nevtot*100))
        return content

    def SaveToJson(self, name):
        content = self._loopOverBins()
        print ('Couplings : ',self.couplings)

        xbinning = []
        ybinning = []
        data = []
        for (mHHLow,mHHUp) in content.keys():
            if mHHLow not in xbinning:
                xbinning.append(mHHLow)
            if mHHUp not in xbinning:
                xbinning.append(mHHUp)
            mHH_entry = {'bin':[mHHLow,mHHUp],'values':[]}
            for (cosThetaLow,cosThetaUp),weight in content[(mHHLow,mHHUp)].items():
                if cosThetaLow not in ybinning:
                    ybinning.append(cosThetaLow)
                if cosThetaUp not in ybinning:
                    ybinning.append(cosThetaUp)
                mHH_entry['values'].append({'bin':[cosThetaLow,cosThetaUp],'value':weight,'error_low':0.,'error_high':0.})
            data.append(mHH_entry)

        json_content = {'dimension': 2, 
                        'variables': ['Eta', 'Pt'], 
                        'binning': {'x': xbinning,
                                    'y': ybinning},
                        'data': data, 
                        'error_type': 'absolute'}                                                           
        with open(f'{name}.json','w') as handle:
            json.dump(json_content,handle,indent=2)
        print (f"Saved to file {name}.json")


    def SetBenchmark(self,idx):
        if idx not in [i for i in range(0,13)]:
            raise RuntimeError("The possible benchmark index ranges from 0 (SM) to 12")
        self.couplings = self.model.getBenchmark(idx)

    def SetCouplings(self,couplings):
        """ ["kl", "kt", "c2", "cg", "c2g"] """
        if len(couplings) != 5:
            raise RuntimeError("There must be 5 couplings")
        self.couplings = couplings


couplings = {
    # Basic BM #
    'SM'  : 0,
    '1'   : 1,
    '2'   : 2,
    '3'   : 3,
    '4'   : 4,
    '5'   : 5,
    '6'   : 6,
    '7'   : 7,
    '8'   : 8,
    '9'   : 9,
    '10'  : 10,
    '11'  : 11,
    '12'  : 12,
    # Additional node for "classic" BM 
    # -> https://link.springer.com/article/10.1007/JHEP09(2018)057
    '8a'  : [1.,1.,0.5,0.8/3,0.],
    # NLO samples
    'cHHH0'    : [0.,1.,0.,0.,0.],
    'cHHH1'    : [1.,1.,0.,0.,0.],
    'cHHH2p45' : [2.45,1.,0.,0.,0.],
    'cHHH5'    : [5.,1.,0.,0.,0.],
    # https://arxiv.org/pdf/1908.08923.pdf
    # -> to be reparameterized !!!
    # https://gitlab.cern.ch/hh/eft-benchmarks
    'cluster1' : [3.94,0.94,-1./3,0.5*1.5,1./3.*(-3.)],
    'cluster2' : [6.84,0.61,1./3,0.*1.5,-1./3*(-3.)],
    'cluster3' : [2.21,1.05,-1./3,0.5*1.5,0.5*(-3.)],
    'cluster4' : [2.79,0.61,1./3,-0.5*1.5,1/6*(-3.)],
    'cluster5' : [3.95,1.17,-1./3,1./6*1.5,-0.5*(-3.)],
    'cluster6' : [5.68,0.83,1./3,-0.5*1.5,1./3*(-3.)],
    'cluster7' : [-0.10,0.94,1.,1./6*1.5,-1./6*(-3.)],
}


if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Produce LO weights')
    parser.add_argument('--path', action='store', required=True, type=str,
                        help='Path to directory with the root files')
    parser.add_argument('--era', action='store', required=True, type=str,
                        help='era')
    args = parser.parse_args()

    for f in sorted(glob.glob(os.path.join(args.path,'*root'))):
        if "__skeleton__" in f:
            continue
        sample = os.path.basename(f)
        print (sample)
        reweight = ReweightingNLO(f,"mHHvsCosThetaStar")

        for BM,coupl in couplings.items():
            if isinstance(coupl,int):
                reweight.SetBenchmark(coupl)
            elif isinstance(coupl,list):
                assert len(coupl) == 5
                reweight.SetCouplings(coupl)
            else:
                raise RuntimeError('Not understood couplings')
            reweight.SaveToJson("weights/{}_to_Benchmark{}_{}".format(sample.replace('.root',''),BM,args.era))
