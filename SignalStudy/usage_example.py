# this script provides an example of how to reweight a certain input 
# sample to a given output benchmark. Two options are viable:
#
# 1. If you know the coefficients used to generate the input sample, and,
#    in case it is a LO Madgraph sample, you have the LO reweightin tool, the 
#    reweight can be derived as:
#    reweight = [XS(mHH,costh|outputBM)/XS(mHH,costh|inputBM)] * [XStot(inputBM)/XStot(outputBM)]         
#
# 2. If you know the 2D distribution in (mHH,costh) of your input sample Nev(mHH,costh), 
#    you can derive the reweight as:
#    reweight = [XS(mHH,costh|outputBM)/Nev(mHH,costh)] * [Nevtot/XStot(outputBM)]
#    NOTE: the binning for Nev(mHH,costh) has to be the same defined in NonResonantModelNLO        
#
# The tools offer also the possibility to reweight only in mHH, 
# neglecting the dependence on costhetaHH.
# 
 
import  NonResonantModelNLO
import ROOT 
mymodel = NonResonantModelNLO.NonResonantModelNLO()
mymodel.ReadCoefficients("pm_pw_NLO_Ais_13TeV_V2.txt")

# event example
event_mHH=516
event_costhetaHH=0.5
event_weight=0.1234

# assume that I want to reweight a given input sample for benchmark XYZ
# to all the NLO benchmarks using the approach 2 as described above
inputevfile = ROOT.TFile("/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SignalLOHistograms_2016/results/GluGluToHHTo2B2VTo2L2Nu_node_10.root")
histo_Nev = inputevfile.Get("mHHvsCosThetaStar") # this is a TH2 histo
Nevtot = histo_Nev.Integral() 

for iBM in range(0,12):
    BMcouplings = mymodel.getBenchmark(iBM)
    print (iBM)
    kl, kt, c2, cg, c2g = BMcouplings[0], BMcouplings[1], BMcouplings[2], BMcouplings[3], BMcouplings[4]
    print ('couplings',BMcouplings)
    print ("totalXS", mymodel.getTotalXS(kl, kt, c2, cg, c2g))
    print ("diffXS(event_mHH,event_costhetaHH)", mymodel.getDifferentialXS2D(event_mHH, event_costhetaHH, kl, kt, c2, cg, c2g))
    print ("diffXS(event_mHH)", mymodel.getDifferentialXSmHH(event_mHH, kl, kt, c2, cg, c2g))
    Nev = histo_Nev.GetBinContent( histo_Nev.FindBin(event_mHH, event_costhetaHH) )
    print (event_mHH, event_costhetaHH, histo_Nev.FindBin(event_mHH, event_costhetaHH))
    XS = mymodel.getDifferentialXS2D(event_mHH, event_costhetaHH, kl, kt, c2, cg, c2g)
    # before using the differential XS, scale it by the bin area to properly compare it with histo_Nev content 
    Noutputev = XS * mymodel.getmHHbinwidth(event_mHH) * mymodel.getcosthetabinwidth(event_costhetaHH) 
    XStot = mymodel.getTotalXS(kl, kt, c2, cg, c2g)
    reweight = event_weight * Noutputev/Nev * Nevtot/XStot
    print (Noutputev/Nev * Nevtot/XStot)
    print ()

