import os
import sys
import yaml
from copy import deepcopy
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

with open('/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/TestLO_to_NLO/plots.yml','r') as handle:
    config = yaml.load(handle)
    

path_NLO = {
    '2b2l2nu cHHH0'     : '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/TestLO_to_NLO/results/GluGluToHHTo2B2VTo2L2Nu_node_cHHH0.root',
    '2b2tau cHHH0'      : '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/TestLO_to_NLO/results/GluGluToHHTo2B2Tau_node_cHHH0.root',
    '2b2l2nu cHHH1'     : '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/TestLO_to_NLO/results/GluGluToHHTo2B2VTo2L2Nu_node_cHHH1.root',
    '2b2tau cHHH1'      : '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/TestLO_to_NLO/results/GluGluToHHTo2B2Tau_node_cHHH1.root',
    '2b2l2nu cHHH2p45'  : '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/TestLO_to_NLO/results/GluGluToHHTo2B2VTo2L2Nu_node_cHHH2p45.root',
    '2b2tau cHHH2p45'   : '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/TestLO_to_NLO/results/GluGluToHHTo2B2Tau_node_cHHH2p45.root',
    '2b2l2nu cHHH5'     : '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/TestLO_to_NLO/results/GluGluToHHTo2B2VTo2L2Nu_node_cHHH5.root',
    '2b2tau cHHH5'      : '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/TestLO_to_NLO/results/GluGluToHHTo2B2Tau_node_cHHH5.root',
}
path_rew = {
    '2b2l2nu cHHH0'     : '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/TestLO_to_NLO/results/GluGluToHHTo2B2VTo2L2Nu_NLOBenchmarkcHHH0.root', 
    '2b2tau cHHH0'      : '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/TestLO_to_NLO/results/GluGluToHHTo2B2Tau_NLOBenchmarkcHHH0.root',
    '2b2l2nu cHHH1'     : '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/TestLO_to_NLO/results/GluGluToHHTo2B2VTo2L2Nu_NLOBenchmarkcHHH1.root', 
    '2b2tau cHHH1'      : '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/TestLO_to_NLO/results/GluGluToHHTo2B2Tau_NLOBenchmarkcHHH1.root',
    '2b2l2nu cHHH2p45'  : '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/TestLO_to_NLO/results/GluGluToHHTo2B2VTo2L2Nu_NLOBenchmarkcHHH2p45.root', 
    '2b2tau cHHH2p45'   : '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/TestLO_to_NLO/results/GluGluToHHTo2B2Tau_NLOBenchmarkcHHH2p45.root',
    '2b2l2nu cHHH5'     : '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/TestLO_to_NLO/results/GluGluToHHTo2B2VTo2L2Nu_NLOBenchmarkcHHH5.root', 
    '2b2tau cHHH5'      : '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/TestLO_to_NLO/results/GluGluToHHTo2B2Tau_NLOBenchmarkcHHH5.root',
}


histName = "ElEl_Has2FakeableElElOSWithTriggersPtCutsPreMllCutOutZTightSelectedTwoAk4JetsExclusiveResolvedOneBtag_combined_mHH"

h_NLO = {}
h_rew = {}

lumi = 59740.565201546

for key in path_NLO.keys():
    print ('-------')
    print (key)
    F_NLO = ROOT.TFile(path_NLO[key],'READ')
    F_rew = ROOT.TFile(path_rew[key],'READ')
    param_NLO = config['files'][os.path.basename(path_NLO[key])]
    param_rew = config['files'][os.path.basename(path_rew[key])]
    xsec = 1.

    h_NLO[key] = deepcopy(F_NLO.Get(histName))
    h_rew[key] = deepcopy(F_rew.Get(histName))

    print (h_NLO[key].Integral(),h_rew[key].Integral())
    h_NLO[key].Scale(lumi*xsec*param_NLO['branching-ratio']/param_NLO['generated-events'])
    h_rew[key].Scale(lumi*xsec*param_rew['branching-ratio']/param_rew['generated-events'])
    print (param_NLO['branching-ratio'],param_rew['branching-ratio'])
    print (param_NLO['generated-events'],param_rew['generated-events'])
    print (h_NLO[key].Integral(),h_rew[key].Integral())
    
    F_NLO.Close()
    F_rew.Close()

C = ROOT.TCanvas()
C.Print("Comparisonv2.pdf[")

for key in h_NLO.keys():
    C.Clear()

    legend = ROOT.TLegend(0.5,0.7,0.88,0.88)
    legend.SetFillStyle(4000)
    legend.SetBorderSize(0)

    h_NLO[key].SetLineWidth(2)
    h_NLO[key].SetLineStyle(1)
    h_NLO[key].SetLineColor(601)
    h_NLO[key].GetXaxis().SetTitle("Reco M_{HH} (e^{\pm}e^{\mp} channel) [GeV]")
    h_NLO[key].GetYaxis().SetTitle("Yield")
    h_NLO[key].SetMaximum(h_NLO[key].GetMaximum()*1.1)
    h_NLO[key].SetMinimum(0.)

    h_rew[key].SetLineWidth(2)
    h_rew[key].SetLineStyle(1)
    h_rew[key].SetLineColor(418)
    
    h_NLO[key].Draw("H")
    h_rew[key].Draw("H same")
    legend.AddEntry(h_NLO[key],key+" (NLO)")
    legend.AddEntry(h_rew[key],key+" (LO reweighted)")

    legend.Draw()
    C.Print("Comparisonv2.pdf")

C.Print("Comparisonv2.pdf]")



