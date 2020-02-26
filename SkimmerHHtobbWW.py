import os
import sys

from bamboo import treefunctions as op

sys.path.append('/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/') # Add scripts in this directory -- TODO : make cleaner
from BaseHHtobbWW import BaseNanoHHtobbWW

#===============================================================================================#
#                                 PlotterHHtobbWW                                               #
#===============================================================================================#
class SkimmerNanoHHtobbWW(BaseNanoHHtobbWW):
    """ Plotter module: HH->bbW(->e/µ nu)W(->e/µ nu) histograms from NanoAOD """
    def __init__(self, args):
        super(SkimmerNanoHHtobbWW, self).__init__(args)

    def prepareTree(self, tree, sample=None, sampleCfg=None):
        return super(SkimmerNanoHHtobbWW,self).prepareTree(tree, sample, sampleCfg)

    def defineSkimSelection(self, t, noSel, sample=None, sampleCfg=None): 
        noSel = super(SkimmerNanoHHtobbWW,self).prepareObjects(t, noSel, sample, sampleCfg)

        # Initialize varsToKeep dict #
        varsToKeep = dict()  

        # Event variables #
        varsToKeep["event"] = t.event
        varsToKeep["run"] = t.run
        varsToKeep["ls"] = t.luminosityBlock


#         varsToKeep = {"nMuon": None, "Muon_eta": None, "Muon_pt": None} ## from input file         
#         varsToKeep["nSelMuons"] = op.static_cast("UInt_t", op.rng_len(muons)) ## TBranch doesn't accept size_t         
#         varsToKeep["selMu_miniPFRelIsoNeu"] = op.map(muons, lambda mu : mu.miniPFRelIso_all - mu.miniPFRelIso_chg)


        return noSel, varsToKeep

