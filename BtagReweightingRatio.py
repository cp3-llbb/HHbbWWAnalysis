import os
import sys
from copy import copy

from itertools import chain

import logging
logger = logging.getLogger(__name__) 

import bamboo
from bamboo.analysismodules import HistogramsModule, DataDrivenBackgroundHistogramsModule

from bamboo import treefunctions as op
from bamboo.plots import CutFlowReport, Plot, EquidistantBinning, SummedPlot

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)))) # Add scripts in this directory
from BaseHHtobbWW import BaseNanoHHtobbWW
from plotDef import *

#===============================================================================================#
#                                       PlotterHHtobbWW                                         #
#===============================================================================================#
class BtagReweightingRatioNano(BaseNanoHHtobbWW,DataDrivenBackgroundHistogramsModule):
    def __init__(self, args):
        super(BtagReweightingRatioNano, self).__init__(args)

    def initialize(self):
        super(BtagReweightingRatioNano, self).initialize(forSkimmer=True)

    def definePlots(self, t, noSel, sample=None, sampleCfg=None): 
        noSel = super(BtagReweightingRatioNano,self).prepareObjects(t, noSel, sample, sampleCfg, 'DL')

        plots = []

        era = sampleCfg['era']


        #----- Ratio reweighting variables (before lepton and jet selection) -----#
        if self.args.BtagReweightingOff or self.args.BtagReweightingOn:
            noSel = self.beforeJetselection(noSel)
            plots.append(objectsNumberPlot(channel="NoChannel",suffix='NoSelection',sel=noSel,objCont=self.ak4Jets,objName='Ak4Jets',Nmax=15,xTitle='N(Ak4 jets)'))
            plots.append(CutFlowReport("BtagReweightingCutFlowReport",noSel,printInLog=True,recursive=True))
        else:
            raise RuntimeError("Either --BtagReweightingOff or --BtagReweightingOn must be used")

        return plots
