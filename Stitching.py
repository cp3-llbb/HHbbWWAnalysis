import os
import sys
from copy import copy

from itertools import chain

import logging
logger = logging.getLogger(__name__) 

import bamboo
from bamboo.analysismodules import HistogramsModule

from bamboo import treefunctions as op
from bamboo.plots import CutFlowReport, Plot, EquidistantBinning, SummedPlot

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)))) # Add scripts in this directory
from BaseHHtobbWW import BaseNanoHHtobbWW
from plotDef import *

#===============================================================================================#
#                                       PlotterHHtobbWW                                         #
#===============================================================================================#
class StitchingNano(BaseNanoHHtobbWW):
    def __init__(self, args):
        super(StitchingNano, self).__init__(args)

    def definePlots(self, t, noSel, sample=None, sampleCfg=None): 
        noSel = super(StitchingNano,self).prepareObjects(t, noSel, sample, sampleCfg, 'DL')
        #---- Parameters -----#

        plots = []

        era = sampleCfg['era']

        if self.args.DYStitchingPlots and sampleCfg['group'] != 'DY':
            raise RuntimeError("Stitching is only done on DY MC samples")
        if self.args.WJetsStitchingPlots and sampleCfg['group'] != 'Wjets':
            raise RuntimeError("Stitching is only done on WJets MC samples")
        plots.extend(makeLHEPlots(noSel,t.LHE))
        plots.append(objectsNumberPlot(channel="NoChannel",suffix='NoSelection',sel=noSel,objCont=self.ak4Jets,objName='Ak4Jets',Nmax=15,xTitle='N(Ak4 jets)'))
        plots.append(CutFlowReport("DYStitchingCutFlowReport",noSel,printInLog=True,recursive=True))
        return plots

