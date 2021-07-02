import os
from bamboo.scalefactors import get_scalefactor, binningVariables_nano
from bamboo.treedecorators import NanoAODDescription
from bamboo.analysismodules import NanoAODModule, NanoAODHistoModule, NanoAODSkimmerModule, DataDrivenBackgroundHistogramsModule
from bamboo.plots import SelectionWithDataDriven
from bamboo.plots import Plot, EquidistantBinning, VariableBinning
from bamboo import treefunctions as op

class SignalReweightingHistogramNano(NanoAODHistoModule,DataDrivenBackgroundHistogramsModule):
    """ Base module for NanoAOD mHH vs cos(theta*) histograms """
    def __init__(self, args):
        super(SignalReweightingHistogramNano, self).__init__(args)

    def addArgs(self,parser):
        super(SignalReweightingHistogramNano, self).addArgs(parser)
        parser.add_argument("-s",
                            "--subset", 
                            type        = str,
                            required    = False,
                            help="Subset of samples to be run over, keys defined in the 'samples' entry of the yaml analysis config")
        parser.add_argument("-r",
                            "--reweighting", 
                            action      = 'store_true',
                            required    = False,
                            default     = False,
                            help="Whether to apply reweighting")


    def initialize(self):
        self._customizeAnalysisCfg(self.analysisConfig)
        super(SignalReweightingHistogramNano, self).initialize()

    def _customizeAnalysisCfg(self,analysisCfg):
        import yaml
        samples = {} 
        if self.args.subset is None:
            return
        reqArgs = self.args.subset.split(',')
        foundArgs = set()
        subsets = []
        for item in analysisCfg['samples']:
            if not 'keys' in item.keys() or not 'config' in item.keys():
                continue
            keys = item['keys']
            if not isinstance(keys,list):
                keys = [keys]
            if all([key in reqArgs for key in keys]) or 'all' in reqArgs:
                foundArgs.update(keys)
                subsets.append(item['config'])
                with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'Yaml',item['config'])) as handle:
                    samples.update(yaml.load(handle,yaml.SafeLoader))
        self.analysisConfig['samples'] = samples
        notFoundArgs = [arg for arg in reqArgs if arg not in foundArgs]
        if len(notFoundArgs)>0:
            raise RuntimeError('The following subsets have not been found in the keys of the analysis yaml file : '+', '.join(notFoundArgs))
        if len(subsets) > 0:
            print ("Imported following yaml subsets :")
            for subset in subsets:
                print ('... {}'.format(subset))


    def prepareTree(self, tree, sample=None, sampleCfg=None):
        era = sampleCfg.get("era") if sampleCfg else None
        isMC = self.isMC(sample)
        tree,noSel,be,lumiArgs = super(SignalReweightingHistogramNano, self).prepareTree(tree, sample=sample, sampleCfg=sampleCfg,description=NanoAODDescription.get("v7", year=era, isMC=isMC))
        return tree,noSel,be,lumiArgs

    def definePlots(self, t, noSel, sample=None, sampleCfg=None): 
        if 'type' not in sampleCfg.keys() or sampleCfg["type"] != "signal": 
            raise RuntimeError("Sample needs to be HH signal LO GGF sample")

        era = sampleCfg.get("era") if sampleCfg else None

        # Select gen level Higgs #
        genh = op.select(t.GenPart,lambda g : op.AND(g.pdgId==25, g.statusFlags & ( 0x1 << 13)))
        HH_p4 = genh[0].p4 + genh[1].p4 
        cm = HH_p4.BoostToCM() 
        boosted_h = op.extMethod("ROOT::Math::VectorUtil::boost", returnType=genh[0].p4._typeName)(genh[0].p4,cm)
        mHH = op.invariant_mass(genh[0].p4,genh[1].p4) 
        cosHH = op.abs(boosted_h.Pz()/boosted_h.P())

        # Apply reweighting #

        benchmarks = [
            'BenchmarkSM',  
            'Benchmark1',  
            'Benchmark2',  
            'Benchmark3',  
            'Benchmark4',  
            'Benchmark5',  
            'Benchmark6',  
            'Benchmark7',  
            'Benchmark8',  
            'Benchmark8a',  
            'Benchmark9',  
            'Benchmark10',  
            'Benchmark11',  
            'Benchmark12',  
            'BenchmarkcHHH0',  
            'BenchmarkcHHH1',  
            'BenchmarkcHHH2p45',  
            'BenchmarkcHHH5',  
            'Benchmarkcluster1',  
            'Benchmarkcluster2',  
            'Benchmarkcluster3',  
            'Benchmarkcluster4',  
            'Benchmarkcluster5',  
            'Benchmarkcluster6',  
            'Benchmarkcluster7',  
        ]
        selections = {'':noSel}
        reweights = {}
        if self.args.reweighting:
            for benchmark in benchmarks:
                json_file = os.path.join(os.path.abspath(os.path.dirname(__file__)),'data','ScaleFactors_GGF_LO','{}_to_{}_{}.json'.format(sample,benchmark,era))
                if os.path.exists(json_file):
                    print("Found file {}".format(json_file))
                    reweightLO = get_scalefactor("lepton", json_file, paramDefs={'Eta': lambda x : mHH, 'Pt': lambda x : cosHH})
                    selections[benchmark] = SelectionWithDataDriven.create(
                                                           parent   = noSel,
                                                           name     = 'noSel'+benchmark,
                                                           ddSuffix = benchmark,
                                                           cut      = op.c_bool(True),
                                                           ddCut    = op.c_bool(True),
                                                           weight   = op.c_float(1.),
                                                           ddWeight = reweightLO(op.c_float(1.)),
                                                           enable   = True)
                    reweights[benchmark] = reweightLO(op.c_float(1.))
                else:
                    print("Could not find file {}".format(json_file))

        # Plots #
        plots = []

        for name,reweight in reweights.items():
            plots.append(Plot.make1D("weight_{}".format(name),
                                     reweight,
                                     noSel,
                                     EquidistantBinning(100,0,5.),
                                     xTitle = 'weight'))
                    
        for selName,sel in selections.items():
            plots.append(Plot.make2D(f"mHHvsCosThetaStar{selName}",
                                     [mHH,cosHH],
                                     sel,
                                     [VariableBinning([250.,270.,290.,310.,330.,
                                                       350.,370.,390.,410.,430.,
                                                       450.,470.,490.,510.,530.,
                                                       550.,570.,590.,610.,630.,
                                                       650.,670.,700.,750.,800.,
                                                       850.,900.,950.,1000.,1100.,
                                                       1200.,1300.,1400.,1500.,1750.,2000.,5000.]),
                                      VariableBinning([ 0.0, 0.4, 0.6, 0.8, 1.0 ])],
                                     xTitle = 'm_{HH}',
                                     yTitle = 'cos(#theta^{*})'))
            plots.append(Plot.make1D(f"mHH{selName}",
                                     mHH,
                                     sel,
                                     VariableBinning([250.,270.,290.,310.,330.,
                                                      350.,370.,390.,410.,430.,
                                                      450.,470.,490.,510.,530.,
                                                      550.,570.,590.,610.,630.,
                                                      650.,670.,700.,750.,800.,
                                                      850.,900.,950.,1000.,1100.,
                                                      1200.,1300.,1400.,1500.,1750.,2000.,5000.]),
                                     xTitle = 'm_{HH}'))
            plots.append(Plot.make1D(f"cosThetaStar{selName}",
                                     cosHH,
                                     sel,
                                     VariableBinning([ 0.0, 0.4, 0.6, 0.8, 1.0 ]),
                                     xTitle = 'cos(#theta^{*})'))
                    
        return plots


