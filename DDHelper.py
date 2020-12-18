from bamboo.analysismodules import DataDrivenContribution

class DataDrivenFake(DataDrivenContribution):
    def usesSample(self, sampleName, sampleConfig):
        return not ('type' in sampleConfig.keys() and sampleConfig["type"]=="signal") 
        # Not to be done on signal samples
    def replacesSample(self, sampleName, sampleConfig):
        return False # No sample to be replaced : new contribution
    def modifiedSampleConfig(self, sampleName, sampleConfig, lumi=None):
        # Fake contribution = reweighted data - reweighted MC (in Fake CR)
        modCfg = dict(sampleConfig)
        if 'type' in sampleConfig.keys() and sampleConfig["type"]=="signal":
            return {}
        if sampleConfig["group"]=="data":
            # Data : add as MC with correct normalization so that lumi*Xsec/generated-events = 1
            modCfg.update({"type": "mc",
                           "generated-events": lumi,
                           "cross-section": 1.,
                           "group": "Fake"})
        else:
            # MC : Change sign of contribution so that in the end Fake = Data-MC
            modCfg.update({"branching-ratio" : -1.,
                           "group": "Fake"})

        return modCfg

class DataDrivenDY(DataDrivenContribution):
    def __init__(self, name, config, pseudodata=False):
        self.pseudodata = pseudodata
        super(DataDrivenDY,self).__init__(name,config)
    def usesSample(self, sampleName, sampleConfig):
        if 'group' not in sampleConfig.keys():
            return False # Signal samples not to be used
        if not self.pseudodata:
            return sampleConfig['group'] == 'data' # Data driven 
        else:
            return sampleConfig['group'] != 'data' # Closure
    def replacesSample(self, sampleName, sampleConfig):
        if 'group' not in sampleConfig.keys():
            return False # Signal samples not to be used
        return sampleConfig['group'] == 'DY' # Replace DY by data evaluation
    def modifiedSampleConfig(self, sampleName, sampleConfig, lumi=None):
        modCfg = dict(sampleConfig)
        if 'type' in sampleConfig.keys() and sampleConfig["type"]=="signal":
            return {}
        else:
            if sampleConfig["group"] == "data" and not self.pseudodata:
                # Data : add as MC with correct normalization so that lumi*Xsec/generated-events = 1
                modCfg.update({"type": "mc",
                               "generated-events": lumi,
                               "cross-section": 1.,
                               "group": "DYEstimation"})
            if sampleConfig["group"] != "data" and self.pseudodata:
                # Closure 
                modCfg.update({"group": "DYEstimation"})

        return modCfg

class DataDrivenPseudoData(DataDrivenContribution):
    def usesSample(self, sampleName, sampleConfig):
        return 'group' in sampleConfig.keys() and sampleConfig['group'] != 'data' 
    def replacesSample(self, sampleName, sampleConfig):
        return 'group' in sampleConfig.keys() and sampleConfig['group'] == 'data' 
    def modifiedSampleConfig(self, sampleName, sampleConfig, lumi=None):
        modCfg = dict(sampleConfig)
        if 'type' in sampleConfig.keys() and sampleConfig["type"]=="signal":
            return {}
        if sampleConfig["group"] != "data": 
            modCfg.update({"group": "pseudodata",
                           "stack-index":1,
                           "type":'mc'})
        return modCfg
