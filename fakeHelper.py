from bamboo.analysismodules import DataDrivenContribution

class DataDrivenFake(DataDrivenContribution):
    def usesSample(self, sampleName, sampleConfig):
        return True # We want all samples to be used
    def replacesSample(self, sampleName, sampleConfig):
        return False # No sample to be replaced 
    def modifiedSampleConfig(self, sampleName, sampleConfig, lumi=None):
        # Fake contribution = reweighted data - reweighted MC (in Fake CR)
        modCfg = dict(sampleConfig)
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


