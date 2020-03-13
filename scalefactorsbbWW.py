import os 
import sys
from bamboo import scalefactors
from bamboo import treefunctions as op

from scalefactorsutils import MakeScaleFactorsDict

class ScaleFactorsbbWW:
    def __init__(self):
        self.binningVariables = {
        "Eta"       : lambda obj : obj.p4.Eta(),
        "ClusEta"   : lambda obj : obj.p4.Eta() + obj.deltaEtaSC,
        "AbsEta"    : lambda obj : op.abs(obj.p4.Eta()),
        "AbsClusEta": lambda obj : op.abs(obj.clusterEta) +op.abs(obj.deltaEtaSC),
        "Pt"        : lambda obj : obj.p4.Pt(),
        }

         # Scale factors dictionary construction #
         # (better to do that here so that it is not buil everytime we ask for scalefactors #
        instance = MakeScaleFactorsDict(efficiencies_path  = os.path.join(os.path.dirname(os.path.abspath(__file__)),'data','ScaleFactors_ttH'),
                                        scalefactors_path  = os.path.join(os.path.dirname(os.path.abspath(__file__)),'data','ScaleFactors_ttH'),
                                        check_path         = True)
           # 2016 legacy:
           # Electron Reco efficiency
           #    https://gitlab.cern.ch/ttH_leptons/doc/blob/master/Legacy/data_to_mc_corrections.md#electron-id-efficiency-scale-factors-for-loose-lepton-id 

        # Electrons Loose #
        instance.addScaleFactor(key_entry       = 'electron_loosereco_2016',
                                base_key        = "electron_loosereco_{pt}",
                                base_str        = "Electron_EGamma_SF2D_RecoEff2016{pt}.json",
                                format_dict     = {'pt':['PtLt20','PtGt20']},
                                lowercase_keys  = True)
        instance.addEfficiency(key_entry   = 'electron_looseid_2016',
                               base_str    = "Electron_EGamma_SF2D_LooseMVA2016.json")
        instance.addEfficiency(key_entry   = 'electron_looseeff_2016',
                               base_str    = "Electron_EGamma_SF2D_LooseEff2016.json")

                                format_dict     = {'algo':["DeepCSV"],'wp':["loose", "medium"],('flav', 'calib'):[("lightjets", "incl"), ("cjets", "lt"), ("bjets","lt")]})
        # Muons Loose #
        instance.addEfficiency(key_entry   = 'muon_loose_2016',
                               base_str    = "Muon_EGamma_SF2D_Loose2016.json")

        # Electrons tight #
        instance.addEfficiency(key_entry   = 'electron_tightMVA_2016',
                               base_str    = "Electron_EGamma_SF2D_TightMVA2016.json")

        # Muons tight #
        instance.addEfficiency(key_entry   = 'muon_tightMVA_2016',
                               base_str    = "Muon_EGamma_SF2D_TightMVA2016.json")
        # Btagging #
        instance.addScaleFactor(key_entry       = 'btag_2016_94X',
                                base_key        = '{algo}_{wp}',
                                base_str        = "BTagging_{wp}_{flav}_{calib}_{algo}.json",
                                format_dict     = {'algo':["DeepCSV", "DeepJet"],'wp':["loose", "medium", "tight"],('flav', 'calib'):[("lightjets", "incl"), ("cjets", "comb"), ("bjets","comb")]})
        instance.addScaleFactor(key_entry       = 'subjet_btag_2016_94X',
                                base_key        = '{algo}_{wp}',
                                base_str        = "BTagging_{wp}_{flav}_{calib}_subjet_{algo}.json",
                                format_dict     = {'algo':["DeepCSV"],'wp':["loose", "medium"],('flav', 'calib'):[("lightjets", "incl"), ("cjets", "lt"), ("bjets","lt")]})

        # Get full dict #
        self.all_scalefactors = instance.getScaleFactorsDict()

        import pprint
        pprint.pprint(self.all_scalefactors)

    def get_scalefactor(self,objType, key, periods=None, combine=None, additionalVariables=dict(), systName=None):

        return scalefactors.get_scalefactor(objType             = objType, 
                                            key                 = key, 
                                            periods             = periods, 
                                            combine             = combine, 
                                            additionalVariables = additionalVariables, 
                                            sfLib               = self.all_scalefactors, 
                                            paramDefs           = self.binningVariables, 
                                            getFlavour          = (lambda j : j.hadronFlavour), 
                                            systName            = systName)
    
