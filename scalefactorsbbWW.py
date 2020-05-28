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
        instance = MakeScaleFactorsDict(paths       = {'ttH_SF' : os.path.join(os.path.dirname(os.path.abspath(__file__)),'data','ScaleFactors_ttH'),
                                                       'DY_SF'  : os.path.join(os.path.dirname(os.path.abspath(__file__)),'data','ScaleFactors_DY'),
                                                       'POG_SF' : os.path.join(os.path.dirname(os.path.abspath(__file__)),'data','ScaleFactors_POG')},
                                        check_path  = True)
        #----- 2016 -----#
        # Electrons Loose #
           # https://gitlab.cern.ch/ttH_leptons/doc/blob/master/Legacy/data_to_mc_corrections.md#electron-id-efficiency-scale-factors-for-loose-lepton-id 
        instance.AddScaleFactorWithWorkingPoint(path_key        = 'ttH_SF',
                                                entry_key       = 'electron_loosereco_2016',
                                                base_key        = "electron_loosereco_{pt}",
                                                base_str        = "Electron_EGamma_SF2D_RecoEff2016{pt}.json",
                                                format_dict     = {'pt':['PtLt20','PtGt20']},
                                                lowercase_keys  = True)
        instance.AddScaleFactor(path_key    = 'ttH_SF',
                                entry_key   = 'electron_looseid_2016',
                                base_str    = "Electron_EGamma_SF2D_LooseMVA2016.json")
        instance.AddScaleFactor(path_key    = 'ttH_SF',
                                entry_key   = 'electron_looseeff_2016',
                                base_str    = "Electron_EGamma_SF2D_LooseEff2016.json")
        # Muons Loose #
            # https://gitlab.cern.ch/ttH_leptons/doc/-/blob/master/Legacy/data_to_mc_corrections.md#muon-id-efficiency-scale-factors-for-loose-lepton-id
        instance.AddScaleFactor(path_key    = 'ttH_SF',
                                entry_key   = 'muon_loose_2016',
                                base_str    = "Muon_EGamma_SF2D_Loose2016.json")

        # Electrons tight #
            # https://gitlab.cern.ch/ttH_leptons/doc/-/blob/master/Legacy/data_to_mc_corrections.md#lepton-id-efficiency-scale-factors-for-tight-lepton-id
        instance.AddScaleFactor(path_key    = 'ttH_SF',
                                entry_key   = 'electron_tightMVA_2016',
                                base_str    = "TTHSF_EGamma_SF2D_ElectronTight2016.json")

        # Muons tight #
            # https://gitlab.cern.ch/ttH_leptons/doc/-/blob/master/Legacy/data_to_mc_corrections.md#lepton-id-efficiency-scale-factors-for-tight-lepton-id
        instance.AddScaleFactor(path_key    = 'ttH_SF',
                                entry_key   = 'muon_tightMVA_2016',
                                base_str    = "TTHSF_EGamma_SF2D_MuonTight2016.json")
        # Btagging #
            # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation2016Legacy#Supported_Algorithms_and_Operati
        instance.AddScaleFactorWithWorkingPoint(path_key    = 'POG_SF',
                                                entry_key   = 'btag_2016',
                                                base_key    = '{algo}_{wp}',
                                                base_str    = "BTagging_{wp}_{flav}_{calib}_{algo}_2016.json",
                                                format_dict = {'algo':["DeepCSV", "DeepJet"],'wp':["loose", "medium", "tight"],('flav', 'calib'):[("lightjets", "incl"), ("cjets", "comb"), ("bjets","comb")]})
            # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation2016Legacy#Subjet_b_tagging
        instance.AddScaleFactorWithWorkingPoint(path_key    = 'POG_SF',
                                                entry_key   = 'subjet_btag_2016',
                                                base_key    = '{algo}_{wp}',
                                                base_str    = "BTagging_{wp}_{flav}_{calib}_subjet_{algo}_2016.json",
                                                format_dict = {'algo':["DeepCSV"],'wp':["loose", "medium"],('flav', 'calib'):[("lightjets", "incl"), ("cjets", "lt"), ("bjets","lt")]})

        #-----  DY weight for 1 and 2 btag -----#
        instance.AddScaleFactorWithWorkingPoint(path_key    = 'DY_SF',
                                                entry_key   = 'DY_2016',
                                                base_key    = '{channel}_{type}_{btag}',
                                                base_str    = 'weight_{channel}_{type}_weight_{btag}_2016.json',
                                                format_dict = {'channel':['ElEl','MuMu'],'type':['data','mc'],'btag':['1b','2b']})

        #----- 2017 -----#
        # Check links of 2016 #
        # Electrons Loose #
        instance.AddScaleFactorWithWorkingPoint(path_key        = 'ttH_SF',
                                                entry_key       = 'electron_loosereco_2017',
                                                base_key        = "electron_loosereco_{pt}",
                                                base_str        = "Electron_EGamma_SF2D_RecoEff2017{pt}.json",
                                                format_dict     = {'pt':['PtLt20','PtGt20']},
                                                lowercase_keys  = True)
        instance.AddScaleFactor(path_key    = 'ttH_SF',
                                entry_key   = 'electron_looseid_2017',
                                base_str    = "Electron_EGamma_SF2D_LooseMVA2017.json")
        instance.AddScaleFactor(path_key    = 'ttH_SF',
                                entry_key   = 'electron_looseeff_2017',
                                base_str    = "Electron_EGamma_SF2D_LooseEff2017.json")
        # Muons Loose #
        instance.AddScaleFactor(path_key    = 'ttH_SF',
                                entry_key   = 'muon_loose_2017',
                                base_str    = "Muon_EGamma_SF2D_Loose2017.json")

        # Electrons tight #
        instance.AddScaleFactor(path_key    = 'ttH_SF',
                                entry_key   = 'electron_tightMVA_2017',
                                base_str    = "TTHSF_EGamma_SF2D_ElectronTight2017.json")

        # Muons tight #
        instance.AddScaleFactor(path_key    = 'ttH_SF',
                                entry_key   = 'muon_tightMVA_2017',
                                base_str    = "TTHSF_EGamma_SF2D_MuonTight2017.json")
        # Btagging #
            # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X#Supported_Algorithms_and_Operati
        instance.AddScaleFactorWithWorkingPoint(path_key    = 'POG_SF',
                                                entry_key   = 'btag_2017',
                                                base_key    = '{algo}_{wp}',
                                                base_str    = "BTagging_{wp}_{flav}_{calib}_{algo}_2017.json",
                                                format_dict = {'algo':["DeepCSV", "DeepJet"],'wp':["loose", "medium", "tight"],('flav', 'calib'):[("lightjets", "incl"), ("cjets", "comb"), ("bjets","comb")]})
            # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X#Boosted_event_topologies
        instance.AddScaleFactorWithWorkingPoint(path_key    = 'POG_SF',
                                                entry_key   = 'subjet_btag_2017',
                                                base_key    = '{algo}_{wp}',
                                                base_str    = "BTagging_{wp}_{flav}_{calib}_subjet_{algo}_2017.json",
                                                format_dict = {'algo':["DeepCSV"],'wp':["loose", "medium"],('flav', 'calib'):[("lightjets", "incl"), ("cjets", "lt"), ("bjets","lt")]})


        #----- 2018 -----#
        # Check links of 2016 #
        # Electrons Loose #
        instance.AddScaleFactor(path_key    = 'ttH_SF',
                                entry_key   = 'electron_loosereco_2018',
                                base_str    = "Electron_EGamma_SF2D_RecoEff2018.json")
        instance.AddScaleFactor(path_key    = 'ttH_SF',
                                entry_key   = 'electron_looseid_2018',
                                base_str    = "Electron_EGamma_SF2D_LooseMVA2018.json")
        instance.AddScaleFactor(path_key    = 'ttH_SF',
                                entry_key   = 'electron_looseeff_2018',
                                base_str    = "Electron_EGamma_SF2D_LooseEff2018.json")
        # Muons Loose #
        instance.AddScaleFactor(path_key    = 'ttH_SF',
                                entry_key   = 'muon_loose_2018',
                                base_str    = "Muon_EGamma_SF2D_Loose2018.json")

        # Electrons tight #
        instance.AddScaleFactor(path_key    = 'ttH_SF',
                                entry_key   = 'electron_tightMVA_2018',
                                base_str    = "TTHSF_EGamma_SF2D_ElectronTight2018.json")

        # Muons tight #
        instance.AddScaleFactor(path_key    = 'ttH_SF',
                                entry_key   = 'muon_tightMVA_2018',
                                base_str    = "TTHSF_EGamma_SF2D_MuonTight2018.json")
        # Btagging #
            # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
        instance.AddScaleFactorWithWorkingPoint(path_key    = 'POG_SF',
                                                entry_key   = 'btag_2018',
                                                base_key    = '{algo}_{wp}',
                                                base_str    = "BTagging_{wp}_{flav}_{calib}_{algo}_2018.json",
                                                format_dict = {'algo':["DeepCSV", "DeepJet"],'wp':["loose", "medium", "tight"],('flav', 'calib'):[("lightjets", "incl"), ("cjets", "comb"), ("bjets","comb")]})
            # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X#Boosted_event_topologies
        instance.AddScaleFactorWithWorkingPoint(path_key    = 'POG_SF',
                                                entry_key   = 'subjet_btag_2018',
                                                base_key    = '{algo}_{wp}',
                                                base_str    = "BTagging_{wp}_{flav}_{calib}_subjet_{algo}_2018.json",
                                                format_dict = {'algo':["DeepCSV"],'wp':["loose", "medium"],('flav', 'calib'):[("lightjets", "incl"), ("cjets", "lt"), ("bjets","lt")]})




        # Get full dict #
        self.all_scalefactors = instance.GetScaleFactorsDict()

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
    
if __name__ == "__main__":
    instance = ScaleFactorsbbWW()
    import pprint
    pprint.pprint(instance.all_scalefactors)

