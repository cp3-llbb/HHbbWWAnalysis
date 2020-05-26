import os 
import sys
#from bamboo import scalefactors
#from bamboo import treefunctions as op

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
        instance = MakeScaleFactorsDict(efficiencies_path  = os.path.join(os.path.dirname(os.path.abspath(__file__)),'data','TriggerEfficienciesStudies'),
                                        scalefactors_path  = os.path.join(os.path.dirname(os.path.abspath(__file__)),'data','ScaleFactors_FullRunII'),
                                        check_path         = True)
           # 2016 legacy:
           # https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaRunIIRecommendations#Fall17v2
           # https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonReferenceEffs2016LegacyRereco#Efficiencies
           # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation2016Legacy
        #-------- scale factors ----------#
        instance.addScaleFactor(key_entry       = 'electron_2016_94X',
                                base_key        = "id_{wp}",
                                base_str        = "Electron_EGamma_SF2D_2016Legacy_{wp}_Fall17V2.json",
                                format_dict     = {"wp":["Loose", "Medium", "Tight", "MVA80","MVA90", "MVA80noiso", "MVA90noiso"]},
                                lowercase_keys  = True)
        instance.addScaleFactor(key_entry       = 'muon_2016_94X',
                                base_key        = 'id_{wp}',
                                base_str        = "Muon_NUM_{wp}ID_DEN_genTracks_eta_pt_{uncer}_2016Run{era}.json",
                                format_dict     = {'wp':["Loose", "Medium", "Tight"],'uncer':['syst','stat'],'era':["BCDEF", "GH"]},
                                sublevel        = {'era':{'BCDEF':('Run2016B','Run2016C','Run2016D','Run2016E','Run2016F'),'GH':('Run2016G', 'Run2016H')}},
                                lowercase_keys  = True)
        instance.addScaleFactor(key_entry       = 'muon_2016_94X',
                                base_key        = 'iso_{isowp}_id_{idwp}{ip}',
                                base_str        = "Muon_NUM_{isowp}RelIso_DEN_{idwp}ID{ip}_eta_pt_{uncer}_2016Run{era}.json",
                                format_dict     = {('isowp','idwp','ip'):[("Loose", "Loose", ""), ("Loose", "Medium", ""), ("Loose", "Tight", "andIPCut"),("Tight", "Medium", ""), ("Tight", "Tight", "andIPCut")],'uncer':['syst','stat'],'era':["BCDEF", "GH"]},
                                sublevel        = {'era':{'BCDEF':('Run2016B','Run2016C','Run2016D','Run2016E','Run2016F'),'GH':('Run2016G', 'Run2016H')}},
                                lowercase_keys  = True)
        instance.addScaleFactor(key_entry       = 'muon_2016_94X',
                                base_key        = 'idtrk_{wp}',
                                base_str        = "Muon_NUM_{wp}ID_DEN_genTracks_eta_pair_newTuneP_probe_pt_{uncer}_2016Run{era}.json",
                                format_dict     = {'wp':['HighPt'],'era':["BCDEF", "GH"],'uncer':["syst", "stat"]},
                                sublevel        = {'era':{'BCDEF':('Run2016B','Run2016C','Run2016D','Run2016E','Run2016F'),'GH':('Run2016G', 'Run2016H')}},
                                lowercase_keys  = True)
        instance.addScaleFactor(key_entry       = 'muon_2016_94X',
                                base_key        = 'isotrk_{isowp}_idtrk_{idwp}',
                                base_str        = "Muon_NUM_{isowp}RelTkIso_DEN_{idwp}_eta_pair_newTuneP_probe_pt_{uncer}_2016Run{era}.json",
                                format_dict     = {('isowp','idwp'):[("Loose", "HighPtIDandIPCut")],'uncer':['syst','stat'],'era':["BCDEF", "GH"]},
                                sublevel        = {'era':{'BCDEF':('Run2016B','Run2016C','Run2016D','Run2016E','Run2016F'),'GH':('Run2016G', 'Run2016H')}},
                                lowercase_keys  = True)
        instance.addScaleFactor(key_entry       = 'btag_2016_94X',
                                base_key        = '{algo}_{wp}',
                                base_str        = "BTagging_{wp}_{flav}_{calib}_{algo}.json",
                                format_dict     = {'algo':["DeepCSV", "DeepJet"],'wp':["loose", "medium", "tight"],('flav', 'calib'):[("lightjets", "incl"), ("cjets", "comb"), ("bjets","comb")]})
        instance.addScaleFactor(key_entry       = 'subjet_btag_2016_94X',
                                base_key        = '{algo}_{wp}',
                                base_str        = "BTagging_{wp}_{flav}_{calib}_subjet_{algo}.json",
                                format_dict     = {'algo':["DeepCSV"],'wp':["loose", "medium"],('flav', 'calib'):[("lightjets", "incl"), ("cjets", "lt"), ("bjets","lt")]})
 
        #------- single muon trigger ----------#
        instance.addEfficiency(key_entry   = 'mutrig_2016_94X',
                               base_str    = '{trig}_PtEtaBins_2016Run{eras}.json',
                               format_dict = {'trig':["IsoMu24_OR_IsoTkMu24","Mu50_OR_TkMu50"],'eras':["BtoF", "GtoH"]})
        
        #------- dilepton trigger ----------#
        instance.addEfficiency(key_entry   = 'doubleEleLeg_HHMoriond17_2016',
                               base_str    = "{wp}.json",
                               format_dict = {'wp':["Electron_IsoEle23Leg", "Electron_IsoEle12Leg", "Electron_IsoEle23Leg", "Electron_IsoEle12Leg"]})
        instance.addEfficiency(key_entry   = 'doubleMuLeg_HHMoriond17_2016',
                               base_str    = "{wp}.json",
                               format_dict = {'wp':["Muon_DoubleIsoMu17Mu8_IsoMu17leg", "Muon_DoubleIsoMu17TkMu8_IsoMu8legORTkMu8leg", "Muon_DoubleIsoMu17Mu8_IsoMu17leg", "Muon_DoubleIsoMu17TkMu8_IsoMu8legORTkMu8leg"]})
        instance.addEfficiency(key_entry   = 'mueleLeg_HHMoriond17_2016',
                               base_str    = "{wp}.json",
                               format_dict = {'wp':["Muon_XPathIsoMu23leg", "Muon_XPathIsoMu8leg", "Electron_IsoEle23Leg", "Electron_IsoEle12Leg"]})
        instance.addEfficiency(key_entry   = 'elemuLeg_HHMoriond17_2016',
                               base_str    = "{wp}.json",
                               format_dict = {'wp':["Electron_IsoEle23Leg", "Electron_IsoEle12Leg", "Muon_XPathIsoMu23leg", "Muon_XPathIsoMu8leg"]})


            # 2017: 
            # https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffs2017
            # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
        #-------- scale factors ----------#
        #instance.addScaleFactor(key_entry       = 'electron_2017_94X',
        #                        base_key        = "id_{wp}",
        #                        base_str        = "Electron_EGamma_SF2D_2017{wp}.json",
        #                        format_dict     = {"wp":["Loose", "Medium", "Tight"]},
        #                        lowercase_keys  = True)
        #instance.addScaleFactor(key_entry       = 'muon_2017_94X',
        #                        base_key        = 'id_{wp}',
        #                        base_str        = "Muon_NUM_{wp}ID_DEN_genTracks_pt_abseta_2017RunBCDEF.json",
        #                        format_dict     = {'wp':["Loose", "Medium", "Tight","Soft", "MediumPrompt"]},
        #                        lowercase_keys  = True)
        #instance.addScaleFactor(key_entry       = 'muon_2017_94X',
        #                        base_key        = 'iso_{isowp}_id_{idwp}',
        #                        base_str        = "Muon_NUM_{isowp}RelIso_DEN_{idwp}ID_pt_abseta_2017RunBCDEF.json",
        #                        format_dict     = {('isowp','idwp'):[("Loose", "Loose"), ("Loose", "Medium"), ("Loose", "TightIDandIPCut"),  ("Tight", "Medium"), ("Tight", "TightIDandIPCut")]},
        #                        lowercase_keys  = True)
        #instance.addScaleFactor(key_entry       = 'btag_2017_94X',
        #                        base_key        = '{algo}_{wp}',
        #                        base_str        = "BTagging_{wp}_{flav}_{calib}_{algo}_2017BtoF.json",
        #                        format_dict     = {'algo':["DeepFlavour","CSVv2", "DeepCSV"],'wp':["loose", "medium", "tight"],('flav', 'calib'):[("lightjets", "incl"), ("cjets", "comb"), ("bjets","comb")]})
        ##------- single muon trigger ----------#
        #instance.addEfficiency(key_entry   = 'mutrig_2017_94X',
        #                       base_str    = "{trig}_PtEtaBins_2017RunBtoF.json",
        #                       format_dict = {'trig':["IsoMu27", "Mu50"]})
       #
       # # Get full dict #
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
    
instance = ScaleFactorsbbWW()
