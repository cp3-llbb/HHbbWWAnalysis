import os 
import sys
from bamboo import scalefactors
from bamboo import treefunctions as op

from scalefactorsutils import MakeScaleFactorsDict

class ScaleFactorsbbWW:
    def __init__(self):
        self.binningVariables = {
        "Eta"       : lambda obj : obj.eta,
        "ClusEta"   : lambda obj : obj.eta + obj.deltaEtaSC,
        "AbsEta"    : lambda obj : op.abs(obj.eta),
        "AbsClusEta": lambda obj : op.abs(obj.clusterEta+obj.deltaEtaSC),
        "Pt"        : lambda obj : obj.pt,
        }

         # Scale factors dictionary construction #
         # (better to do that here so that it is not buil everytime we ask for scalefactors #
        instance = MakeScaleFactorsDict(paths       = {'ttH_SF' : os.path.join(os.path.dirname(os.path.abspath(__file__)),'data','ScaleFactors_ttH'),
                                                       'DY_SF'  : os.path.join(os.path.dirname(os.path.abspath(__file__)),'data','ScaleFactors_DY'),
                                                       'POG_SF' : os.path.join(os.path.dirname(os.path.abspath(__file__)),'data','ScaleFactors_POG'),
                                                       'Btag_SF': os.path.join(os.path.dirname(os.path.abspath(__file__)),'data','ScaleFactors_Btag'),
                                                       'FR'     : os.path.join(os.path.dirname(os.path.abspath(__file__)),'data','FakeRates')},
                                        check_path  = True)
        #----- 2016 -----#
        # Single Trigger SF # (Double triggers are single numbers and are in Base)
        instance.AddScaleFactor(path_key    = 'ttH_SF',
                                entry_key   = 'singleTrigger_electron_2016',
                                base_str    = "TTH_trigger_SingleElectron_2016.json")
        instance.AddScaleFactor(path_key    = 'ttH_SF',
                                entry_key   = 'singleTrigger_muon_2016',
                                base_str    = "TTH_trigger_SingleMuon_2016.json")

        ### ttH mva ###
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
        
        ### POG ID ###
            # https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffs2016LegacyRereco
            #   -> SF root in "Systematic uncertainties" : ID + ISO /!\ needs to be lumi reweighted
        instance.AddScaleFactor(path_key    = 'POG_SF',
                                entry_key   = 'muon_POGSF_ID_2016',
                                base_str    = 'Muon_NUM_TightID_DEN_genTracks_pt_eta_statPlusSyst_2016_RunBCDEFGH_ID.json')
        instance.AddScaleFactor(path_key    = 'POG_SF',
                                entry_key   = 'muon_POGSF_ISO_2016',
                                base_str    = 'Muon_NUM_TightRelIso_DEN_TightIDandIPCut_pt_eta_statPlusSyst_2016_RunBCDEFGH_ISO.json')
        instance.AddScaleFactor(path_key    = 'POG_SF',
                                entry_key   = 'electron_POGSF_2016',
                                base_str    = 'Electron_EGamma_SF2D_2016.json')


#        # DY weight  # 
#        instance.AddScaleFactorWithWorkingPoint(path_key    = 'DY_SF',
#                                                entry_key   = 'DY_resolved_2016',
#                                                base_key    = '{channel}_{variable}_{type}_{btag}',
#                                                base_str    = 'weight_{variable}_{channel}_{type}_1D_weight_{btag}_2016.json',
#                                                format_dict = {'channel':['ElEl','MuMu','SSDL'],'type':['data','mc'],'btag':['1b','2b'],'variable':['leadjetPt']})
#        instance.AddScaleFactorWithWorkingPoint(path_key    = 'DY_SF',
#                                                entry_key   = 'DY_boosted_2016',
#                                                base_key    = '{channel}_{variable}_{type}_{btag}',
#                                                base_str    = 'weight_{variable}_{channel}_{type}_1D_weight_{btag}_2016.json',
#                                                format_dict = {'channel':['ElEl','MuMu','SSDL'],'type':['data','mc'],'btag':['1b'],'variable':['fatjetsoftDropmass']})


        #  Fake rates #
        instance.AddScaleFactorWithWorkingPoint(path_key    = 'FR',
                                                entry_key   = 'electron_fakerates_2016',
                                                base_key    = '{wp}_{syst}_syst',
                                                base_str    = 'TTHFakeRates_{wp}MVA_Electron_2016_{syst}Syst.json',
                                                format_dict = {'wp':['Loose','Tight'],'syst':['pt','barrel','norm']})
        instance.AddScaleFactorWithWorkingPoint(path_key    = 'FR',
                                                entry_key   = 'muon_fakerates_2016',
                                                base_key    = '{wp}_{syst}_syst',
                                                base_str    = 'TTHFakeRates_{wp}MVA_Muon_2016_{syst}Syst.json',
                                                format_dict = {'wp':['Loose','Tight'],'syst':['pt','barrel','norm']})
        # PU ID SF #
        instance.AddScaleFactorWithWorkingPoint(path_key    = 'Btag_SF',
                                                entry_key   = 'jet_puid',
                                                base_key    = '{eom}_{mcsf}_{era}_{wp}',
                                                base_str    = 'PUID_80X_{eom}_{mcsf}_{era}_{wp}.json',
                                                format_dict = {'eom':["eff", "mistag"],'mcsf':["mc", "sf"],'wp':['L','M','T'],'era':['2016','2017','2018']})
                                                
        # Ak8 efficiency btagging #
        instance.AddScaleFactorWithWorkingPoint(path_key    = 'Btag_SF',
                                                entry_key   = 'ak8btag_eff_2016',
                                                base_key    = 'eff_{flav}',
                                                base_str    = 'BtagAk8_{flav}_2016.json',
                                                format_dict = {'flav':['bjets','cjets','lightjets']})

        #----- 2017 -----#
        # Check links of 2016 #
        # Single Trigger SF # (Double triggers are single numbers and are in Base)
        instance.AddScaleFactor(path_key    = 'ttH_SF',
                                entry_key   = 'singleTrigger_electron_2017',
                                base_str    = "TTH_trigger_SingleElectron_2017.json")
        instance.AddScaleFactor(path_key    = 'ttH_SF',
                                entry_key   = 'singleTrigger_muon_2017',
                                base_str    = "TTH_trigger_SingleMuon_2017.json")
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

        ### POG ID ###
            # https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffs2017
            #   -> SF root in "Systematic uncertainties" 
        instance.AddScaleFactor(path_key    = 'POG_SF',
                                entry_key   = 'muon_POGSF_ID_2017',
                                base_str    = 'Muon_NUM_TightID_DEN_genTracks_abseta_pt_statPlusSyst_2017_RunBCDEF_ID.json')
        instance.AddScaleFactor(path_key    = 'POG_SF',
                                entry_key   = 'muon_POGSF_ISO_2017',
                                base_str    = 'Muon_NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt_statPlusSyst_2017_RunBCDEF_ISO.json')
        instance.AddScaleFactor(path_key    = 'POG_SF',
                                entry_key   = 'electron_POGSF_2017',
                                base_str    = 'Electron_EGamma_SF2D_2017.json')


        #  Fake rates #
        instance.AddScaleFactorWithWorkingPoint(path_key    = 'FR',
                                                entry_key   = 'electron_fakerates_2017',
                                                base_key    = '{wp}_{syst}_syst',
                                                base_str    = 'TTHFakeRates_{wp}MVA_Electron_2017_{syst}Syst.json',
                                                format_dict = {'wp':['Loose','Tight'],'syst':['pt','barrel','norm']})
        instance.AddScaleFactorWithWorkingPoint(path_key    = 'FR',
                                                entry_key   = 'muon_fakerates_2017',
                                                base_key    = '{wp}_{syst}_syst',
                                                base_str    = 'TTHFakeRates_{wp}MVA_Muon_2017_{syst}Syst.json',
                                                format_dict = {'wp':['Loose','Tight'],'syst':['pt','barrel','norm']})

#        # DY weight  # 
#        instance.AddScaleFactorWithWorkingPoint(path_key    = 'DY_SF',
#                                                entry_key   = 'DY_resolved_2017',
#                                                base_key    = '{channel}_{variable}_{type}_{btag}',
#                                                base_str    = 'weight_{variable}_{channel}_{type}_1D_weight_{btag}_2017.json',
#                                                format_dict = {'channel':['ElEl','MuMu','SSDL'],'type':['mc'],'btag':['1b','2b'],'variable':['leadjetPt']})
#        instance.AddScaleFactorWithWorkingPoint(path_key    = 'DY_SF',
#                                                entry_key   = 'DY_boosted_2017',
#                                                base_key    = '{channel}_{variable}_{type}_{btag}',
#                                                base_str    = 'weight_{variable}_{channel}_{type}_1D_weight_{btag}_2017.json',
#                                                format_dict = {'channel':['ElEl','MuMu','SSDL'],'type':['mc'],'btag':['1b'],'variable':['fatjetsoftDropmass']})
#
        #----- 2018 -----#
        # Check links of 2016 #
        # Single Trigger SF # (Double triggers are single numbers and are in Base)
        instance.AddScaleFactor(path_key    = 'ttH_SF',
                                entry_key   = 'singleTrigger_electron_2018',
                                base_str    = "TTH_trigger_SingleElectron_2018.json")
        instance.AddScaleFactor(path_key    = 'ttH_SF',
                                entry_key   = 'singleTrigger_muon_2018',
                                base_str    = "TTH_trigger_SingleMuon_2018.json")
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

        ### POG ID ###
            # https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffs2017
            #   -> SF root in "Systematic uncertainties" 
        instance.AddScaleFactor(path_key    = 'POG_SF',
                                entry_key   = 'muon_POGSF_ID_2018',
                                base_str    = 'Muon_NUM_TightID_DEN_TrackerMuons_abseta_pt_statPlusSyst_2018_RunABCD_ID.json')
        instance.AddScaleFactor(path_key    = 'POG_SF',
                                entry_key   = 'muon_POGSF_ISO_2018',
                                base_str    = 'Muon_NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt_statPlusSyst_2018_RunABCD_ISO.json')
        instance.AddScaleFactor(path_key    = 'POG_SF',
                                entry_key   = 'electron_POGSF_2018',
                                base_str    = 'Electron_EGamma_SF2D_2018.json')

        #  Fake rates #
        instance.AddScaleFactorWithWorkingPoint(path_key    = 'FR',
                                                entry_key   = 'electron_fakerates_2018',
                                                base_key    = '{wp}_{syst}_syst',
                                                base_str    = 'TTHFakeRates_{wp}MVA_Electron_2018_{syst}Syst.json',
                                                format_dict = {'wp':['Loose','Tight'],'syst':['pt','barrel','norm']})
        instance.AddScaleFactorWithWorkingPoint(path_key    = 'FR',
                                                entry_key   = 'muon_fakerates_2018',
                                                base_key    = '{wp}_{syst}_syst',
                                                base_str    = 'TTHFakeRates_{wp}MVA_Muon_2018_{syst}Syst.json',
                                                format_dict = {'wp':['Loose','Tight'],'syst':['pt','barrel','norm']})

#        # DY weight  # 
#        instance.AddScaleFactorWithWorkingPoint(path_key    = 'DY_SF',
#                                                entry_key   = 'DY_resolved_2018',
#                                                base_key    = '{channel}_{variable}_{type}_{btag}',
#                                                base_str    = 'weight_{variable}_{channel}_{type}_1D_weight_{btag}_2018.json',
#                                                format_dict = {'channel':['ElEl','MuMu','SSDL'],'type':['mc'],'btag':['1b','2b'],'variable':['leadjetPt']})
#        instance.AddScaleFactorWithWorkingPoint(path_key    = 'DY_SF',
#                                                entry_key   = 'DY_boosted_2018',
#                                                base_key    = '{channel}_{variable}_{type}_{btag}',
#                                                base_str    = 'weight_{variable}_{channel}_{type}_1D_weight_{btag}_2018.json',
#                                                format_dict = {'channel':['ElEl','MuMu','SSDL'],'type':['mc'],'btag':['1b'],'variable':['fatjetsoftDropmass']})

        # Get full dict #
        self.all_scalefactors = instance.GetScaleFactorsDict()

    def get_scalefactor(self,objType, key, periods=None, combine=None, additionalVariables=dict(), systName=None, defineOnFirstUse=True):
        paramDefs = self.binningVariables
        if additionalVariables is not None:
            paramDefs.update(additionalVariables)

        return scalefactors.get_scalefactor(objType             = objType, 
                                            key                 = key, 
                                            periods             = periods, 
                                            combine             = combine, 
                                            sfLib               = self.all_scalefactors, 
                                            paramDefs           = paramDefs,
                                            getFlavour          = (lambda j : j.hadronFlavour), 
                                            systName            = systName,
                                            defineOnFirstUse    = defineOnFirstUse)
    
if __name__ == "__main__":
    instance = ScaleFactorsbbWW()
    import pprint
    pprint.pprint(instance.all_scalefactors)

