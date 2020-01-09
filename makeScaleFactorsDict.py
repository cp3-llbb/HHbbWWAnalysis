import os
import sys
import re
import itertools
import copy

from pprint import pprint
class MakeScaleFactorsDict():
    def __init__(self,trigger_path,scalefactors_path):
        self.trigger_path = trigger_path
        self.scalefactors_path = scalefactors_path
        self.full_dict = {}

    def getScaleFactorsDict(self):
        return self.full_dict

    def addEfficiency(self,key,base_str,format_dict):
        """
        Add an efficiency entry to the dict -> {'key':tuple(jsonx4)}
        base_str is a string containing {}(s) for formatting (with keys)
        format_dict is a dict whose keys correspond to the ones in base_str and values are list of str to be formatted
        Note : if key already present, the content will be appended to the tuple
        """
        N_format = base_str.count('{')
        if len(format_dict.keys())>1: # More than one formatting
            assert N_format == len(format_dict.keys())
            l = []
            comb = itertools.product(*list(format_dict.values()))
            for c in comb:
                d = dict((k,v) for k,v in zip(format_dict.keys(),c))
                l.append(os.path.join(self.trigger_path,base_str.format(**d)))
            aTup = tuple(l)
        else:
            assert N_format == 1
            l = []
            for k,val in format_dict.items():
                for v in val:
                    d = {k:v}
                    l.append(os.path.join(self.trigger_path,base_str.format(**d)))
            aTup = tuple(l)

        # Check if key already present #
        if key in self.full_dict.keys():
            prevTup = self.full_dict[key]
            fullTup = tuple(list(prevTup)+list(aTup))
            self.full_dict[key] = fullTup
        else:
            self.full_dict[key] = aTup

    def addScaleFactor(self,key,base_key,base_str,format_dict):
        N_format_key = base_key.count('{') 
        N_format_str = base_str.count('{')
        print ('------------------')
        print (key)
        pprint(format_dict)
        if len(format_dict.keys())>1: # More than one formatting
            assert max(N_format_key,N_format_str) == len(format_dict.keys())
            aDict = {}
            # Find the format keys in base_key and base_str #
            base_key_formatkeys = re.findall(r"\{([A-Za-z0-9_]+)\}", base_key)
            base_str_formatkeys = re.findall(r"\{([A-Za-z0-9_]+)\}", base_str)
            # get format keys in base_str but not in base_key #
            onlybase_str_formatkeys = [b for b in base_str_formatkeys if b not in base_key_formatkeys]
            # Generate base_key combinations #
            al = [format_dict[k] for k in base_key_formatkeys]
            comb_key = itertools.product(*al)
            for comb in comb_key:
                dk = dict((k,v) for k,v in zip(base_key_formatkeys,comb))
                aList = []
                for k in onlybase_str_formatkeys:
                    for v in format_dict[k]:
                        dv = copy.deepcopy(dk)
                        dv[k] = v
                        aList.append(os.path.join(self.scalefactors_path,base_str.format(**dv)))
                aDict[base_key.format(**dk)] = tuple(aList)
        else:
            assert N_format_key == 1 and N_format_str == 1  # Simple scalefactor
            aDict = {}
            for k,val in format_dict.items():
                for v in val:
                    dk = {k:v.lower()}
                    dv = {k:v}
                    aDict[base_key.format(**dk)] = os.path.join(self.scalefactors_path,base_str.format(**dv))

        # Check if key already present #
        if key in self.full_dict.keys():
            self.full_dict[key].update(aDict)
        else:
            self.full_dict[key] = aDict


        

inst = MakeScaleFactorsDict(os.path.join(os.path.abspath('.'),'data','TriggerEfficienciesStudies'),os.path.join(os.path.abspath('.'),'data','ScaleFactors_FullRunII'))
#inst.addEfficiency('doubleEleLeg_HHMoriond17_2016','{wp}.json',{'wp':["Electron_IsoEle23Leg", "Electron_IsoEle12Leg", "Electron_IsoEle23Leg", "Electron_IsoEle12Leg"]})
#inst.addEfficiency('mutrig_2016_94X','{trig}_PtEtaBins_2016Run{eras}.json',{'trig':["IsoMu24_OR_IsoTkMu24","Mu50_OR_TkMu50"],'eras':["BtoF", "GtoH"]})
#inst.addScaleFactor('electron_2016_94X',"id_{wp}","Electron_EGamma_SF2D_2016Legacy_{wp}_Fall17V2.json",{"wp":["Loose", "Medium", "Tight", "MVA80","MVA90", "MVA80noiso", "MVA90noiso"]})
#inst.addScaleFactor('btag_2016_94X','{algo}_{wp}',"BTagging_{wp}_{flav}_{algo}.json",{'algo':["DeepCSV", "DeepJet"],'wp':["loose", "medium", "tight"],'flav':["lightjets_incl","cjets_comb","bjets_comb"]})
inst.addScaleFactor('btag_2016_94X','{algo}_{wp}',"BTagging_{wp}_{flav}_{calib}_{algo}_2017BtoF.json",{'algo':["DeepCSV", "DeepJet"],'wp':["loose", "medium", "tight"],('flav','calib'):[("lightjets", "incl"),("cjets","comb"),("bjets","comb")]})
#inst.addScaleFactor('muon_2017_94X','id_{wp}',"Muon_NUM_{wp}ID_DEN_genTracks_pt_abseta_2017RunBCDEF.json",{'wp':["Loose", "Medium", "Tight", "Soft", "MediumPrompt"]})
#inst.addScaleFactor('muon_2017_94X_supp','iso_{isowp}_id_{idwp}',"Muon_NUM_{isowp}RelIso_DEN_{idwp}ID_pt_abseta_2017RunBCDEF.json",{'wp':["Loose", "Medium", "Tight", "Soft", "MediumPrompt"]})
dico = inst.getScaleFactorsDict()
pprint (dico)

