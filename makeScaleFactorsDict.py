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

    def flattenTuple(self,nestedTup):
        l = []
        for t in nestedTup:
            if isinstance(t,tuple):    
                l += list(t)
            else:
                l.append(t)
        return tuple(l)


    def addEfficiency(self,key_entry,base_str,format_dict):
        """
        Add an efficiency entry to the dict -> {'key':tuple(jsonx4)}
        base_str is a string containing {}(s) for formatting (with keys)
        format_dict is a dict whose keys correspond to the ones in base_str and values are list of str to be formatted
        Note : if key already present, the content will be appended to the tuple
        """
        N_format = base_str.count('{')
        if len(format_dict.keys())>1: # More than one formatting
            # Flatten the nested keys if need be #
            keys = self.flattenTuple(tuple(format_dict.keys()))
            aList = []
            comb = itertools.product(*list(format_dict.values()))
            for c in comb:
                # Flatten the nested comb if need be #
                t = self.flattenTuple(c)
                d = dict((k,v) for k,v in zip(keys,t))
                aList.append(os.path.join(self.trigger_path,base_str.format(**d)))
            aTup = tuple(aList)
        else:
            l = []
            for k,val in format_dict.items():
                if not isinstance(k,tuple):
                    for v in val:
                        d = {k:v}
                        l.append(os.path.join(self.trigger_path,base_str.format(**d)))
                else:
                    for v in val:
                        d = dict((nk,nv) for nk,nv in zip(k,v))
                        l.append(os.path.join(self.trigger_path,base_str.format(**d)))
            aTup = tuple(l)

        # Check if key already present #
        if key_entry in self.full_dict.keys():
            prevTup = self.full_dict[key_entry]
            fullTup = tuple(list(prevTup)+list(aTup))
            self.full_dict[key_entry] = fullTup
        else:
            self.full_dict[key_entry] = aTup

    def addScaleFactor(self,key_entry,base_key,base_str,format_dict,sublevel=None,lowercase_keys=False):
        N_format_key = base_key.count('{') 
        N_format_str = base_str.count('{')
        if len(format_dict.keys())>1: # More than one formatting
            aDict = {}
            # Find the format keys in base_key and base_str #
            base_key_formatkeys = re.findall(r"\{([A-Za-z0-9_]+)\}", base_key)
            base_str_formatkeys = re.findall(r"\{([A-Za-z0-9_]+)\}", base_str)
            # Check whether the format keys in base_key are tuples in format_dict keys #
            tuple_in_base_key = False
            for key_format in base_key_formatkeys:
                for key_dict in format_dict.keys():
                    if isinstance(key_dict,tuple) and key_format in key_dict:
                        tuple_in_base_key = True
            # Check whether the format keys in base_str are tuples in format_dict keys #
            tuple_in_base_str = False
            for key_format in base_str_formatkeys:
                for key_dict in format_dict.keys():
                    if isinstance(key_dict,tuple) and key_format in key_dict:
                        tuple_in_base_str = True
            # get format keys in base_str but not in base_key #
            onlybase_str_formatkeys = [b for b in base_str_formatkeys if b not in base_key_formatkeys]
            # Generate base_key combinations #
            # Check if tuple in the key of format_dict #
            if tuple_in_base_key:
                al = []
                for key in format_dict.keys():
                    if isinstance(key,tuple):
                        for k in key: # We assume that if one of the val in tuple is in onlybase_str_formatkeys, the others are as well
                            if k in base_key_formatkeys:
                                al.append(format_dict[key])
                                break
                    else:
                        if key in base_key_formatkeys:
                            al.append(format_dict[key])
            else:
                al = [format_dict[k] for k in base_key_formatkeys]
            if (len(al)==1): # Otherwise will do product on tuple
                comb_key = itertools.product(al[0])
            elif (len(al)>1):
                comb_key = itertools.product(*al)
            # Loop over key combinations #
            for combk in comb_key:
                if lowercase_keys:
                    dk = dict((k,v.lower()) for k,v in zip(base_key_formatkeys,self.flattenTuple(combk)))
                else:
                    dk = dict((k,v) for k,v in zip(base_key_formatkeys,self.flattenTuple(combk)))
                aList = []
                # Generate values combination for each key combination #
                # Check if tuple, must not be inter product #
                if tuple_in_base_str:
                    bl = []
                    for key in format_dict.keys():
                        if isinstance(key,tuple):
                            for k in key: # We assume that if one of the val in tuple is in onlybase_str_formatkeys, the others are as well
                                if k in onlybase_str_formatkeys: 
                                    bl.append(key)
                                    break
                        else:
                            if key in onlybase_str_formatkeys:
                                bl.append(key) 
                    comb_val = itertools.product(*[format_dict[k] for k in bl])
                else:
                    comb_val = itertools.product(*[format_dict[k] for k in onlybase_str_formatkeys])
                # Generate list of values comb and add that to tmp dict #
                for combv in comb_val:
                    print (combv)
                    nd = dict((k,v) for k,v in zip(onlybase_str_formatkeys,self.flattenTuple(combv)))
                    dv = {**dk,**nd}
                    print (dv)
                    print (os.path.join(self.scalefactors_path,base_str.format(**dv)))
                    if sublevel is not None: 
                        for key in dv.keys():# Find the associated key in the format_dict
                            if key in sublevel.keys(): 
                                tup = sublevel[key][dv[key]] # Recover tuple corresponding to entry
                                aList.append(tuple((tup,os.path.join(self.scalefactors_path,base_str.format(**dv)))))
                    else:
                        aList.append(os.path.join(self.scalefactors_path,base_str.format(**dv)))
                        
                aDict[base_key.format(**dk)] = tuple(aList)
        else:# Simple scalefactor
            aDict = {}
            key = list(format_dict.keys())[0] # Can be str or tuple
            val = list(format_dict.values())[0]    # Must be a list
            for v in val:
                if isinstance(key,tuple):
                    d = dict((k,av) for k,av in zip(key,v))
                else:
                    d = {key:v}
                if lowercase_keys:
                    aDict[base_key.format(**d).lower()] = os.path.join(self.scalefactors_path,base_str.format(**d))
                else:
                    aDict[base_key.format(**d)] = os.path.join(self.scalefactors_path,base_str.format(**d))

        # Check if key already present #
        if key_entry in self.full_dict.keys():
            self.full_dict[key_entry].update(aDict)
        else:
            self.full_dict[key_entry] = aDict


        

inst = MakeScaleFactorsDict(os.path.join(os.path.abspath('.'),'data','TriggerEfficienciesStudies'),os.path.join(os.path.abspath('.'),'data','ScaleFactors_FullRunII'))
#inst.addEfficiency('doubleEleLeg_HHMoriond17_2016','{wp}.json',{'wp':["Electron_IsoEle23Leg", "Electron_IsoEle12Leg", "Electron_IsoEle23Leg", "Electron_IsoEle12Leg"]})
#inst.addEfficiency('doubleEleLeg_HHMoriond17_2016','{wp}_{calib}.json',{('wp','calib'):[("Electron_IsoEle23Leg",'testA'), ("Electron_IsoEle12Leg",'testB'), ("Electron_IsoEle23Leg",'testC'), ("Electron_IsoEle12Leg",'testA')]})
#inst.addEfficiency('mutrig_2016_94X','{trig}_PtEtaBins_2016Run{eras}.json',{'trig':["IsoMu24_OR_IsoTkMu24","Mu50_OR_TkMu50"],'eras':["BtoF", "GtoH"]})
#inst.addEfficiency('doubleEleLeg_HHMoriond17_2016','{type}_{other}_{wp}_{calib}.json',{'type':['TypeI','TypeII'],'other':['otherA','otherB'],('wp','calib'):[("Electron_IsoEle23Leg",'testA'), ("Electron_IsoEle12Leg",'testB'), ("Electron_IsoEle23Leg",'testC'), ("Electron_IsoEle12Leg",'testA')]})
#inst.addScaleFactor('electron_2016_94X',"id_{wp}","Electron_EGamma_SF2D_2016Legacy_{wp}_Fall17V2.json",{"wp":["Loose", "Medium", "Tight", "MVA80","MVA90", "MVA80noiso", "MVA90noiso"]})
#inst.addScaleFactor('btag_2016_94X','{algo}_{wp}',"BTagging_{wp}_{flav}_{algo}.json",{'algo':["DeepCSV", "DeepJet"],'wp':["loose", "medium", "tight"],'flav':["lightjets_incl","cjets_comb","bjets_comb"]})
#inst.addScaleFactor('muon_2017_94X','id_{wp}',"Muon_NUM_{wp}ID_DEN_genTracks_pt_abseta_2017RunBCDEF.json",{'wp':["Loose", "Medium", "Tight", "Soft", "MediumPrompt"]},True)
#inst.addScaleFactor('muon_2017_94X_supp','iso_{isowp}_id_{idwp}',"Muon_NUM_{isowp}RelIso_DEN_{idwp}ID_pt_abseta_2017RunBCDEF.json",{'wp':["Loose", "Medium", "Tight", "Soft", "MediumPrompt"]})
inst.addScaleFactor('muon_2016_94X','id_{wp}',"Muon_NUM_{wp}ID_DEN_genTracks_eta_pt_{uncer}_2016Run{era}.json",
                    {'wp':["Loose", "Medium", "Tight"],
                     'uncer':['sys','stat'],
                     'era':["BCDEF", "GH"]},
                    sublevel = {'era':{'BCDEF':('Run2016B','Run2016C','Run2016D','Run2016E','Run2016F'),'GH':('Run2016G', 'Run2016H')}},
                    lowercase_keys=True)
dico = inst.getScaleFactorsDict()
pprint (dico)

