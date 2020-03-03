import os
import sys
import re
import itertools
import copy
from pprint import pprint

class MakeScaleFactorsDict():
    def __init__(self,efficiencies_path,scalefactors_path,check_path=False):
        self.efficiencies_path       = efficiencies_path
        self.scalefactors_path  = scalefactors_path
        self.check_path         = check_path 
        self.full_dict          = {}

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

    def makeFullEfficiencyPath(self,name):
        full_path = os.path.join(self.efficiencies_path,name)
        if os.path.isfile(full_path) or not self.check_path:
            return full_path
        else:
            print ("[ERROR] Incorrect path to efficiency json file : %s"%full_path)

    def makeFullScaleFactorPath(self,name):
        full_path = os.path.join(self.scalefactors_path,name)
        if os.path.isfile(full_path) or not self.check_path:
            return full_path
        else:
            print ("[ERROR] Incorrect path to scalefactor json file : %s"%full_path)



    def addEfficiency(self,key_entry,base_str,format_dict,sublevel=None):
        """
        key_entry   : new efficiency entry of the dict 
        base_str    : string to be formatted with {keyword} (not empty, can be several)
        format_dict : dict whose keys correspond to the ones in base_str and values are list of str to be formatted
        sublevel    : dict where keys are in format_dict and values are run numbers
        Note : if key already present, the content will be appended to the tuple
        Note : if one doesn't want all combinations, the keys and related values can be expressed as tuples in format_dict
        Will generate one entry in the dict as a *tuple*
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
                if sublevel is not None: 
                    for key in d.keys(): # Find the associated key 
                        if key in sublevel.keys(): 
                            tup = sublevel[key][d[key]] # Recover tuple corresponding to entry
                            aList.append(tuple((tup,self.makeFullEfficiencyPath(base_str.format(**d)))))
                else:
                    aList.append(self.makeFullEfficiencyPath(base_str.format(**d)))
            aTup = tuple(aList)
        else:
            l = []
            for k,val in format_dict.items():
                if not isinstance(k,tuple):
                    for v in val:
                        d = {k:v}
                        l.append(self.makeFullEfficiencyPath(base_str.format(**d)))
                else:
                    for v in val:
                        d = dict((nk,nv) for nk,nv in zip(k,v))
                        l.append(self.makeFullEfficiencyPath(base_str.format(**d)))
            aTup = tuple(l)

        # Check if key already present #
        if key_entry in self.full_dict.keys():
            prevTup = self.full_dict[key_entry]
            fullTup = tuple(list(prevTup)+list(aTup))
            self.full_dict[key_entry] = fullTup
        else:
            self.full_dict[key_entry] = aTup

    def addScaleFactor(self,key_entry,base_key,base_str,format_dict,sublevel=None,lowercase_keys=False):
        """
        key_entry   : new efficiency entry of the dict 
        base_key    : string to be formatted with {keyword} (not empty, can be several) and used as nested dict key
        base_str    : string to be formatted with {keyword} (not empty, can be several) and used as nested dict value
        format_dict : dict whose keys correspond to the ones in base_str and base_key and values are list of str to be formatted
        sublevel    : dict where keys are in format_dict and values are run numbers
        lowercase_keys : lower case the keys in the entry dict
        Note : if key already present, the content will be appended to the tuple
        Note : if one doesn't want all combinations, the keys and related values can be expressed as tuples in format_dict
        Will generate one entry in the dict as a *dict* whose keys are working points
        """
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
                dk_forvalues = dict((k,v) for k,v in zip(base_key_formatkeys,self.flattenTuple(combk))) # Save for values in case lowercase keys acked
                if lowercase_keys:
                    dk = dict((k,v.lower()) for k,v in dk_forvalues.items())
                else:
                    dk = dict((k,v) for k,v in dk_forvalues.items())
                    
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
                    nd = dict((k,v) for k,v in zip(onlybase_str_formatkeys,self.flattenTuple(combv)))
                    dv = {**dk_forvalues,**nd}
                    if sublevel is not None: 
                        for key in dv.keys():# Find the associated key in the format_dict
                            if key in sublevel.keys(): 
                                tup = sublevel[key][dv[key]] # Recover tuple corresponding to entry
                                aList.append(tuple((tup,self.makeFullScaleFactorPath(base_str.format(**dv)))))
                    else:
                        aList.append(self.makeFullScaleFactorPath(base_str.format(**dv)))
                        
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
                    aDict[base_key.format(**d).lower()] = self.makeFullScaleFactorPath(base_str.format(**d))
                else:
                    aDict[base_key.format(**d)] = self.makeFullScaleFactorPath(base_str.format(**d))

        # Check if key already present #
        if key_entry in self.full_dict.keys():
            self.full_dict[key_entry].update(aDict)
        else:
            self.full_dict[key_entry] = aDict


