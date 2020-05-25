import os
import sys
import re
import itertools
import copy
from pprint import pprint

class MakeScaleFactorsDict():
    def __init__(self,paths,check_path=False):
        """
        paths       : dict containing paths to different directories (kays will be used when registering SF)
        check_path  : (bool) whether to check the json paths or not
        """
        self.paths      = paths 
        self.check_path = check_path 
        self.full_dict  = {}

    def GetScaleFactorsDict(self):
        """ Returns the full SF dictionary path """
        return self.full_dict

    def _flattenTuple(self,nestedTup):
        """ Flatten nested tuple into a simple tuple """
        l = []
        for t in nestedTup:
            if isinstance(t,tuple):    
                l += list(t)
            else:
                l.append(t)
        return tuple(l)

    def _makeFullSFPath(self,name,path_key):
        """ find the paths associated to the path_key and produced the absoluet path to the json file """
        full_path = os.path.join(self.paths[path_key],name)
        if os.path.isfile(full_path) or not self.check_path:
            return full_path
        else:
            raise KeyError("Could not find file %s"%full_path) 

    def AddScaleFactor(self,path_key,entry_key,base_str,format_dict={},sublevel=None):
        """
        path_key    : select which one of the paths to search for the json
        entry_key   : new efficiency entry of the dict 
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
            keys = self._flattenTuple(tuple(format_dict.keys()))
            aList = []
            comb = itertools.product(*list(format_dict.values()))
            for c in comb:
                # Flatten the nested comb if need be #
                t = self._flattenTuple(c)
                d = dict((k,v) for k,v in zip(keys,t))
                if sublevel is not None: 
                    for key in d.keys(): # Find the associated key 
                        if key in sublevel.keys(): 
                            tup = sublevel[key][d[key]] # Recover tuple corresponding to entry
                            aList.append(tuple((tup,self._makeFullSFPath(base_str.format(**d)),path_key)))
                else:
                    aList.append(self._makeFullSFPath(base_str.format(**d),path_key))
            aTup = tuple(aList)
        elif len(format_dict.keys()) == 1: # Simple case
            l = []
            for k,val in format_dict.items():
                if not isinstance(k,tuple):
                    for v in val:
                        d = {k:v}
                        l.append(self._makeFullSFPath(base_str.format(**d),path_key))
                else:
                    for v in val:
                        d = dict((nk,nv) for nk,nv in zip(k,v))
                        l.append(self._makeFullSFPath(base_str.format(**d),path_key))
            aTup = tuple(l)
        else:   # No formatting required
            aTup = self._makeFullSFPath(base_str,path_key)

        # Check if key already present #
        if entry_key in self.full_dict.keys():
            prevTup = self.full_dict[entry_key]
            fullTup = tuple(list(prevTup)+list(aTup))
            self.full_dict[entry_key] = fullTup
        else:
            self.full_dict[entry_key] = aTup

    def AddScaleFactorWithWorkingPoint(self,path_key,entry_key,base_key,base_str,format_dict={},sublevel=None,lowercase_keys=False):
        """
        path_key    : select which one of the paths to search for the json
        entry_key   : new efficiency entry of the dict 
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
                dk_forvalues = dict((k,v) for k,v in zip(base_key_formatkeys,self._flattenTuple(combk))) # Save for values in case lowercase keys acked
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
                    nd = dict((k,v) for k,v in zip(onlybase_str_formatkeys,self._flattenTuple(combv)))
                    dv = {**dk_forvalues,**nd}
                    if sublevel is not None: 
                        for key in dv.keys():# Find the associated key in the format_dict
                            if key in sublevel.keys(): 
                                tup = sublevel[key][dv[key]] # Recover tuple corresponding to entry
                                aList.append(tuple((tup,self._makeFullSFPath(base_str.format(**dv),path_key))))
                    else:
                        aList.append(self._makeFullSFPath(base_str.format(**dv),path_key))
                        
                aDict[base_key.format(**dk)] = tuple(aList) if len(aList)>1 else aList[0]
        elif len(format_dict.keys()) == 1: # Simple scalefactor
            aDict = {}
            key = list(format_dict.keys())[0] # Can be str or tuple
            val = list(format_dict.values())[0]    # Must be a list
            for v in val:
                if isinstance(key,tuple):
                    d = dict((k,av) for k,av in zip(key,v))
                else:
                    d = {key:v}
                if lowercase_keys:
                    aDict[base_key.format(**d).lower()] = self._makeFullSFPath(base_str.format(**d),path_key)
                else:
                    aDict[base_key.format(**d)] = self._makeFullSFPath(base_str.format(**d),path_key)
        else:   # No formatting required
            aDict[base_key] = self._makeFullSFPath(base_str,path_key)

        # Check if key already present #
        if entry_key in self.full_dict.keys():
            self.full_dict[entry_key].update(aDict)
        else:
            self.full_dict[entry_key] = aDict


