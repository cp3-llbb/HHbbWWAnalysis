import os
import sys
from itertools import chain
import json
import yaml
import subprocess
import logging
from root_numpy import tree2array
from cppyy import gbl
from collections import defaultdict


triggers_per_era = {}
triggers_per_era["2016"] = ["IsoMu22","IsoTkMu22","IsoMu22_eta2p1","IsoTkMu22_eta2p1","IsoMu24","IsoTkMu24",
                            "Ele27_WPTight_Gsf","Ele25_eta2p1_WPTight_Gsf","Ele27_eta2p1_WPLoose_Gsf",
                            "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL","Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ","Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL","Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",
                            "Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                            "Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL","Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ","Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL","Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ"]
triggers_per_era["2017"] = ["IsoMu24","IsoMu27",
                            "Ele35_WPTight_Gsf","Ele32_WPTight_Gsf",
                            "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8","Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",
                            "Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
                            "Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ","Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ","Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL","Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"]
triggers_per_era["2018"] = ["IsoMu24","IsoMu27",
                            "Ele32_WPTight_Gsf",
                            "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
                            "Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
                            "Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ","Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ","Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL"]

analysis_files = {}
analysis_files["2016"] = "/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/analysis2016_v7.yml"
analysis_files["2017"] = "/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/analysis2017_v7.yml"
analysis_files["2018"] = "/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/analysis2018_v7.yml"
                 

storageBase = "/storage/data/cms"

golden = {
    "2016": os.path.join(os.path.dirname(__file__), "Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"),
    "2017": os.path.join(os.path.dirname(__file__), "Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt"),
    "2018": os.path.join(os.path.dirname(__file__), "Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt")
    }

#years = ["2016","2017","2018"]
years = ["2017"]

cacheFile = "cacheTriggers.json"

missingTriggers = {}
if not os.path.exists(cacheFile):
    print ('Generating cache file')
    for year in years:
        with open(golden[year]) as jsf:
            goldenLumis = json.load(jsf)

        with open(analysis_files[year],"r") as f:
            samples = yaml.load(f)['samples']
        daspaths = defaultdict(list)
        for sample,sampleCfg in samples.items():
            if 'group' not in sampleCfg.keys() or sampleCfg['group']!='data':
                continue
            if isinstance(sampleCfg['db'],list):
                daspaths[sample].extend([db.replace('das:','') for db in sampleCfg['db']])
            else:
                daspaths[sample].append(sampleCfg['db'].replace('das:',''))

        triggers = ['HLT_'+trig for trig in triggers_per_era[year]]

        missingTriggers.update({k:{} for k in daspaths.keys()})

        for sample,daspath in daspaths.items():
            for das in daspath:
                dasout = subprocess.check_output(["dasgoclient", "-query", "file dataset={0}".format(das)])
                lfnList = [storageBase + ln.strip() for ln in dasout.decode().split("\n") if ln.strip() ]
                for pfn in lfnList:
                    if not os.path.exists(pfn):
                        raise RuntimeError('Could not find file %s'%pfn)
                    f = gbl.TFile.Open(pfn, "READ")
                    tEvt = f.Get("Events")
                    lvs = set(lv.GetName() for lv in tEvt.GetListOfLeaves())
                    trig_missing = [ trig for trig in triggers if trig not in lvs ]
                    if trig_missing:
                        tRuns = f.Get("Runs")
                        runs_in_file = [str(run) for run in list(tree2array(tRuns,'run'))]
                        print("File {} is missing {} branches".format(pfn,len(trig_missing)))
                        if not any(run in goldenLumis for run in runs_in_file):
                            print("Contains events from runs {"+', '.join(runs_in_file)+"} none of which on golden, so the file can be skipped")
                            missingTriggers[sample][pfn] = ['not golden lumi']
                        else:
                            missingTriggers[sample][pfn] = trig_missing
                    else:
                        missingTriggers[sample][pfn] = []

    with open(cacheFile,'w') as f:
        json.dump(missingTriggers,f,indent=4)
    print ('Cache file generated at %s'%missingTriggers)
        
print ('Loading cache file')
with open(cacheFile,'r') as f:
    missingTriggers = json.load(f)

groups = {}
for sample, mistrigs_per_file in missingTriggers.items():
    groups[sample] = {}
    combinations = {}
    for f,mistrig in mistrigs_per_file.items():
        mistrig = tuple(mistrig)
        if not mistrig in combinations.keys():
            combinations[mistrig] = 1   
            groups[sample][mistrig] = [f]
        else:
            combinations[mistrig] += 1
            groups[sample][mistrig].append(f)
    print ("In sample %s, found missing triggers combinations below :"%sample)
    for comb,rep in combinations.items():
        print ('...',comb)
        print ('\t-> %d occurences / %d files'%(rep,len(mistrigs_per_file)))

path_out = '/home/users/f/b/fbury/bamboodev/HHbbWWAnalysis/utils/DataRuns'
if not os.path.exists(path_out):
    os.makedirs(path_out)

for sample in groups.keys():
    idx = 0
    for mistrig,files in groups[sample].items():
        if mistrig == ('not golden lumi',):
            name = os.path.join(path_out,'{}_{}.txt'.format(sample,'NotGolden'))
        else:
            name = os.path.join(path_out,'{}_{}.txt'.format(sample,idx))
            idx += 1
        with open(name,'w') as handle:
            for f in files:
                handle.write(f+'\n')
    print ('Produced file {}'.format(name))
        
    
    
     

        


