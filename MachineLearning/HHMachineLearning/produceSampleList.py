import os
import sys
import yaml
from glob import glob
from pprint import pprint
from collections import defaultdict 
from itertools import chain

"""
Usage :
produceSampleList.py name_of_yaml.yml
"""

dict_dir = {
'SkimDL{era}_Background_Resolved1Btag_ElEl' : 'resolved_1b_ElEl_{node}',
'SkimDL{era}_Background_Resolved1Btag_ElMu' : 'resolved_1b_ElMu_{node}',
'SkimDL{era}_Background_Resolved1Btag_MuMu' : 'resolved_1b_MuMu_{node}',
'SkimDL{era}_Signal_Resolved1Btag_ElEl' : 'resolved_1b_ElEl_HH',
'SkimDL{era}_Signal_Resolved1Btag_ElMu' : 'resolved_1b_ElMu_HH',
'SkimDL{era}_Signal_Resolved1Btag_MuMu' : 'resolved_1b_MuMu_HH',
'SkimDL{era}_Background_Resolved2Btag_ElEl' : 'resolved_2b_ElEl_{node}',
'SkimDL{era}_Background_Resolved2Btag_ElMu' : 'resolved_2b_ElMu_{node}',
'SkimDL{era}_Background_Resolved2Btag_MuMu' : 'resolved_2b_MuMu_{node}',
'SkimDL{era}_Signal_Resolved2Btag_ElEl' : 'resolved_2b_ElEl_HH',
'SkimDL{era}_Signal_Resolved2Btag_ElMu' : 'resolved_2b_ElMu_HH',
'SkimDL{era}_Signal_Resolved2Btag_MuMu' : 'resolved_2b_MuMu_HH',
'SkimDL{era}_Background_Boosted_ElEl' : 'boosted_ElEl_{node}',
'SkimDL{era}_Background_Boosted_ElMu' : 'boosted_ElMu_{node}',
'SkimDL{era}_Background_Boosted_MuMu' : 'boosted_MuMu_{node}',
'SkimDL{era}_Signal_Boosted_ElEl' : 'boosted_ElEl_HH',
'SkimDL{era}_Signal_Boosted_ElMu' : 'boosted_ElMu_HH',
'SkimDL{era}_Signal_Boosted_MuMu' : 'boosted_MuMu_HH',
}


base_path = '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/'

sampleDict = {} 
for era in ['2016','2017','2018']:
    back_yaml = '/home/users/f/b/fbury/bamboodev/HHbbWWAnalysis/analysis{era}_v6.yml'.format(era=era)
    sig_yaml = '/home/users/f/b/fbury/bamboodev/HHbbWWAnalysis/analysis{era}_v6_signal.yml'.format(era=era)
    
    with open(back_yaml,'r') as f:
        yamlBack = yaml.load(f)
    with open(sig_yaml,'r') as f:
        yamlSig = yaml.load(f)
    
    samplesBack = {k:v['group'] for k,v in yamlBack['samples'].items()}
    samplesSig= [k for k in yamlSig['samples'].keys()]
    conv_dict = {'DY':'DY','ttbar':'TT','ST':'ST','VVV':'VVV','ttVX':'TTVX','SM':'H','Rares':'Rare','Wjets':'Rares'}
    
    nodeDict = {}
    dirDict = defaultdict(list) 
    for directory, name in dict_dir.items():
        directory = directory.format(era=era)
        path = os.path.join(base_path,directory,'results','*.root')
        keys = []
        for f in glob(path):
            if 'skeleton' in f:
                continue
            sample = os.path.basename(f).replace('.root','')
            sample_path = os.path.join(directory,'results',sample+'.root')
            if '{node}' in name:
                if sample not in samplesBack.keys():
                    print ('[WARNING] Could not find sample %s in yaml file %s'%(sample,back_yaml))
                else:
                    group = samplesBack[sample]
                    node = conv_dict[group]
                    dirDict[name.format(node=node)].append(sample_path)
            else:
                if sample in samplesSig:
                    dirDict[name].append(sample_path)
                else:
                    print ('[WARNING] Could not find sample %s in yaml file'%(sample,sig_yaml))
    sampleDict[era] = dict(dirDict) 
    # Check if all samples are present #
    all_samples = list(set(os.path.basename(s).replace('.root','') for s in chain.from_iterable(dirDict.values())))
    for sample in list(yamlBack['samples'].keys())+list(yamlSig['samples'].keys()):
        if sample not in all_samples:
            print ('[WARNING] sample %s in yaml file of era %s has not been produced and registered'%(sample,era))


yamlDict = {'sampleDir':base_path,'sampleDict':sampleDict}
with open(sys.argv[1],'w') as f:
    yaml.dump(yamlDict,f)

    
    


