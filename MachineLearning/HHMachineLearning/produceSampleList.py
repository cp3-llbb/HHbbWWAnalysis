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
    'Skim{era}_SLBackground_BoostedHbb2Wj_El' : 'boostedHbb2Wj_El_{node}',
    'Skim{era}_SLBackground_BoostedHbb2Wj_Mu' : 'boostedHbb2Wj_Mu_{node}',
    'Skim{era}_SLSignal_BoostedHbb2Wj_El' : 'boostedHbb2Wj_El_HH',
    'Skim{era}_SLSignal_BoostedHbb2Wj_Mu' : 'boostedHbb2Wj_Mu_HH',
    'Skim{era}_SLBackground_BoostedHbb1Wj_El' : 'boostedHbb1Wj_El_{node}',
    'Skim{era}_SLBackground_BoostedHbb1Wj_Mu' : 'boostedHbb1Wj_Mu_{node}',
    'Skim{era}_SLSignal_BoostedHbb1Wj_El' : 'boostedHbb1Wj_El_HH',
    'Skim{era}_SLSignal_BoostedHbb1Wj_Mu' : 'boostedHbb1Wj_Mu_HH',
    'Skim{era}_SLBackground_BoostedHbb0Wj_El' : 'boostedHbb0Wj_El_{node}',
    'Skim{era}_SLBackground_BoostedHbb0Wj_Mu' : 'boostedHbb0Wj_Mu_{node}',
    'Skim{era}_SLSignal_BoostedHbb0Wj_El' : 'boostedHbb0Wj_El_HH',
    'Skim{era}_SLSignal_BoostedHbb0Wj_Mu' : 'boostedHbb0Wj_Mu_HH',
}


base_path = '/nfs/scratch/fynu/gsaha/SkimmedNtuplesJPA/'

sampleDict = {} 
for era in ['2016','2017','2018']:
    back_yaml = '/home/users/g/s/gsaha/bamboodev/HHbbWWAnalysis/analysis{era}_v7_MC.yml'.format(era=era)
    sig_yaml = '/home/users/g/s/gsaha/bamboodev/HHbbWWAnalysis/analysis{era}_v7_signal_SL.yml'.format(era=era)
    
    with open(back_yaml,'r') as f:
        yamlBack = yaml.load(f)
    with open(sig_yaml,'r') as f:
        yamlSig = yaml.load(f)
    
    samplesBack = {k:v['group'] for k,v in yamlBack['samples'].items()}
    samplesSig= [k for k in yamlSig['samples'].keys()]
    #conv_dict = {'DY':'DY','ttbar':'TT','ST':'ST','VVV':'VVV','ttVX':'TTVX','SM':'H','Rares':'Rare','Wjets':'Rares'} # DL
    conv_dict = {'DY':'DY','ttbar':'TT','ST':'ST','VVV':'Rare','ttVX':'Rare','SM':'H','Rares':'Rare','Wjets':'WJets'} # SL
    
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

    
    


