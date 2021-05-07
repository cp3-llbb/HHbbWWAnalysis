import os
import sys
import re
import subprocess
import shutil
import argparse
import time
import ROOT
from utils.finalize import finalize

parser = argparse.ArgumentParser(description='Launch Skim production')
parser.add_argument('--remove', action='store_true', required=False, default=False,
                    help='Clears the directories')
parser.add_argument('--send', action='store_true', required=False, default=False,
                    help='Sends the jobs')
parser.add_argument('--resubmit', action='store_true', required=False, default=False,
                    help='Resubmit the jobs')
parser.add_argument('--debug', action='store_true', required=False, default=False,
                    help='Print the resubmit commands without sending them')
parser.add_argument('--id', action='store_true', required=False, default=False,
                    help='Return ID of the sent jobs')
parser.add_argument('--era', action='store', required=True, type=str,
                    help='Era')
args = parser.parse_args()

era = args.era

cmd = [['bambooRun','--distributed=driver','-m','SkimmerHHtobbWWDL.py:SkimmerNanoHHtobbWWDL',f'Yaml/analysis{era}_v7.yml','-s=mc,bbww2l,graviton,radion','--Resolved1Btag','--Channel','ElEl','-o',f'/nfs/scratch/fynu/gsaha/SkimmedNtuples_Resonant_DL/Skim{era}_Resolved1Btag_ElEl'],
       ['bambooRun','--distributed=driver','-m','SkimmerHHtobbWWDL.py:SkimmerNanoHHtobbWWDL',f'Yaml/analysis{era}_v7.yml','-s=mc,bbww2l,graviton,radion','--Resolved1Btag','--Channel','ElMu','-o',f'/nfs/scratch/fynu/gsaha/SkimmedNtuples_Resonant_DL/Skim{era}_Resolved1Btag_ElMu'],
       ['bambooRun','--distributed=driver','-m','SkimmerHHtobbWWDL.py:SkimmerNanoHHtobbWWDL',f'Yaml/analysis{era}_v7.yml','-s=mc,bbww2l,graviton,radion','--Resolved1Btag','--Channel','MuMu','-o',f'/nfs/scratch/fynu/gsaha/SkimmedNtuples_Resonant_DL/Skim{era}_Resolved1Btag_MuMu'],
       ['bambooRun','--distributed=driver','-m','SkimmerHHtobbWWDL.py:SkimmerNanoHHtobbWWDL',f'Yaml/analysis{era}_v7.yml','-s=mc,bbww2l,graviton,radion','--Resolved2Btag','--Channel','ElEl','-o',f'/nfs/scratch/fynu/gsaha/SkimmedNtuples_Resonant_DL/Skim{era}_Resolved2Btag_ElEl'],
       ['bambooRun','--distributed=driver','-m','SkimmerHHtobbWWDL.py:SkimmerNanoHHtobbWWDL',f'Yaml/analysis{era}_v7.yml','-s=mc,bbww2l,graviton,radion','--Resolved2Btag','--Channel','ElMu','-o',f'/nfs/scratch/fynu/gsaha/SkimmedNtuples_Resonant_DL/Skim{era}_Resolved2Btag_ElMu'],
       ['bambooRun','--distributed=driver','-m','SkimmerHHtobbWWDL.py:SkimmerNanoHHtobbWWDL',f'Yaml/analysis{era}_v7.yml','-s=mc,bbww2l,graviton,radion','--Resolved2Btag','--Channel','MuMu','-o',f'/nfs/scratch/fynu/gsaha/SkimmedNtuples_Resonant_DL/Skim{era}_Resolved2Btag_MuMu'],
       ['bambooRun','--distributed=driver','-m','SkimmerHHtobbWWDL.py:SkimmerNanoHHtobbWWDL',f'Yaml/analysis{era}_v7.yml','-s=mc,bbww2l,graviton,radion','--Boosted1Btag','--Channel','ElEl','-o',f'/nfs/scratch/fynu/gsaha/SkimmedNtuples_Resonant_DL/Skim{era}_Boosted1Btag_ElEl'],
       ['bambooRun','--distributed=driver','-m','SkimmerHHtobbWWDL.py:SkimmerNanoHHtobbWWDL',f'Yaml/analysis{era}_v7.yml','-s=mc,bbww2l,graviton,radion','--Boosted1Btag','--Channel','ElMu','-o',f'/nfs/scratch/fynu/gsaha/SkimmedNtuples_Resonant_DL/Skim{era}_Boosted1Btag_ElMu'],
       ['bambooRun','--distributed=driver','-m','SkimmerHHtobbWWDL.py:SkimmerNanoHHtobbWWDL',f'Yaml/analysis{era}_v7.yml','-s=mc,bbww2l,graviton,radion','--Boosted1Btag','--Channel','MuMu','-o',f'/nfs/scratch/fynu/gsaha/SkimmedNtuples_Resonant_DL/Skim{era}_Boosted1Btag_MuMu']]

def remove():
    for c in cmd:
        path = c[-1]
        if os.path.exists(path):
            shutil.rmtree(path)
            print ("Removed %s"%path)

def send():
    for c in cmd:
        path = c[-1]
        if os.path.exists(path):
            print ('Path exists, will not launch jobs : %s'%path)
        else:
            subprocess.Popen(c)
        #print (' '.join(c))

def getSlurmID():
    cid = []
    for c in cmd:
        path = c[-1]
        c_id_path = os.path.join(path,'batch','input','cluster_id') 
        with open(c_id_path,'r') as f:
            cid.append(f.readline())
    print ('ID of jobs')
    print (' '.join(cid))

def resubmit(debug):
    res_cmd = []
    for c in cmd:
        path = c[-1]
        res = finalize(path,False,False,False)
        if res == 0:
            print ('Path %s : complete'%path)
        else:
            res_cmd.append(res)
    print ('='*80)
    print ('Commands to resubmit')
    cid = []
    for c in res_cmd:
        if debug:
            print (c)
        else:
            p = subprocess.Popen([a for a in c.split(' ') if len(a)!=0],stdout=subprocess.PIPE)
            out, err = p.communicate()
            cid.append(int(re.findall(r'\d+',str(out))[0]))
    if debug:
        with open('skim_resubmit.txt','w') as f:
            f.write('\n'.join(res_cmd))
        print ('Commands saved in skim_resubmit.txt')
    else:
        print ('List of resubmit IDs')
        print (' '.join([str(c) for c in cid]))

if args.remove:
    remove()
if args.send:
    send()
if args.id:
    getSlurmID()
if args.resubmit:
    resubmit(args.debug)
