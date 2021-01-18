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

cmd = [['bambooRun','--TTHIDLoose','--distributed=driver','-m','SkimmerHHtobbWWSL.py:SkimmerNanoHHtobbWWSL',f'analysis{era}_v7_MC.yml','--Res2b2Wj','--Channel','El','-o',f'/nfs/scratch/fynu/gsaha/SkimmedNtuplesJPA/Skim{era}_SLBackground_Resolved2b2Wj_El'],
       ['bambooRun','--TTHIDLoose','--distributed=driver','-m','SkimmerHHtobbWWSL.py:SkimmerNanoHHtobbWWSL',f'analysis{era}_v7_MC.yml','--Res2b2Wj','--Channel','Mu','-o',f'/nfs/scratch/fynu/gsaha/SkimmedNtuplesJPA/Skim{era}_SLBackground_Resolved2b2Wj_Mu'],
       ['bambooRun','--TTHIDLoose','--distributed=driver','-m','SkimmerHHtobbWWSL.py:SkimmerNanoHHtobbWWSL',f'analysis{era}_v7_MC.yml','--Res2b1Wj','--Channel','El','-o',f'/nfs/scratch/fynu/gsaha/SkimmedNtuplesJPA/Skim{era}_SLBackground_Resolved2b1Wj_El'],
       ['bambooRun','--TTHIDLoose','--distributed=driver','-m','SkimmerHHtobbWWSL.py:SkimmerNanoHHtobbWWSL',f'analysis{era}_v7_MC.yml','--Res2b1Wj','--Channel','Mu','-o',f'/nfs/scratch/fynu/gsaha/SkimmedNtuplesJPA/Skim{era}_SLBackground_Resolved2b1Wj_Mu'],
       ['bambooRun','--TTHIDLoose','--distributed=driver','-m','SkimmerHHtobbWWSL.py:SkimmerNanoHHtobbWWSL',f'analysis{era}_v7_MC.yml','--Res2b0Wj','--Channel','El','-o',f'/nfs/scratch/fynu/gsaha/SkimmedNtuplesJPA/Skim{era}_SLBackground_Resolved2b0Wj_El'],
       ['bambooRun','--TTHIDLoose','--distributed=driver','-m','SkimmerHHtobbWWSL.py:SkimmerNanoHHtobbWWSL',f'analysis{era}_v7_MC.yml','--Res2b0Wj','--Channel','Mu','-o',f'/nfs/scratch/fynu/gsaha/SkimmedNtuplesJPA/Skim{era}_SLBackground_Resolved2b0Wj_Mu'],
       ['bambooRun','--TTHIDLoose','--distributed=driver','-m','SkimmerHHtobbWWSL.py:SkimmerNanoHHtobbWWSL',f'analysis{era}_v7_MC.yml','--Res1b2Wj','--Channel','El','-o',f'/nfs/scratch/fynu/gsaha/SkimmedNtuplesJPA/Skim{era}_SLBackground_Resolved1b2Wj_El'],
       ['bambooRun','--TTHIDLoose','--distributed=driver','-m','SkimmerHHtobbWWSL.py:SkimmerNanoHHtobbWWSL',f'analysis{era}_v7_MC.yml','--Res1b2Wj','--Channel','Mu','-o',f'/nfs/scratch/fynu/gsaha/SkimmedNtuplesJPA/Skim{era}_SLBackground_Resolved1b2Wj_Mu'],
       ['bambooRun','--TTHIDLoose','--distributed=driver','-m','SkimmerHHtobbWWSL.py:SkimmerNanoHHtobbWWSL',f'analysis{era}_v7_MC.yml','--Res1b1Wj','--Channel','El','-o',f'/nfs/scratch/fynu/gsaha/SkimmedNtuplesJPA/Skim{era}_SLBackground_Resolved1b1Wj_El'],
       ['bambooRun','--TTHIDLoose','--distributed=driver','-m','SkimmerHHtobbWWSL.py:SkimmerNanoHHtobbWWSL',f'analysis{era}_v7_MC.yml','--Res1b1Wj','--Channel','Mu','-o',f'/nfs/scratch/fynu/gsaha/SkimmedNtuplesJPA/Skim{era}_SLBackground_Resolved1b1Wj_Mu'],
       ['bambooRun','--TTHIDLoose','--distributed=driver','-m','SkimmerHHtobbWWSL.py:SkimmerNanoHHtobbWWSL',f'analysis{era}_v7_MC.yml','--Res1b0Wj','--Channel','El','-o',f'/nfs/scratch/fynu/gsaha/SkimmedNtuplesJPA/Skim{era}_SLBackground_Resolved1b0Wj_El'],
       ['bambooRun','--TTHIDLoose','--distributed=driver','-m','SkimmerHHtobbWWSL.py:SkimmerNanoHHtobbWWSL',f'analysis{era}_v7_MC.yml','--Res1b0Wj','--Channel','Mu','-o',f'/nfs/scratch/fynu/gsaha/SkimmedNtuplesJPA/Skim{era}_SLBackground_Resolved1b0Wj_Mu'],
       ['bambooRun','--TTHIDLoose','--distributed=driver','-m','SkimmerHHtobbWWSL.py:SkimmerNanoHHtobbWWSL',f'analysis{era}_v7_MC.yml','--Res0b','--Channel','El','-o',f'/nfs/scratch/fynu/gsaha/SkimmedNtuplesJPA/Skim{era}_SLBackground_Resolved0b_El'],
       ['bambooRun','--TTHIDLoose','--distributed=driver','-m','SkimmerHHtobbWWSL.py:SkimmerNanoHHtobbWWSL',f'analysis{era}_v7_MC.yml','--Res0b','--Channel','Mu','-o',f'/nfs/scratch/fynu/gsaha/SkimmedNtuplesJPA/Skim{era}_SLBackground_Resolved0b_Mu'],
       ['bambooRun','--TTHIDLoose','--distributed=driver','-m','SkimmerHHtobbWWSL.py:SkimmerNanoHHtobbWWSL',f'analysis{era}_v7_signal_SL.yml','--Res2b2Wj','--Channel','El','-o',f'/nfs/scratch/fynu/gsaha/SkimmedNtuplesJPA/Skim{era}_SLSignal_Resolved2b2Wj_El'],
       ['bambooRun','--TTHIDLoose','--distributed=driver','-m','SkimmerHHtobbWWSL.py:SkimmerNanoHHtobbWWSL',f'analysis{era}_v7_signal_SL.yml','--Res2b2Wj','--Channel','Mu','-o',f'/nfs/scratch/fynu/gsaha/SkimmedNtuplesJPA/Skim{era}_SLSignal_Resolved2b2Wj_Mu'],
       ['bambooRun','--TTHIDLoose','--distributed=driver','-m','SkimmerHHtobbWWSL.py:SkimmerNanoHHtobbWWSL',f'analysis{era}_v7_signal_SL.yml','--Res2b1Wj','--Channel','El','-o',f'/nfs/scratch/fynu/gsaha/SkimmedNtuplesJPA/Skim{era}_SLSignal_Resolved2b1Wj_El'],
       ['bambooRun','--TTHIDLoose','--distributed=driver','-m','SkimmerHHtobbWWSL.py:SkimmerNanoHHtobbWWSL',f'analysis{era}_v7_signal_SL.yml','--Res2b1Wj','--Channel','Mu','-o',f'/nfs/scratch/fynu/gsaha/SkimmedNtuplesJPA/Skim{era}_SLSignal_Resolved2b1Wj_Mu'],
       ['bambooRun','--TTHIDLoose','--distributed=driver','-m','SkimmerHHtobbWWSL.py:SkimmerNanoHHtobbWWSL',f'analysis{era}_v7_signal_SL.yml','--Res2b0Wj','--Channel','El','-o',f'/nfs/scratch/fynu/gsaha/SkimmedNtuplesJPA/Skim{era}_SLSignal_Resolved2b0Wj_El'],
       ['bambooRun','--TTHIDLoose','--distributed=driver','-m','SkimmerHHtobbWWSL.py:SkimmerNanoHHtobbWWSL',f'analysis{era}_v7_signal_SL.yml','--Res2b0Wj','--Channel','Mu','-o',f'/nfs/scratch/fynu/gsaha/SkimmedNtuplesJPA/Skim{era}_SLSignal_Resolved2b0Wj_Mu'],
       ['bambooRun','--TTHIDLoose','--distributed=driver','-m','SkimmerHHtobbWWSL.py:SkimmerNanoHHtobbWWSL',f'analysis{era}_v7_signal_SL.yml','--Res1b2Wj','--Channel','El','-o',f'/nfs/scratch/fynu/gsaha/SkimmedNtuplesJPA/Skim{era}_SLSignal_Resolved1b2Wj_El'],
       ['bambooRun','--TTHIDLoose','--distributed=driver','-m','SkimmerHHtobbWWSL.py:SkimmerNanoHHtobbWWSL',f'analysis{era}_v7_signal_SL.yml','--Res1b2Wj','--Channel','Mu','-o',f'/nfs/scratch/fynu/gsaha/SkimmedNtuplesJPA/Skim{era}_SLSignal_Resolved1b2Wj_Mu'],
       ['bambooRun','--TTHIDLoose','--distributed=driver','-m','SkimmerHHtobbWWSL.py:SkimmerNanoHHtobbWWSL',f'analysis{era}_v7_signal_SL.yml','--Res1b1Wj','--Channel','El','-o',f'/nfs/scratch/fynu/gsaha/SkimmedNtuplesJPA/Skim{era}_SLSignal_Resolved1b1Wj_El'],
       ['bambooRun','--TTHIDLoose','--distributed=driver','-m','SkimmerHHtobbWWSL.py:SkimmerNanoHHtobbWWSL',f'analysis{era}_v7_signal_SL.yml','--Res1b1Wj','--Channel','Mu','-o',f'/nfs/scratch/fynu/gsaha/SkimmedNtuplesJPA/Skim{era}_SLSignal_Resolved1b1Wj_Mu'],
       ['bambooRun','--TTHIDLoose','--distributed=driver','-m','SkimmerHHtobbWWSL.py:SkimmerNanoHHtobbWWSL',f'analysis{era}_v7_signal_SL.yml','--Res1b0Wj','--Channel','El','-o',f'/nfs/scratch/fynu/gsaha/SkimmedNtuplesJPA/Skim{era}_SLSignal_Resolved1b0Wj_El'],
       ['bambooRun','--TTHIDLoose','--distributed=driver','-m','SkimmerHHtobbWWSL.py:SkimmerNanoHHtobbWWSL',f'analysis{era}_v7_signal_SL.yml','--Res1b0Wj','--Channel','Mu','-o',f'/nfs/scratch/fynu/gsaha/SkimmedNtuplesJPA/Skim{era}_SLSignal_Resolved1b0Wj_Mu'],
       ['bambooRun','--TTHIDLoose','--distributed=driver','-m','SkimmerHHtobbWWSL.py:SkimmerNanoHHtobbWWSL',f'analysis{era}_v7_signal_SL.yml','--Res0b','--Channel','El','-o',f'/nfs/scratch/fynu/gsaha/SkimmedNtuplesJPA/Skim{era}_SLSignal_Resolved0b_El'],
       ['bambooRun','--TTHIDLoose','--distributed=driver','-m','SkimmerHHtobbWWSL.py:SkimmerNanoHHtobbWWSL',f'analysis{era}_v7_signal_SL.yml','--Res0b','--Channel','Mu','-o',f'/nfs/scratch/fynu/gsaha/SkimmedNtuplesJPA/Skim{era}_SLSignal_Resolved0b_Mu']]

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
