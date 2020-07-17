import os
import sys
import subprocess
import shutil
import argparse
import time
import ROOT


base_cmd = ['bambooRun','--distributed=driver','-m','SkimmerHHtobbWWDL.py:SkimmerNanoHHtobbWWDL','analysis{era}_synchro_HH.yml','--Synchronization']

add_args = [['--outputTreeName' , 'syncTree'],
            ['--outputTreeName' , 'syncTree_hhbb2l_SR' , '--Tight' , '--Resolved1Btag' , '--Channel', 'ElEl'],
            ['--outputTreeName' , 'syncTree_hhbb2l_SR' , '--Tight' , '--Resolved1Btag' , '--Channel', 'MuMu'],
            ['--outputTreeName' , 'syncTree_hhbb2l_SR' , '--Tight' , '--Resolved1Btag' , '--Channel', 'ElMu'],
            ['--outputTreeName' , 'syncTree_hhbb2l_SR' , '--Tight' , '--Resolved2Btag' , '--Channel', 'ElEl'],
            ['--outputTreeName' , 'syncTree_hhbb2l_SR' , '--Tight' , '--Resolved2Btag' , '--Channel', 'MuMu'],
            ['--outputTreeName' , 'syncTree_hhbb2l_SR' , '--Tight' , '--Resolved2Btag' , '--Channel', 'ElMu'],
            ['--outputTreeName' , 'syncTree_hhbb2l_SR' , '--Tight' , '--Boosted' , '--Channel', 'ElEl'],
            ['--outputTreeName' , 'syncTree_hhbb2l_SR' , '--Tight' , '--Boosted' , '--Channel', 'MuMu'],
            ['--outputTreeName' , 'syncTree_hhbb2l_SR' , '--Tight' , '--Boosted' , '--Channel', 'ElMu'],
            ['--outputTreeName' , 'syncTree_hhbb2l_Fake' , '--FakeExtrapolation' , '--Resolved1Btag' , '--Channel', 'ElEl'],
            ['--outputTreeName' , 'syncTree_hhbb2l_Fake' , '--FakeExtrapolation' , '--Resolved1Btag' , '--Channel', 'MuMu'],
            ['--outputTreeName' , 'syncTree_hhbb2l_Fake' , '--FakeExtrapolation' , '--Resolved1Btag' , '--Channel', 'ElMu'],
            ['--outputTreeName' , 'syncTree_hhbb2l_Fake' , '--FakeExtrapolation' , '--Resolved2Btag' , '--Channel', 'ElEl'],
            ['--outputTreeName' , 'syncTree_hhbb2l_Fake' , '--FakeExtrapolation' , '--Resolved2Btag' , '--Channel', 'MuMu'],
            ['--outputTreeName' , 'syncTree_hhbb2l_Fake' , '--FakeExtrapolation' , '--Resolved2Btag' , '--Channel', 'ElMu'],
            ['--outputTreeName' , 'syncTree_hhbb2l_Fake' , '--FakeExtrapolation' , '--Boosted' , '--Channel', 'ElEl'],
            ['--outputTreeName' , 'syncTree_hhbb2l_Fake' , '--FakeExtrapolation' , '--Boosted' , '--Channel', 'MuMu'],
            ['--outputTreeName' , 'syncTree_hhbb2l_Fake' , '--FakeExtrapolation' , '--Boosted' , '--Channel', 'ElMu']]


path_args = [['-o', '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SynchronizationHH_{era}_noSel'],
             ['-o', '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SynchronizationHH_{era}_hhbb2l_SR_Res1b_ElEl'],
             ['-o', '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SynchronizationHH_{era}_hhbb2l_SR_Res1b_MuMu'],
             ['-o', '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SynchronizationHH_{era}_hhbb2l_SR_Res1b_ElMu'],
             ['-o', '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SynchronizationHH_{era}_hhbb2l_SR_Res2b_ElEl'],
             ['-o', '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SynchronizationHH_{era}_hhbb2l_SR_Res2b_MuMu'],
             ['-o', '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SynchronizationHH_{era}_hhbb2l_SR_Res2b_ElMu'],
             ['-o', '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SynchronizationHH_{era}_hhbb2l_SR_Boosted_ElEl'],
             ['-o', '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SynchronizationHH_{era}_hhbb2l_SR_Boosted_MuMu'],
             ['-o', '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SynchronizationHH_{era}_hhbb2l_SR_Boosted_ElMu'],
             ['-o', '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SynchronizationHH_{era}_hhbb2l_Fake_Res1b_ElEl'],
             ['-o', '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SynchronizationHH_{era}_hhbb2l_Fake_Res1b_MuMu'],
             ['-o', '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SynchronizationHH_{era}_hhbb2l_Fake_Res1b_ElMu'],
             ['-o', '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SynchronizationHH_{era}_hhbb2l_Fake_Res2b_ElEl'],
             ['-o', '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SynchronizationHH_{era}_hhbb2l_Fake_Res2b_MuMu'],
             ['-o', '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SynchronizationHH_{era}_hhbb2l_Fake_Res2b_ElMu'],
             ['-o', '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SynchronizationHH_{era}_hhbb2l_Fake_Boosted_ElEl'],
             ['-o', '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SynchronizationHH_{era}_hhbb2l_Fake_Boosted_MuMu'],
             ['-o', '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/SynchronizationHH_{era}_hhbb2l_Fake_Boosted_ElMu']]

def remove():
    for _,path in path_args:
        if os.path.exists(path):
            shutil.rmtree(path)
            print ("Removed %s"%path)


def send():
    for add_arg,path_arg in zip(add_args,path_args):
        cmd = base_cmd + path_arg + add_arg
        subprocess.Popen(cmd)
        #print (' '.join(cmd))
        


def check():
    for idx in range(len(path_args)):
        full_path = os.path.join(path_args[idx][1],'results','GluGluToRadionToHHTo2B2VTo2L2Nu_M-750.root')
        cmd = ''
        if os.path.exists(full_path):
            f = ROOT.TFile(full_path)
            if not f.GetListOfKeys().Contains(add_args[idx][1]):
                print ("File %s does not contain %s tree"%(path_args[idx][1],add_args[idx][1]))
                cmd = base_cmd + path_args[idx] + add_args[idx]
        else:
            print ('Path %s is missing'%full_path)
            cmd = base_cmd + path_args[idx] + add_args[idx]
        if cmd != '':
            print (' '.join(cmd))
            #subprocess.Popen(cmd)
            #time.sleep(0.5)

def merge(era):
    path_merge = os.path.join(os.path.dirname(__file__),'merge_GluGluToRadionToHHTo2B2VTo2L2Nu_M-750_{}.root'.format(era))
    path_tmp = os.path.join(os.path.dirname(__file__),'tmp_GluGluToRadionToHHTo2B2VTo2L2Nu_M-750_{}.root'.format(era))
    if os.path.exists(path_merge):
        os.remove(path_merge)
    if os.path.exists(path_tmp):
        os.remove(path_tmp)
    f = ROOT.TFile(path_merge,"recreate")
    f.Close()

    for add_arg,path_arg in zip(add_args,path_args):
        treeName = add_arg[1]
        path = os.path.join(path_arg[1],'results','GluGluToRadionToHHTo2B2VTo2L2Nu_M-750.root')
        if not os.path.exists(path):
            continue
        cmd = ['hadd',path_tmp,path_merge,path]
        subprocess.run(cmd)
        shutil.move(path_tmp,path_merge)

    print ("Merged file:",path_merge)

def concatenate(era):
    path_merge = os.path.join(os.path.dirname(__file__),'concatenate_GluGluToRadionToHHTo2B2VTo2L2Nu_M-750_{}.root'.format(era))
    f = ROOT.TFile(path_merge,"recreate")

    for add_arg,path_arg in zip(add_args,path_args):
        treeName = add_arg[1]
        if 'SR' in treeName or 'Fake' in treeName:
            newTreeName = 'syncTree' + os.path.basename(path_arg[1]).replace("SynchronizationHH",'')
        else:
            newTreeName = treeName
        path = os.path.join(path_arg[1],'results','GluGluToRadionToHHTo2B2VTo2L2Nu_M-750.root')
        cmd = ['rootcp', '{}:{}'.format(path,treeName) , '{}:{}'.format(path_merge,newTreeName)]
        print (' '.join(cmd))
        #subprocess.run(cmd)
        #subprocess.Popen(['rootls',path_merge])
    print ('Concatenated file:',path_merge)
    f.Close()

def selectEra(args,era):
    l = []
    for arg in args:
        if isinstance(arg,list):
            sl = []
            for a in arg:
                if '{era}' in a:
                    a = a.format(era=era)
                sl.append(a)
            l.append(sl)
        else:
            if '{era}' in arg:
                arg = arg.format(era=era)
            l.append(arg)
    return l
    
    


parser = argparse.ArgumentParser(description='Integration with Cuba')
parser.add_argument('--remove', action='store_true', required=False, default=False,
                    help='Clears the sync directories')
parser.add_argument('--send', action='store_true', required=False, default=False,
                    help='Sends the jobs')
parser.add_argument('--check', action='store_true', required=False, default=False,
                    help='Checks if all jobs have finished, if output not present will resubmit') 
parser.add_argument('--merge', action='store_true', required=False, default=False,
                    help='Merge channel and jet selections in single root file')
parser.add_argument('--concatenate', action='store_true', required=False, default=False,
                    help='Concatenate channel and jet selections in single root file withtout merging')
parser.add_argument('--era', action='store', required=True, type=str,
                    help='Era')
args = parser.parse_args()

base_cmd = selectEra(base_cmd,args.era)
add_args = selectEra(add_args,args.era)
path_args = selectEra(path_args,args.era)

if args.remove:
    remove()
if args.check:
    check()
if args.send:
    send()
if args.merge:
    merge(args.era)
if args.concatenate:
    concatenate(args.era)
