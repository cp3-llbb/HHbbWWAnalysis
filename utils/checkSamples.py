import os
import sys
import subprocess

def getNanoAODList(daspaths,version):
    nanoAODfiles = []
    for daspath in daspaths:
        cmd = ['dasgoclient', '-query', 'child dataset=%s'%daspath]
        try:
            results = subprocess.check_output(cmd).decode().split()
        except subprocess.CalledProcessError:
            print ("Error in DAS for sample %s"%daspath)
            return []
        if len(results)==0:
            print ("No NanoAOD child of sample %s"%daspath)
        addnano = [result for result in results if version in result]
        if len(addnano) > 1:
            print ('MiniAOD sample %s has several children'%daspath)
            for nano in addnano:
                print ('...',nano)
        if len(addnano)==0:
            print ('Could not find %s child of MiniAOD sample %s'%(version,daspath))
            print ('Available NanoAOD versions')
            for result in results:
                print ('...',result)
        else:
            nanoAODfiles.extend(addnano)
    return nanoAODfiles

def readTextFile(path):
    with open(path, "r") as f:
        list_mini = [line.strip() for line in f if len(line)!=0]
    list_mini = [l for l in list_mini if len(l)!=0]
    return list_mini

def checkSample(sample):
    cmd = ['dasgoclient', '-query', 'file dataset=%s'%sample]
    files = subprocess.check_output(cmd).decode().split()
    missing_files = [f for f in files if not os.path.exists('/storage/data/cms'+f)]
    print ("Sample : %s"%sample)
    print ("... Local : %d / %d"%(len(files)-len(missing_files),len(files)))
    if len(missing_files)>0:
        for f in missing_files:
            print ('... Missing file :',f)
    return len(missing_files) == 0

def saveToText(list_samples,filename):
    with open(filename,'w') as f:
        for sample in list_samples:
            f.write(sample+'\n')
    print('Saved text file as %s'%filename)
     
def processAndSave(miniaod_file=None,nanoaod_file=None,version=None): 
    if miniaod_file is not None:
        list_mini = readTextFile(miniaod_file)
        list_nano = getNanoAODList(list_mini,version)
    else:
        list_nano = readTextFile(nanoaod_file)
    list_missing_nano = [nano for nano in list_nano if not checkSample(nano)]
    rucio_list = ['cms:'+nano for nano in list_nano if nano not in list_missing_nano]
    rucio_cmd = 'rucio add-rule '+' '.join(rucio_list) + ' 1 T2_BE_UCL'
    print ('\nRucio command below')
    print ('Note that you should do the commands before :')
    print ('source /cvmfs/cms.cern.ch/cmsset_default.sh')
    print ('source /cvmfs/cms.cern.ch/rucio/setup.sh')
    print ('export RUCIO_ACCOUNT=your_username')
    print ()
    print (rucio_cmd)
    print ()
    if miniaod_file is not None:
        saveToText(list_nano,nanoaod_file)
    saveToText(list_missing_nano,nanoaod_file.replace('.txt','_missing.txt'))

#processAndSave(miniaod_file='signal_miniaod_samples_2016.txt',nanoaod_file='signal_nanoaod_samples_2016_v6.txt',version='NanoAODv6')
#processAndSave(miniaod_file='signal_miniaod_samples_2017.txt',nanoaod_file='signal_nanoaod_samples_2017_v6.txt',version='NanoAODv6')
#processAndSave(miniaod_file='signal_miniaod_samples_2018.txt',nanoaod_file='signal_nanoaod_samples_2018_v6.txt',version='NanoAODv6')
#processAndSave(miniaod_file='background_miniaod_samples_2016.txt',nanoaod_file='background_nanoaod_samples_2016_v6.txt',version='NanoAODv6')
#processAndSave(miniaod_file='background_miniaod_samples_2017.txt',nanoaod_file='background_nanoaod_samples_2017_v6.txt',version='NanoAODv6')
#processAndSave(miniaod_file='background_miniaod_samples_2018.txt',nanoaod_file='background_nanoaod_samples_2018_v6.txt',version='NanoAODv6')
processAndSave(miniaod_file='all_miniaod_mc.txt',nanoaod_file='all_nanoaod_mc.txt',version='NanoAODv7')
#processAndSave(miniaod_file='all_miniaod_data.txt',nanoaod_file='all_nanoaod_data.txt',version='Nano25Oct2019')
#processAndSave(miniaod_file='data_miniaod_samples_2016.txt',nanoaod_file='data_nanoaod_samples_2016_Nano25Oct2019.txt',version='Nano25Oct2019')
#processAndSave(miniaod_file='data_miniaod_samples_2017.txt',nanoaod_file='data_nanoaod_samples_2017_Nano25Oct2019.txt',version='Nano25Oct2019')
#processAndSave(miniaod_file='data_miniaod_samples_2018.txt',nanoaod_file='data_nanoaod_samples_2018_Nano25Oct2019.txt',version='Nano25Oct2019')
#processAndSave(miniaod_file='signal_miniaod_samples_2016.txt',nanoaod_file='signal_nanoaod_samples_2016_v7.txt',version='NanoAODv7')
#processAndSave(miniaod_file='signal_miniaod_samples_2017.txt',nanoaod_file='signal_nanoaod_samples_2017_v7.txt',version='NanoAODv7')
#processAndSave(miniaod_file='signal_miniaod_samples_2018.txt',nanoaod_file='signal_nanoaod_samples_2018_v7.txt',version='NanoAODv7')
#processAndSave(miniaod_file='background_miniaod_samples_2016.txt',nanoaod_file='background_nanoaod_samples_2016_v7.txt',version='NanoAODv7')
#processAndSave(miniaod_file='background_miniaod_samples_2017.txt',nanoaod_file='background_nanoaod_samples_2017_v7.txt',version='NanoAODv7')
#processAndSave(miniaod_file='background_miniaod_samples_2018.txt',nanoaod_file='background_nanoaod_samples_2018_v7.txt',version='NanoAODv7')
#processAndSave(nanoaod_file='data_nanoaod_samples_02Apr2020.txt')
