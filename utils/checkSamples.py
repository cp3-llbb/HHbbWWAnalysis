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
    n_total = len(subprocess.check_output(cmd).decode().split())
    cmd = ['dasgoclient', '-query', 'file dataset=%s site=T2_BE_UCL'%sample]
    n_local = len(subprocess.check_output(cmd).decode().split())
    print ("Sample : %s"%sample)
    print ("... Local : %d / %d"%(n_local,n_total))
    if n_total<n_local:
        raise RuntimeError("More files locally than globally, how is that possible ?!")
    return n_total == n_local

def saveToText(list_samples,filename):
    with open(filename,'w') as f:
        for sample in list_samples:
            f.write(sample+'\n')
    print('Saved text file as %s'%filename)
     
def processAndSave(txtfile,filename,version): 
    list_mini = readTextFile(txtfile)
    list_nano = getNanoAODList(list_mini,version)
    list_missing_nano = [nano for nano in list_nano if not checkSample(nano)]
    saveToText(list_nano,filename)
    saveToText(list_missing_nano,filename.replace('.txt','_missing.txt'))

processAndSave('signal_miniaod_samples_2016.txt','signal_nanoaod_samples_2016.txt','NanoAODv6')
processAndSave('signal_miniaod_samples_2017.txt','signal_nanoaod_samples_2017.txt','NanoAODv6')
processAndSave('signal_miniaod_samples_2018.txt','signal_nanoaod_samples_2018.txt','NanoAODv6')
processAndSave('background_miniaod_samples_2016.txt','background_nanoaod_samples_2016.txt','NanoAODv6')
processAndSave('background_miniaod_samples_2017.txt','background_nanoaod_samples_2017.txt','NanoAODv6')
processAndSave('background_miniaod_samples_2018.txt','background_nanoaod_samples_2018.txt','NanoAODv6')
processAndSave('data_miniaod_samples_2016.txt','data_nanoaod_samples_2016.txt','Nano25Oct2019')
processAndSave('data_miniaod_samples_2017.txt','data_nanoaod_samples_2017.txt','Nano25Oct2019')
processAndSave('data_miniaod_samples_2018.txt','data_nanoaod_samples_2018.txt','Nano25Oct2019')
