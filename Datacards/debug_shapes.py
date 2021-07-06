import os
import sys
import glob
import subprocess
import ROOT
from multiprocessing import Process, current_process

paths = [
    'dataset_fit_TTHIDLoose_DNN10_2018_Boosted1B_em_datadriven_syst',
    'dataset_fit_TTHIDLoose_DNN10_2018_Resolved1B_em_datadriven_syst',
    'dataset_fit_TTHIDLoose_DNN10_2018_Resolved2B_em_datadriven_syst',
    'dataset_fit_TTHIDLoose_DNN10_2018_Boosted1B_ee_mm_datadriven_syst',
    'dataset_fit_TTHIDLoose_DNN10_2018_Resolved1B_ee_mm_datadriven_syst',
    'dataset_fit_TTHIDLoose_DNN10_2018_Resolved2B_ee_mm_datadriven_syst',
]   

for path in paths:
    if not os.path.exists(path):
        raise RuntimeError('Directory {} does not exist'.format(path))

outdir = 'debug_shapes'
if not os.path.exists(outdir):
    os.makedirs(outdir)

def callSystematicPlot(path):
    rootFiles = glob.glob(os.path.join(path,'plotit','root','*root'))
    for idx,rootFile in enumerate(rootFiles):
        sample = os.path.basename(rootFile).replace('.root','')
        if 'data' in sample:
            continue
        print ('Process {} - Progress {}/{} [{:3.1f}%]- Root file : {}'.format(current_process(),idx,len(rootFiles),idx*100/len(rootFiles),rootFile))
        F = ROOT.TFile(rootFile)
        for key in F.GetListOfKeys():   
            if '__' in key.GetName():
                continue
            cmd = ['python','SystematicPlot.py','--file',rootFile,'--nominal',key.GetName(),'--output','{}/{}_{}.pdf'.format(outdir,sample,key.GetName()),'--xlabel',key.GetName(),'--rebin','20']
            subprocess.check_output(cmd)
        F.Close()

processes = [Process(target=callSystematicPlot,args=(path,)) for path in paths]

for p in processes:
    p.start()

for p in processes:
    p.join()


