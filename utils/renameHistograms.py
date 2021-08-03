import os
import sys
import glob 
import argparse
import ROOT
import enlighten
import multiprocessing as mp

from IPython import embed

ROOT.gROOT.ProcessLine('#include "{}"'.format(os.path.join(os.path.dirname(__file__),'rename_header.h')))

def rename_in_file(f,input,output,debug=False):
    ROOT.rename(f,input,output,debug)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Renaming root files histograms')
    parser.add_argument('input', action='store',  type=str,
                        help='Input string (to be replaced)')
    parser.add_argument('output', action='store', type=str,
                        help='Output string (used as replacement)')
    parser.add_argument('regex', action='store',  type=str,
                        help='Path to root files (can be a regex, then use quotes)')
    parser.add_argument('--debug', action='store_true', required=False, default=False,
                        help='Print without changing the names')
    parser.add_argument('-j', dest='jobs',action='store', required=False, default=None, type=int,
                        help='Number of processes to spawn')
    args = parser.parse_args()


    rootFiles = glob.glob(args.regex)
    if len(rootFiles) == 0:
        raise RuntimeError('Regex returned no file')
    print (f'Looking at {len(rootFiles)} root files')
    pbar = enlighten.Counter(total=len(rootFiles), desc='Progress', unit='files')
    if args.jobs is None:
        for rootFile in rootFiles:
            rename_in_file(rootFile,args.input,args.output)
            pbar.update()
    else:
        with mp.Pool(processes=args.jobs) as pool:
            plotIt_configs = pool.starmap(rename_in_file,[(f,args.input,args.output,args.debug) 
                                                        for f in rootFiles])




