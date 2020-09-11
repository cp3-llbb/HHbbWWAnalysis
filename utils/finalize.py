import os
import sys
import subprocess
import glob
import argparse
from pprint import pprint

def finalize(path,force=False,verbose=False,onlyFailed=False):
    path = os.path.abspath(path)

    with open(os.path.join(path,'batch','input','cluster_id'),'r') as f:
        cluster_id = str(int(f.readline()))
    p = subprocess.Popen(['sacct', '-j', cluster_id, '--format=jobid%25,State%25', '-X', '--noheader'], stdout=subprocess.PIPE)
    out, _ = p.communicate()
    status = {}
    for line in out.decode("utf-8").splitlines():
        n, s = line.split()[:2]
        n = n.replace(cluster_id+"_","")
        if '-' in n: # Sometimes sacct returns range [..-..]
            bounds = n.replace('[','').replace(']','').split('-')
            for i in range(int(bounds[0]),int(bounds[1])+1):
                status[str(i)] = s
        else:
            status[n] = s

    state = ["COMPLETED","FAILED","TIMEOUT","OUT_OF_MEMORY","CANCELLED"]
    for s in state:
        print (("   Status %s"%s).ljust(30,' ')+"%4d jobs / %d"%(sum(v==s for v in status.values()),len(status.keys())))

    if not onlyFailed:
        content = {}
        for num in status.keys():
            p = os.path.join(path,'batch','output',num)
            if verbose:
                print ('Looking at %s'%p)
            if not os.path.exists(p):
                content[p] = None
            else:
                #content[p] = os.listdir(p)[0] if len(os.listdir(p))!=0 else ''
                content[p] = os.listdir(p) if len(os.listdir(p))!=0 else None
            if verbose:
                print ('... Found',content[p],'[%s/%s]'%(num,len(status.keys())))
    else:
        content = {os.path.basename(k):None for k,v in status.items() if v in  ['FAILED','TIMEOUT','OUT_OF_MEMORY']}
        print ('Following command will only be for failed jobs')

    if None in content.values():
        print ("Some outputs are missing")
        sbatch_cmd = "sbatch --array="
        list_num = [os.path.basename(k) for k,v in content.items() if v is None]
        sbatch_cmd += ','.join(list_num)+" "

        with open(os.path.expanduser("~/.config/bamboorc"),'r') as f:
            for line in f:
                if line.startswith('sbatch_partition'):
                    sbatch_cmd += "--partition=%s "%line.split()[-1]
                if line.startswith('sbatch_qos'):
                    sbatch_cmd += "--qos=%s "%line.split()[-1]
                if line.startswith('sbatch_time'):
                    sbatch_cmd += "--time=%s "%line.split()[-1]
                if line.startswith('sbatch_additionalOptions'):
                    sbatch_cmd += ' '.join(line.replace('sbatch_additionalOptions =','').split(', ')).replace('\n','')+' '
                if line.startswith('sbatch_mem'):
                    sbatch_cmd += "--mem-per-cpu=%s "%line.split()[-1]
        sbatch_cmd += " "+os.path.join(path,'batch','slurmSubmission.sh')
        print ("Command to resubmit them is below")
        print (sbatch_cmd)
        return sbatch_cmd
    else:
        print ("All the outputs are present, will hadd them now") 
        samples = sorted(list(set([item for sublist in content.values() for item in sublist])))
        if force:
            print ("Careful ! Force hadding the output")
            content  = {k:v for k,v in content.items() if v is not None}
        else:
            samples = [sample for sample in samples if sample not in [os.path.basename(f) for f in glob.glob(os.path.join(path,'results','*.root'))]]
            if len(samples) == 0:
                print ('All root files are in results dir, if you want to force hadd, use --force')
                return 0
        for sample in samples:
            if sample is '':
                continue
            list_sample = [os.path.join(path,f) for path,files in content.items() for f in files if f==sample]
            cmd = ['hadd','-f',os.path.join(path,'results',sample)]+list_sample
            if verbose:
                print ("Command : ",' '.join(cmd))
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            out, _ = p.communicate()
            print ("   Produced file %s"%os.path.join(path,'results',sample))
        return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Utility finalize bamboo')
    parser.add_argument('--path', action='store', required=True, type=str, 
                        help='Path to bamboo output')
    parser.add_argument('--failed', action='store_true', required=False, default=False,
                        help='Will only resubmit the jobs that have failed (script, time out or memory exhaustion)')
    parser.add_argument('--force', action='store_true', required=False, default=False,
                        help='Force hadding, dangerous')
    parser.add_argument('--verbose', action='store_true', required=False, default=False,
                        help='Verbose mode')
    args = parser.parse_args()

    finalize(args.path,args.force,args.verbose,args.failed)
