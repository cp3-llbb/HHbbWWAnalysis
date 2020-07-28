import os
import sys
import subprocess
import argparse


def finalize(path,force=True,verbose=False,onlyFailed=False):
    path = os.path.abspath(path)

    with open(os.path.join(path,'batch','input','cluster_id'),'r') as f:
        cluster_id = f.readline()

    p = subprocess.Popen(['sacct', '-j', cluster_id, '--format=jobid%25,State', '-X', '--noheader'], stdout=subprocess.PIPE)
    out, _ = p.communicate()
    status = {}
    for line in out.decode("utf-8").splitlines():
        n, s = line.split()
        n = n.replace(cluster_id+"_","")
        status[n] = s

    state = ["COMPLETED","FAILED","TIMEOUT","OUT_OF_MEMORY","CANCELLED"]
    for s in state:
        print (("   Status %s"%s).ljust(30,' ')+"%4d jobs / %d"%(sum(v==s for v in status.values()),len(status.keys())))

    if not onlyFailed:
        content = {}
        for num in status.keys():
            p = os.path.join(path,'batch','output',num)
            if not os.path.exists(p):
                content[p] = None
            else:
                #content[p] = os.listdir(p)[0] if len(os.listdir(p))!=0 else ''
                content[p] = os.listdir(p)[0] if len(os.listdir(p))!=0 else None
    else:
        content = {os.path.basename(k):None for k,v in status.items() if v in  ['FAILED','TIMEOUT','OUT_OF_MEMORY']}
        print ('Following command will only be for failed jobs')
    #import pprint
    #pprint.pprint(content)

    if None in content.values() and not force:
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
                if line.startswith('sbatch_mem'):
                    sbatch_cmd += "--mem-per-cpu=%s "%line.split()[-1]
                if line.startswith('sbatch_additionalOptions'):
                    sbatch_cmd += ' '.join(line.split()[2:]).replace(',','')+" "
        sbatch_cmd += " "+os.path.join(path,'batch','slurmSubmission.sh')
        print ("Command to resubmit them is below")
        print (sbatch_cmd)
    else:
        print ("All the outputs are present, will hadd them now") 
        if force:
            print ("Careful ! Force hadding the output")
            content  = {k:v for k,v in content.items() if v is not None}
        samples = sorted(list(set(content.values())))
        for sample in samples:
            if sample is '':
                continue
            list_sample = [os.path.join(k,v) for k,v in content.items() if v==sample]
            cmd = ['hadd','-f',os.path.join(path,'results',sample)]+list_sample
            if verbose:
                print ("Command : ",' '.join(cmd))
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            out, _ = p.communicate()
            print ("   Produced file %s"%os.path.join(path,'results',sample))


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
