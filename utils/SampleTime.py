import os
import sys
import glob
from collections import defaultdict
from pprint import pprint
import subprocess
from prettytable import PrettyTable

if len(sys.argv) != 2:
    print ('You need to provide the path to a bamboo output dir as argument')
    sys.exit(1)

path = os.path.abspath(os.path.join(sys.argv[1],'batch','logs'))

class InfoJob:
    def __init__(self,slurmid,arrayid,state,time):
        self.slurmid = slurmid
        self.arrayid = arrayid
        self.state   = state
        self.time    = time
        self.totalid = self.slurmid+'_'+self.arrayid

    def convert_time(self):
        h,m,s = self.time.split(':')
        return int(h)*3600+int(m)*60+int(s)
    def __repr__(self):
        return 'Job ID : {} - array id : {:4s} - State : {:15s} - Elapsed time : {} = {:6d} s'.format(self.slurmid,self.arrayid,self.state,self.time,self.convert_time())

def convert_time_to_str(s):
    s = int(s)
    h = s // 3600
    s -= h * 3600
    m = s // 60
    s -= m* 60
    return "{:02d}:{:02d}:{:02d}".format(h,m,s)


def getSlurmIdPerSample(path):
    samplesID = defaultdict(list)
    logs = glob.glob(os.path.join(path,"*.out"))
    for log in logs:
        logbase = os.path.basename(log)
        slurm_id = logbase.replace('slurm-','').replace('.out','')
        sample = None
        with open(log,'r') as f:
            for line in f:
                if '--sample' in line:
                    for l in line.split(' '):
                        if '--sample' in l:
                            sample = (l.split('=')[1]).replace('\n','')
        if sample is None:
            raise RuntimeError('Could not figure out sample from log {}'.format(logbase))
    
        samplesID[sample].append(slurm_id)
    return samplesID

def getSacctInfo(slurmid):
    p = subprocess.Popen(['sacct', '-j', slurmid, '--format=jobid%25,State%25,Elapsed', '-X', '--noheader'], stdout=subprocess.PIPE)
    out, _ = p.communicate()
    list_job = []
    for line in out.decode("utf-8").splitlines():
        linesplit = line.split()
        sid = linesplit[0]
        state = linesplit[1]
        time = linesplit[-1]
        list_job.append(InfoJob(sid.split('_')[0],sid.split('_')[1],state,time))

    return list_job


jobid_per_sample = getSlurmIdPerSample(path)

slurmIds = list(set([job.split('_')[0] for jobid in jobid_per_sample.values() for job in jobid]))
sacctSummary = {}
for slurmId in slurmIds:
    sacctSummary.update({jobdid.totalid:jobdid for jobdid in getSacctInfo(slurmId)})

samples = sorted(list(jobid_per_sample.keys()))
sample_time = {}
for sample in samples:
    jobids = sorted(jobid_per_sample[sample])
    time_per_job = {jobid.split('_')[1]:[] for jobid in jobids}
    for jobid in jobids:
        summary = sacctSummary[jobid]
        if summary.state == 'COMPLETED':
            time_per_job[str(summary.arrayid)].append(summary.convert_time())

    invalid_jobs = [job for job,time in time_per_job.items() if len(time)==0]
    if len(invalid_jobs)>0:
        print ("Warning : no completed job for sample {} with array id : ".format(sample)+','.join(invalid_jobs))
    time_per_job = {k:max(t) for k,t in time_per_job.items() if len(t)>0}
    times = [t for j,t in time_per_job.items() if j not in invalid_jobs]
    if len(times)==0:
        times = [0]
    sample_time[sample] = [len(jobids),len(jobids)-len(invalid_jobs),sum(times),sum(times)/max(1,len(times)),min(times),max(times)]
    
table = PrettyTable(["Sample","Sent jobs","Completed jobs","Total time","Mean time","Min time","Max time"])
for sample in samples:
    times = sample_time[sample]
    times_str = [convert_time_to_str(t) for t in times[2:]]
    table.add_row([sample]+times[0:2]+times_str)
table.align = "r"
table.align['Sample'] = "l"
print (table)
