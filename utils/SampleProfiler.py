import os
import re
import sys
import glob
import subprocess
from prettytable import PrettyTable

if len(sys.argv) != 2:
    raise RuntimeError('You need to provide the path to a bamboo output directory as argument')

main_path = os.path.abspath(sys.argv[1])
if not os.path.exists(main_path):
    raise RuntimeError('The path your provided does not exist')

class InfoJob:
    def __init__(self,slurmid,arrayid,state,time):
        self.slurmid = slurmid
        self.arrayid = arrayid
        self.state   = state
        self.time    = time
        self.totalid = self.slurmid+'_'+self.arrayid

    def convert_time(self):
        assert len(self.time) == 8
        h,m,s = self.time.split(':')
        return int(h)*3600+int(m)*60+int(s)
    def __repr__(self):
        return 'Job ID : {} - array id : {:4s} - State : {:15s} - Elapsed time : {} = {:6d} s'.format(self.slurmid,self.arrayid,self.state,self.time,self.convert_time())

class InfoLog:
    def __init__(self,sample,exitcode,jobid,arrayid,def_time,def_mem,fill_time,fill_mem):
        self.sample    = sample
        self.exitcode  = exitcode
        self.jobid     = jobid
        self.arrayid   = arrayid
        self.def_time  = def_time
        self.def_mem   = def_mem
        self.fill_time = fill_time
        self.fill_mem  = def_mem

def convert_time_to_str(s):
    s = int(s)
    h = s // 3600
    s -= h * 3600
    m = s // 60
    s -= m* 60
    return "{:02d}:{:02d}:{:02d}".format(h,m,s)


def inspectLogs(path):
    logs = glob.glob(os.path.join(path,"*.out"))
    def findTime(s):
        return float(re.findall(r'\d+\.\d+s',s)[0].replace('s',''))
    def findMem(s):
        return float(re.findall(r'\d+\.\d+MB',s)[0].replace('MB',''))
    infoLogs = []
    for log in logs:
        logbase = os.path.basename(log)
        slurm_id = logbase.replace('slurm-','').replace('.out','')
        assert '_' in slurm_id
        sample    = None
        def_time  = None
        def_mem   = None
        fill_time = None
        fill_mem  = None
        exitcode  = None
        with open(log,'r') as f:
            for line in f:
                if line.startswith('taskcmd = bambooRun'):
                    for l in line.split(' '):
                        if '--sample' in l:
                            sample = (l.split('=')[1]).replace('\n','')
                if line.startswith('INFO:bamboo.analysismodules'):
                    if 'plots defined' in line:
                        def_time = findTime(line)
                        def_mem = findMem(line)
                    if 'Plots finished' in line:
                        fill_time = findTime(line)
                        fill_mem = findMem(line)
                if line.startswith('Payload exit code: '):
                    exitcode = int(re.findall('\d+',line)[0])
                    
        if sample is None:
            print(f'[WARNING] Could not figure out sample from log {logbase}')
        if exitcode is None:
            print(f'[WARNING] Exit code not found from log {logbase}')
        infoLogs.append(InfoLog(sample,exitcode,slurm_id.split('_')[0],slurm_id.split('_')[1],def_time,def_mem,fill_time,fill_mem))

    return infoLogs

def getSacctInfo(slurmids):
    p = subprocess.Popen(['sacct', '-j', ','.join(slurmids), '--format=jobid%25,State%25,Elapsed', '-X', '--noheader'], stdout=subprocess.PIPE)
    out, _ = p.communicate()
    list_job = []
    for line in out.decode("utf-8").splitlines():
        linesplit = line.split()
        sid = linesplit[0]
        state = linesplit[1]
        time = linesplit[-1]
        list_job.append(InfoJob(sid.split('_')[0],sid.split('_')[1],state,time))
    return list_job

def inspectInFiles(path):
    sampleFiles = dict()
    sampleJobs  = dict()
    for infile in glob.glob(os.path.join(path,'*txt')):
        sample = os.path.basename(infile).split('_in_')[0]
        if sample not in sampleFiles.keys():
            sampleFiles[sample] = 0
            sampleJobs[sample] = 1
        else:
            sampleJobs[sample] += 1
        with open(infile,'r') as handle:
            sampleFiles[sample] += sum([1 for line in handle])

    return sampleFiles,sampleJobs

def inspectSubmissionScript(path):
    sampleIds = dict()
    arrayid = 1
    with open(path,'r') as handle:
        for line in handle:
            if 'bambooRun' in line:
                for l in line.split(' '):
                    if '--sample' in l:
                        sample = (l.split('=')[1]).replace('\n','').replace('"','')
                        if sample not in sampleIds.keys():
                            sampleIds[sample] = [str(arrayid)]
                        else:
                            sampleIds[sample].append(str(arrayid))
                        arrayid += 1
    return sampleIds

sampleFiles,sampleJobs = inspectInFiles(os.path.join(main_path,'infiles'))
sampleIds = inspectSubmissionScript(os.path.join(main_path,'batch','slurmSubmission.sh'))

logSummary = {}
slurmIds = set()
for log in inspectLogs(os.path.join(main_path,'batch','logs')):
    if log.arrayid not in logSummary.keys():
        logSummary[log.arrayid] = [log]
    else:
        logSummary[log.arrayid].append(log)
    slurmIds.add(log.jobid)

print ('Find following slurm IDs in the logs :')
for slurmId in slurmIds:
    print (f'... {slurmId}')

        
sacctSummary = {}
for job in getSacctInfo(slurmIds):
    if job.arrayid not in sacctSummary.keys():
        sacctSummary[job.arrayid] = [job]
    else:
        sacctSummary[job.arrayid].append(job)
    
    
samples = sorted(list(sampleFiles.keys()))
table = PrettyTable(["Sample","Files","Jobs","Run time [hh:mm:ss]","Define time [hh:mm:ss]","Fill time [hh:mm:ss]","Define memory [MB]","Fill memory [MB]"])
table.add_row(["","","fin / tot","Min < Mean < Max","Min < Mean < Max","Min < Mean < Max","Min < Mean < Max","Min < Mean < Max",])
invalid_arrayids = []
for sample in samples:
    arrayids = sampleIds[sample]
    values_per_job = {'run_time'  : {},
                      'def_time'  : {},
                      'fill_time' : {},
                      'def_mem'  : {},
                      'fill_mem' : {}}

    for arrayid in arrayids:
        for key in values_per_job.keys():
            values_per_job[key][arrayid] = -1
        for job in sacctSummary[arrayid]:
            if job.state == 'COMPLETED':
                values_per_job['run_time'][arrayid] = max(values_per_job['run_time'][arrayid],job.convert_time())
        if arrayid in logSummary.keys():
            for log in logSummary[arrayid]:
                if log.exitcode == 0:
                    for key in values_per_job.keys():
                        if key == 'run_time':
                            continue
                        values_per_job[key][arrayid] = max(values_per_job[key][arrayid],getattr(log,key))

    invalid_jobs = [arrayid for arrayid,time in values_per_job['run_time'].items() if time==-1]
    invalid_arrayids.extend(invalid_jobs)

    values = {}
    for key in values_per_job.keys():
        values[key] = [val for arrayid,val in values_per_job[key].items() if arrayid not in invalid_jobs]
        if len(values[key]) == 0:
            values[key] = [0]
    N_files = sampleFiles[sample]
    N_jobs  = sampleJobs[sample]
    N_success = N_jobs-len(invalid_jobs)
    line = [sample,N_files,f"{N_success:3d} / {N_jobs:3d}"]
    for key,val in values.items():
        if 'time' in key:
            line.append("{} < {} < {}".format(convert_time_to_str(min(val)),convert_time_to_str(sum(val)/max(1,len(val))),convert_time_to_str(max(val))))
        else:
            line.append("{:8.2f} < {:8.2f} < {:8.2f}".format(min(val),sum(val)/max(1,len(val)),max(val)))
    table.add_row(line)

table.align = "c"
table.align['Sample'] = "l"
print (table)
if len(invalid_arrayids) > 0:
    print ('Resubmit jobs with')
    cmd = f"sbatch --licenses=cms_storage:3 --array {','.join(invalid_arrayids)} {os.path.join(main_path,'batch','slurmSubmission.sh')}"
    print (cmd)
