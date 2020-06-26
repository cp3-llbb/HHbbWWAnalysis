import os
import sys
import subprocess
import argparse
import pandas as pd
import time
import curses


parser = argparse.ArgumentParser(description='Utility to monirot jobs')
parser.add_argument('-u','--user', action='store', required=True, type=str, 
                    help='User name')
parser.add_argument('-j','--jobid', action='store', required=False, nargs='+', type=str,
                    help='Job ID (can be several)')
parser.add_argument('-l','--loop', action='store', required=False, type=int,
                    help='Looping time in seconds')
parser.add_argument('--active', action='store_true', required=False, default=False,
                    help='Wether to only display active jobs (with RUNNING or PENDING jobs)')

args = parser.parse_args()

def report_progress(stdscr,printout):
    stdscr.clear()
    for i,line in enumerate(printout.splitlines()):
        try:
            stdscr.addstr(i,0,line)
        except curses.error:
            pass
    stdscr.refresh()

def parse_sacct(args):
    cmd = ['sacct', '-u', args.user, '--format=jobid%25,State,Partition', '-X', '--noheader'] 
    if args.jobid:
        cmd += ['-j',','.join(args.jobid)]

    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    out, _ = p.communicate()

    dict_jobs = None
    for line in out.decode("utf-8").splitlines():
        j, s, p = line.split()
        if '_' in j:
            j, a = j.split('_')
        else:
            a = 1
        if dict_jobs is None:
            dict_jobs = {'jobid':[j],'array':[a],'status':[s],'partition':[p]}
        else:
            dict_jobs['jobid'].append(j)
            dict_jobs['array'].append(a)
            dict_jobs['status'].append(s)
            dict_jobs['partition'].append(p)

    df_jobs = pd.DataFrame.from_dict(dict_jobs)

    status_list = ['PENDING','RUNNING','COMPLETED']

    printout = ""
    for jobid in pd.unique(df_jobs['jobid']):
        jobiddf = df_jobs[df_jobs['jobid']==jobid]
        if args.active and 'PENDING' not in pd.unique(jobiddf['status']) and 'RUNNING' not in pd.unique(jobiddf['status']):
            continue
        printout += "Job ID : %s\n"%jobid
        for partition in pd.unique(jobiddf['partition']):
            partdf = jobiddf[jobiddf['partition']==partition]
            statuses = [l for l in status_list if l in pd.unique(partdf['status'])]+sorted([l for l in pd.unique(partdf['status']) if l not in status_list])
            printout += "\t %s\n"%partition
            N = partdf.shape[0]
            for status in statuses:
                num = partdf[partdf['status']==status].shape[0]
                if status == "PENDING": # Need to develop array line
                    statdf = partdf[partdf['status']==status]
                    arrays = pd.unique(statdf['array'])
                    for array in arrays:
                        if "[" in array: 
                            ps = subprocess.Popen(['squeue','-j',jobid,'-r','--noheader','--long'], stdout=subprocess.PIPE)
                            outs,_ = ps.communicate()
                            sup = len([o for o in outs.decode("utf-8").splitlines() if 'PENDING' in o and partition in o])
                            num += sup -1
                            N += sup -1
                try:
                    printout += "\t\t"+status.ljust(20,' ')+"{:5d}  [{:6.2f}%]\n".format(num,num/N*100)
                except ZeroDivisionError:
                    printout += "\t\t"+status.ljust(20,' ')+"{:5d}\n".format(num)
            printout += "\t\t"+'TOTAL'.ljust(20,' ')+"%5d\n"%N
    return printout


def main(args):
    # Single printout #
    if not args.loop:
        printout = parse_sacct(args)
        print (printout)
    # Loop display #
    else:
        try:
            # Init #
            stdscr = curses.initscr()
            curses.noecho()
            curses.cbreak()
            stdscr.keypad(True)
            # Loop #
            while True:
                printout = parse_sacct(args)
                if len(printout) != 0:
                    report_progress(stdscr,printout)
                time.sleep(args.loop)
            else:
                stdscr.clear()
                stdscr.refresh()
        except KeyboardInterrupt:
            print ("User stop")
            curses.nocbreak()
            stdscr.keypad(False)
            curses.echo()
            curses.endwin()

if __name__ == '__main__':
    main(args)
