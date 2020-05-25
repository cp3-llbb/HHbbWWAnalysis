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
parser.add_argument('-j','--jobid', action='store', required=False, type=str, 
                    help='Job ID')
parser.add_argument('-l','--loop', action='store', required=False, type=int, 
                    help='Looping time in seconds')

args = parser.parse_args()

def report_progress(stdscr,printout):
    for i,line in enumerate(printout.splitlines()):
        stdscr.addstr(i,0,line)
    stdscr.refresh()

def parse_sacct(args):
    cmd = ['sacct', '-u', args.user, '--format=jobid%25,State,Partition', '-X', '--noheader'] 
    if args.jobid:
        cmd += ['-j',args.jobid]

    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    out, _ = p.communicate()

    dict_jobs = None
    for line in out.decode("utf-8").splitlines():
        j, s, p = line.split()
        j, a = j.split('_')
        if dict_jobs is None:
            dict_jobs = {'jobid':[j],'array':[a],'status':[s],'partition':[p]}
        else:
            dict_jobs['jobid'].append(j)
            dict_jobs['array'].append(a)
            dict_jobs['status'].append(s)
            dict_jobs['partition'].append(p)

    df_jobs = pd.DataFrame.from_dict(dict_jobs)

    printout = ""
    for jobid in pd.unique(df_jobs['jobid']):
        printout += "Job ID : %s\n"%jobid
        jobdf = df_jobs[df_jobs['jobid']==jobid]
        for partition in pd.unique(jobdf['partition']):
            printout += "\t %s\n"%partition
            partdf = df_jobs[df_jobs['jobid']==jobid]
            statuses = pd.unique(partdf['status'])
            N = partdf.shape[0]
            for status in sorted(statuses):
                num = partdf[partdf['status']==status].shape[0]
                if status == "PENDING": # Need to develop array line
                    statdf = partdf[partdf['status']==status]
                    arrays = pd.unique(statdf['array'])
                    for array in arrays:
                        if "[" in array: 
                            ps = subprocess.Popen(['squeue','-j',jobid,'-r','--noheader','--long'], stdout=subprocess.PIPE)
                            outs,_ = ps.communicate()
                            sup = len([o for o in outs.decode("utf-8").splitlines() if 'PENDING' in o])
                            num += sup -1
                            N += sup -1
                printout += "\t\t"+status.ljust(20,' ')+"{:5d}  [{:6.2f}%]\n".format(num,num/N*100)
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

            # Clear screen #
            stdscr.clear()

            # Loop #
            while True:
                printout = parse_sacct(args)
                report_progress(stdscr,printout)
                time.sleep(args.loop)
        except KeyboardInterrupt:
            print ("User stop")
            curses.nocbreak()
            stdscr.keypad(False)
            curses.echo()
            curses.endwin()

if __name__ == '__main__':
    main(args)
