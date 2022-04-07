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
parser.add_argument('-s','--sum', action='store_true', required=False, default=False,
                    help='Only displays the sum of all jobs status from the different job ids')
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

    dict_jobs = {'jobid':[],'array':[],'status':[],'partition':[]}
    for line in out.decode("utf-8").splitlines():
        j, s, p = line.split()
        arrs = []
        if '_' in j:
            j, a = j.split('_')
            if s == 'PENDING': # Need to develop array line
                ps = subprocess.Popen(['squeue','-j',j,'-r','-t','pending','--noheader'], stdout=subprocess.PIPE)
                outs,_ = ps.communicate()
                for pline in outs.decode("utf-8").splitlines():
                    arrs.append(pline.split()[0].split('_')[1])
            else:
                arrs.append(a)
        else:
            arrs = [1]
        for a in arrs:
            dict_jobs['jobid'].append(j)
            dict_jobs['array'].append(a)
            dict_jobs['status'].append(s)
            dict_jobs['partition'].append(p)

    df_jobs = pd.DataFrame.from_dict(dict_jobs)

    status_list = ['PENDING','RUNNING','COMPLETED']

    if args.sum:
        df_jobs['jobid'] = 'total'

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
