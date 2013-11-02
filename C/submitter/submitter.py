#!/usr/bin/env python2

""" 
Logs in to cplab computers and launches jobs in parallel
"""

from multiprocessing import Process
from subprocess import call
from progressbar import Bar, ETA, Percentage, ProgressBar
from time import sleep
from sys import argv
from glob import glob
import os

devnull = open("/dev/null", "w")
errlog = open("errlog", "w")

def setup_remote(host, runs, folder, dt=0.0005, samples=7):
    call(["ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no "
          "s0905879@%(host)s \' cd /dev/shm; rm -rf %(folder)s ;"
          "mkdir %(folder)s; cp ~/integrator %(folder)s; cd %(folder)s; "
          "./integrator %(samples)s %(dt)s 10000 %(runs)d\'" 
          % dict(host=host, runs=runs, folder=folder, dt=dt, samples=samples)], 
          shell=True, stdout=devnull, stderr=errlog)

    call(["scp -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no "
          "s0905879@%(host)s:/dev/shm/%(folder)s/output %(host)s-%(folder)s.out" 
          % dict(host=host, folder=folder)], 
            shell=True, stderr=errlog, stdout=devnull)

def run_set(dt, samples):    
    hosts = int(argv[1])
    runs = int(argv[2])

    processes = []
    print("Starting dt = %s, samples = %s" % (dt, samples))

    # Spawning two threads per host
    for i in range(hosts):
        processes.append(Process(target=setup_remote, 
            args=["cplab%03d" % (i,), runs, 'turb1', dt, samples]))
        processes.append(Process(target=setup_remote, 
            args=["cplab%03d" % (i,), runs, 'turb2', dt, samples]))

    for process in processes:
        process.start()

    counter = 0
    widgets = ['Computing in CPLB:', Percentage(), ' ', 
            Bar(marker='#', left='[', right=']'), ' ', ETA(), ' ']
    pbar = ProgressBar(widgets=widgets, maxval=len(processes))
    while counter < len(processes):
        counter = 0
        sleep(1)
        for process in processes:
            if not process.is_alive():
                counter += 1
        pbar.update(counter)
    pbar.finish()

if __name__ == '__main__':
    s_dt = 0.0005
    for dt in [2 * s_dt, s_dt, s_dt / 2, s_dt / 4]:
        for samples in range(7, 11):
            run_set(dt, samples)
            with open('output_' + str(samples) + '_' + str(dt), 'w') as f:
                for filename in glob('*.out'):
                    with open(filename, 'r') as input_f:
                        f.write(input_f.read())
                    os.remove(filename)


