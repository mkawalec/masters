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

processes = []
hosts = int(argv[1])
runs = int(argv[2])

devnull = open("/dev/null", "w")
errlog = open("errlog", "w")

def setup_remote(host, runs, folder, dt=0.0005, samples=7, directory='.', tmax=2000):
    call(["ssh -tt -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no "
          "s0905879@%(host)s \' cd /dev/shm; rm -rf %(folder)s ;"
          "mkdir %(folder)s; cp ~/integrator %(folder)s; cd %(folder)s; "
          "./integrator %(samples)s %(dt)s %(tmax)s %(runs)d\'" 
          % dict(host=host, runs=runs, folder=folder, dt=dt, samples=samples,
              tmax=tmax)], 
          shell=True, stdout=devnull, stderr=errlog)

    call(["scp -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no "
          "s0905879@%(host)s:/dev/shm/%(folder)s/output %(directory)s/%(host)s-%(folder)s.out" 
          % dict(host=host, folder=folder, directory=directory)], 
            shell=True, stderr=errlog, stdout=devnull)

def kill(host):
    call(["ssh -tt -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no "
          "s0905879@%(host)s \' pkill -u s0905879 \'" % dict(host=host)],
          shell=True, stdout=devnull, stderr=errlog)

def kill_all():
    procs = []
    for i in range(hosts):
        procs.append(Process(target=kill,
            args=["cplab%03d" % (i,)]))

    for process in procs:
        process.start()

    counter = 0
    while counter < len(procs):
        sleep(1)
        counter = 0
        for process in procs:
            if not process.is_alive():
                counter += 1

def run_set(dt, samples, directory, tmax):    
    kill_all()
    processes = []

    print("Starting dt = %s, samples = %s" % (dt, samples))

    # Spawning two threads per host
    for i in range(hosts):
        processes.append(Process(target=setup_remote, 
            args=["cplab%03d" % (i,), runs, 'turb1', dt, samples, directory, tmax]))
        processes.append(Process(target=setup_remote, 
            args=["cplab%03d" % (i,), runs, 'turb2', dt, samples, directory, tmax]))

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


    kill_all()
    pbar.finish()

if __name__ == '__main__':
    maxt = 7000
    s_dt = 0.0005
    for n,dt in enumerate([s_dt, s_dt/2, s_dt/4, s_dt/8]):
        for samples in [7, 8, 9, 10]:
            if samples == 10 and (n == 0 or n == 1):
                continue

            directory = 'new_' + str(samples) + '_' + str(dt)
            os.mkdir(directory)

            run_set(dt, samples, directory, maxt)
            with open(directory + '/output', 'w') as f:
                for filename in glob(directory + '/*.out'):
                    with open(filename, 'r') as input_f:
                        f.write(input_f.read())


