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

import numpy as np
from scipy.optimize import curve_fit

import signal
import sys


processes = []
hosts = int(argv[1])
runs = int(argv[2])
current_dir = ""

devnull = open("/dev/null", "w")
errlog = open("errlog", "w")

def setup_remote(host, runs, folder, dt=0.0005, 
        samples=7, directory='.', tmax=2000, R=1.0):

    call(["ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no "
          "s0905879@%(host)s \' cd /dev/shm; rm -rf %(folder)s ;"
          "mkdir %(folder)s; cp ~/integrator %(folder)s; cd %(folder)s; "
          "./integrator --samples %(samples)s --dt %(dt)s "
          "--end-time %(tmax)s --runs %(runs)s -R %(R)s\'" 
          % dict(host=host, runs=runs, folder=folder, dt=dt, samples=samples,
              tmax=tmax, R=R)], 
          shell=True, stdout=devnull, stderr=errlog)

    call(["scp -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no "
          "s0905879@%(host)s:/dev/shm/%(folder)s/output %(directory)s/%(host)s-%(folder)s.out" 
          % dict(host=host, folder=folder, directory=directory)], 
            shell=True, stderr=errlog, stdout=devnull)

def kill(host):
    call(["ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no "
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

def run_set(dt, samples, directory, tmax, R):    
    kill_all()
    processes = []

    print("Starting dt = %s, samples = %s" % (dt, samples))

    # Spawning two threads per host
    for i in range(hosts):
        processes.append(Process(target=setup_remote, 
            args=["cplab%03d" % (i,), runs, 'turb1', dt, samples, directory, tmax, R]))
        processes.append(Process(target=setup_remote, 
            args=["cplab%03d" % (i,), runs, 'turb2', dt, samples, directory, tmax, R]))

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

        if counter > 0.9 * len(processes):
            print("Finishing")
            break

    kill_all()
    pbar.finish()

def fit_func(x, a, c):
    ''' Example function used for curve fitting '''
    return a * np.exp(c * x)

def fit(values, func):
    ''' Fit the values to the function '''

    # Decay probabilities
    y = np.empty([len(values)])

    # Points at which we have decay measurements
    x = np.empty([len(values)])
    for i, time in enumerate(values):
        y[i] = (len(values) - i) / len(values)
        x[i] = time

    popt, pcov = curve_fit(func, x, y)
    print("Parameter values are", popt)

def finalize(directory):
    ''' Finalize computation inside a dir '''
    with open(directory + '/output', 'w') as f:
        values = []

        for filename in glob(directory + '/*.out'):
            with open(filename, 'r') as input_f:
                lines = input_f.read()
                values.extend(map(lambda x: int(x), lines.split('\n')))
                f.write(lines)

            os.remove(filename)
        fit(values, fit_func)

def signal_handler(signal, frame):
    ''' Handle ctrl-c '''
    kill_all()
    finalize(current_dir) if len(current_dir) > 0 else None
    sys.exit(0)


if __name__ == '__main__':
    signal.signal(signal.SIGINT, signal_handler)

    maxt = 7000
    s_dt = 0.0005
    dt = 0.0005
    samples = 7
    R = 1.0

    while R < 1.08:
        print "starting at R =", R
        directory = 'R_' + str(R)
        current_dir = directory
        os.mkdir(directory)

        run_set(dt, samples, directory, maxt, R)
        finalize(directory)

        R += 0.1


