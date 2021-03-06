#!/usr/bin/env python2

""" 
Logs in to cplab computers and launches jobs in parallel
"""

from __future__ import division
from multiprocessing import Process
from subprocess import call
from progressbar import Bar, ETA, Percentage, ProgressBar
from time import sleep
from sys import argv
import sys
from glob import glob
import os

import numpy as np
from scipy.optimize import curve_fit

import signal
import sys


processes = []
try:
    hosts = int(argv[1])
except ValueError:
    hosts = argv[1]

if len(argv) > 2:
    runs = int(argv[2])

current_dir = ""

devnull = open("/dev/null", "w")
errlog = open("errlog", "w")

def setup_remote(host, runs, folder, dt=0.0005, 
        samples=7, directory='.', tmax=2000, R=1.0):

    call(["ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no "
          "s0905879@%(host)s \' cd /dev/shm; rm -rf %(folder)s ;"
          "mkdir %(folder)s; cp ~/integrator %(folder)s; cd %(folder)s; "
          "./integrator -c decay-path --samples %(samples)s --dt %(dt)s "
          "--end-time %(tmax)s --runs %(runs)s -R %(R)s --fast-threshold 2000\'" 
          % dict(host=host, runs=runs, folder=folder, dt=dt, samples=samples,
              tmax=tmax, R=R)], 
          shell=True, stdout=devnull, stderr=errlog)

    call(["scp -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no "
          "s0905879@%(host)s:/dev/shm/%(folder)s/stationary %(directory)s/%(host)s-%(folder)s.out" 
          % dict(host=host, folder=folder, directory=directory)], 
            shell=True, stderr=errlog, stdout=devnull)

def kill(host):
    call(["ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no "
          "s0905879@%(host)s \' pkill -u s0905879 \'" % dict(host=host)],
          shell=True, stdout=devnull, stderr=errlog)

def kill_all():
    procs = []
    for i in range(1, hosts):
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

    # Spawning two threads per host
    for i in range(1, hosts):
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
            break

    kill_all()
    pbar.finish()

def fit_func(x, a, c):
    ''' Example function used for curve fitting '''
    return a * np.exp(-c * x)

def fit(values, func):
    ''' Fit the values to the function '''

    # Decay probabilities
    y = np.empty([len(values)])

    # Points at which we have decay measurements
    x = np.empty([len(values)])
    values = sorted(values)

    for i, time in enumerate(values):
        y[i] = (len(values) - i) / len(values)
        x[i] = time

    start = 4 * len(values) / 5
    popt, pcov = curve_fit(func, x[start:], y[start:], [1, 0.005])
    print("Parameter values are", popt)
    print("Computed with %d values\n" % (len(values)))
    return popt

def finalize(directory):
    ''' Finalize computation inside a dir '''
    print("Finishing")
    with open(directory + '/output', 'w') as f:
        values = []

        for filename in glob(directory + '/*.out'):
            with open(filename, 'r') as input_f:
                lines = input_f.read()
                values.extend(map(lambda x: float(x), 
                    filter(lambda x: len(x) > 0, lines.split('\n'))))
                f.write(lines)

            os.remove(filename)
        fit(values, fit_func)

def gen_signal_hdl():
    handled = [False]

    def signal_handler(signal, frame):
        ''' Handle ctrl-c '''
        if handled[0]:
            return

        handled[0] = True
        kill_all()
        finalize(current_dir) if len(current_dir) > 0 else None
        sys.exit(0)

    return signal_handler

def fit_R():
    fit_output = open('fit_results', 'w')

    for filename in glob('R*/output'):
        with open(filename, 'r') as f:
            print("Fitting %s" % filename)
            lines = f.read().split('\n')
            parsed = map(lambda x: float(x), 
                    filter(lambda x: len(x) > 0, lines))
            fit_results = fit(parsed, fit_func)
            fit_output.write("%s\t%06e\t%06e\n" % 
                    (filename, fit_results[0], fit_results[1]))

    fit_output.close()


if __name__ == '__main__':
    if hosts == 'fit':
        fit_R()
        sys.exit()

    signal.signal(signal.SIGINT, gen_signal_hdl())

    maxt = 10000
    s_dt = 0.0005
    dt = 0.0005
    samples = 7
    R = 1.04

    print("starting at R =", R)
    directory = 'R_' + str(R)
    current_dir = directory
    os.mkdir(directory)

    run_set(dt, samples, directory, maxt, R)
    #finalize(directory)


