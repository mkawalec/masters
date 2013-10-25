#!/usr/bin/env python2

""" 
Logs in to cplab computers and launches jobs in parallel
"""

from multiprocessing import Process
from subprocess import call
from progressbar import Bar, ETA, Percentage, ProgressBar
from time import sleep
from sys import argv

devnull = open("/dev/null", "w")
errlog = open("errlog", "w")

def setup_remote(host, runs):
    call(["scp -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no "
         "integrator s0905879@%(host)s:/dev/shm" % dict(host=host)], 
         shell=True, stdout=devnull, stderr=errlog)
    call(["ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no "
          "s0905879@%(host)s \' cd /dev/shm; rm -rf turb ;"
          "mkdir turb; mv integrator turb; cd turb; "
          "./integrator 7 0.0005 2000 %(runs)d\'" % dict(host=host, runs=runs)], 
          shell=True, stdout=devnull, stderr=errlo)

    call(["scp -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no "
          "s0905879@%(host)s:/dev/shm/turb/output %(host)s.out" % \
            dict(host=host)], shell=True, stderr=errlog, stdout=devnull)

if __name__ == '__main__':
    hosts = int(argv[1])
    runs = int(argv[2])

    processes = []

    # Spawning two threads per host
    for i in range(hosts):
        processes.append(Process(target=setup_remote, args=["cplab%03d" % (i,), runs]))
        processes.append(Process(target=setup_remote, args=["cplab%03d" % (i,), runs]))

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
