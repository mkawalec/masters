#!/usr/bin/env python2

""" 
Logs in to cplab computers and launches jobs in parallel
"""

from multiprocessing import Process
from subprocess import call
from progressbar import Bar, ETA, Percentage, ProgressBar
from time import sleep

devnull = open("/dev/null", "w")

def setup_remote(host):
    call(["ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no "
          "s0905879@%(host)s \' cd /dev/shm; rm -rf turb ;"
          "git clone -b ssh_bomb https://github.com/mkawalec/masters turb; "
          "cd turb/C; mkdir build; cd build; cmake ..; make -j3; "
          "./integrator 7 0.0005 1000 10\'" % dict(host=host)], shell=True,
          stdout=devnull, stderr=devnull)
    call(["scp -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no "
          "s0905879@%(host)s:/dev/shm/turb/C/build/output %(host)s.out" % \
            dict(host=host)], shell=True, stderr=devnull, stdout=devnull)


if __name__ == '__main__':
    processes = []
    for i in range(20):
        processes.append(Process(target=setup_remote, args=["cplab%03d" % (i,)]))

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
