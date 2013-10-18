#!/usr/bin/env python2

""" 
Logs in to cplab computers and launches jobs in parallel
"""

from multiprocessing import Process
from subprocess import call
from progressbar import Bar, ETA, Percentage, ProgressBar
from time import sleep

def setup_remote(host):
    call(["ssh s0905879@%(host)s \' cd /dev/shm; "
         "git clone https://github.com/mkawalec/masters turb; "
         "cd turb; mkdir build; cd build; cmake ..; make -j3; "
         "./integrator 7 0.0005 1000 10" % (host=host)], shell=True)
    call(["scp s0905879@%(host)s:/dev/shm/turb/build/output %(host)s.out" % \
            (host=host)], shell=True)


if __name__ == '__main__':
    for i in xrange(10):
        setup_remote("cplab%03d" % (i,))
