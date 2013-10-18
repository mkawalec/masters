#!/usr/bin/env python2

""" 
Logs in to cplab computers and launches jobs in parallel
"""

from multiprocessing import Process
from subprocess import call
#from progressbar import Bar, ETA, Percentage, ProgressBar
#from time import sleep

def setup_remote(host):
    call(["ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no "
          "s0905879@%(host)s \' cd /dev/shm; rm -rf turb ;"
          "git clone https://github.com/mkawalec/masters turb; "
          "cd turb/C; git checkout ssh_bomb; "
          "mkdir build; cd build; cmake ..; make -j3; "
          "./integrator 7 0.0005 1000 10\'" % dict(host=host)], shell=True)
    call(["scp -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no "
          "s0905879@%(host)s:/dev/shm/turb/C/build/output %(host)s.out" % \
            dict(host=host)], shell=True)


if __name__ == '__main__':
    processes = []
    for i in range(10):
        processes.append(Process(target=setup_remote, args=["cplab%03d" % (i,)]))

    for process in processes:
        process.start()
