#!/usr/bin/env python2

from math import pi
from sys import argv
import numpy as np


def create_submitters():
    with open(argv[1], 'r') as f:
        template = f.read()

    if argv[2] == 'mult':
        create_input_mult(template)
    else:
        create_diff_domains(template)


def create_input_mult(template):
    output = open("submit-input-mult", "w")
    output.write(template)

    for mult in np.logspace(0, 2, 20):
        output.write("aprun -n 2400 -N 24 /work/d59/d59/s0905879/integrator "
                     "-c decay-mult -r 40 -f 1 --end-time 7000 --dt 0.00005 "
                     "--start-mult %f --find-zeros 0 --use-output 0 --prefix '%s'\n" %
                     (mult, 'input-mult-' + str(mult) + '-'))

    output.close()


def create_diff_domains(template):
    # Domain size from 24 to 2 pi
    for domain in range(24, 0, -2):
        with open("submit-" + str(domain) + 'pi.sh', 'w') as f:
            f.write(template)

            # R range from 0.9 to 1.09
            for R in [x * 0.01 for x in range(90, 110)]:
                f.write("aprun -n 2400 -N 24 /work/d59/d59/s0905879/integrator "
                        "-c decay-mult -r 40 -f 1 --end-time 7000 -d %s "
                        "--find-zeros 0 --use-output 0 -R %s --prefix '%s'\n"
                        % (str(domain * pi), str(R), str(domain) + '-' + str(R) + '-'))


if __name__ == '__main__':
    create_submitters()
