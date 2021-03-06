#!/usr/bin/env python2

from glob import glob
import re
import matplotlib.pyplot as plt
import numpy as np
from sys import argv


def get_a1(pattern):
    a1 = {}
    for fit_file in glob(pattern):
        with open(fit_file) as f:
            line = f.readline()
            coeffs = line.split(' ')
            fit_params = fit_file.split('-')

            if fit_params[0] not in a1:
                a1[fit_params[0]] = []

            a1[fit_params[0]].append((float(fit_params[1]), float(coeffs[1])))

    # Sort and remove the soring hints
    for key in a1.keys():
        a1[key] = sorted(a1[key], key=lambda x: x[0])
        a1[key] = dict(y=map(lambda x: float(x[1]), a1[key]),
                       x=map(lambda x: float(x[0]), a1[key]))
    return a1

def plot_a1():
    a1 = get_a1(argv[1])
    fig, ax = plt.subplots()

    for domain in sorted(a1.keys(), key=lambda x: float(x)):
        ax.plot(a1[domain]['x'], a1[domain]['y'],
                label='%s pi' % (domain))

    ax.legend(loc=0)
    fig.savefig('a1.png', dpi=300)
    plt.show()


if __name__ == '__main__':
    plot_a1()
