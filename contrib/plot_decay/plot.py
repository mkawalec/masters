#!/usr/bin/env python2

from glob import glob
import re


def plot_decay():
    decay_rates = []
    for fit_file in glob('*fit'):
        with open(fit_file, 'r') as f:
            line = f.readline()
            numbers = line.split(' ')
            decay_rates.append(dict(c0=float(numbers[0]),
                                    c1=float(numbers[1]),
                                    R=float(fit_file.split('-')[0])))

    print decay_rates



if __name__ == '__main__':
    plot_decay()
