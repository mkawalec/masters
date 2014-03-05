#!/usr/bin/env python2

from glob import glob
import re
import matplotlib.pyplot as plt
import numpy as np


def plot_decay():
    decay_rates = []
    for fit_file in glob('*fit'):
        with open(fit_file, 'r') as f:
            line = f.readline()
            numbers = line.split(' ')
            decay_rates.append(dict(c0=float(numbers[0]),
                                    c1=float(numbers[1]),
                                    R=float(fit_file.split('-')[0])))

    decay_rates = sorted(decay_rates, key=lambda x: x['R'])
    x = np.linspace(0, 1500, 1500)
    fig, ax = plt.subplots()

    for rate in decay_rates:
        y = rate['c0'] * np.exp(- rate['c1'] * x)
        ax.plot(x, y, label='R = %s' % (rate['R']))

    ax.legend(loc=0)
    plt.show()



if __name__ == '__main__':
    plot_decay()