#!/usr/bin/env python2

from sys import argv
import scipy.optimize as opt
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr as pear

from plot2 import get_a1
import params


def model_function(x, a0, a1, a2):
    return a0 * np.exp(-a1 * np.exp(a2 * x))

def execute_fit():
    a1 = get_a1(argv[1])
    fig, ax = plt.subplots(2)

    coefficients = {}

    for domain in sorted(a1.keys(), key=lambda x: float(x)):
        try:
            popt, pcov = opt.curve_fit(model_function,
                                   np.array(a1[domain]['x']), np.array(a1[domain]['y']),
                                   (0.5324, -0.0090589, 6.2374), maxfev=100000)
            coefficients[domain] = popt
        except RuntimeError:
            continue

        ax[0].plot(a1[domain]['x'], a1[domain]['y'],
                label='data %s' % (domain))

        ax[0].plot(a1[domain]['x'],
                [model_function(x, *popt) for x in a1[domain]['x']],
                label='fit %s' % (domain))

        if domain == "24":
            values = [model_function(x, *popt) for x in a1[domain]['x']]
            residuals = map(lambda x: abs(x[0] - x[1]) / x[0], zip(a1[domain]['y'], values))
            ax[1].plot(a1[domain]['x'], residuals, marker='o')
            ax[1].set_yscale('log')
            print pear(values, a1[domain]['y'])


    ax[0].legend(loc=0)
    plot_coeffs(coefficients)
    plt.show()


def plot_coeffs(coeffs):
    fig, ax = plt.subplots()
    toplot = dict(x=[], y1=[], y2=[])

    for domain in sorted(coeffs.keys(), key=lambda x: float(x)):
        toplot['x'].append(float(domain))
        toplot['y1'].append(coeffs[domain][1])
        toplot['y2'].append(coeffs[domain][2])

    ax.plot(toplot['x'], toplot['y1'], marker='o')
    ax.plot(toplot['x'], toplot['y2'], marker='o')
    ax.set_yscale('log')


if __name__ == '__main__':
    execute_fit()


