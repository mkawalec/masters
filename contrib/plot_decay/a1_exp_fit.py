#!/usr/bin/env python2

from sys import argv
import scipy.optimize as opt
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr as pear

from plot2 import get_a1
from surv_prob_plot import gen_fits
from colorsys import hsv_to_rgb

import matplotlib.gridspec as gridspec

import params


def model_function(x, a0, a1, a2):
    return a0 * np.exp(-a1 * np.exp(a2 * x))

def execute_fit():
    fits = gen_fits(0.3, 0.1)

    gs = gridspec.GridSpec(5, 1)
    gs.update(wspace=0)
    ax1 = plt.subplot(gs[:4, 0])
    ax2 = plt.subplot(gs[4, 0])

    coefficients = {}

    for n, domain in enumerate(sorted(fits.keys(), key=lambda x: float(x))):
        curr_keys = sorted(fits[domain].keys(), key=lambda x: float(x))
        x = map(lambda x: float(x), curr_keys)
        y = map(lambda x: fits[domain][x]['avg'][1], curr_keys)

        try:
            popt, pcov = opt.curve_fit(model_function, np.array(x), np.array(y),
                                       (0.5324, -0.0090589, 6.2374),
                                       maxfev=100000)
            coefficients[domain] = popt
        except RuntimeError:
            continue

        color1 = hsv_to_rgb(float(2 * n) / (2 * len(fits.keys())), 0.7, 0.9)
        color2 = hsv_to_rgb(float(2 * n + 1) / (2 * len(fits.keys())), 0.7, 0.9)

        if int(domain)%4 == 0:
            ax1.plot(x, y,
                    label='data %s' % (domain),
                    color=color1)

            ax1.plot(x,
                    [model_function(x_val, *popt) for x_val in x],
                    label='fit %s' % (domain), color=color2)

        if domain == "24":
            values = [model_function(x_val, *popt) for x_val in x]
            residuals = map(lambda x: abs(x[0] - x[1]) / x[0], zip(y, values))
            ax2.plot(x, residuals, marker='o')
            ax2.set_yscale('log')

    ax1.legend(loc=0)
    plot_coeffs(coefficients)
    plt.tight_layout()
    plt.show()


def plot_coeffs(coeffs):
    fig, ax = plt.subplots()
    toplot = dict(x=[], y1=[], y2=[])

    for domain in sorted(coeffs.keys(), key=lambda x: float(x)):
        toplot['x'].append(float(domain))
        toplot['y1'].append(coeffs[domain][1])
        toplot['y2'].append(coeffs[domain][2])

    ax.plot(toplot['x'], toplot['y1'], marker='o',
            color=hsv_to_rgb(0, 0.7, 0.9), label='b1')
    ax.plot(toplot['x'], toplot['y2'], marker='o',
            color=hsv_to_rgb(0.5, 0.7, 0.9), label='b2')

    ax.legend(loc=0)
    ax.set_yscale('log')


if __name__ == '__main__':
    execute_fit()


