#!/usr/bin/env python2

from sys import argv
import scipy.optimize as opt
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr as pear
from pylab import savefig

from plot2 import get_a1
from surv_prob_plot import gen_fits
from colorsys import hsv_to_rgb
from matplotlib import rcParams

import matplotlib.gridspec as gridspec
from get_dims import get_dims

import params


def model_function(x, a0, a1, a2):
    return a0 * np.exp(-a1 * np.exp(a2 * x))

def execute_fit():
    fits = gen_fits(0.3, 0.1)

    rcParams['figure.figsize'] = 4.981320049813201, 4
    gs = gridspec.GridSpec(5, 1)
    gs.update(wspace=0)
    ax1 = plt.subplot(gs[:4, 0])
    ax2 = plt.subplot(gs[4, 0], sharex=ax1)

    coefficients = {}

    for n, domain in enumerate(sorted(fits.keys(), key=lambda x: float(x))):
        if domain == 10:
            continue

        curr_keys = sorted(fits[domain].keys(), key=lambda x: float(x))
        x = map(lambda x: float(x), curr_keys)
        y = map(lambda x: fits[domain][x]['avg'][1], curr_keys)

        with open('domain_' + str(domain), 'w') as f:
            for i,R in enumerate(x):
                f.write("%f %f\n" % (R, y[i]))

        try:
            popt, pcov = opt.curve_fit(model_function, np.array(x), np.array(y),
                                       (0.5324, -0.0090589, 6.2374),
                                       maxfev=100000)
            coefficients[domain] = popt
        except RuntimeError:
            continue

        color1 = hsv_to_rgb(float(n) / len(fits.keys()), 0.7, 0.9)
        color2 = hsv_to_rgb(float(n + 1) / len(fits.keys()), 0.7, 0.9)

        if int(domain)%4 == 0 and int(domain) != 20:
            ax1.plot(x, y,
                    label='Data at %s$\pi$' % (domain),
                    color=color1)

            ax1.plot(x,
                    [model_function(x_val, *popt) for x_val in x],
                    label='Fit at %s$\pi$' % (domain), color=color2)

        if domain == 24:
            values = [model_function(x_val, *popt) for x_val in x]
            residuals = map(lambda x: abs(x[0] - x[1]) / x[0], zip(y, values))
            ax2.plot(x, residuals, marker='o')
            ax2.set_yscale('log')

    ax1.legend(loc=0)
    ax1.xaxis.set_visible(False)
    ax2.set_xlabel('R')
    ax1.set_ylabel('Value of $a_1$')
    ax2.set_ylabel('Residue\nat 24$\pi$')
    ax2.yaxis.tick_right()
    plt.subplots_adjust(bottom=0.15)

    savefig('double_exp_fit.png', dpi=600)
    plot_coeffs(coefficients)
    #plt.show()


def plot_coeffs(coeffs):
    rcParams['figure.figsize'] = get_dims()
    fig, ax1 = plt.subplots()

    toplot = dict(x=[], y1=[], y2=[])

    for domain in sorted(coeffs.keys(), key=lambda x: float(x)):
        toplot['x'].append(float(domain))
        toplot['y1'].append(coeffs[domain][1])
        toplot['y2'].append(coeffs[domain][2])

    ax1.plot(toplot['x'], toplot['y1'], marker='o',
            color=hsv_to_rgb(0, 0.7, 0.9), label='$b_1$')

    plt.subplots_adjust(bottom=0.3, left=0.20)
    ax1.set_yscale('log')
    ax1.set_xlabel('Domain size')
    savefig('b1_vs_domain.png', dpi=600)

    fig, ax2 = plt.subplots()
    ax2.plot(toplot['x'], toplot['y2'], marker='o',
            color=hsv_to_rgb(0.5, 0.7, 0.9), label='$b_2$')

    ax2.set_xlabel('Domain size')
    plt.subplots_adjust(bottom=0.3)
    savefig('b2_vs_domain.png', dpi=600)


if __name__ == '__main__':
    execute_fit()


