#!/usr/bin/env python2

import numpy as np
import matplotlib.pyplot as plt
from pylab import savefig

from colorsys import hsv_to_rgb
from surv_prob_plot import gen_fits
import params


def plot_mult():
    fits = gen_fits(0.3, 0.1, prefix='mult/')
    fig, ax = plt.subplots()

    x = sorted(fits.keys())
    y = map(lambda x: fits[x]['avg'][1], x)
    yerr = map(lambda x: fits[x]['variance'][1], x)
    color2 = hsv_to_rgb(0.5, 0.7, 0.9)

    ax.plot(x, y, marker='o', color=color2)
    ax.set_yscale('log')
    ax.set_ylabel('$a_1$')

    ax.set_xlabel('Input conditions multiplier', size=11)
    ax.set_xscale('log')
    plt.subplots_adjust(bottom=0.15)
    savefig('a1_mult.png', dpi=600)


if __name__ == '__main__':
    plot_mult()
    #plt.show()
