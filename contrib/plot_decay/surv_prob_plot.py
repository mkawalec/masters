#!/usr/bin/env python2

import matplotlib.pyplot as plt
import scipy.optimize as opt
from sys import argv
import numpy as np
from glob import glob
from colorsys import hsv_to_rgb
from pylab import savefig

import redis
import pickle

from os import urandom
import random as rnd

import params


r = redis.StrictRedis(host='localhost', port=6379, db=0)

def model_function(x, a0, a1):
    return a0 * np.exp(-a1 * x)


def gen_color(curr, total):
    return hsv_to_rgb(float(curr) / total, 0.7, 0.9)


def random_selection(points, prop):
    rnd.seed(urandom(32))
    selected = {}
    new_points = []

    for i in range(int(len(points) * prop)):
        while True:
            current = rnd.randrange(len(points))
            if current not in selected:
                selected[current] = True
                new_points.append(points[current])
                break

    return sorted(new_points, key=lambda x: x[0])

def gen_sum(acc, el):
    acc[0] += el[0]
    acc[1] += el[1]
    return acc

def get_free_params(filename):
    params = []
    for value in filename.split('-'):
        try:
            value = float(value)
            params.append(value)
        except ValueError:
            pass

    return params


def set_value(params, container, value):
    if len(params) == 1:
        container[params[0]] = value
        return

    if params[0] not in container:
        container[params[0]] = {}

    set_value(params[1:], container[params[0]], value)


def gen_fits(fit_prop, prop=1, prefix='', cache=True):
    ''' Generates the fits to the experimental data,
        and caches them on redis. Please use this function
        to access the fits and not generate them yourself
    '''

    key = 'fits-%f-%f-%s' % (fit_prop, prop, prefix)
    fits = r.get(key)

    if fits and cache:
        return pickle.loads(fits)
    else:
        fits = {}

    for dist in glob(str(prefix) + '*-survival_probability'):
        params = get_free_params(dist)
        print(dist)

        with open(dist) as f:
            points = map(lambda x: (float(x.split(' ')[0]),
                         float(x.split(' ')[1])),
                         f.readlines())


        random_popt = []

        for i in range(int(1/prop)):
            new_points = random_selection(points, prop)

            fit_from = int(len(new_points) * (1 - fit_prop))
            popt, pcov = opt.curve_fit(model_function,
                                np.array(map(lambda x: x[0], new_points[fit_from:])),
                                np.array(map(lambda x: x[1], new_points[fit_from:])),
                                (1, 0), maxfev=100000)
            random_popt.append(popt)

        # Getting the average
        avg = map(lambda x: x / len(random_popt), reduce(gen_sum, random_popt, [0, 0]))
        variance = [np.var(map(lambda x: x[i], random_popt)) for i in range(2)]

        set_value(params, fits, dict(variance=variance, avg=avg))

    r.set(key, pickle.dumps(fits))
    return fits


def plot_fits(fit_prop, prop=1, cache=True):
    fits = gen_fits(fit_prop, prop, cache=cache)

    fig, ax = plt.subplots()
    #plt.title('Fit to tail of %.1f, with %.2f in each sample' % (fit_prop, prop))
    for n, domain in enumerate(sorted(fits.keys(), key=lambda x: float(x))):
        if domain == 10 or domain == 12:
            continue

        x = sorted(fits[domain].keys(), key=lambda x: float(x))
        y = map(lambda x: fits[domain][x]['avg'][1], x)
        err = map(lambda x: fits[domain][x]['variance'][1], x)

        ax.plot(x, y, label='%s$\pi$' % (domain),
                color=hsv_to_rgb(float(n-2) / (len(fits.keys())-2), 0.7, 0.9))

    ax.legend(loc=0)
    ax.set_ylabel('Value of $a_1$')
    ax.set_xlabel('R')
    plt.subplots_adjust(bottom=0.15)
    savefig('a1_vs_R_vs_domain_size_fits.png', dpi=600)


def gen_surv_plot():
    with open('24-1.05-survival_probability') as f:
        points = map(lambda x: (float(x.split(' ')[0]),
                                float(x.split(' ')[1])),
                     f.readlines())

    fits = gen_fits(0.1, 0.1)
    a0, a1= fits[24][1.05]['avg']
    color1 = hsv_to_rgb(0, 0.7, 0.9)
    color2 = hsv_to_rgb(0.5, 0.7, 0.9)

    fig, ax = plt.subplots()
    ax.plot(map(lambda x: x[0], points), map(lambda x: x[1], points),
            label="Data", marker=',', color=color1)
    ax.plot(map(lambda x: x[0], points),
            map(lambda x: model_function(x, a0, a1),
                map(lambda x: x[0], points)),
            label="Fit", color=color2, alpha=0.7)

    #plt.title('Data and fit at R=1.06')
    ax.legend(loc=0)
    ax.set_ylabel('$S(t)$')
    ax.set_xlabel('$t$')
    ax.set_xlim([0, 7000])
    plt.hlines(0.1, 0, 7000, linestyles='dashed')
    plt.text(5000, 0.11, 'Fit threshold')
    plt.subplots_adjust(bottom=0.15)

    savefig('surv_prob.png', dpi=600)

if __name__ == "__main__":
    plot_fits(0.3, 0.010)
    gen_surv_plot()
    #plt.show()



