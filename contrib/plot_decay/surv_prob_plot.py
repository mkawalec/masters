#!/usr/bin/env python2

import matplotlib.pyplot as plt
import scipy.optimize as opt
from sys import argv
import numpy as np
from glob import glob
from colorsys import hsv_to_rgb

import redis
import pickle

from os import urandom
import random as rnd

import params

r = redis.StrictRedis(host='localhost', port=6379, db=0)

def model_function(x, a0, a1):
    return a0 * np.exp(-a1 * x)


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

    return new_points

def gen_sum(acc, el):
    acc[0] += el[0]
    acc[1] += el[1]
    return acc


def gen_fits(fit_prop, prop=1, cache=True):
    key = 'fits-%f-%f' % (fit_prop, prop)
    fits = r.get(key)

    if fits and cache:
        return pickle.loads(fits)
    else:
        fits = {}

    for dist in glob('*-survival_probability'):
        domain, R = dist.split('-')[:2]

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
                                (1, 0))
            random_popt.append(popt)

        # Getting the average
        avg = map(lambda x: x / len(random_popt), reduce(gen_sum, random_popt))
        variance = [np.var(map(lambda x: x[i], random_popt)) for i in range(len(popt))]

        if domain not in fits:
            fits[domain] = {}

        fits[domain][R] = dict(variance=variance, avg=avg)

    r.set(key, pickle.dumps(fits))
    return fits


def plot_fits(fit_prop):
    fits = gen_fits(fit_prop)

    fig, ax = plt.subplots()
    plt.title(fit_prop)
    for n, domain in enumerate(sorted(fits.keys(), key=lambda x: float(x))):
        x = sorted(fits[domain].keys(), key=lambda x: float(x))
        y = map(lambda x: fits[domain][x][1], x)
        ax.plot(x, y, label='%s pi' % (domain),
                color=hsv_to_rgb(float(n) / len(fits.keys()), 0.7, 0.9))

    ax.legend(loc=0)


def gen_surv_plot():
    with open('24-1.08-survival_probability') as f:
        points = map(lambda x: (float(x.split(' ')[0]),
                                float(x.split(' ')[1])),
                     f.readlines())

    with open('24-1.08-fit') as f:
        a0, a1 = map(lambda x: float(x), f.read().split(' ')[:2])


    popt, pcov = opt.curve_fit(model_function,
                               np.array(map(lambda x: x[0], points[len(points)*3/9:])),
                               np.array(map(lambda x: x[1], points[len(points)*3/9:])),
                               (a0, a1))

    fig, ax = plt.subplots()
    ax.plot(map(lambda x: x[0], points), map(lambda x: x[1], points),
            label="data")
    ax.plot(map(lambda x: x[0], points),
            map(lambda x: model_function(x, a0, a1),
                map(lambda x: x[0], points)),
            label="fit")

    ax.plot(map(lambda x: x[0], points),
            map(lambda x: model_function(x, popt[0], popt[1]),
                map(lambda x: x[0], points)),
            label="taily fit")

    ax.legend(loc=0)
    plt.show()

if __name__ == "__main__":
    #gen_surv_plot()
    #plot_fits(0.3)
    #plot_fits(0.2)
    #plot_fits(0.1)
    #plot_fits(0.05)
    #plt.show()
    print gen_fits(0.3, 0.1)



