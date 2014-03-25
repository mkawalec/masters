#!/usr/bin/env python2


from glob import glob
import numpy as np

import matplotlib.pyplot as plt
from pylab import savefig
import pylab

from colorsys import hsv_to_rgb
import params


def get_avg(prefix='../../build/', maxlen=-1):
    u = []
    for filename in glob(prefix + 'process0-output*'):
        with open(filename) as f:
            lines = f.readlines()
            u.append(map(lambda x: float(x.split(' ')[1]),
                         filter(lambda x: len(x.split(' ')) > 2, lines)))

    u = filter(lambda x: len(x) > 0, u)
    for points in u:
        if len(points) > maxlen:
            maxlen = len(points)
    for points in u:
        for i in range(maxlen - len(points)):
            points.append(0)


    transformed = map(lambda x: np.fft.fft(np.array(x)), u)
    trans_len = transformed[0].shape[0]
    avg = [0] * trans_len

    for i in range(trans_len):
        for j in range(len(transformed)):
                avg[i] += transformed[j][i] * np.conjugate(transformed[j][i])

    for i in range(len(avg)):
        avg[i] /= len(transformed)

    return dict(x=map(lambda x: 1 / 0.05 * float(x) / len(avg), range(len(avg))), y=avg)

def plot_avg():
    maxlen = -1
    for filename in glob('/tmp/e-02/process0-output*'):
        with open(filename) as f:
            length = len(f.readlines())
            if length > maxlen:
                maxlen = length

    for filename in glob('/tmp/e-05/process0-output*'):
        with open(filename) as f:
            length = len(f.readlines())
            if length > maxlen:
                maxlen = length

    second = get_avg('/tmp/e-02/', maxlen)
    fifth = get_avg('/tmp/e-05/', maxlen)


    color1 = hsv_to_rgb(0, 0.7, 0.9)
    color2 = hsv_to_rgb(0.5, 0.7, 0.9)
    fig, ax = plt.subplots()
    ax.plot(fifth['x'], fifth['y'], label='1e-05', color=color2)
    ax.plot(second['x'], second['y'], alpha=0.8, label='1e-02', color=color1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc=0)
    plt.subplots_adjust(bottom=0.15)
    ax.set_xlabel('frequency')
    savefig('fourier.png', dpi=600)

def plot_portrait():
    with open('/tmp/e-02/process0-output16') as f:
        lines = f.readlines()
        pos02 = map(lambda y: map(lambda x: float(x), y.split(' ')[:3]),
                                  filter(lambda z: len(z.split(' ')) > 2,
                                         lines))

    with open('/tmp/e-05/process0-output02') as f:
        lines = f.readlines()
        pos05 = map(lambda y: map(lambda x: float(x), y.split(' ')[:3]),
                                  filter(lambda z: len(z.split(' ')) > 2,
                                         lines))

    fig, ax = plt.subplots()
    color1 = hsv_to_rgb(0, 0.7, 0.9)
    color2 = hsv_to_rgb(0.5, 0.7, 0.9)

    ax.plot(map(lambda x: x[1], pos02), map(lambda x: x[2], pos02),
            label='$10^{-2}$', color=color1)
    ax.plot(map(lambda x: x[1], pos05), map(lambda x: x[2], pos05),
            label='$10^{-5}$', color=color2, alpha=0.8)
    ax.legend(loc=0)
    savefig('phase.png', dpi=600)


if __name__ == '__main__':
    plot_avg()
    plot_portrait()
    plt.show()

