#!/usr/bin/env python2


from glob import glob
import numpy as np

import matplotlib.pyplot as plt
from pylab import savefig
import pylab

from colorsys import hsv_to_rgb
import params


def get_avg(prefix='../../build/'):
    u = []
    for filename in glob(prefix + 'process0-output*'):
        with open(filename) as f:
            lines = f.readlines()
            u.append(map(lambda x: float(x.split(' ')[1]),
                         filter(lambda x: len(x.split(' ')) > 2, lines)))

    u = filter(lambda x: len(x) > 0, u)
    transformed = map(lambda x: np.fft.fft(np.array(x)), u)
    avg = []
    maxlen = -1
    for fft in transformed:
        if fft.shape[0] > maxlen:
            maxlen = fft.shape[0]
    for i in range(maxlen):
        avg.append(0)

    for i in range(maxlen):
        for j in range(len(transformed)):
            if np.shape(transformed[j])[0] > i:
                location = i * int(maxlen / np.shape(transformed[j])[0])
                avg[location] += transformed[j][i] * np.conjugate(transformed[j][i])

    for i in range(len(avg)):
        avg[i] /= len(transformed)

    return dict(x=map(lambda x: 1 / 0.0005 * float(x) / len(avg), range(len(avg))), y=avg)

def plot_avg():
    second = get_avg('/tmp/e-02/')
    fifth = get_avg('/tmp/e-05/')


    color1 = hsv_to_rgb(0, 0.7, 0.9)
    color2 = hsv_to_rgb(0.5, 0.7, 0.9)
    fig, ax = plt.subplots()
    ax.plot(fifth['x'], fifth['y'], label='1e-05', color=color1)
    ax.plot(second['x'], second['y'], alpha=0.8, label='1e-02', color=color2)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(loc=0)
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

