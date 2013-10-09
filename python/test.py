from math import pi, cos, pow, sqrt
from functools import reduce

def l2_norm(values):
    return sqrt(reduce(lambda acc, x: acc + pow(x, 2), values, 0))

def print_norm(samples):
    values = []
    for i in range(samples):
        x = i/samples * 24 * pi
        values.append(2.0 * cos(x) + 0.03 * cos(11.0 * x / 12.0))

    return l2_norm(values)

for i, norm in zip(range(60, 257), map(print_norm, range(60, 257))):
    print("%s %s" % (str(i), str(norm)))
