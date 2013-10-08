from math import pi, cos, pow, sqrt

def l2_norm(values):
    norm = 0.0
    for value in sorted(values):
        norm += pow(value, 2)

    return sqrt(norm)

def print_norm(samples):
    values = []
    for i in range(samples):
        x = i/samples * 24 * pi
        values.append(2.0 * cos(x) + 0.03 * cos(11.0 * x / 12.0))

    return l2_norm(values)

for i in range(60, 257):
    print("%s %s" % (str(i), str(print_norm(i))))
