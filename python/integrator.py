from __future__ import division
from numpy import fft
from math import cos, pi, pow, sqrt, floor
from datetime import datetime
from functools import reduce

def frange(stop, step):
    ''' A generator accepting any 'sane' values,
        as both range and xrange only accept integers '''
    r = 0
    while r < stop:
        yield r
        r += step

def l2_norm(values):
    return sqrt(reduce(lambda init, value: init + pow(value, 2), values, 0))

class Integrator(object):
    e = -0.1
    a = 0.125
    b = 0.004
    R = 1
    D = 40

    @staticmethod
    def initial_function(x):
        ''' The initial function, to be overloaded.

            Returns a tuple, with values (u, v) inside
        '''
        return (2.0 * cos(x) + 0.03 * cos(11.0/12.0 * x), 1.0)

    def linear_operators(self, x):
        ''' The linear operators, they would also like
            to be overloaded '''
        
        k = 2 * pi * x
        Lu = - pow(k, 4) - pow(k, 2) - (1 - self.e)
        Lv = self.D * pow(k, 2) - 1
        print(k)

        return [(1 + 0.5 * self.dt * Lu) / (1 - 0.5 * self.dt * Lu),
                (1 + 0.5 * self.dt * Lv) / (1 - 0.5 * self.dt * Lv)]
       

    def compute_stuff(self, what, step_scale=1):
        ''' The loop returning a full list of points
            for the initial (real) conditions '''

        values = map(what, frange(self.max_x, 
            step_scale*self.max_x / self.resolution))
        return reduce(lambda init, vals: (init[0] + [vals[0]], init[1] + [vals[1]]),
                values, ([], []))

    def compute_derivative(self, init, val):
        return init + [complex(-val[1] * 2 * pi / self.max_x * val[0].imag,
                                val[1] * 2 * pi / self.max_x * val[0].real)]

    def compute_step_value(self, init, val):
        ''' Computes the numerical value of the current 
            step in the nonlinear computation and returns
            the partially reduced result 
            
            The structure of val is:
                val = (u, v, du)
            '''
        # Compute the u component
        init[0].append(self.dt * ((1 - self.e) * \
                (self.a * val[1] + self.b * pow(val[1], 2)) * \
                val[0] + val[0] * val[2]))

        # Compute the v component
        init[1].append(self.dt * self.R * pow(val[0], 2))
        return init

    def nonlinear_step(self):
        derivative = reduce(self.compute_derivative, 
                            zip(self.current_values[0], 
                                range(len(self.current_values[0]))), 
                            [])

        trans_u = fft.irfft(self.current_values[0])
        trans_v = fft.irfft(self.current_values[1])
        trans_du = fft.irfft(derivative)

        computed = reduce(self.compute_step_value, 
                zip(trans_u, trans_v, trans_du), 
                [[], []])
        transformed = (fft.rfft(computed[0]), fft.rfft(computed[1]))

        # Nullify the padding
        for i in range(int(2 / 3 * len(transformed[0])), len(transformed[0])):
            transformed[0][i] = 0
            transformed[1][i] = 0

        return transformed

    def linear_step(self):
        ''' Computes the step of linear simulation '''
        for i in range(len(self.operators[0])):
            self.current_values[0][i] *= self.operators[0][i]
            self.current_values[1][i] *= self.operators[1][i]

    def __init__(self, resolution=1024, max_x=40, dt=0.001):
        ''' Please note that the resolution MUST be even, and
            should be a power of two
        '''

        self.resolution = resolution
        self.max_x = max_x
        self.dt = dt

        self.initialized = self.compute_stuff(self.initial_function)
        self.operators = self.compute_stuff(self.linear_operators, 
                self.resolution / (self.resolution * 3 / 4))
        self.current_values = [fft.rfft(self.initialized[0], int(len(self.initialized[0]) * 3 / 2)),
                               fft.rfft(self.initialized[1], int(len(self.initialized[1]) * 3 / 2))]

    def run(self, end_time=1):
        current_time = 0

        with open('output_v', 'w') as f:
            while current_time < end_time:
                current_time += self.dt

                #self.nonlinear_step()
                self.linear_step()
                f.write('%s %s\n' % (current_time, 
                    l2_norm(self.current_values[1])))

if __name__ == '__main__':
    integrator = Integrator()
    integrator.run(2)
