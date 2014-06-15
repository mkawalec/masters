#!/usr/bin/env python2
from surv_prob_plot import gen_fits


def dump_keys():
    fits = gen_fits(0.3)
    for domain,Rs in fits.iteritems():
        with open('domains/' + str(domain), 'w') as f:
            for R in sorted(Rs.keys()):
                avg = Rs[R]['avg']
                f.write('%s\t%s\t%s\n' % (R, avg[0], avg[1]))

    

if __name__ == '__main__':
    dump_keys()
