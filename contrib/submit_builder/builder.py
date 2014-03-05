#!/usr/bin/env python2

from math import pi

def create_submitters():
    with open('scaffold.tpl', 'r') as f:
        template = f.read()

    # Domain size from 24 to 2 pi
    for domain in range(24, 0, -2):
        with open("submit-" + str(domain) + 'pi.sh', 'w') as f:
            f.write(template)

            # R range from 0.9 to 1.09
            for R in [x * 0.01 for x in range(90, 110)]:
                f.write("aprun -n 256 -N 32 /work/d54/d54/s0905879/integrator "
                        "-c decay-mult -r 40 -f 1 --end-time 7000 -d %s "
                        "--find-zeros 0 --use-output 0 -R %s --prefix '%s'\n"
                        % (str(domain * pi), str(R), str(domain) + '-' + str(R) + '-'))

if __name__ == '__main__':
    create_submitters()
