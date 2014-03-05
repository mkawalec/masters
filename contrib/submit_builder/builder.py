#!/usr/bin/env python2


def create_submitters():
    with open('scaffold.tpl', 'r') as f:
        template = f.read()

    # R range from 0.9 to 1.09
    for R in [x * 0.01 for x in range(90, 110)]:
        with open("submit-" + str(R) + '.sh', 'w') as f:
            f.write(template)
            f.write("aprun -n 256 -N 32 /work/d54/d54/s0905879/integrator "
                    "-c decay-mult -r 40 -f 1 --fast-threshold 7000 "
                    "--find-zeros 0 --use-output 0 -R %s --prefix '%s'"
                    % (str(R), str(R) + '-'))

if __name__ == '__main__':
    create_submitters()
