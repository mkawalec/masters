#!/bin/bash
#$ -V
#$ -cwd

mpirun -np $NSLOTS ./integrator -c decay-path -r 1 -f 0 --fast-threshold 1000 --static-interval 1000
