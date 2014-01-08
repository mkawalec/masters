#!/bin/bash
#PBS -N Integrator
#PBS -l mppwidth=1
#PBS -l mppnppn=1
#PBS -l walltime=00:05:00
#PBS -A d54

aprun -n 1 -N 1 /home/d54/d54/s0905879/masters/build/integrator -c decay-path -r 1 -f 0 --fast-threshold 1000 --static-interval 1000
