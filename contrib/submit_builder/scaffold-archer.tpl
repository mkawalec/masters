#!/bin/bash
#PBS -N turbulence
#PBS -l select=10
#PBS -l walltime=02:30:00
#PBS -A d59

export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)
cd $PBS_O_WORKDIR

