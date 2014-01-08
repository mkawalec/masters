#!/bin/bash

module load git
module load cmake

module load fftw/3.3.0.3
module swap PrgEnv-cray PrgEnv-gnu
module load boost

export CMAKE_PREFIX_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
export CXX=g++
export CC=gcc
