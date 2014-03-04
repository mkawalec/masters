#!/bin/bash

current=`pwd`
cd /tmp
wget http://sourceforge.net/projects/arma/files/armadillo-4.100.1.tar.gz
tar -xvf armadillo-4.100.1.tar.gz
cd armadillo-4.100.1
./configure
make -j3

cp libarmadillo.so.4.100.1 $current
cp -r include arma_include

cd $current
ln -s libarmadillo.so.4.100.1 libarmadillo.so
