#!/bin/bash

# Based on https://aur.archlinux.org/packages/alglib/

BUILD_DIR=/tmp/turb-build
VERSION=3.7.0
SOURCE=http://www.alglib.net/translator/re/alglib-${VERSION}.cpp.zip

DOCS_DIR=/usr/share/doc/alglib
HDR_DIR=/usr/include/libalglib
LIB_DIR=/usr/lib

mkdir -p $BUILD_DIR
cd $BUILD_DIR
wget $SOURCE
unzip alglib-${VERSION}.cpp.zip
cd cpp/src

# make static lib
gcc -I. -c *.cpp || return 1
ar rcs libalglib.a *.o || return 1

# make shared lib
rm -f *.o
gcc -I. -fPIC -c *.cpp || return 1
gcc -shared -Wl,-soname,libalglib.so.2 -o libalglib.so.${VERSION} *.o


# install docs
sudo install -d $DOCS_DIR
sudo install ../manual.cpp.html $DOCS_DIR

# install headers
sudo install -d $HDR_DIR
sudo install *.h $HDR_DIR

# install library
sudo install -d $LIB_DIR
sudo install libalglib.a $LIB_DIR
sudo install libalglib.so.${VERSION} $LIB_DIR
sudo ln -s libalglib.so.${VERSION} ${LIB_DIR}/libalglib.so.2
sudo ln -s libalglib.so.${VERSION} ${LIB_DIR}/libalglib.so
