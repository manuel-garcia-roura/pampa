#!/bin/bash

set -e

source setenv.sh

rm -rf lib bin
mkdir lib bin

cd src
make clean
make
cd ..

cd src/cxx
make clean
make
cd ../..

cd src/c
make clean
make
cd ../..

cd src/f90
make clean
make
cd ../..

cd bin
ln -sf pampa-f90 pampa
cd ..
